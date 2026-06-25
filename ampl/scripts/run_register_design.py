import subprocess

import matplotlib.pyplot as plt
import numpy as np
from amplpy import AMPL, modules
from pulser import Register
from pulser.devices import DigitalAnalogDevice
from scipy.optimize import differential_evolution
from scipy.spatial.distance import pdist, squareform


def check_ipopt():
    try:
        ipopt_path = modules.find("ipopt")
        subprocess.run(
            [ipopt_path, "-v"],
            check=True,
            capture_output=True,
            text=True,
        )
    except Exception as e:
        print("IPOPT cannot be executed in this environment.")
        executable = ipopt_path if "ipopt_path" in locals() else "not found"
        print(f"Executable found: {executable}")
        print(f"Error: {e}")
        print(
            "On Apple Silicon Macs, this usually happens when the installed IPOPT "
            "binary is x86_64 and Rosetta 2 is not installed."
        )
        return False

    return True


def check_solver(solver):
    try:
        modules.find(solver)
    except Exception as e:
        print(f"{solver} cannot be found in the installed AMPL modules.")
        print(f"Error: {e}")
        return False

    return True


def evaluate_mapping(new_coords, Q):
    """
    Asymmetric Cost Function: Strict on large weights, forgiving on zeros.

    Regularisation: C6 / (d^6 + epsilon), matching ampl/register_design.mod.
    """
    new_coords = np.reshape(new_coords, (len(Q), 2))
    raw_dists = pdist(new_coords)
    U = squareform(DigitalAnalogDevice.interaction_coeff / (raw_dists**6 + 1e-6))

    Q_off = np.copy(Q)
    np.fill_diagonal(Q_off, 0.0)

    mask_active = Q_off > 0
    mask_zero = (Q_off == 0) & (~np.eye(len(Q), dtype=bool))

    # Standard error for active connections
    error_active = np.sum((U[mask_active] - Q_off[mask_active]) ** 2)

    # Asymmetric zero penalty: only penalize interactions greater than 0.5
    threshold = 0.5
    error_zero = np.sum(np.maximum(U[mask_zero] - threshold, 0) ** 2)

    return np.sqrt(error_active + error_zero)


def run_diff_evolution(Q, upper_bound=35.0):
    """
    Executes a high-effort global optimization using Differential Evolution
    to minimize the established evaluate_mapping cost function.
    """
    num_atoms = len(Q)

    # Define the physical bounds for the search area.
    bounds = [(0.0, upper_bound)] * (num_atoms * 2)

    print(f"Running high-effort Differential Evolution for {num_atoms} atoms...")
    print("This will take longer as we increased population size and tolerances.")

    # Deep global search configuration
    res = differential_evolution(
        evaluate_mapping,
        bounds=bounds,
        args=(Q,),
        strategy="best2bin",  # Uses two difference vectors for better diversity
        popsize=80,  # Population size.
        maxiter=100,  # Maximum number of generations allowed
        tol=1e-6,  # Tighter relative tolerance
        atol=1e-6,  # Tighter absolute tolerance
        mutation=(0.5, 1.2),
        recombination=0.8,
        polish=True,  # Ensures a local optimizer (L-BFGS-B) refines the final point
        seed=999,
        disp=True,  # Monitors the progress of each generation.
    )

    coords = np.reshape(res.x, (num_atoms, 2))
    error = evaluate_mapping(res.x, Q)

    print("\nHigh-effort Optimization finished.")
    print("Final error:", error)

    return coords


def plot_register(coords):
    # Creates the Pulser register with the optimized coordinates.
    qubits = {f"q{i}": coord for (i, coord) in enumerate(coords)}
    reg = Register(qubits)

    # Draws the physical register layout.
    reg.draw(
        blockade_radius=DigitalAnalogDevice.rydberg_blockade_radius(1.0),
        draw_graph=False,
        draw_half_radius=True,
    )

    # Renders the plot to the screen.
    plt.show()


def get_q_matrix(ampl, N):
    Q = np.zeros((N, N))
    Q_vals = ampl.get_parameter("Q").get_values().to_dict()
    for i in range(1, N + 1):
        for j in range(1, N + 1):
            Q[i - 1, j - 1] = Q_vals[i, j]

    return Q


def run_ampl_model(
    mod_file="ampl/register_design.mod",
    dat_file="data/set_partitioning/ampl/modeling_test_qubo.dat",
):
    print("Initializing the AMPL environment...")

    # Instantiate AMPL with the local modules available, including ipopt.
    ampl = AMPL()

    # 1. Loads the model and data.
    try:
        ampl.read(mod_file)
        ampl.read_data(dat_file)
    except Exception as e:
        print(f"Error while reading the files: {e}")
        return None, None

    # Overrides the hardcoded C6 in the .mod file with the authoritative value
    # from the Pulser library, keeping Python and AMPL consistent.
    ampl.get_parameter("C6").set(DigitalAnalogDevice.interaction_coeff)

    # 2. Configures the solver.
    solver = "gurobi"
    ampl.set_option("solver", solver)
    ampl.set_option("gurobi_options", "outlev=1")
    # ampl.set_option("ipopt_options", "tol=1e-6 max_iter=3000")

    if not check_solver(solver):
        return None, None

    # 3. Random initialization to avoid singularities.
    # Reads N and L from the data file.
    N = int(ampl.get_parameter("N").value())
    L = ampl.get_parameter("L").value()
    Q = get_q_matrix(ampl, N)

    # Captures the x and y decision variables from the model.
    x = ampl.get_variable("x")
    y = ampl.get_variable("y")

    # Generates an initial solution with Differential Evolution,
    # then pass it to the solver.
    try:
        initial_coords = run_diff_evolution(Q, upper_bound=L)
    except Exception as e:
        print(
            f"Error while generating the Differential Evolution initial solution: {e}"
        )
        initial_coords = None

    # Injects an initial guess for each atom.
    np.random.seed(42)  # Seed for reproducibility.
    for i in range(1, N + 1):
        if initial_coords is None:
            x[i].set_value(np.random.uniform(10.0, L - 10.0))
            y[i].set_value(np.random.uniform(10.0, L - 10.0))
        else:
            x[i].set_value(initial_coords[i - 1, 0])
            y[i].set_value(initial_coords[i - 1, 1])

    print(f"Solving the model for {N} atoms with {solver.upper()}...")

    # 4. Runs the optimization.
    try:
        ampl.solve()
    except Exception as e:
        print(f"Error while solving the model: {e}")
        return None, None

    # 5. Extracts and analyzes the results.
    solve_result = ampl.get_value("solve_result")

    if solve_result == "solved":
        print("\n--- Optimization Completed Successfully ---")
    else:
        print(f"\n--- Optimization Finished (Status: {solve_result}) ---")

    # Extracts the final residual error.
    final_error = ampl.get_objective("Frobenius_Error").value()
    print(f"Residual Error (Frobenius Norm): {final_error:.6f}")

    # Extracts coordinates into a NumPy array ready for Pulser.
    coords = np.zeros((N, 2))
    x_vals = x.get_values().to_dict()
    y_vals = y.get_values().to_dict()

    print("\nFinal Coordinates:")
    for i in range(1, N + 1):
        coords[i - 1, 0] = x_vals[i]
        coords[i - 1, 1] = y_vals[i]
        print(f"Atom {i}: x = {x_vals[i]:.4f}, y = {y_vals[i]:.4f}")

    mapping_error = evaluate_mapping(coords, Q)
    print(f"Evaluate Mapping Error: {mapping_error:.6f}")

    return coords, mapping_error


if __name__ == "__main__":
    final_coords, final_error = run_ampl_model()
    if final_coords is not None:
        plot_register(final_coords)

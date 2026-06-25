"""
Implementation of algorithms for mapping a QUBO matrix to physical coordinates.

The goal is to find coordinates that minimize the difference between the desired
interaction matrix Q and the physical interaction matrix induced by the register.
"""

import itertools
import math

import cma
import matplotlib.pyplot as plt
import numpy as np
from pulser import Register
from pulser.devices import DigitalAnalogDevice
from scipy.optimize import (
    basinhopping,
    differential_evolution,
    direct,
    dual_annealing,
    minimize,
    shgo,
)
from scipy.spatial.distance import pdist, squareform


def computing_interaction_matrix(coords):
    """Computes the interaction matrix for a given set of coordinates."""
    new_Q = squareform(DigitalAnalogDevice.interaction_coeff / pdist(coords) ** 6)
    return new_Q


def evaluate_mapping2(new_coords, Q):
    """
    Cost function with thresholding for zero-interaction terms.
    Reduces the geometric frustration caused by mathematically exact zeros.
    """
    new_coords = np.reshape(new_coords, (len(Q), 2))

    # Computes the physical interaction matrix U = C6 / r^6
    new_Q = squareform(DigitalAnalogDevice.interaction_coeff / pdist(new_coords) ** 6)

    # Creates a copy of Q and sets its diagonal to zero
    Q_off_diag = np.copy(Q)
    np.fill_diagonal(Q_off_diag, 0.0)

    # Creates boolean masks to separate the matrix entries
    mask_non_zero = Q_off_diag > 0
    # The mask for zeros must exclude the diagonal (which is naturally zero)
    mask_zero = (Q_off_diag == 0) & (~np.eye(len(Q), dtype=bool))

    # 1. Exact penalty for expected interactions
    error_non_zero = new_Q[mask_non_zero] - Q_off_diag[mask_non_zero]
    cost_non_zero = np.sum(error_non_zero**2)

    # 2. Soft penalty for zero interactions
    # Epsilon defines the acceptable threshold for physical interaction
    epsilon = 1.0

    # Only penalizes if the physical interaction exceeds epsilon
    # np.maximum guarantees that if interaction is below epsilon, error is 0
    error_zero = np.maximum(new_Q[mask_zero] - epsilon, 0.0)
    cost_zero = np.sum(error_zero**2)

    # Returns the square root of the total squared errors (similar to Frobenius norm)
    return np.sqrt(cost_non_zero + cost_zero)


def evaluate_mapping_baseline(new_coords, Q):
    """Cost function to minimize. Ideally, the pairwise distances are conserved."""
    new_coords = np.reshape(new_coords, (len(Q), 2))
    # Computes the matrix of the distances between all coordinate pairs
    new_Q = squareform(DigitalAnalogDevice.interaction_coeff / pdist(new_coords) ** 6)

    # Creates a copy of Q and sets its diagonal to zero
    # This prevents the solver from optimizing linear terms via coordinates.
    Q_off_diag = np.copy(Q)
    np.fill_diagonal(Q_off_diag, 0.0)

    return np.linalg.norm(new_Q - Q_off_diag)


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


def plot_register(coords):
    # Creates the Pulser register with the optimized coordinates
    qubits = {f"q{i}": coord for (i, coord) in enumerate(coords)}
    reg = Register(qubits)

    # Draws the physical register layout
    reg.draw(
        blockade_radius=DigitalAnalogDevice.rydberg_blockade_radius(1.0),
        draw_graph=False,
        draw_half_radius=True,
    )

    # Renders the plot to the screen
    plt.show()


def print_results(coords, Q):
    """Prints the results of the optimization."""

    error = evaluate_mapping_baseline(coords, Q)

    # Displays the raw coordinate results
    new_Q = computing_interaction_matrix(coords)
    print("Optimized coordinates:")
    print(coords)
    print("\nOriginal QUBO matrix Q:")
    print(Q)
    print("\nComputed interaction matrix from optimized coordinates:")
    print(new_Q)
    print("Error:")
    print(error)


def load_q_matrix(path):
    return np.loadtxt(path, delimiter=",")


def run_baseline(Q):
    np.random.seed(0)

    # Generates an initial random guess for the 2D coordinates
    x0 = np.random.random(len(Q) * 2)

    # Runs the Nelder-Mead optimization algorithm
    res = minimize(
        evaluate_mapping_baseline,
        x0,
        args=(Q,),
        method="Nelder-Mead",
        tol=1e-9,
        options={"maxiter": 200000, "maxfev": None},
    )

    # Reshapes the flat result array back into pairs of (x, y) coordinates
    coords = np.reshape(res.x, (len(Q), 2))

    return coords


def get_bounds_arrays(bounds):
    lower_bounds = np.array([lower for lower, _ in bounds])
    upper_bounds = np.array([upper for _, upper in bounds])
    return lower_bounds, upper_bounds


def get_coordinate_limits(bounds):
    lower_bounds, upper_bounds = get_bounds_arrays(bounds)
    return float(np.min(lower_bounds)), float(np.max(upper_bounds))


def run_cma_es(Q, bounds):
    """
    Executes CMA-ES, the state-of-the-art evolutionary algorithm.
    """
    num_atoms = len(Q)
    lower_bounds, upper_bounds = get_bounds_arrays(bounds)

    # CMA-ES needs a starting point, even if random, and a sigma.
    x0 = np.random.uniform(lower_bounds, upper_bounds)

    # Initial sigma (standard deviation): how far it should search at the start.
    sigma0 = float(np.max(upper_bounds - lower_bounds))

    print(f"Running CMA-ES for {num_atoms} atoms...")

    # Physical bounds and seed configuration.
    opts = cma.CMAOptions()
    opts["bounds"] = [lower_bounds.tolist(), upper_bounds.tolist()]
    opts["seed"] = 42
    opts["popsize"] = 500  # Population size.
    # opts['maxiter'] = 5000

    # Runs the algorithm using the unchanged evaluate_mapping function.
    res = cma.fmin(lambda x: evaluate_mapping(x, Q), x0, sigma0, options=opts)

    # res[0] contains the best coordinate vector found.
    coords = np.reshape(res[0], (num_atoms, 2))
    error = evaluate_mapping(res[0], Q)

    print("\nCMA-ES Optimization finished.")
    print("Final error:", error)

    return coords


def run_dual_annealing(Q, bounds):
    """
    Executes a global optimization using Dual Annealing.
    Highly effective for physics-inspired landscapes with deep local minima.
    """
    num_atoms = len(Q)

    print(f"Running Dual Annealing for {num_atoms} atoms...")

    # maxiter in dual_annealing is the number of global search iterations.
    # It runs a local search (like L-BFGS-B) internally automatically.
    res = dual_annealing(
        evaluate_mapping,
        bounds=bounds,
        args=(Q,),
        # maxiter=50000,      # Number of global search iterations.
        seed=24,
    )

    coords = np.reshape(res.x, (num_atoms, 2))
    error = evaluate_mapping(res.x, Q)

    print("\nDual Annealing Optimization finished.")
    print("Final error:", error)

    return coords


def run_diff_evolution(Q, bounds):
    """
    Executes a high-effort global optimization using Differential Evolution
    to minimize the established evaluate_mapping cost function.
    """
    num_atoms = len(Q)

    print(f"Running high-effort Differential Evolution for {num_atoms} atoms...")
    print("This will take longer as we increased population size and tolerances.")

    # Deep global search configuration
    res = differential_evolution(
        evaluate_mapping,
        bounds=bounds,
        args=(Q,),
        strategy="best2bin",  # Uses two difference vectors for better diversity
        popsize=25,  # Population size.
        maxiter=500,  # Maximum number of generations allowed
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


def run_basin_hopping(Q, bounds):
    """
    Executes Basin-Hopping: a continuous space meta-heuristic that
    mimics the local-escape philosophy of Tabu Search.
    """
    num_atoms = len(Q)
    lower_bounds, upper_bounds = get_bounds_arrays(bounds)

    # Basin-Hopping needs an initial starting point.
    x0 = np.random.uniform(lower_bounds, upper_bounds)

    # Options passed to the local minimizer at the bottom of each basin.
    minimizer_kwargs = {"method": "L-BFGS-B", "bounds": bounds, "args": (Q,)}

    print(f"Running Basin-Hopping for {num_atoms} atoms...")

    # maxiter controls how many global jumps the algorithm will take.
    # T is the acceptance temperature: higher values temporarily accept worse basins.
    # stepsize is the jump size in coordinate space.
    res = basinhopping(
        evaluate_mapping,
        x0,
        minimizer_kwargs=minimizer_kwargs,
        niter=100,  # Number of global jumps.
        T=1.0,  # Acceptance temperature.
        stepsize=15.0,  # Jump up to 15 micrometers away to escape the local basin.
        seed=42,
        disp=False,
    )

    coords = np.reshape(res.x, (num_atoms, 2))
    error = evaluate_mapping(res.x, Q)

    print("\nBasin-Hopping Optimization finished.")
    print("Final error:", error)

    return coords


def run_shgo(Q, bounds):
    """
    Executes Simplicial Homology Global Optimization (SHGO).
    Uses mathematical topology to guarantee finding local minima within the bounds.
    """
    num_atoms = len(Q)

    print(f"Running SHGO (Topological Search) for {num_atoms} atoms...")
    print("Note: SHGO can be computationally heavy for high dimensions.")

    # n: Sampling points generated per dimension in the initial complex.
    # iters: Number of mesh refinement iterations.
    res = shgo(
        evaluate_mapping,
        bounds=bounds,
        args=(Q,),
        # Initial mesh points.
        n=1000,
        iters=3,  # How many times the mesh is subdivided.
        options={"disp": True},  # Show progress in the terminal.
    )

    coords = np.reshape(res.x, (num_atoms, 2))
    error = evaluate_mapping(res.x, Q)

    print("\nSHGO Optimization finished.")
    print("Final error:", error)

    return coords


def run_direct(Q, bounds):
    """
    Executes the DIRECT algorithm.
    A deterministic global search that recursively divides the search space.
    """
    num_atoms = len(Q)

    print(f"Running DIRECT algorithm for {num_atoms} atoms...")
    print(
        "Note: DIRECT is highly exhaustive and relies on maxfun/maxiter to terminate."
    )

    # maxiter: Limit on rectangle divisions.
    # vol_tol: Volume tolerance; stop slicing once a rectangle becomes tiny.
    res = direct(
        evaluate_mapping,
        bounds=bounds,
        args=(Q,),
        # maxiter=1000,     # Maximum number of algorithm iterations.
        # maxfun=100000,    # Maximum calls to evaluate_mapping for runtime safety.
        # vol_tol=1e-100     # Minimum rectangle size before it is ignored.
    )

    coords = np.reshape(res.x, (num_atoms, 2))
    error = evaluate_mapping(res.x, Q)

    print("\nDIRECT Optimization finished.")
    print("Final error:", error)
    print("Status message:", res.message)

    return coords


def run_grid_multistart(Q, bounds, grid_size=3):
    N = len(Q)
    lower_bound, upper_bound = get_coordinate_limits(bounds)
    # Defines the grid points.
    grid_points = np.linspace(lower_bound, upper_bound, grid_size)
    # Generates all grid positions (x, y) for a single atom.
    single_atom_positions = list(itertools.product(grid_points, repeat=2))
    total_combinations = math.comb(len(single_atom_positions), N)
    print(f"Total possible combinations: {total_combinations}", flush=True)

    best_overall_error = float("inf")
    best_overall_coords = None
    progress_step_percent = 10
    next_progress_count = math.ceil(total_combinations * progress_step_percent / 100)

    # Iterates over initial combinations.
    for combination_index, start_config in enumerate(
        itertools.combinations(single_atom_positions, N),
        start=1,
    ):
        while progress_step_percent <= 100 and combination_index >= next_progress_count:
            print(f"{progress_step_percent}% of combinations reached.", flush=True)
            progress_step_percent += 10
            next_progress_count = math.ceil(
                total_combinations * progress_step_percent / 100
            )

        x0 = np.array(start_config).flatten()

        # Local polishing with L-BFGS-B.
        res = minimize(
            fun=lambda x: evaluate_mapping(x, Q),
            x0=x0,
            method="L-BFGS-B",
            bounds=bounds,
        )

        if res.fun < best_overall_error:
            print(f"New best error found: {res.fun}")
            best_overall_error = res.fun
            best_overall_coords = np.reshape(res.x, (N, 2))

    return best_overall_coords


def run_random_grid_multistart(Q, bounds, grid_size=3, num_starts=5000, seed=42):
    N = len(Q)
    rng = np.random.default_rng(seed)
    lower_bound, upper_bound = get_coordinate_limits(bounds)

    grid_points = np.linspace(lower_bound, upper_bound, grid_size)
    single_atom_positions = np.array(list(itertools.product(grid_points, repeat=2)))

    total_possible = math.comb(len(single_atom_positions), N)
    print(f"Total possible combinations: {total_possible}", flush=True)
    print(f"Randomly testing {num_starts} initial configurations.", flush=True)

    best_overall_error = float("inf")
    best_overall_coords = None

    for start_idx in range(1, num_starts + 1):
        selected_indices = rng.choice(
            len(single_atom_positions),
            size=N,
            replace=False,
        )
        x0 = single_atom_positions[selected_indices].flatten()

        res = minimize(
            fun=lambda x: evaluate_mapping(x, Q),
            x0=x0,
            method="L-BFGS-B",
            bounds=bounds,
        )

        if res.fun < best_overall_error:
            print(f"New best error found: {res.fun}")
            best_overall_error = res.fun
            best_overall_coords = np.reshape(res.x, (N, 2))

        if start_idx % max(1, num_starts // 10) == 0:
            progress = 100 * start_idx / num_starts
            print(f"{progress:.1f}% of random starts completed.", flush=True)

    return best_overall_coords


def main():
    q_matrix_path = "data/set_partitioning/csv/modeling_test_qubo.csv"
    # q_matrix_path = "data/set_partitioning/csv/16x6_qubo.csv"
    Q = load_q_matrix(q_matrix_path)
    bounds = [(0.0, 35.0)] * (len(Q) * 2)

    # coords = run_cma_es(Q, bounds)
    # coords = run_dual_annealing(Q, bounds)
    # coords = run_diff_evolution(Q, bounds)
    # coords = run_baseline(Q)
    # coords = run_basin_hopping(Q, bounds)
    # coords = run_shgo(Q, bounds) # Does not work
    # coords = run_direct(Q, bounds)
    # coords = run_grid_multistart(Q, bounds, grid_size=3)
    coords = run_random_grid_multistart(Q, bounds, grid_size=6, num_starts=1000)
    print_results(coords, Q)
    plot_register(coords)


if __name__ == "__main__":
    main()

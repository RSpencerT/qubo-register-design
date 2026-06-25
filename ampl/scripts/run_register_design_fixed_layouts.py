import matplotlib.pyplot as plt
import numpy as np
from amplpy import AMPL, modules
from pulser import Register
from pulser.devices import DigitalAnalogDevice


def check_solver(solver):
    try:
        modules.find(solver)
    except Exception as e:
        print(f"{solver} cannot be found in the installed AMPL modules.")
        print(f"Error: {e}")
        return False

    return True


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


def run_register_design_fixed_layouts_model(
    mod_file="ampl/register_design_fixed_layouts.mod",
    dat_file="data/set_partitioning/ampl/modeling_test_qubo-fixed.dat",
    solver="gurobi",
):
    print("Initializing the AMPL environment...")

    ampl = AMPL()

    try:
        ampl.read(mod_file)
        ampl.read_data(dat_file)
    except Exception as e:
        print(f"Error while reading the files: {e}")
        return None, None

    if not check_solver(solver):
        return None, None

    ampl.set_option("solver", solver)
    if solver == "gurobi":
        ampl.set_option("gurobi_options", "outlev=1")

    print(f"Solving the fixed-layout register design model with {solver.upper()}...")

    try:
        ampl.solve()
    except Exception as e:
        print(f"Error while solving the model: {e}")
        return None, None

    solve_result = ampl.get_value("solve_result")
    print(f"\n--- Optimization Finished (Status: {solve_result}) ---")

    if solve_result != "solved":
        print("No feasible solution was found.")
        return None, None

    squared_error = ampl.get_objective("Squared_Frobenius_Error").value()
    print(f"Squared Frobenius Error: {squared_error:g}")

    use_layout_values = ampl.get_variable("use_layout").get_values().to_dict()
    selected_layouts = [
        layout for layout, value in use_layout_values.items() if round(value) == 1
    ]

    selected_layout = selected_layouts[0] if selected_layouts else None
    print(f"\nSelected Layout: {selected_layout}")

    assign_values = ampl.get_variable("assign").get_values().to_dict()
    assignments = []
    for (layout, atom, site), value in assign_values.items():
        if round(value) == 1:
            assignments.append((atom, layout, site))

    assignments.sort(key=lambda item: item[0])

    print("\nAtom Assignments:")
    for atom, layout, site in assignments:
        print(f"Atom {atom} -> Layout {layout}, Site {site}")

    site_x = ampl.get_parameter("Site_X").get_values().to_dict()
    site_y = ampl.get_parameter("Site_Y").get_values().to_dict()
    coords = np.zeros((len(assignments), 2))

    for atom, layout, site in assignments:
        atom_index = int(atom) - 1
        coords[atom_index, 0] = site_x[layout, site]
        coords[atom_index, 1] = site_y[layout, site]

    plot_register(coords)

    return assignments, squared_error


if __name__ == "__main__":
    run_register_design_fixed_layouts_model()

from amplpy import AMPL, modules


def check_solver(solver):
    try:
        modules.find(solver)
    except Exception as e:
        print(f"{solver} cannot be found in the installed AMPL modules.")
        print(f"Error: {e}")
        return False

    return True


def run_set_partitioning_model(
    mod_file="ampl/set_partitioning.mod",
    # dat_file="data/set_partitioning/ampl/set_partitioning.dat",
    # dat_file="data/set_partitioning/ampl/6x6.dat",
    # dat_file="data/set_partitioning/ampl/16x6.dat",
    dat_file="data/set_partitioning/ampl/modeling_test_recover.dat",
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

    print(f"Solving the set partitioning model with {solver.upper()}...")

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

    total_cost = ampl.get_objective("Total_Cost").value()
    print(f"Total Cost: {total_cost:g}")

    x_values = ampl.get_variable("x").get_values().to_dict()
    selected_subsets = [
        subset for subset, value in x_values.items() if round(value) == 1
    ]

    print("\nSelected Subsets:")
    if selected_subsets:
        for subset in selected_subsets:
            print(f"{subset}: x = 1")
    else:
        print("No subsets selected.")

    return selected_subsets, total_cost


if __name__ == "__main__":
    run_set_partitioning_model()

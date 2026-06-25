from amplpy import AMPL, modules


def check_solver(solver):
    try:
        modules.find(solver)
    except Exception as e:
        print(f"{solver} cannot be found in the installed AMPL modules.")
        print(f"Error: {e}")
        return False

    return True


def run_qubo_model(
    mod_file="ampl/qubo.mod",
    dat_file="data/set_partitioning/ampl/modeling_test_qubo.dat",
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

    print(f"Solving the QUBO model with {solver.upper()}...")

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

    objective_value = ampl.get_objective("QUBO_Objective").value()
    print(f"QUBO Objective: {objective_value:g}")

    x_values = ampl.get_variable("x").get_values().to_dict()
    selected_vars = [index for index, value in x_values.items() if round(value) == 1]

    print("\nSolution:")
    for index in sorted(x_values):
        print(f"x[{index}] = {round(x_values[index])}")

    print("\nSelected Variables:")
    if selected_vars:
        print(", ".join([f"x[{index}]" for index in selected_vars]))
    else:
        print("No variables selected.")

    return selected_vars, objective_value


if __name__ == "__main__":
    run_qubo_model()

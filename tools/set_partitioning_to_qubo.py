import argparse

import numpy as np


def parse_costs(lines, subset_names):
    costs = {subset_name: 1.0 for subset_name in subset_names}
    in_cost = False

    for line in lines:
        line = line.strip()
        if not line:
            continue

        if line.startswith("param Cost default"):
            parts = line.replace(";", "").split()
            default_cost = float(parts[-1])
            costs = {subset_name: default_cost for subset_name in subset_names}
            continue

        if line.startswith("param Cost :="):
            in_cost = True
            continue

        if in_cost:
            if line == ";" or line.startswith(";"):
                break
            parts = line.replace(";", "").split()
            if len(parts) >= 2:
                costs[parts[0]] = float(parts[1])

    return np.array([costs[subset_name] for subset_name in subset_names])


def format_number(value):
    if float(value).is_integer():
        return f"{int(value):2d}"
    return f"{value:g}"


def set_partitioning_dat_to_qubo(
    dat_filename="data/set_partitioning/ampl/set_partitioning.dat",
    penalty_factor=10,
    output_filename="data/set_partitioning/ampl/modeling_test_qubo.dat",
):
    """
    Reads an AMPL .dat incidence matrix and writes the corresponding QUBO .dat.
    """
    try:
        with open(dat_filename, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: file '{dat_filename}' was not found.")
        return

    matrix_lines = []
    subset_names = []
    in_matrix = False

    # 1. Parses the AMPL .dat file.
    for line in lines:
        line = line.strip()
        if line.startswith("param A :"):
            subset_names = line.replace(":=", "").split()[3:]
            in_matrix = True
            continue

        if in_matrix:
            if line == ";" or line.startswith(";"):
                break
            if not line:
                continue

            parts = line.split()[1:]
            matrix_lines.append([int(p) for p in parts])

    A = np.array(matrix_lines, dtype=int)
    num_items, num_subsets = A.shape
    costs = parse_costs(lines, subset_names)

    print(f"1. Incidence matrix A ({num_items}x{num_subsets}) read successfully.")

    # 2. QUBO algebra.
    Q = np.dot(A.T, A) * penalty_factor
    diagonal = costs - penalty_factor * A.sum(axis=0)

    # 3. Reconstructs the linear terms on the diagonal.
    np.fill_diagonal(Q, diagonal)

    print("\n2. Reconstructed QUBO matrix:")
    print(Q)

    # 4. Exports to the .dat format expected by register_design.mod.
    with open(output_filename, "w") as f:
        f.write(f"param N := {num_subsets};\n\n")

        headers = "  ".join([str(i + 1) for i in range(num_subsets)])
        f.write(f"param Q :  {headers} :=\n")

        for i in range(num_subsets):
            row_str = " ".join([format_number(val) for val in Q[i]])
            f.write(f"       {i + 1}   {row_str}\n")
        f.write(";\n")

    print(
        f"\n3. Success! File '{output_filename}' was generated and is ready for AMPL."
    )
    return Q


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert an AMPL Set Partitioning .dat instance into a QUBO .dat."
    )
    parser.add_argument(
        "--input",
        default="data/set_partitioning/ampl/modeling_test_recover.dat",
        help="Input Set Partitioning .dat file.",
    )
    parser.add_argument(
        "--output",
        default="data/set_partitioning/ampl/modeling_test_qubo.dat",
        help="Output QUBO .dat file.",
    )
    parser.add_argument(
        "--penalty-factor",
        type=float,
        default=10.0,
        help="Penalty factor for equality constraints.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    set_partitioning_dat_to_qubo(
        dat_filename=args.input,
        penalty_factor=args.penalty_factor,
        output_filename=args.output,
    )

import argparse
import time

import numpy as np


def write_set_partitioning_dat(A, filename, costs=None):
    num_items, num_subsets = A.shape
    if costs is None:
        costs = np.ones(num_subsets, dtype=float)
    costs = np.array(costs, dtype=float)

    if len(costs) != num_subsets:
        raise ValueError("costs must have one value per subset.")

    with open(filename, "w") as f:
        f.write("# =========================================================\n")
        f.write("# Reconstructed Set Partitioning Instance\n")
        f.write("# =========================================================\n\n")

        items_list = " ".join([f"I{i + 1}" for i in range(num_items)])
        subsets_list = " ".join([f"S{j + 1}" for j in range(num_subsets)])

        f.write(f"set ITEMS := {items_list};\n")
        f.write(f"set SUBSETS := {subsets_list};\n\n")

        f.write("param Cost :=\n")
        for j in range(num_subsets):
            f.write(f"     S{j + 1}    {costs[j]:g}\n")
        f.write(";\n\n")

        headers = "  ".join([f"S{j + 1}" for j in range(num_subsets)])
        f.write(f"param A : {headers} :=\n")

        for i in range(num_items):
            row_values = "  ".join(map(str, A[i]))
            f.write(f"     I{i + 1}    {row_values}\n")

        f.write(";\n")


def find_pairwise_exact_cover(C):
    """
    Finds a feasible exact cover for the pairwise reconstruction, if one exists.

    The reconstruction creates one item for each positive pair (i, j), so every
    row enforces x[i] + x[j] = 1. This is feasible exactly when the graph induced
    by positive off-diagonal entries is bipartite.
    """
    N = len(C)
    colors = np.full(N, -1, dtype=int)

    for start in range(N):
        if colors[start] != -1:
            continue

        colors[start] = 0
        stack = [start]

        while stack:
            i = stack.pop()
            neighbors = np.flatnonzero(C[i] > 0)
            for j in neighbors:
                if i == j:
                    continue
                if colors[j] == -1:
                    colors[j] = 1 - colors[i]
                    stack.append(j)
                elif colors[j] == colors[i]:
                    return None

    return colors


def coverage_error(A, x):
    coverage = np.dot(A, x)
    return int(np.sum((coverage - 1) ** 2))


def simulated_annealing_feasible_incidence(
    A0,
    x0=None,
    max_iterations=50000,
    initial_temperature=5.0,
    cooling_rate=0.9995,
    seed=42,
):
    """
    Repairs a binary incidence matrix until it has a feasible exact-cover vector.

    The state is the pair (A, x), where A is the incidence matrix and x is a
    binary subset-selection vector. The search flips entries in A and x until it
    reaches A @ x == 1 for every item.
    """
    rng = np.random.default_rng(seed)
    A = np.array(A0, dtype=int).copy()

    if A.ndim != 2:
        raise ValueError("A0 must be a 2-dimensional incidence matrix.")
    if not np.all((A == 0) | (A == 1)):
        raise ValueError("A0 must be binary.")

    num_items, num_subsets = A.shape
    if x0 is None:
        x = rng.integers(0, 2, size=num_subsets, dtype=int)
        if not np.any(x):
            x[rng.integers(0, num_subsets)] = 1
    else:
        x = np.array(x0, dtype=int).copy()
        if x.shape != (num_subsets,):
            raise ValueError("x0 must have one entry per subset.")
        if not np.all((x == 0) | (x == 1)):
            raise ValueError("x0 must be binary.")

    current_error = coverage_error(A, x)
    best_A = A.copy()
    best_x = x.copy()
    best_error = current_error
    temperature = initial_temperature

    for iteration in range(1, max_iterations + 1):
        if current_error == 0:
            print(f"Feasible incidence matrix found at iteration {iteration}.")
            return A, x

        candidate_A = A.copy()
        candidate_x = x.copy()

        if rng.random() < 0.8:
            i = rng.integers(0, num_items)
            j = rng.integers(0, num_subsets)
            candidate_A[i, j] = 1 - candidate_A[i, j]
        else:
            j = rng.integers(0, num_subsets)
            candidate_x[j] = 1 - candidate_x[j]
            if not np.any(candidate_x):
                candidate_x[j] = 1

        candidate_error = coverage_error(candidate_A, candidate_x)
        delta = candidate_error - current_error

        if delta <= 0 or rng.random() < np.exp(-delta / max(temperature, 1e-12)):
            A = candidate_A
            x = candidate_x
            current_error = candidate_error

            if current_error < best_error:
                best_A = A.copy()
                best_x = x.copy()
                best_error = current_error

        temperature *= cooling_rate

    if best_error == 0:
        return best_A, best_x

    raise ValueError(
        "Simulated annealing did not find a feasible set-partitioning instance. "
        f"Best coverage error: {best_error}."
    )


def binary_matrix_to_feasible_set_partitioning_dat(
    A0,
    filename="data/set_partitioning/ampl/annealed_set_partitioning.dat",
    costs=None,
    x0=None,
    max_iterations=50000,
    initial_temperature=5.0,
    cooling_rate=0.9995,
    seed=42,
):
    """
    Repairs a binary incidence matrix with simulated annealing and writes a
    feasible set-partitioning AMPL .dat file.
    """
    A, x = simulated_annealing_feasible_incidence(
        A0,
        x0=x0,
        max_iterations=max_iterations,
        initial_temperature=initial_temperature,
        cooling_rate=cooling_rate,
        seed=seed,
    )
    write_set_partitioning_dat(A, filename, costs=costs)

    selected = [f"S{i + 1}" for i, value in enumerate(x) if value == 1]
    print(f"Feasible exact-cover assignment: {selected}")
    print(f"File '{filename}' generated with {len(A)} set-partitioning items.")
    return A, x


def validate_qubo_matrix(Q_matrix, penalty_factor):
    Q = np.array(Q_matrix, dtype=float)
    N = len(Q)

    if Q.shape != (N, N):
        raise ValueError("Q_matrix must be square.")
    if not np.allclose(Q, Q.T):
        raise ValueError("Q_matrix must be symmetric.")

    C = np.rint(Q / penalty_factor).astype(int)
    np.fill_diagonal(C, 0)

    off_diagonal_mask = ~np.eye(N, dtype=bool)
    if not np.allclose(Q[off_diagonal_mask], (C * penalty_factor)[off_diagonal_mask]):
        raise ValueError(
            "Q_matrix off-diagonal entries must be multiples of penalty_factor."
        )
    if np.any(C < 0):
        raise ValueError("Q_matrix off-diagonal entries must be nonnegative.")

    return Q, C


def pairwise_incidence_from_intersections(C):
    rows = []
    N = len(C)

    for i in range(N):
        for j in range(i + 1, N):
            for _ in range(C[i, j]):
                row = np.zeros(N, dtype=int)
                row[i] = 1
                row[j] = 1
                rows.append(row)

    A = np.array(rows, dtype=int)
    if A.size == 0:
        A = np.zeros((0, N), dtype=int)

    return A


def reconstructed_qubo_from_incidence(A, diagonal, penalty_factor):
    Q = np.dot(A.T, A) * penalty_factor
    np.fill_diagonal(Q, diagonal)
    return Q


def off_diagonal_overlap_error(A, target_intersections):
    overlaps = np.dot(A.T, A)
    mask = ~np.eye(len(target_intersections), dtype=bool)
    return float(np.sum((overlaps[mask] - target_intersections[mask]) ** 2))


def anneal_feasible_incidence_near_intersections(
    A0,
    target_intersections,
    x0=None,
    max_iterations=200000,
    initial_temperature=25.0,
    cooling_rate=0.99995,
    feasibility_weight=10000.0,
    seed=42,
):
    """
    Repairs A0 while keeping A.T @ A close to target_intersections.
    """
    rng = np.random.default_rng(seed)
    A = np.array(A0, dtype=int).copy()
    target_intersections = np.array(target_intersections, dtype=int)

    if A.ndim != 2:
        raise ValueError("A0 must be a 2-dimensional incidence matrix.")
    if not np.all((A == 0) | (A == 1)):
        raise ValueError("A0 must be binary.")

    num_items, num_subsets = A.shape
    if x0 is None:
        x = rng.integers(0, 2, size=num_subsets, dtype=int)
        if not np.any(x):
            x[rng.integers(0, num_subsets)] = 1
    else:
        x = np.array(x0, dtype=int).copy()

    def score(candidate_A, candidate_x):
        feasible_part = coverage_error(candidate_A, candidate_x)
        q_part = off_diagonal_overlap_error(candidate_A, target_intersections)
        return feasibility_weight * feasible_part + q_part

    current_score = score(A, x)
    best_A = None
    best_x = None
    best_q_error = float("inf")
    temperature = initial_temperature

    for iteration in range(1, max_iterations + 1):
        candidate_A = A.copy()
        candidate_x = x.copy()

        if rng.random() < 0.9:
            i = rng.integers(0, num_items)
            j = rng.integers(0, num_subsets)
            candidate_A[i, j] = 1 - candidate_A[i, j]
        else:
            j = rng.integers(0, num_subsets)
            candidate_x[j] = 1 - candidate_x[j]
            if not np.any(candidate_x):
                candidate_x[j] = 1

        candidate_score = score(candidate_A, candidate_x)
        delta = candidate_score - current_score

        if delta <= 0 or rng.random() < np.exp(-delta / max(temperature, 1e-12)):
            A = candidate_A
            x = candidate_x
            current_score = candidate_score

            if coverage_error(A, x) == 0:
                q_error = off_diagonal_overlap_error(A, target_intersections)
                if q_error < best_q_error:
                    best_A = A.copy()
                    best_x = x.copy()
                    best_q_error = q_error
                    print(
                        f"Iteration {iteration}: feasible instance found "
                        f"with overlap error {best_q_error:g}."
                    )
                    if best_q_error == 0:
                        return best_A, best_x

        temperature *= cooling_rate

    if best_A is None:
        raise ValueError(
            "Simulated annealing did not find a feasible set-partitioning instance."
        )

    return best_A, best_x


def timed_multistart_feasible_incidence_search(
    A0,
    target_intersections,
    time_limit_seconds=300,
    max_iterations_per_start=200000,
    initial_temperature=25.0,
    cooling_rate=0.99995,
    feasibility_weight=10000.0,
    seed=42,
):
    """
    Runs repeated annealing attempts and returns the best feasible incidence found.
    """
    rng = np.random.default_rng(seed)
    deadline = time.monotonic() + time_limit_seconds
    best_A = None
    best_x = None
    best_error = float("inf")
    start_count = 0

    while time.monotonic() < deadline:
        start_count += 1
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            break

        run_seed = int(rng.integers(0, np.iinfo(np.int32).max))
        try:
            A, x = anneal_feasible_incidence_near_intersections(
                A0,
                target_intersections,
                max_iterations=max_iterations_per_start,
                initial_temperature=initial_temperature,
                cooling_rate=cooling_rate,
                feasibility_weight=feasibility_weight,
                seed=run_seed,
            )
        except ValueError:
            continue

        error = off_diagonal_overlap_error(A, target_intersections)
        if error < best_error:
            best_A = A
            best_x = x
            best_error = error
            elapsed = time_limit_seconds - max(0, deadline - time.monotonic())
            print(
                f"Multi-start {start_count}: best feasible overlap error "
                f"{best_error:g} after {elapsed:.1f}s."
            )
            if best_error == 0:
                break

    if best_A is None:
        raise ValueError(
            "No feasible set-partitioning instance found within the time limit."
        )

    return best_A, best_x


def write_qubo_dat(Q, output_filename):
    num_subsets = len(Q)

    with open(output_filename, "w") as f:
        f.write(f"param N := {num_subsets};\n\n")

        headers = "  ".join([str(i + 1) for i in range(num_subsets)])
        f.write(f"param Q :  {headers} :=\n")

        for i in range(num_subsets):
            row_str = " ".join([format_number(val) for val in Q[i]])
            f.write(f"       {i + 1}   {row_str}\n")
        f.write(";\n")


def format_number(value):
    if float(value).is_integer():
        return f"{int(value):2d}"
    return f"{value:g}"


def qubo_to_feasible_set_partitioning_dat(
    Q_matrix,
    penalty_factor=10,
    filename="data/set_partitioning/ampl/set_partitioning_from_qubo.dat",
    reconstructed_qubo_filename="data/set_partitioning/ampl/modeling_test_qubo.dat",
    x0=None,
    max_iterations=200000,
    initial_temperature=25.0,
    cooling_rate=0.99995,
    feasibility_weight=10000.0,
    time_limit_seconds=None,
    seed=42,
):
    """
    Builds a feasible set-partitioning instance from a QUBO matrix.

    This starts from the exact pairwise reconstruction and then repairs the
    incidence matrix with simulated annealing until an exact-cover assignment is
    found. If the original QUBO is not exactly representable by a feasible
    pairwise reconstruction, the repaired instance may no longer reproduce the
    original off-diagonal Q entries exactly.
    """
    Q, C = validate_qubo_matrix(Q_matrix, penalty_factor)
    A0 = pairwise_incidence_from_intersections(C)
    if time_limit_seconds is None:
        A, x = anneal_feasible_incidence_near_intersections(
            A0,
            C,
            x0=x0,
            max_iterations=max_iterations,
            initial_temperature=initial_temperature,
            cooling_rate=cooling_rate,
            feasibility_weight=feasibility_weight,
            seed=seed,
        )
    else:
        A, x = timed_multistart_feasible_incidence_search(
            A0,
            C,
            time_limit_seconds=time_limit_seconds,
            max_iterations_per_start=max_iterations,
            initial_temperature=initial_temperature,
            cooling_rate=cooling_rate,
            feasibility_weight=feasibility_weight,
            seed=seed,
        )

    costs = np.diag(Q) + penalty_factor * A.sum(axis=0)
    write_set_partitioning_dat(A, filename, costs=costs)

    repaired_Q = reconstructed_qubo_from_incidence(A, np.diag(Q), penalty_factor)
    if reconstructed_qubo_filename is not None:
        write_qubo_dat(repaired_Q, reconstructed_qubo_filename)

    off_diagonal_mask = ~np.eye(len(Q), dtype=bool)
    max_offdiag_diff = np.max(
        np.abs(repaired_Q[off_diagonal_mask] - Q[off_diagonal_mask])
    )
    frobenius_diff = np.linalg.norm(repaired_Q - Q)

    selected = [f"S{i + 1}" for i, value in enumerate(x) if value == 1]
    print(f"Feasible exact-cover assignment: {selected}")
    print(f"File '{filename}' generated with {len(A)} set-partitioning items.")
    print(f"Max off-diagonal Q difference after repair: {max_offdiag_diff:g}")
    print(f"Frobenius Q difference after repair: {frobenius_diff:g}")
    return A, x, repaired_Q


def qubo_to_set_partitioning_dat(
    Q_matrix,
    penalty_factor=10,
    filename="data/set_partitioning/ampl/set_partitioning.dat",
    require_feasible=True,
):
    """
    Reconstructs a set-partitioning incidence matrix from a QUBO matrix and
    writes it as an AMPL .dat file.

    This is exact for the matrix generated by set_partitioning_dat_to_qubo:

        Q_offdiag = penalty_factor * (A.T @ A)_offdiag
        Q_diag[j] = Cost[j] - penalty_factor * column_sum[j]

    The diagonal is represented through the subset Cost parameter.

    Important: this reverse direction is not unique. The pairwise reconstruction
    used here is feasible as a set-partitioning instance only when the graph of
    positive QUBO interactions is bipartite. Set require_feasible=False to write
    the Q-preserving reconstruction even when it has no exact-cover solution.
    """
    Q, C = validate_qubo_matrix(Q_matrix, penalty_factor)

    feasible_solution = find_pairwise_exact_cover(C)
    if feasible_solution is None and require_feasible:
        raise ValueError(
            "This QUBO cannot be converted into a feasible set-partitioning "
            "instance with the pairwise reconstruction. The graph induced by "
            "positive off-diagonal Q entries is not bipartite. Use a source "
            "set-partitioning .dat file or set require_feasible=False if you "
            "only need a Q-preserving, possibly infeasible reconstruction."
        )

    A = pairwise_incidence_from_intersections(C)

    column_sums = A.sum(axis=0)
    costs = np.diag(Q) + penalty_factor * column_sums

    write_set_partitioning_dat(A, filename, costs=costs)
    print(f"File '{filename}' generated with {len(A)} set-partitioning items.")
    if feasible_solution is not None:
        selected = [
            f"S{i + 1}" for i, value in enumerate(feasible_solution) if value == 1
        ]
        print(f"Feasible exact-cover assignment for reconstructed instance: {selected}")
    return A


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a QUBO CSV into a feasible Set Partitioning .dat."
    )
    parser.add_argument(
        "--input",
        default="data/set_partitioning/csv/modeling_test_qubo.csv",
        help="Input QUBO matrix in CSV format.",
    )
    parser.add_argument(
        "--output",
        default="data/set_partitioning/ampl/set_partitioning_from_qubo.dat",
        help="Output Set Partitioning .dat file.",
    )
    parser.add_argument(
        "--output-qubo",
        default="data/set_partitioning/ampl/modeling_test_qubo.dat",
        help="Output repaired QUBO .dat file.",
    )
    parser.add_argument(
        "--penalty-factor",
        type=float,
        default=10.0,
        help="Penalty factor for equality constraints.",
    )
    parser.add_argument(
        "--time-limit",
        type=int,
        default=1600,
        help="Time limit for simulated annealing in seconds.",
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=20000000,
        help="Max iterations for simulated annealing.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    Q_input = np.loadtxt(args.input, delimiter=",")

    qubo_to_feasible_set_partitioning_dat(
        Q_input,
        penalty_factor=args.penalty_factor,
        filename=args.output,
        reconstructed_qubo_filename=args.output_qubo,
        time_limit_seconds=args.time_limit,
        max_iterations=args.max_iterations,
        seed=args.seed,
    )

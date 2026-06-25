import argparse
import math
import random
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class FixedLayoutsInstance:
    Q: np.ndarray
    layouts: list
    sites: list
    available: dict
    site_x: dict
    site_y: dict
    interaction: dict


@dataclass
class HeuristicSolution:
    layout: str
    assignment: dict
    squared_error: float


def _read_block(text, pattern, name):
    match = re.search(pattern, text, flags=re.DOTALL)
    if not match:
        raise ValueError(f"Could not find {name} in the data file.")
    return match.group(1)


def _read_groups(text, pattern, name):
    match = re.search(pattern, text, flags=re.DOTALL)
    if not match:
        raise ValueError(f"Could not find {name} in the data file.")
    return match.groups()


def _parse_set(text, name):
    block = _read_block(
        text,
        rf"set\s+{re.escape(name)}\s*:=\s*(.*?)\s*;",
        f"set {name}",
    )
    return block.split()


def _parse_table_param(text, name):
    header_text, body_text = _read_groups(
        text,
        rf"param\s+{re.escape(name)}\s*:\s*(.*?)\s*:=\s*(.*?)\s*;",
        f"parameter {name}",
    )
    columns = header_text.split()
    table = {}

    for raw_line in body_text.strip().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        parts = line.split()
        row = parts[0]
        values = parts[1:]
        if len(values) != len(columns):
            raise ValueError(
                f"Parameter {name} row {row} has {len(values)} values, "
                f"but {len(columns)} were expected."
            )

        for column, value in zip(columns, values):
            table[row, column] = float(value)

    return columns, table


def _parse_interaction(text):
    body = _read_block(
        text,
        r"param\s+Interaction\s*:=\s*(.*?)\s*;",
        "parameter Interaction",
    )
    interaction = {}

    for raw_line in body.strip().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        layout, site_i, site_j, value = line.split()
        interaction[layout, site_i, site_j] = float(value)

    return interaction


def load_fixed_layouts_instance(dat_file):
    text = Path(dat_file).read_text()

    n_match = re.search(r"param\s+N\s*:=\s*(\d+)\s*;", text)
    if not n_match:
        raise ValueError("Could not find parameter N in the data file.")
    num_atoms = int(n_match.group(1))

    layouts = _parse_set(text, "LAYOUTS")
    sites = _parse_set(text, "SITES")

    q_columns, q_table = _parse_table_param(text, "Q")
    if len(q_columns) != num_atoms:
        raise ValueError("The Q table does not match parameter N.")

    Q = np.zeros((num_atoms, num_atoms))
    for i in range(1, num_atoms + 1):
        for j in range(1, num_atoms + 1):
            Q[i - 1, j - 1] = q_table[str(i), str(j)]

    _, available = _parse_table_param(text, "Available")
    _, site_x = _parse_table_param(text, "Site_X")
    _, site_y = _parse_table_param(text, "Site_Y")
    interaction = _parse_interaction(text)

    return FixedLayoutsInstance(
        Q=Q,
        layouts=layouts,
        sites=sites,
        available=available,
        site_x=site_x,
        site_y=site_y,
        interaction=interaction,
    )


def available_sites(instance, layout):
    return [
        site
        for site in instance.sites
        if round(instance.available.get((layout, site), 0.0)) == 1
    ]


def evaluate_assignment(instance, layout, assignment):
    Q = instance.Q
    num_atoms = len(Q)
    squared_error = 0.0

    for i in range(num_atoms):
        atom_i = i + 1
        site_i = assignment[atom_i]
        for j in range(i + 1, num_atoms):
            atom_j = j + 1
            site_j = assignment[atom_j]
            induced = instance.interaction.get((layout, site_i, site_j), 0.0)
            squared_error += (induced - Q[i, j]) ** 2

    return squared_error


def random_assignment(num_atoms, sites, rng):
    selected_sites = rng.sample(sites, num_atoms)
    return {atom: site for atom, site in zip(range(1, num_atoms + 1), selected_sites)}


def greedy_assignment(instance, layout, sites, rng):
    Q = instance.Q
    num_atoms = len(Q)
    atom_scores = []

    for atom in range(1, num_atoms + 1):
        row = np.copy(Q[atom - 1])
        row[atom - 1] = 0.0
        atom_scores.append((float(np.sum(np.abs(row))), atom))

    ordered_atoms = [atom for _, atom in sorted(atom_scores, reverse=True)]
    remaining_sites = list(sites)
    assignment = {}

    first_atom = ordered_atoms[0]
    assignment[first_atom] = rng.choice(remaining_sites)
    remaining_sites.remove(assignment[first_atom])

    for atom in ordered_atoms[1:]:
        best_site = None
        best_partial_error = math.inf

        for candidate_site in remaining_sites:
            partial_error = 0.0
            for assigned_atom, assigned_site in assignment.items():
                i = min(atom, assigned_atom) - 1
                j = max(atom, assigned_atom) - 1
                induced = instance.interaction.get(
                    (layout, candidate_site, assigned_site),
                    instance.interaction.get(
                        (layout, assigned_site, candidate_site),
                        0.0,
                    ),
                )
                partial_error += (induced - Q[i, j]) ** 2

            if partial_error < best_partial_error:
                best_partial_error = partial_error
                best_site = candidate_site

        assignment[atom] = best_site
        remaining_sites.remove(best_site)

    return assignment


def propose_neighbor(assignment, all_sites, rng):
    candidate = dict(assignment)
    atoms = list(candidate)
    used_sites = set(candidate.values())
    unused_sites = [site for site in all_sites if site not in used_sites]

    if unused_sites and rng.random() < 0.35:
        atom = rng.choice(atoms)
        candidate[atom] = rng.choice(unused_sites)
        return candidate

    atom_i, atom_j = rng.sample(atoms, 2)
    candidate[atom_i], candidate[atom_j] = candidate[atom_j], candidate[atom_i]
    return candidate


def run_simulated_annealing(
    instance,
    layout,
    initial_assignment,
    iterations,
    initial_temperature,
    cooling_rate,
    rng,
):
    current_assignment = dict(initial_assignment)
    current_value = evaluate_assignment(instance, layout, current_assignment)
    best_assignment = dict(current_assignment)
    best_value = current_value
    sites = available_sites(instance, layout)
    temperature = initial_temperature

    for _ in range(iterations):
        candidate_assignment = propose_neighbor(current_assignment, sites, rng)
        candidate_value = evaluate_assignment(instance, layout, candidate_assignment)
        delta = candidate_value - current_value

        if delta < 0 or rng.random() < math.exp(-delta / max(temperature, 1e-12)):
            current_assignment = candidate_assignment
            current_value = candidate_value

            if current_value < best_value:
                best_assignment = dict(current_assignment)
                best_value = current_value

        temperature *= cooling_rate

    return best_assignment, best_value


def run_fixed_layouts_heuristic(
    instance,
    restarts=50,
    iterations=5000,
    initial_temperature=100.0,
    cooling_rate=0.995,
    seed=42,
):
    rng = random.Random(seed)
    num_atoms = len(instance.Q)
    best_solution = None
    improvement_tolerance = 1e-9

    print("Running the fixed-layout heuristic...")
    print(f"Atoms: {num_atoms}")
    print(f"Layouts: {', '.join(instance.layouts)}")
    print(f"Restarts per layout: {restarts}")
    print(f"Iterations per restart: {iterations}")

    for layout in instance.layouts:
        sites = available_sites(instance, layout)

        if len(sites) < num_atoms:
            print(f"Skipping layout {layout}: only {len(sites)} available sites.")
            continue

        print(f"\nLayout: {layout}")
        layout_best = math.inf

        for restart in range(1, restarts + 1):
            if restart == 1:
                initial_assignment = greedy_assignment(instance, layout, sites, rng)
            else:
                initial_assignment = random_assignment(num_atoms, sites, rng)

            assignment, squared_error = run_simulated_annealing(
                instance=instance,
                layout=layout,
                initial_assignment=initial_assignment,
                iterations=iterations,
                initial_temperature=initial_temperature,
                cooling_rate=cooling_rate,
                rng=rng,
            )

            if squared_error < layout_best - improvement_tolerance:
                layout_best = squared_error
                print(
                    f"  Restart {restart:>4}: "
                    f"new layout best squared error = {squared_error:.6f}"
                )

            if (
                best_solution is None
                or squared_error < best_solution.squared_error - improvement_tolerance
            ):
                best_solution = HeuristicSolution(
                    layout=layout,
                    assignment=assignment,
                    squared_error=squared_error,
                )
                print(
                    f"  Restart {restart:>4}: "
                    f"new global best squared error = {squared_error:.6f}"
                )

    if best_solution is None:
        raise ValueError("No feasible layout assignment was found.")

    return best_solution


def solution_coordinates(instance, solution):
    num_atoms = len(instance.Q)
    coords = np.zeros((num_atoms, 2))

    for atom, site in solution.assignment.items():
        coords[atom - 1, 0] = instance.site_x[solution.layout, site]
        coords[atom - 1, 1] = instance.site_y[solution.layout, site]

    return coords


def induced_interaction_matrix(instance, solution):
    num_atoms = len(instance.Q)
    matrix = np.zeros((num_atoms, num_atoms))

    for i in range(num_atoms):
        atom_i = i + 1
        site_i = solution.assignment[atom_i]
        for j in range(i + 1, num_atoms):
            atom_j = j + 1
            site_j = solution.assignment[atom_j]
            value = instance.interaction.get((solution.layout, site_i, site_j), 0.0)
            matrix[i, j] = value
            matrix[j, i] = value

    return matrix


def print_solution(instance, solution):
    coords = solution_coordinates(instance, solution)
    induced = induced_interaction_matrix(instance, solution)

    print("\n--- Best Heuristic Solution ---")
    print(f"Selected layout: {solution.layout}")
    print(f"Squared Frobenius error: {solution.squared_error:.6f}")
    print(f"Frobenius error: {math.sqrt(solution.squared_error):.6f}")

    print("\nAtom assignments:")
    for atom in sorted(solution.assignment):
        site = solution.assignment[atom]
        x, y = coords[atom - 1]
        print(f"Atom {atom} -> Site {site} ({x:.6g}, {y:.6g})")

    print("\nTarget QUBO matrix:")
    print(instance.Q)

    print("\nInduced interaction matrix:")
    print(induced)


def plot_register(coords):
    import matplotlib.pyplot as plt
    from pulser import Register
    from pulser.devices import DigitalAnalogDevice

    qubits = {f"q{i}": coord for (i, coord) in enumerate(coords)}
    register = Register(qubits)
    register.draw(
        blockade_radius=DigitalAnalogDevice.rydberg_blockade_radius(1.0),
        draw_graph=False,
        draw_half_radius=True,
    )
    plt.show()


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run a simulated-annealing heuristic for Fixed Layout Register Design."
        )
    )
    parser.add_argument(
        "--dat-file",
        default="data/set_partitioning/ampl/modeling_test_qubo-fixed.dat",
        help="AMPL .dat instance for register_design_fixed_layouts.mod.",
    )
    parser.add_argument(
        "--restarts",
        type=int,
        default=50,
        help="Number of initial assignments per layout.",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=5000,
        help="Number of simulated-annealing iterations per restart.",
    )
    parser.add_argument(
        "--initial-temperature",
        type=float,
        default=100.0,
        help="Initial simulated-annealing temperature.",
    )
    parser.add_argument(
        "--cooling-rate",
        type=float,
        default=0.995,
        help="Multiplicative temperature cooling rate.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Disable the Pulser/Matplotlib register plot.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    instance = load_fixed_layouts_instance(args.dat_file)
    solution = run_fixed_layouts_heuristic(
        instance=instance,
        restarts=args.restarts,
        iterations=args.iterations,
        initial_temperature=args.initial_temperature,
        cooling_rate=args.cooling_rate,
        seed=args.seed,
    )
    print_solution(instance, solution)

    if not args.no_plot:
        plot_register(solution_coordinates(instance, solution))


if __name__ == "__main__":
    main()

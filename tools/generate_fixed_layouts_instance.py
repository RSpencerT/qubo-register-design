import argparse
import math
from pathlib import Path

import numpy as np
from pulser.devices import DigitalAnalogDevice

# Authoritative Pulser constant to ensure all Python components agree on C6.
DEFAULT_C6 = DigitalAnalogDevice.interaction_coeff


def format_number(value):
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.10g}"


def center_coordinates(coords):
    coords = np.array(coords, dtype=float)
    return coords - coords.mean(axis=0)


def line_layout(num_sites, spacing):
    coords = [(site * spacing, 0.0) for site in range(num_sites)]
    return center_coordinates(coords)


def grid_layout(num_sites, spacing):
    cols = math.ceil(math.sqrt(num_sites))
    coords = []

    for site in range(num_sites):
        row = site // cols
        col = site % cols
        coords.append((col * spacing, row * spacing))

    return center_coordinates(coords)


def circle_layout(num_sites, spacing):
    if num_sites == 1:
        return np.array([[0.0, 0.0]])

    radius = spacing / (2 * math.sin(math.pi / num_sites))
    coords = []

    for site in range(num_sites):
        angle = 2 * math.pi * site / num_sites
        coords.append((radius * math.cos(angle), radius * math.sin(angle)))

    return center_coordinates(coords)


def generate_layouts(num_sites, spacing):
    return {
        "line": line_layout(num_sites, spacing),
        "grid": grid_layout(num_sites, spacing),
        "circle": circle_layout(num_sites, spacing),
    }


def interaction_matrix(coords, c6=DEFAULT_C6):
    coords = np.array(coords, dtype=float)
    num_sites = len(coords)
    interactions = np.zeros((num_sites, num_sites))

    for i in range(num_sites):
        for j in range(num_sites):
            if i == j:
                continue

            distance = np.linalg.norm(coords[i] - coords[j])
            interactions[i, j] = c6 / distance**6

    return interactions


def load_q_matrix(csv_filename):
    Q = np.loadtxt(csv_filename, delimiter=",")
    if Q.ndim != 2 or Q.shape[0] != Q.shape[1]:
        raise ValueError("The QUBO CSV must contain a square matrix.")
    if not np.allclose(Q, Q.T):
        raise ValueError("The QUBO matrix must be symmetric.")
    return Q


def write_fixed_layouts_dat(Q, layouts, output_filename, c6=DEFAULT_C6):
    output_path = Path(output_filename)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    num_atoms = len(Q)
    num_sites = len(next(iter(layouts.values())))
    site_names = [f"S{i + 1}" for i in range(num_sites)]
    layout_names = list(layouts.keys())

    with output_path.open("w") as f:
        f.write("# =========================================================\n")
        f.write("# Fixed Layout Register Design Instance\n")
        f.write("# =========================================================\n\n")

        f.write(f"param N := {num_atoms};\n\n")
        f.write(f"set LAYOUTS := {' '.join(layout_names)};\n")
        f.write(f"set SITES := {' '.join(site_names)};\n\n")

        headers = "  ".join(str(i + 1) for i in range(num_atoms))
        f.write(f"param Q :  {headers} :=\n")
        for i in range(num_atoms):
            row = " ".join(format_number(value) for value in Q[i])
            f.write(f"       {i + 1}   {row}\n")
        f.write(";\n\n")

        f.write(f"# C6 used to precompute Interaction: {format_number(c6)}\n\n")

        f.write(f"param Available : {'  '.join(site_names)} :=\n")
        for layout_name in layout_names:
            availability = "  ".join("1" for _ in site_names)
            f.write(f"    {layout_name}  {availability}\n")
        f.write(";\n\n")

        f.write(f"param Site_X : {'  '.join(site_names)} :=\n")
        for layout_name, coords in layouts.items():
            values = "  ".join(format_number(x) for x, _ in coords)
            f.write(f"    {layout_name}  {values}\n")
        f.write(";\n\n")

        f.write(f"param Site_Y : {'  '.join(site_names)} :=\n")
        for layout_name, coords in layouts.items():
            values = "  ".join(format_number(y) for _, y in coords)
            f.write(f"    {layout_name}  {values}\n")
        f.write(";\n\n")

        f.write("param Interaction :=\n")
        for layout_name, coords in layouts.items():
            interactions = interaction_matrix(coords, c6=c6)
            for i, site_i in enumerate(site_names):
                for j, site_j in enumerate(site_names):
                    if i == j:
                        continue
                    value = interactions[i, j]
                    f.write(
                        f"    {layout_name} {site_i} {site_j} {format_number(value)}\n"
                    )
        f.write(";\n")


def generate_fixed_layouts_instance(
    qubo_csv_filename,
    output_filename,
    num_sites=None,
    spacing=7.5,
    c6=DEFAULT_C6,
):
    Q = load_q_matrix(qubo_csv_filename)
    num_atoms = len(Q)

    if num_sites is None:
        num_sites = num_atoms
    if num_sites < num_atoms:
        raise ValueError(
            "num_sites must be greater than or equal to the number of atoms."
        )

    layouts = generate_layouts(num_sites=num_sites, spacing=spacing)
    write_fixed_layouts_dat(Q, layouts, output_filename, c6=c6)

    print(f"Generated fixed-layout instance: {output_filename}")
    print(f"Atoms: {num_atoms}")
    print(f"Sites per layout: {num_sites}")
    print(f"Layouts: {', '.join(layouts)}")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Generate an AMPL .dat instance for register_design_fixed_layouts.mod."
        )
    )
    parser.add_argument(
        "--input",
        default="data/set_partitioning/csv/modeling_test_qubo.csv",
        help="Input QUBO matrix in CSV format.",
    )
    parser.add_argument(
        "--output",
        default="data/set_partitioning/ampl/modeling_test_qubo-fixed.dat",
        help="Output AMPL .dat file.",
    )
    parser.add_argument(
        "--num-sites",
        type=int,
        default=None,
        help="Number of calibrated sites per layout. Defaults to the QUBO size.",
    )
    parser.add_argument(
        "--spacing",
        type=float,
        default=7.5,
        help="Nominal spacing between neighboring trapping sites.",
    )
    parser.add_argument(
        "--c6",
        type=float,
        default=DEFAULT_C6,
        help="Rydberg interaction coefficient.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    generate_fixed_layouts_instance(
        qubo_csv_filename=args.input,
        output_filename=args.output,
        num_sites=args.num_sites,
        spacing=args.spacing,
        c6=args.c6,
    )

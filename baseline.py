"""Baseline code adapted from the official Pulser tutorial.

Source:
https://github.com/pasqal-io/Pulser/blob/develop/tutorials/applications/QAA%20to%20solve%20a%20QUBO%20problem.ipynb
https://github.com/pasqal-io/Pulser/blob/develop/tutorials/retired/QAOA%20and%20QAA%20to%20solve%20a%20QUBO%20problem.ipynb

This file is kept as a reference baseline for comparison purposes.
"""

import matplotlib.pyplot as plt
import numpy as np
from pulser import Register
from pulser.devices import DigitalAnalogDevice
from scipy.optimize import minimize
from scipy.spatial.distance import pdist, squareform


def computing_interaction_matrix(coords):
    """Computes the interaction matrix for a given set of coordinates."""
    new_Q = squareform(DigitalAnalogDevice.interaction_coeff / pdist(coords) ** 6)
    return new_Q


def evaluate_mapping(new_coords, Q):
    """Cost function to minimize. Ideally, the pairwise distances are conserved."""
    new_coords = np.reshape(new_coords, (len(Q), 2))
    # Computes the matrix of the distances between all coordinate pairs
    new_Q = squareform(DigitalAnalogDevice.interaction_coeff / pdist(new_coords) ** 6)

    # Creates a copy of Q and sets its diagonal to zero
    # This prevents the solver from optimizing linear terms via coordinates.
    Q_off_diag = np.copy(Q)
    np.fill_diagonal(Q_off_diag, 0.0)

    return np.linalg.norm(new_Q - Q_off_diag)


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

    error = evaluate_mapping(coords, Q)

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


def run_baseline(Q):
    np.random.seed(0)

    # Generates an initial random guess for the 2D coordinates
    x0 = np.random.random(len(Q) * 2)

    # Runs the Nelder-Mead optimization algorithm
    res = minimize(
        evaluate_mapping,
        x0,
        args=(Q,),
        method="Nelder-Mead",
        tol=1e-9,
        options={"maxiter": 200000, "maxfev": None},
    )

    # Reshapes the flat result array back into pairs of (x, y) coordinates
    coords = np.reshape(res.x, (len(Q), 2))

    return coords


def main():
    # Defines the target QUBO matrix Q from the problem definition image
    Q = np.array(
        [
            [-17, 10, 10, 10, 0, 20],
            [10, -18, 10, 10, 10, 20],
            [10, 10, -29, 10, 20, 20],
            [10, 10, 10, -19, 10, 10],
            [0, 10, 20, 10, -17, 10],
            [20, 20, 20, 10, 10, -28],
        ]
    )
    # Q = np.array(
    #     [
    #         [-10.0, 19.7365809, 19.7365809, 5.42015853, 5.42015853],
    #         [19.7365809, -10.0, 20.67626392, 0.17675796, 0.85604541],
    #         [19.7365809, 20.67626392, -10.0, 0.85604541, 0.17675796],
    #         [5.42015853, 0.17675796, 0.85604541, -10.0, 0.32306662],
    #         [5.42015853, 0.85604541, 0.17675796, 0.32306662, -10.0],
    #     ]
    # )
    coords = run_baseline(Q)
    print_results(coords, Q)
    plot_register(coords)


if __name__ == "__main__":
    main()

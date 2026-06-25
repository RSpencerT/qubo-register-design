# 6 - Provide a heuristic to solve the chosen version of the problem.

The chosen version is the **Fixed Layout Register Design** problem. The exact MIQP model is implemented in:

* [ampl/register_design_fixed_layouts.mod](../ampl/register_design_fixed_layouts.mod)

The fixed-layout instance generator is available in:

* [tools/generate_fixed_layouts_instance.py](../tools/generate_fixed_layouts_instance.py)

The exact AMPL runner is available in:

* [ampl/scripts/run_register_design_fixed_layouts.py](../ampl/scripts/run_register_design_fixed_layouts.py)

The heuristic implementation proposed in this answer is available in:

* [heuristics/fixed_layout_register_design.py](../heuristics/fixed_layout_register_design.py)

It can be executed with:

```bash
.venv/bin/python heuristics/fixed_layout_register_design.py
```

For a quick terminal-only run without plotting:

```bash
.venv/bin/python heuristics/fixed_layout_register_design.py --no-plot
```

Because the fixed-layout version contains binary layout-selection and atom-to-site assignment decisions, a natural heuristic is a **multi-start layout-aware local search**. The heuristic searches over the finite catalog of layouts and, for each layout, tries to find a good assignment of atoms to available trapping sites.

## 1. Objective Used by the Heuristic

For a selected layout $\ell$, let $\pi(i)$ be the site assigned to atom $i$. The physical interaction induced by the assignment is:

$$
\hat Q_{ij}
=
I_{\ell,\pi(i),\pi(j)}
\quad
\forall i<j
$$

where $I_{\ell st}$ is the precomputed interaction between sites $s$ and $t$ in layout $\ell$.

The heuristic minimizes the same objective as the exact model:

$$
\min
\sum_{i<j}
\left(
I_{\ell,\pi(i),\pi(j)}
-Q_{ij}
\right)^2
$$

This is the squared Frobenius mismatch between the target QUBO interaction matrix and the interaction matrix produced by the selected fixed layout.

## 2. Main Idea

The heuristic has three levels:

1. Try each feasible layout.
2. Generate several initial atom-to-site assignments for that layout.
3. Improve each assignment using local search or simulated annealing.

The best solution found over all layouts and restarts is returned.

## 3. Initial Solution Construction

For each layout $\ell$, first check whether it has at least $N$ available sites. If not, skip it.

Then create one or more initial assignments using either:

* **Random initialization**: randomly assign atoms to distinct available sites.
* **Greedy initialization**: match the largest target interactions in $Q$ with the largest physical interactions in the layout.

The greedy initialization can be described as follows:

1. Sort atom pairs $(i,j)$ by decreasing $|Q_{ij}|$.
2. Sort available site pairs $(s,t)$ by decreasing $I_{\ell st}$.
3. Assign atoms appearing in high-priority pairs to sites appearing in high-interaction site pairs whenever this does not violate the one-atom-per-site rule.
4. Assign any remaining atoms randomly to unused available sites.

This produces a structured starting point that already tries to align strong QUBO interactions with strong physical interactions.

## 4. Local Search Neighborhood

Given a current assignment $\pi$, the heuristic explores neighboring assignments using simple moves:

* **Swap move**: exchange the sites of two atoms $i$ and $j$.
* **Relocation move**: move one atom $i$ to an unused available site, if the layout has more sites than atoms.
* **Layout move**: restart the assignment search on another layout.

For each move, compute the change in the objective. The current implementation recomputes the full objective after each move, which costs $O(N^2)$. A more optimized implementation can use delta evaluation: a swap or relocation only affects the interactions involving the moved atoms, so the objective change can be computed in $O(N)$.

For example, if atom $i$ moves from site $s$ to site $s'$, only terms involving $i$ change:

$$
\Delta
=
\sum_{j\ne i}
\left[
\left(
I_{\ell,s',\pi(j)}-Q_{ij}
\right)^2
-
\left(
I_{\ell,s,\pi(j)}-Q_{ij}
\right)^2
\right]
$$

The move is accepted if it improves the objective.

## 5. Simulated Annealing Variant

To avoid getting stuck in poor local minima, the local search can be replaced by simulated annealing.

At each iteration:

1. Generate a random swap or relocation move.
2. Compute the objective change $\Delta$.
3. If $\Delta < 0$, accept the move.
4. If $\Delta \ge 0$, accept the move with probability:

$$
P(\text{accept})
=
\exp\left(-\frac{\Delta}{T}\right)
$$

where $T$ is the current temperature.

The temperature is gradually reduced:

$$
T \leftarrow \alpha T
$$

with $0<\alpha<1$, for example $\alpha=0.995$.

This allows the heuristic to occasionally accept worse assignments early in the search, while becoming more selective as the search progresses.

## 6. Pseudocode

```text
best_solution = None
best_value = infinity

for each layout l in L:
    available_sites = sites available in layout l

    if number of available_sites < number of atoms:
        continue

    for restart in 1,...,R:
        assignment = build_initial_assignment(Q, l, available_sites)
        value = evaluate_assignment(Q, l, assignment)

        T = initial_temperature

        for iteration in 1,...,K:
            candidate = random_swap_or_relocation(assignment)
            candidate_value = evaluate_assignment(Q, l, candidate)
            delta = candidate_value - value

            if delta < 0:
                assignment = candidate
                value = candidate_value
            else:
                with probability exp(-delta / T):
                    assignment = candidate
                    value = candidate_value

            T = alpha * T

            if value < best_value:
                best_solution = (l, assignment)
                best_value = value

return best_solution, best_value
```

## 7. Complexity

Let:

* $L$ be the number of layouts;
* $R$ be the number of restarts per layout;
* $K$ be the number of local-search iterations per restart;
* $N$ be the number of atoms.

The current implementation recomputes the objective from scratch after each move, so its complexity is:

$$
O(LRK N^2)
$$

If delta evaluation is added later for swap and relocation moves, each iteration would cost $O(N)$, giving:

$$
O(LR(N^2 + KN))
$$

The $N^2$ term comes from evaluating the initial assignment.

## 8. Use as a Warm Start

The heuristic can also be used as a warm start for the exact AMPL model. It returns:

* the selected layout $\ell$;
* the atom-to-site assignment $\pi(i)$;
* the induced interaction matrix $\hat Q$;
* the objective value.

These values can initialize the binary variables $u_\ell$ and $a_{\ell i s}$ in the MIQP model before calling a solver such as Gurobi.

This is useful because the exact model can be expensive for large instances, while a good heuristic solution gives the solver a strong incumbent early in the search.

## 9. Related Code in the Project

The project already contains related heuristic ideas for the free-space version in:

* [main.py](../main.py)

In particular, the following functions are useful references:

* `run_random_grid_multistart`
* `run_diff_evolution`
* `run_basin_hopping`
* `evaluate_mapping`

For the fixed-layout version, the same evaluation idea can be reused, but the search space changes from continuous atom coordinates to discrete assignments of atoms to calibrated trapping sites.

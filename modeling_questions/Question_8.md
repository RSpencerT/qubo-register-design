# 8 - Discuss the complexity of solving the proposed heuristic.

The proposed heuristic is the **multi-start simulated annealing heuristic** for the **Fixed Layout Register Design** problem.

The implementation is available in:

* [heuristics/fixed_layout_register_design.py](../heuristics/fixed_layout_register_design.py)

The exact AMPL model used for comparison is:

* [ampl/register_design_fixed_layouts.mod](../ampl/register_design_fixed_layouts.mod)

## 1. Notation

Let:

* $N$ be the number of atoms;
* $L$ be the number of fixed layouts;
* $S_\ell$ be the number of available sites in layout $\ell$;
* $S=\max_{\ell} S_\ell$;
* $R$ be the number of restarts per layout;
* $K$ be the number of simulated annealing iterations per restart.

The implemented heuristic searches over feasible atom-to-site assignments for each layout.

For one layout $\ell$, a feasible assignment is an injective mapping:

$$
\pi:
\{1,\dots,N\}
\rightarrow
\mathcal S_\ell
$$

where no two atoms use the same site.

## 2. Size of the Search Space

For a fixed layout $\ell$, the number of possible assignments is:

$$
P(S_\ell,N)
=
\frac{S_\ell!}{(S_\ell-N)!}
$$

If $S_\ell=N$, this becomes:

$$
N!
$$

Across all layouts, the total search space is:

$$
\sum_{\ell=1}^{L}
\frac{S_\ell!}{(S_\ell-N)!}
$$

Therefore, exact enumeration is factorial in the number of atoms. This is the reason a heuristic is useful: even for moderate $N$, trying every assignment quickly becomes impractical.

## 3. Objective Evaluation Cost

For a selected layout $\ell$ and assignment $\pi$, the objective is:

$$
F(\ell,\pi)
=
\sum_{i<j}
\left(
I_{\ell,\pi(i),\pi(j)}
-Q_{ij}
\right)^2
$$

There are:

$$
\frac{N(N-1)}{2}
$$

atom pairs. Therefore, evaluating one complete assignment from scratch costs:

$$
O(N^2)
$$

In the implementation, this is done by:

* `evaluate_assignment`

## 4. Initial Solution Complexity

The heuristic uses two types of initial solutions.

### Random Initialization

Random initialization samples $N$ distinct sites from the available sites. Its complexity is:

$$
O(N)
$$

This is implemented by:

* `random_assignment`

### Greedy Initialization

The greedy initialization first computes an atom score based on the absolute interaction values in $Q$. This requires scanning the matrix:

$$
O(N^2)
$$

Then, for each atom, it tests candidate sites against the atoms already assigned. In the worst case, this costs:

$$
O(NS^2)
$$

If the number of sites is close to the number of atoms, $S=O(N)$, this becomes:

$$
O(N^3)
$$

This is implemented by:

* `greedy_assignment`

Since the greedy initialization is only used once per layout in the current implementation, its cost is usually dominated by the simulated annealing phase for large $K$.

## 5. Neighborhood Move Complexity

The implemented neighborhood contains two moves:

* **swap move**: exchange the sites assigned to two atoms;
* **relocation move**: move one atom to an unused available site.

Generating a neighbor requires copying the current assignment dictionary and modifying one or two entries.

With the current implementation, this costs:

$$
O(N)
$$

because the assignment is copied.

This is implemented by:

* `propose_neighbor`

## 6. Simulated Annealing Complexity

Each simulated annealing iteration:

1. proposes a neighboring assignment;
2. evaluates the candidate assignment;
3. computes the objective difference;
4. accepts or rejects the move;
5. updates the temperature.

Because the current implementation evaluates the candidate assignment from scratch, each iteration costs:

$$
O(N^2)
$$

Therefore, one simulated annealing run with $K$ iterations costs:

$$
O(KN^2)
$$

This is implemented by:

* `run_simulated_annealing`

## 7. Total Time Complexity of the Implemented Heuristic

For each feasible layout, the heuristic performs $R$ restarts. Each restart runs simulated annealing for $K$ iterations.

Ignoring lower-order initialization costs, the total running time is:

$$
O(LRKN^2)
$$

Including initialization, a more detailed bound is:

$$
O\left(
L
\left(
N^3
+
RKN^2
\right)
\right)
$$

when $S=O(N)$.

If the number of available sites can be larger than $N$, a more general expression is:

$$
O\left(
L
\left(
NS^2
+
RKN^2
\right)
\right)
$$

The $NS^2$ term comes from the greedy initialization, and the $RKN^2$ term comes from repeated objective evaluations during simulated annealing.

## 8. Possible Improvement with Delta Evaluation

The theoretical heuristic described in Question 6 can be improved using delta evaluation.

For a swap or relocation move, only the terms involving the moved atoms change. Therefore, instead of recomputing the full objective in:

$$
O(N^2)
$$

the objective difference can be computed in:

$$
O(N)
$$

With delta evaluation, one simulated annealing run would cost:

$$
O(KN)
$$

and the full heuristic would cost:

$$
O(LRKN)
$$

up to initialization costs.

The current implementation favors clarity and correctness, so it recomputes the full objective. For larger instances, implementing delta evaluation would be the most direct performance improvement.

## 9. Space Complexity

The main memory requirements are:

* the QUBO matrix $Q$: $O(N^2)$;
* the interaction matrices $I_{\ell st}$: $O(LS^2)$;
* the availability matrix: $O(LS)$;
* the current and candidate assignments: $O(N)$;
* the induced interaction matrix printed at the end: $O(N^2)$.

Therefore, the total space complexity is:

$$
O(N^2 + LS^2)
$$

If $S=O(N)$, this becomes:

$$
O(LN^2)
$$

## 10. Comparison with Exact Optimization

The exact Fixed Layout Register Design model contains binary variables for:

* layout selection;
* atom-to-site assignment;
* pair assignment linearization.

The number of pair-assignment variables scales approximately as:

$$
O(LN^2S^2)
$$

This makes the exact MIQP model much harder than the heuristic for large instances.

The heuristic does not guarantee global optimality, but it has a predictable polynomial running time once $R$ and $K$ are fixed:

$$
O(LRKN^2)
$$

in the current implementation.

This makes it suitable for:

* quickly generating feasible solutions;
* testing larger fixed-layout instances;
* producing warm starts for the exact AMPL model;
* comparing different layout catalogs before running a full MIQP solve.

## 11. Practical Interpretation

The main parameters controlling runtime are:

* `--restarts`, corresponding to $R$;
* `--iterations`, corresponding to $K$;
* the number of layouts in the `.dat` file, corresponding to $L$;
* the number of atoms, corresponding to $N$.

For example, doubling the number of iterations roughly doubles the runtime. Doubling the number of atoms can increase runtime by approximately a factor of four in the current implementation, because each objective evaluation is $O(N^2)$.

Thus, the heuristic is scalable in the search-control parameters $R$ and $K$, but the atom count $N$ is the most important structural driver of computational cost.

# 9 - Briefly discuss how you would apply the Genetic algorithm and the Tabu Search method to solve the problem. Which one do you think would be more appropriate to each version of the problem?

There are two main versions of the register design problem in this project:

* the **Free Space Register Design** problem, where atom coordinates are continuous decision variables;
* the **Fixed Layout Register Design** problem, where one layout is selected and atoms are assigned to calibrated trapping sites.

The related implementations are:

* [main.py](../main.py) for the free-space heuristic experiments;
* [ampl/register_design.mod](../ampl/register_design.mod) for the free-space AMPL model;
* [ampl/register_design_fixed_layouts.mod](../ampl/register_design_fixed_layouts.mod) for the fixed-layout AMPL model;
* [heuristics/fixed_layout_register_design.py](../heuristics/fixed_layout_register_design.py) for the fixed-layout simulated annealing heuristic.

## 1. Genetic Algorithm

A Genetic Algorithm is a population-based metaheuristic. It keeps several candidate solutions at the same time and evolves them using selection, crossover, and mutation.

### Free Space Version

For the free-space version, a chromosome can encode the continuous coordinates of all atoms:

$$
x =
(x_1,y_1,x_2,y_2,\dots,x_N,y_N)
$$

The fitness function should be the same mapping error already used in the project:

$$
f(x)
=
\left\|
U(x)-Q
\right\|
$$

where \(U(x)\) is the physical interaction matrix induced by the coordinates.

A Genetic Algorithm for the free-space version would use:

* **selection**: keep coordinate vectors with smaller mapping error;
* **crossover**: combine coordinates from two parent registers;
* **mutation**: perturb atom coordinates with random noise;
* **repair**: enforce bounds and minimum-distance constraints if needed.

This is similar in spirit to the global optimization routines already tested in [main.py](../main.py), such as differential evolution and basin hopping.

### Fixed Layout Version

For the fixed-layout version, a chromosome can encode:

$$
(\ell,\pi)
$$

where \(\ell\) is the selected layout and \(\pi\) is a permutation-like assignment of atoms to sites.

For example:

```text
layout = grid
assignment = [S4, S2, S6, S1, S3, S5]
```

The fitness function is:

$$
f(\ell,\pi)
=
\sum_{i<j}
\left(
I_{\ell,\pi(i),\pi(j)}
-Q_{ij}
\right)^2
$$

The Genetic Algorithm would use:

* **selection**: prefer assignments with lower squared Frobenius error;
* **crossover**: combine parts of two atom-to-site assignments;
* **mutation**: swap two assigned sites, relocate one atom, or change the selected layout;
* **repair**: remove duplicate site assignments and restore feasibility.

The main difficulty is that assignments must remain injective: two atoms cannot occupy the same trapping site. Therefore, crossover and mutation must be designed carefully, using permutation-based operators such as order crossover, partially matched crossover, or a repair step.

## 2. Tabu Search

Tabu Search is a local-search metaheuristic. It starts from one feasible solution and repeatedly moves to a neighboring solution, while maintaining a tabu list of recently used moves or assignments to avoid cycling.

### Free Space Version

For the free-space version, a solution is again a continuous coordinate vector:

$$
x =
(x_1,y_1,\dots,x_N,y_N)
$$

Possible moves include:

* move one atom by a small vector;
* swap the positions of two atoms;
* perturb a small group of atoms;
* apply a local continuous optimizer after each move.

However, Tabu Search is less natural for the free-space version because the neighborhood is continuous and very large. One must discretize coordinate changes or define a finite set of perturbation moves. This can work, but it introduces several tuning choices: step size, move directions, tabu tenure, and local refinement frequency.

### Fixed Layout Version

For the fixed-layout version, Tabu Search is very natural because the search space is discrete.

A solution is:

$$
(\ell,\pi)
$$

and the neighborhood can be exactly the same as in the heuristic from Question 6:

* swap the sites of two atoms;
* move one atom to an unused site;
* switch to another layout and rebuild the assignment.

The tabu list can store:

* recently swapped atom pairs;
* recently used atom-site assignments;
* recently visited complete assignments;
* recently selected layouts.

At each iteration, the algorithm evaluates candidate moves and chooses the best non-tabu move. A tabu move may still be accepted if it satisfies an aspiration criterion, for example if it gives the best solution found so far.

The objective is the same:

$$
\sum_{i<j}
\left(
I_{\ell,\pi(i),\pi(j)}
-Q_{ij}
\right)^2
$$

## 3. Which Method Is More Appropriate?

### Free Space Register Design

For the **free-space version**, a Genetic Algorithm is more appropriate than Tabu Search.

The reason is that the variables are continuous coordinates. Genetic Algorithms can naturally handle real-valued chromosomes and global exploration over a continuous domain. They can also be combined with local optimization, producing a hybrid memetic algorithm:

1. use GA for global exploration;
2. refine the best individuals with a local optimizer;
3. keep the best coordinate configuration.

Tabu Search can still be used, but it is less direct because a continuous neighborhood must be discretized.

### Fixed Layout Register Design

For the **fixed-layout version**, Tabu Search is more appropriate than a Genetic Algorithm.

The reason is that the decision space is combinatorial:

* choose one layout;
* assign atoms to distinct sites.

This structure is ideal for local-search moves such as swaps and relocations. Tabu Search can exploit these moves very efficiently, and it naturally avoids cycling between equivalent assignments.

A Genetic Algorithm can also be used for this version, but it requires special permutation crossover and repair mechanisms to preserve feasibility.

## 4. Summary

| Problem version | Best suited method | Reason |
|---|---|---|
| Free Space Register Design | Genetic Algorithm | Coordinates are continuous, and GA supports real-valued global search. |
| Fixed Layout Register Design | Tabu Search | Assignments are discrete and naturally handled by swap/relocation neighborhoods. |

Therefore:

* use **Genetic Algorithms** mainly for the continuous free-space version;
* use **Tabu Search** mainly for the discrete fixed-layout version;
* for the fixed-layout version, Tabu Search is also a natural alternative to the simulated annealing heuristic implemented in [heuristics/fixed_layout_register_design.py](../heuristics/fixed_layout_register_design.py).

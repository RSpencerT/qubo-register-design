# 5 - Provide valid inequalities for the proposed formulation

The proposed formulation is the **Fixed Layout Register Design** MIQP from Question 4. The full AMPL implementation is:

* [ampl/register_design_fixed_layouts.mod](../ampl/register_design_fixed_layouts.mod)

The main binary variables are:

* $u_\ell$: whether layout $\ell$ is selected;
* $a_{\ell i s}$: whether atom $i$ is assigned to site $s$ in layout $\ell$;
* $p_{\ell i j s t}$: whether atom pair $(i,j)$ is assigned to site pair $(s,t)$ in layout $\ell$.

The following valid inequalities strengthen the formulation without removing any feasible solution.

## 1. Layout Capacity Inequality

Each selected layout must contain enough available sites to place all atoms. Let:

$$
M_\ell = \sum_{s\in\mathcal S} B_{\ell s}
$$

where $B_{\ell s}$ indicates whether site $s$ is available in layout $\ell$. Then:

$$
N u_\ell
\le
\sum_{s\in\mathcal S} B_{\ell s}
\quad
\forall \ell\in\mathcal L
$$

Equivalently:

$$
u_\ell = 0
\quad
\text{if}
\quad
M_\ell < N
$$

This inequality removes layouts that cannot physically host all atoms.

## 2. Assignment Activation Inequality

The current model already enforces:

$$
a_{\ell i s}
\le
B_{\ell s}u_\ell
$$

A useful aggregate strengthening is:

$$
\sum_{s\in\mathcal S}
a_{\ell i s}
\le
u_\ell
\quad
\forall \ell\in\mathcal L,\ i\in A
$$

This explicitly states that atom $i$ cannot be assigned to any site of layout $\ell$ unless layout $\ell$ is selected.

## 3. Exact Assignment Inside the Selected Layout

Since exactly one layout is selected, the assignment of each atom can be linked more tightly to each layout:

$$
\sum_{s\in\mathcal S}
a_{\ell i s}
=
u_\ell
\quad
\forall \ell\in\mathcal L,\ i\in A
$$

This is stronger than:

$$
\sum_{\ell\in\mathcal L}
\sum_{s\in\mathcal S}
a_{\ell i s}
=1
$$

together with activation constraints. It forces every atom to have exactly one assignment in the selected layout and zero assignments in non-selected layouts.

## 4. Site Capacity Strengthening

The model already contains:

$$
\sum_{i\in A}
a_{\ell i s}
\le
B_{\ell s}u_\ell
$$

An aggregate version over all sites gives:

$$
\sum_{i\in A}
\sum_{s\in\mathcal S}
a_{\ell i s}
=
N u_\ell
\quad
\forall \ell\in\mathcal L
$$

This ensures that the selected layout receives exactly $N$ atom-site assignments, while every non-selected layout receives none.

## 5. Pair-Assignment Activation Inequality

For all $\ell\in\mathcal L$, $i<j$, and $s\ne t$:

$$
p_{\ell i j s t}
\le
u_\ell
$$

This is implied by:

$$
p_{\ell i j s t}\le a_{\ell i s}
$$

and:

$$
a_{\ell i s}\le B_{\ell s}u_\ell
$$

but adding it can tighten the LP relaxation because it directly links pair variables to layout selection.

## 6. Pair-Assignment Availability Inequality

For all $\ell\in\mathcal L$, $i<j$, and $s\ne t$:

$$
p_{\ell i j s t}
\le
B_{\ell s}
$$

$$
p_{\ell i j s t}
\le
B_{\ell t}
$$

If either site is unavailable, the pair assignment must be zero. These are valid because atoms cannot be assigned to unavailable sites.

## 7. Unique Site-Pair Assignment for Each Atom Pair

For every atom pair $i<j$, exactly one ordered site pair is induced by the selected layout and assignments:

$$
\sum_{\ell\in\mathcal L}
\sum_{\substack{s,t\in\mathcal S\\s\ne t}}
p_{\ell i j s t}
=
1
\quad
\forall i<j
$$

This is implied by the assignment constraints and the linearization, but adding it explicitly strengthens the formulation.

## 8. Symmetry Reduction for Identical Layouts

If two layouts $\ell_1$ and $\ell_2$ have identical interaction matrices:

$$
I_{\ell_1 st}=I_{\ell_2 st}
\quad
\forall s,t
$$

then they are interchangeable. A valid symmetry-breaking inequality is:

$$
u_{\ell_1}
\ge
u_{\ell_2}
$$

for an arbitrary ordering $\ell_1 < \ell_2$. This prevents the solver from exploring equivalent solutions multiple times.

## 9. Atom Ordering Symmetry Breaking

If two atoms $i$ and $j$ have identical interaction profiles in $Q$, meaning:

$$
Q_{ik}=Q_{jk}
\quad
\forall k\ne i,j
$$

then atoms $i$ and $j$ are interchangeable. A symmetry-breaking rule can impose an ordering on their assigned site indices.

Let $idx(s)$ be the index of site $s$. Then:

$$
\sum_{\ell\in\mathcal L}
\sum_{s\in\mathcal S}
idx(s)a_{\ell i s}
\le
\sum_{\ell\in\mathcal L}
\sum_{s\in\mathcal S}
idx(s)a_{\ell j s}
$$

This reduces equivalent permutations of identical atoms.

## 10. Interaction Error Bounds

Let:

$$
I_{\min} =
\min_{\ell,s,t:s\ne t} I_{\ell st}
$$

and:

$$
I_{\max} =
\max_{\ell,s,t:s\ne t} I_{\ell st}
$$

Since the induced interaction for any atom pair must lie between $I_{\min}$ and $I_{\max}$, we can bound the error variable:

$$
I_{\min}-Q_{ij}
\le
e_{ij}
\le
I_{\max}-Q_{ij}
\quad
\forall i<j
$$

These bounds help the solver handle the quadratic objective more tightly.

## 11. Zero-Target Interaction Cut

If $Q_{ij}=0$, the model should prefer site pairs with low physical interaction. A valid inequality using a tolerance threshold $\tau$ can eliminate site-pair assignments that are guaranteed to exceed a user-defined acceptable interaction:

$$
p_{\ell i j s t}=0
\quad
\text{if}
\quad
Q_{ij}=0
\ \text{and}\
I_{\ell st}>\tau
$$

This is not valid for the unconstrained least-squares model unless $\tau$ is declared as a hard modeling requirement. If the problem definition allows a maximum acceptable interaction for zero entries, then this cut is valid and very strong.

## 12. Large-Target Interaction Candidate Cut

Similarly, if $Q_{ij}$ is large, site pairs with interactions far below the target can be excluded under a tolerance parameter $\delta$:

$$
p_{\ell i j s t}=0
\quad
\text{if}
\quad
Q_{ij}>0
\ \text{and}\
|I_{\ell st}-Q_{ij}|>\delta
$$

This is valid only when $\delta$ is imposed as a hard approximation tolerance. It transforms the problem from pure least-squares minimization into a constrained approximation problem.

## 13. Practical AMPL Additions

The following inequalities are the safest to add directly to the current AMPL formulation because they do not change the problem definition:

```ampl
subject to Assign_Each_Atom_In_Selected_Layout {l in LAYOUTS, i in ATOMS}:
    sum {s in SITES} assign[l,i,s] = use_layout[l];

subject to Layout_Assignment_Count {l in LAYOUTS}:
    sum {i in ATOMS, s in SITES} assign[l,i,s] = N * use_layout[l];

subject to Pair_Assign_Layout_Activation {
    l in LAYOUTS,
    i in ATOMS,
    j in ATOMS,
    s in SITES,
    t in SITES:
        i < j and s != t
}:
    pair_assign[l,i,j,s,t] <= use_layout[l];

subject to Unique_Site_Pair_Per_Atom_Pair {i in ATOMS, j in ATOMS: i < j}:
    sum {l in LAYOUTS, s in SITES, t in SITES: s != t}
        pair_assign[l,i,j,s,t] = 1;
```

These inequalities strengthen the relaxation and reduce unnecessary branching without removing feasible solutions.

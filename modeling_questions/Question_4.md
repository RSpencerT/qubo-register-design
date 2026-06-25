# 4 - Provide a complete mathematical formulation that solves the chosen version of problem and optimizes the cost function proposed for question 2.

The chosen version is the **Fixed Layout Register Design** problem. In this variant, atoms cannot be placed anywhere in the continuous plane. Instead, a register must be selected from a finite catalog of calibrated layouts, and each atom must be assigned to one available trapping site in the selected layout.

The AMPL implementation is available in:

* [ampl/register_design_fixed_layouts.mod](../ampl/register_design_fixed_layouts.mod)

For comparison, the other related AMPL models are:

* [ampl/register_design.mod](../ampl/register_design.mod) for the free-space continuous version
* [ampl/qubo.mod](../ampl/qubo.mod) for the classical QUBO formulation
* [ampl/set_partitioning.mod](../ampl/set_partitioning.mod) for the source set partitioning model

## 1. Sets

Let:

* `A={1,...,N}` be the set of atoms.
* `L` be the set of pre-calibrated layouts.
* `S` be the set of trapping sites.

## 2. Parameters

The input data are:

* `Q[i,j]`: target QUBO interaction between atoms `i` and `j`.
* `B[l,s] in {0,1}`: availability of site `s` in layout `l`.
* `X[l,s]`, `Y[l,s]`: coordinates of site `s` in layout `l`.
* `I[l,s,t]`: precomputed interaction between sites `s` and `t` in layout `l`.

The physical interaction is precomputed as:

$$
I_{\ell st}
=
\frac{C_6}
{
\left(
(X_{\ell s}-X_{\ell t})^2
+
(Y_{\ell s}-Y_{\ell t})^2
\right)^3
}
\quad s\ne t
$$

and:

$$
I_{\ell ss}=0
$$

The data instance can be generated with:

* [tools/generate_fixed_layouts_instance.py](../tools/generate_fixed_layouts_instance.py)

## 3. Decision Variables

### Layout Selection

Let:

$$
u_\ell =
\begin{cases}
1, & \text{if layout } \ell \text{ is selected} \\
0, & \text{otherwise}
\end{cases}
$$

### Atom-to-Site Assignment

Let:

$$
a_{\ell i s} =
\begin{cases}
1, & \text{if atom } i \text{ is assigned to site } s \text{ in layout } \ell \\
0, & \text{otherwise}
\end{cases}
$$

### Pair Assignment Linearization

The interaction induced by two assigned atoms contains a product:

$$
a_{\ell i s}a_{\ell j t}
$$

To linearize this product, introduce:

$$
p_{\ell i j s t}
=
a_{\ell i s}a_{\ell j t}
$$

for `i<j` and `s!= t`.

### Interaction Error

Let:

$$
e_{ij}
$$

be the difference between the physical interaction induced by the selected layout and the target QUBO interaction `Q[i,j]`.

## 4. Mathematical Formulation

### Select Exactly One Layout

$$
\sum_{\ell\in\mathcal L} u_\ell = 1
$$

### Assign Each Atom Exactly Once

$$
\sum_{\ell\in\mathcal L}
\sum_{s\in\mathcal S}
a_{\ell i s}
=1
\quad
\forall i\in A
$$

### Use Only Available Sites in the Selected Layout

$$
a_{\ell i s}
\le
B_{\ell s}u_\ell
\quad
\forall \ell\in\mathcal L,\ i\in A,\ s\in\mathcal S
$$

### At Most One Atom Per Site

$$
\sum_{i\in A}
a_{\ell i s}
\le
B_{\ell s}u_\ell
\quad
\forall \ell\in\mathcal L,\ s\in\mathcal S
$$

### Linearization Constraints

For all layouts `l in L`, atom pairs `i<j`, and site pairs `s!=t`:

$$
p_{\ell i j s t}
\le
a_{\ell i s}
$$

$$
p_{\ell i j s t}
\le
a_{\ell j t}
$$

$$
p_{\ell i j s t}
\ge
a_{\ell i s}
+
a_{\ell j t}
-1
$$

Together, these constraints enforce:

$$
p_{\ell i j s t}
=
a_{\ell i s}a_{\ell j t}
$$

for binary variables.

### Interaction Error Definition

For each atom pair `i<j`:

$$
e_{ij}
=
\sum_{\ell\in\mathcal L}
\sum_{\substack{s,t\in\mathcal S\\s\ne t}}
I_{\ell st}
p_{\ell i j s t}
-Q_{ij}
$$

## 5. Objective Function

The objective is to minimize the squared Frobenius mismatch between the target QUBO matrix and the interaction matrix induced by the fixed layout:

$$
\min
\sum_{i<j}
e_{ij}^2
$$

Equivalently:

$$
\min
\sum_{i<j}
\left(
\sum_{\ell\in\mathcal L}
\sum_{\substack{s,t\in\mathcal S\\s\ne t}}
I_{\ell st}
p_{\ell i j s t}
-Q_{ij}
\right)^2
$$

This corresponds to the AMPL objective:

```ampl
minimize Squared_Frobenius_Error:
    sum {i in ATOMS, j in ATOMS: i < j} interaction_error[i,j]^2;
```

## 6. Complete Model Type

The resulting formulation is a **Mixed-Integer Quadratic Program (MIQP)**:

* binary variables select layouts and atom-site assignments;
* linear constraints enforce feasibility;
* the objective is quadratic in the continuous error variables.

This is substantially more solver-friendly than a direct nonlinear formulation with products such as:

$$
u_\ell a_{\ell i s} a_{\ell j t}
$$

because the bilinear assignment products are explicitly linearized through `p[l,i,j,s,t]`.

## 7. Free-Space Reference Formulation

The free-space version is implemented in:

* [ampl/register_design.mod](../ampl/register_design.mod)

Its decision variables are continuous coordinates:

$$
x_i,y_i\in[0,L]
$$

and the objective is:

$$
\min
\sqrt{
\sum_{i\ne j}
\left(
\frac{C_6}
{
((x_i-x_j)^2+(y_i-y_j)^2)^3+\varepsilon
}
-Q_{ij}
\right)^2
}
$$

This formulation directly optimizes atom positions in the continuous plane, while the fixed-layout formulation restricts the search to a finite calibrated set of possible trapping configurations.

# 2 – Provide a cost function to be optimized.

The central optimization goal is to make the physical Rydberg interaction matrix induced by a register layout approximate a target QUBO interaction matrix `Q`. For two atoms `i` and `j` placed at coordinates `r[i]=(x[i],y[i])` and `r[j]=(x[j],y[j])`, the physical interaction is modeled as:

```text
U_{ij}(r) = \frac{C_6}{\|r_i-r_j\|^6}
```

Therefore, a natural cost function is the Frobenius mismatch between the target QUBO matrix and the physical interaction matrix:

```text
\min_{r_1,\dots,r_N}
\left(
\sum_{i \ne j}
\left(
\frac{C_6}{\|r_i-r_j\|^6+\varepsilon}
- Q_{ij}
\right)^2
\right)^{1/2}
```

where:

* `Q[i,j]` is the target interaction between QUBO variables `i` and `j`;
* `C_6` is the Rydberg interaction coefficient;
* `epsilon` is a small safety term to avoid division by zero;
* `r[i]` is the physical position of atom `i`.

## Free-Space Register Design

In the free-space version, atom coordinates are continuous decision variables:

```text
x_i,y_i \in [0,L]
```

The AMPL cost function is:

```text
\min
\sqrt{
\sum_{i \ne j}
\left(
\frac{C_6}{((x_i-x_j)^2+(y_i-y_j)^2)^3+\varepsilon}
-Q_{ij}
\right)^2
}
```

This is implemented in:

* `ampl/register_design.mod`
  * objective: `Frobenius_Error`

The same idea also appears in the Python heuristic code:

* `main.py`
  * function: `evaluate_mapping`
  * function: `evaluate_mapping_baseline`
* `baseline.py`
  * function: `evaluate_mapping`
* `ampl/scripts/run_register_design.py`
  * function: `evaluate_mapping`

The Python version uses:

```text
\sqrt{
\sum_{(i,j):Q_{ij}>0}(U_{ij}-Q_{ij})^2
+
\sum_{(i,j):Q_{ij}=0}\max(U_{ij}-\tau,0)^2
}
```

where `tau` is a small tolerance threshold for zero-interaction entries.

## Fixed-Layout Register Design

In the fixed-layout version, atoms cannot be placed freely. Instead, the model chooses one calibrated layout and assigns atoms to available sites in that layout.

Let:

* `L` be the set of candidate layouts;
* `S` be the set of calibrated sites;
* `I[l,s,t]` be the precomputed physical interaction between sites `s` and `t` in layout `l`;
* `a[l,i,s] in {0,1}` indicate whether atom `i` is assigned to site `s` in layout `l`.

The induced interaction between atoms `i` and `j` is:

```text
\hat Q_{ij}
=
\sum_{\ell \in L}
\sum_{s,t \in S,\,s\ne t}
I_{\ell st}\,a_{\ell i s}\,a_{\ell j t}
```

The cost function is:

```text
\min
\sum_{i<j}
(\hat Q_{ij}-Q_{ij})^2
```

This is implemented in:

* `ampl/register_design_fixed_layouts.mod`
  * objective: `Squared_Frobenius_Error`
  * error variable: `interaction_error[i,j]`

The data instance for this model can be generated with:

* `tools/generate_fixed_layouts_instance.py`
  * function: `generate_fixed_layouts_instance`

## Classical QUBO Objective

For the purely classical QUBO problem, the cost function is:

```text
\min_{x\in\{0,1\}^N} x^TQx
```

This is implemented in:

* `ampl/qubo.mod`
  * objective: `QUBO_Objective`
* `ampl/scripts/run_qubo.py`
  * reads and solves the AMPL QUBO model

## Set Partitioning Objective and QUBO Transformation

The original set partitioning objective is:

```text
\min \sum_j c_jx_j
```

subject to exact-cover constraints:

```text
\sum_j A_{ij}x_j = 1
```

This is implemented in:

* `ampl/set_partitioning.mod`
  * objective: `Total_Cost`

The transformation from set partitioning to QUBO uses the penalty objective:

```text
f(x)
=
\sum_j c_jx_j
+
\lambda
\sum_i
\left(
\sum_j A_{ij}x_j - 1
\right)^2
```

The corresponding implementation is in:

* `tools/set_partitioning_to_qubo.py`
  * function: `set_partitioning_dat_to_qubo`
* `tools/README.md`
  * derivation of the set partitioning to QUBO transformation

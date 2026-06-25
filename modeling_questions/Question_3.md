# 3 - Provide an upper-bound and a lower-bound for the cost function given any arbitrary instance of the problem.

The bounds depend on which version of the problem is considered. The most relevant objective in this project is the register-design mismatch between a target QUBO matrix `Q` and the physical interaction matrix induced by an atom register.

## 1. Free-Space Register Design

The free-space objective is:

$$
F(r)
=
\left(
\sum_{i \ne j}
\left(
\frac{C_6}{\|r_i-r_j\|^6+\varepsilon}
- Q_{ij}
\right)^2
\right)^{1/2}
$$

where `r[i]=(x[i],y[i])`, `x[i], y[i] in [0,L]`, `C_6>0`, and `epsilon>0`.

### Lower Bound

Since the objective is a square root of a sum of squared terms:

$$
F(r) \ge 0
$$

Thus, for any arbitrary instance:

$$
\boxed{LB = 0}
$$

This lower bound is tight if there exists a register geometry whose physical interactions match the target matrix exactly:

$$
\frac{C_6}{\|r_i-r_j\|^6+\varepsilon}=Q_{ij}
\quad \forall i\ne j
$$

In practice, this exact equality is often impossible due to geometric constraints and frustration, but `0` remains the universal lower bound.

### Upper Bound

Because of the safety term `epsilon`, the physical interaction is bounded above:

$$
0 \le
\frac{C_6}{\|r_i-r_j\|^6+\varepsilon}
\le
\frac{C_6}{\varepsilon}
$$

Therefore, for every pair `(i,j)`:

$$
\left|
\frac{C_6}{\|r_i-r_j\|^6+\varepsilon}
-Q_{ij}
\right|
\le
\max
\left\{
|Q_{ij}|,
\left|\frac{C_6}{\varepsilon}-Q_{ij}\right|
\right\}
$$

A valid instance-independent upper bound is:

$$
\boxed{
UB_{\text{free}}
=
\left(
\sum_{i \ne j}
\max
\left\{
|Q_{ij}|^2,
\left(\frac{C_6}{\varepsilon}-Q_{ij}\right)^2
\right\}
\right)^{1/2}
}
$$

This bound is conservative but valid for any arbitrary `Q`.

This objective appears in:

* `ampl/register_design.mod`
  * objective: `Frobenius_Error`
* `ampl/scripts/run_register_design.py`
  * objective value extracted as `Frobenius_Error`

## 2. Fixed-Layout Register Design

In the fixed-layout variant, the set of possible interactions is finite because the layouts and trapping sites are fixed. The objective is:

$$
F_{\text{fixed}}(a)
=
\sum_{i<j}
(\hat Q_{ij}-Q_{ij})^2
$$

where:

$$
\hat Q_{ij}
=
\sum_{\ell \in L}
\sum_{s,t \in S,\,s\ne t}
I_{\ell st}a_{\ell i s}a_{\ell j t}
$$

Here, `I[l,s,t]` is the precomputed interaction between sites `s` and `t` in layout `l`.

### Lower Bound

Again, this is a sum of squared terms, so:

$$
\boxed{LB_{\text{fixed}} = 0}
$$

### Upper Bound

Let:

$$
I_{\max} = \max_{\ell,s,t:s\ne t} I_{\ell st}
$$

Since the selected physical interaction for any atom pair must be one of the available site-pair interactions:

$$
0 \le \hat Q_{ij} \le I_{\max}
$$

Therefore:

$$
\left|\hat Q_{ij}-Q_{ij}\right|
\le
\max
\left\{
|Q_{ij}|,
|I_{\max}-Q_{ij}|
\right\}
$$

A valid upper bound is:

$$
\boxed{
UB_{\text{fixed}}
=
\sum_{i<j}
\max
\left\{
|Q_{ij}|^2,
(I_{\max}-Q_{ij})^2
\right\}
}
$$

Since the fixed-layout problem has a finite feasible set, a tighter upper bound can also be obtained by evaluating any feasible layout assignment. For example, assigning atoms greedily to the first available sites in any available layout gives a constructive feasible solution and therefore a valid upper bound:

$$
UB_{\text{constructive}}
=
F_{\text{fixed}}(a^{\text{greedy}})
$$

This objective appears in:

* `ampl/register_design_fixed_layouts.mod`
  * objective: `Squared_Frobenius_Error`
  * error variable: `interaction_error[i,j]`
* `ampl/scripts/run_register_design_fixed_layouts.py`
  * objective value extracted as `Squared_Frobenius_Error`
* `tools/generate_fixed_layouts_instance.py`
  * generates `Interaction[l,s,t]`, `Site_X`, and `Site_Y`

## 3. Python Heuristic Cost

The main Python heuristic uses an asymmetric cost:

$$
F_{\text{heur}}(r)
=
\left(
\sum_{(i,j):Q_{ij}>0}
(U_{ij}-Q_{ij})^2
+
\sum_{(i,j):Q_{ij}=0}
\max(U_{ij}-\tau,0)^2
\right)^{1/2}
$$

where:

$$
U_{ij}=\frac{C_6}{\|r_i-r_j\|^6}
$$

and `tau` is a tolerance threshold for zero target interactions.

The lower bound is:

$$
\boxed{LB_{\text{heur}} = 0}
$$

If a minimum physical distance `d_min > 0` is enforced, then:

$$
U_{ij} \le \frac{C_6}{d_{\min}^6}
$$

and a valid upper bound is:

$$
\boxed{
UB_{\text{heur}}
=
\left(
\sum_{(i,j):Q_{ij}>0}
\max
\left\{
|Q_{ij}|^2,
\left(\frac{C_6}{d_{\min}^6}-Q_{ij}\right)^2
\right\}
+
\sum_{(i,j):Q_{ij}=0}
\max
\left(
\frac{C_6}{d_{\min}^6}-\tau,
0
\right)^2
\right)^{1/2}
}
$$

This cost appears in:

* `main.py`
  * function: `evaluate_mapping`
* `ampl/scripts/run_register_design.py`
  * function: `evaluate_mapping`

## 4. Classical QUBO Objective

For the binary QUBO problem:

$$
f(x)=x^TQx,
\quad x\in\{0,1\}^N
$$

Because the feasible set is finite, simple universal bounds can be written in terms of the signs of the matrix entries.

Let:

$$
Q_{ij}^{+}=\max(Q_{ij},0),
\quad
Q_{ij}^{-}=\min(Q_{ij},0)
$$

Then:

$$
\boxed{
LB_{\text{QUBO}}
=
\sum_{i,j} Q_{ij}^{-}
}
$$

and:

$$
\boxed{
UB_{\text{QUBO}}
=
\sum_{i,j} Q_{ij}^{+}
}
$$

These bounds are conservative because not every subset of terms can be activated independently, but they are valid for any arbitrary instance.

A tighter constructive upper bound is obtained by evaluating any feasible binary vector, for example `x=0`:

$$
f(0)=0
$$

so:

$$
\boxed{UB_{\text{QUBO}} = 0}
$$

is a valid constructive upper bound when the goal is minimization and `x=0` is allowed.

This objective appears in:

* `ampl/qubo.mod`
  * objective: `QUBO_Objective`
* `ampl/scripts/run_qubo.py`
  * objective value extracted as `QUBO_Objective`

## Summary

For the register-design objectives, the most important universal bound is:

$$
\boxed{0 \le F}
$$

because all proposed costs are sums of squared residuals. Upper bounds are obtained either analytically from the maximum possible physical interaction or constructively by evaluating any feasible register assignment.

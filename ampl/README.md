# AMPL Models

This directory contains the exact mathematical optimization models used in the project, formulated in **AMPL (A Mathematical Programming Language)**.

## What is AMPL?

[AMPL](https://ampl.com) is an algebraic modeling language used for describing and solving high-complexity problems for large-scale mathematical computing (i.e., large-scale optimization and scheduling-type problems).

In this project, AMPL is heavily utilized for **rapid prototyping** and exact formulation of the register design and subset partitioning problems. Its main advantages are:

1. **Algebraic Syntax:** The syntax closely mirrors traditional mathematical notation (e.g., using $\sum$ equivalents, natural indexing). This makes it extremely fast to translate equations from a whiteboard directly into code.
2. **Separation of Model and Data:** The logic of the problem (variables, constraints, objective) is written in `.mod` files, completely separate from the instance data (`.dat` files). This allows testing multiple data configurations without changing the model.
3. **Solver Independence:** An AMPL model can be sent to dozens of different solvers (like Gurobi, CPLEX, IPOPT, Highs, or local heuristical solvers) simply by changing a single configuration string, without rewriting any logic.

## Directory Structure

* **`*.mod` files:** These are the core algebraic models containing the parameters, variables, constraints, and objective functions.
  * `set_partitioning.mod`: Original MILP formulation for subset partitioning.
  * `qubo.mod`: Binary quadratic model representation.
  * `register_design.mod`: Free-space continuous formulation.
  * `register_design_fixed_layouts.mod`: Fixed-layout (discrete assignment) formulation.
* **`scripts/`:** Contains Python wrappers using `amplpy` (the AMPL Python API) to load the models, inject data, trigger solvers (e.g., Gurobi), and extract the results back into Python structures for plotting or analysis.

## References & Further Reading

If you are unfamiliar with AMPL or its Python API, here are some useful resources:

* **Official AMPL Website:** [ampl.com](https://ampl.com/)
* **The AMPL Book:** [A comprehensive guide to AMPL modeling](https://ampl.com/learn/ampl-book/)
* **amplpy Documentation:** [Python API for AMPL](https://amplpy.ampl.com/)
* **AMPL Model Library:** [Examples of classic optimization problems](https://ampl.com/resources/models/)

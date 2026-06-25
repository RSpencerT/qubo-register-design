# Tools

This directory contains utility scripts to generate instances and convert between mathematical models.

## Available Scripts

### `set_partitioning_to_qubo.py`
Converts an AMPL Set Partitioning `.dat` instance into a dense QUBO matrix saved as a CSV file. This is a deterministic mathematical transformation that replaces exact-cover constraints with quadratic penalty terms.

**Usage:**
```bash
.venv/bin/python tools/set_partitioning_to_qubo.py \
    --input data/set_partitioning/ampl/original.dat \
    --output data/set_partitioning/ampl/qubo_result.dat \
    --penalty-factor 10.0
```

### `qubo_to_set_partitioning.py`
Attempts to convert a generic QUBO matrix back into a Set Partitioning `.dat` instance. **Note that this conversion is non-deterministic.** Because the Set Partitioning constraints (binary incidence matrix) are heavily structured, not every QUBO matrix can be factored back into an exact-cover incidence matrix. Finding a valid decomposition requires solving a complex optimization problem itself, which is what this script attempts to do.

**Usage:**
```bash
.venv/bin/python tools/qubo_to_set_partitioning.py \
    --input data/set_partitioning/csv/matrix.csv \
    --output data/set_partitioning/ampl/new_set.dat \
    --output-qubo data/set_partitioning/ampl/repaired_qubo.dat
```

### `generate_fixed_layouts_instance.py`
Generates an AMPL `.dat` file containing fixed layout configurations (e.g., lines, grids, circles) and precomputes the physical Rydberg interactions ($I_{\ell st} = C_6/r^6$) between all trapping sites. This data is used by the `register_design_fixed_layouts.mod` exact solver to find the optimal atom-to-site assignments.

**Usage:**
```bash
.venv/bin/python tools/generate_fixed_layouts_instance.py \
    --input data/set_partitioning/csv/matrix.csv \
    --output data/set_partitioning/ampl/layouts.dat \
    --spacing 7.5
```

---

## How the Conversion Works

The project bridges classical optimization with quantum neutral-atom hardware by using the QUBO formulation as an intermediate language.

### 1. Set Partitioning to QUBO (Deterministic)

The Set Partitioning problem is defined as:

$$\min \sum_{j=1}^{n} c_j x_j \quad \text{s.t.} \quad \sum_{j=1}^{n} A_{ij} x_j = 1, \quad \forall i \in \{1, \dots, m\}$$

To convert this to a QUBO ($\mathbf{x}^T Q \mathbf{x}$), the equality constraints are folded into the objective function using a penalty factor $\lambda$:

$$f(x) = \sum_{j=1}^{n} c_j x_j + \lambda \sum_{i=1}^{m} \left( \sum_{j=1}^{n} A_{ij} x_j - 1 \right)^2$$

Expanding this gives the QUBO matrix $Q$:
* **Off-Diagonal (Overlap):** $Q_{jk} = \lambda (A^T A)_{jk}$ for $j \neq k$. These determine the required interaction strengths ($C_6/r^6$) between atoms.
* **Diagonal (Cost + Degree):** $Q_{jj} = c_j - \lambda \sum_{i=1}^{m} A_{ij}$. These are applied as local detuning ($\Delta_j$) in the quantum hardware.

This transformation ensures that the ground state of the physical Hamiltonian corresponds to the optimal solution of the original Set Partitioning problem.

### 2. QUBO to Set Partitioning (Non-Deterministic)

Going backwards—from a generic QUBO matrix $Q$ to an incidence matrix $A$—is much harder. We must find an incidence matrix $A$ of size $m \times n$ (with binary entries) and a penalty $\lambda$ such that the off-diagonal terms satisfy:

$$\lambda (A^T A)_{jk} = Q_{jk}$$

This is essentially an integer matrix factorization problem. Unlike the forward transformation, this is **non-deterministic** and may not have an exact solution for an arbitrary $Q$. The script formulates this factorization as an optimization problem and attempts to find a matrix $A$ that minimizes the reconstruction error.
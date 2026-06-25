# QUBO Register Design

This project studies how a target QUBO interaction matrix can be represented by neutral-atom register geometries. It combines:

* **Python heuristics** for free-space and fixed-layout register design;
* **AMPL** models for exact optimization formulations;
* **Pulser** for Rydberg interaction constants and register visualization;
* **conversion tools** between Set Partitioning and QUBO instances.

The central modeling question is:

> Given a QUBO matrix $Q$, can we place atoms or assign them to calibrated trapping sites so that the physical interaction matrix $C_6/r^6$ approximates $Q$?

## Project Structure

```text
.
├── ampl/
│   ├── qubo.mod
│   ├── register_design.mod
│   ├── register_design_fixed_layouts.mod
│   ├── set_partitioning.mod
│   └── scripts/
│       ├── run_qubo.py
│       ├── run_register_design.py
│       ├── run_register_design_fixed_layouts.py
│       └── run_set_partitioning.py
├── data/
│   └── set_partitioning/
│       ├── ampl/
│       └── csv/
├── heuristics/
│   └── fixed_layout_register_design.py
├── modeling_questions/
│   └── Question_1.md ... Question_9.md
├── tools/
│   ├── generate_fixed_layouts_instance.py
│   ├── qubo_to_set_partitioning.py
│   ├── set_partitioning_to_qubo.py
│   └── README.md
├── main.py
├── baseline.py
└── requirements.txt
```

## Main Components

### Free-Space Register Design

The free-space version optimizes continuous atom coordinates directly.

Relevant files:

* `main.py`
* `baseline.py`
* `ampl/register_design.mod`
* `ampl/scripts/run_register_design.py`

### Fixed Layout Register Design

The fixed-layout version restricts atoms to a finite catalog of calibrated layouts and assigns each atom to one available trapping site.

Relevant files:

* `ampl/register_design_fixed_layouts.mod`
* `tools/generate_fixed_layouts_instance.py`
* `ampl/scripts/run_register_design_fixed_layouts.py`
* `heuristics/fixed_layout_register_design.py`
* `data/set_partitioning/ampl/modeling_test_qubo-fixed.dat`

### Set Partitioning and QUBO

The project also includes a Set Partitioning model, a classical QUBO model, and transformation tools between both representations.

Relevant files:

* `ampl/set_partitioning.mod`
* `ampl/qubo.mod`
* `ampl/scripts/run_set_partitioning.py`
* `ampl/scripts/run_qubo.py`
* `tools/set_partitioning_to_qubo.py`
* `tools/qubo_to_set_partitioning.py`
* `tools/README.md`

## Setup

Create the virtual environment and install dependencies with:

```bash
make setup
```

This runs:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/pip install -r requirements.txt
```

Manual setup is also possible:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

## Running the Workflows

The project includes a `Makefile` with shortcut targets to run the primary modeling workflows:

Run the continuous free-space model (`main.py`):

```bash
make run-main
```
*(Equivalent to: `.venv/bin/python main.py`)*

Run the multi-start simulated annealing heuristic for the Fixed Layout problem (Question 6):

```bash
make run-heuristic
```
*(Equivalent to: `.venv/bin/python heuristics/fixed_layout_register_design.py`)*

Run the exact AMPL model (MIQP) for the Fixed Layout problem:

```bash
make run-fixed-exact
```
*(Equivalent to: `.venv/bin/python ampl/scripts/run_register_design_fixed_layouts.py`)*

## Running AMPL Models

The AMPL runners are located in `ampl/scripts`.

Run the classical QUBO model:

```bash
.venv/bin/python ampl/scripts/run_qubo.py
```

Run the Set Partitioning model:

```bash
.venv/bin/python ampl/scripts/run_set_partitioning.py
```

Run the free-space AMPL register-design model:

```bash
.venv/bin/python ampl/scripts/run_register_design.py
```

Run the exact fixed-layout AMPL model:

```bash
.venv/bin/python ampl/scripts/run_register_design_fixed_layouts.py
```

These scripts use AMPL through `amplpy`. The default solver is currently `gurobi`.

## Fixed Layout Instance Generation

Generate an AMPL `.dat` instance for the fixed-layout register design model:

```bash
.venv/bin/python tools/generate_fixed_layouts_instance.py
```

By default, this reads:

```text
data/set_partitioning/csv/modeling_test_qubo.csv
```

and writes:

```text
data/set_partitioning/ampl/modeling_test_qubo-fixed.dat
```

You can customize the input QUBO, output file, number of calibrated sites, and site spacing:

```bash
.venv/bin/python tools/generate_fixed_layouts_instance.py \
  --input data/set_partitioning/csv/modeling_test_qubo.csv \
  --output data/set_partitioning/ampl/modeling_test_qubo-fixed.dat \
  --num-sites 6 \
  --spacing 7.5
```

The generated file is intended for:

```text
ampl/register_design_fixed_layouts.mod
```

## Fixed Layout Heuristic

The fixed-layout heuristic uses multi-start simulated annealing over layout choices and atom-to-site assignments.

Run it with:

```bash
.venv/bin/python heuristics/fixed_layout_register_design.py
```

For a terminal-only run without plotting:

```bash
.venv/bin/python heuristics/fixed_layout_register_design.py --no-plot
```

For a stronger run:

```bash
.venv/bin/python heuristics/fixed_layout_register_design.py \
  --restarts 100 \
  --iterations 20000 \
  --no-plot
```

## Modeling Questions

The folder `modeling_questions/` contains the written modeling discussion for the project:

* complexity reduction;
* cost functions and bounds;
* AMPL formulation;
* valid inequalities;
* heuristic design and implementation;
* UML architecture;
* heuristic complexity;
* Genetic Algorithm and Tabu Search discussion.

## Notes

* `make clean` removes the virtual environment and local Python cache folders.
* The AMPL scripts assume that the required AMPL modules are installed through `requirements.txt`.
* Some AMPL runs require a valid local AMPL/Gurobi setup.
* Use `ampl_license.example.txt` as a template for local AMPL activation. The real `ampl_license.txt` file is ignored by Git and should not be committed.

## License

This project is released under the MIT License. See `LICENSE` for details.

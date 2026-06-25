.PHONY: setup run-main run-heuristic run-fixed-exact lint format clean

setup:
	python3 -m venv .venv
	.venv/bin/python -m pip install --upgrade pip
	.venv/bin/pip install -r requirements.txt

run-main:
	.venv/bin/python main.py

run-heuristic:
	.venv/bin/python heuristics/fixed_layout_register_design.py

run-fixed-exact:
	.venv/bin/python ampl/scripts/run_register_design_fixed_layouts.py

lint:
	.venv/bin/ruff check .

format:
	.venv/bin/ruff format .

clean:
	rm -rf .venv __pycache__ .pytest_cache .mypy_cache .ruff_cache
	find . -name ".DS_Store" -delete

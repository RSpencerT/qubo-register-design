# ========================================================================
# Classical QUBO model (Full Symmetric Matrices)
# ========================================================================
#
# This model solves a Quadratic Unconstrained Binary Optimization problem:
#
#     minimize x' Q x
#
# Unlike the standard formulation that strictly requires an upper-triangular 
# matrix, this model iterates through all combinations of 'i' and 'j'. 

param N > 0 integer;
set VARS := 1..N;

# QUBO matrix
param Q {VARS, VARS} default 0;

# Binary decision variables.
var x {VARS} binary;

# Quadratic objective.
minimize QUBO_Objective:
    sum {i in VARS, j in VARS} Q[i,j] * x[i] * x[j];
# 1. Parameters and sets
param N > 0 integer;
set ATOMS := 1..N;

# Target interaction matrix (QUBO)
param Q {ATOMS, ATOMS} default 0.0;

# Physical hardware constants
param C6 := 5420158.53;
param L := 35.0;

# Safety parameter to avoid strict division by zero
param epsilon := 1e-6; 

# 2. Decision variables (continuous coordinates limited by the hardware)
var x {ATOMS} >= 0, <= L;
var y {ATOMS} >= 0, <= L;

# 3. Objective function (Frobenius norm)
# The distance d_ij^6 is exactly ((x_i - x_j)^2 + (y_i - y_j)^2)^3
minimize Frobenius_Error:
    sqrt(
        sum {i in ATOMS, j in ATOMS: i != j} 
        (C6 / (((x[i] - x[j])^2 + (y[i] - y[j])^2)^3 + epsilon) - Q[i,j])^2
    );

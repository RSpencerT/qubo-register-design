# Register design with fixed pre-calibrated layouts
#
# Unlike the free-space model, atoms are not placed at arbitrary continuous
# coordinates. A layout is selected from a finite catalog, and each atom is
# assigned to exactly one calibrated trapping site in that selected layout.

# 1. Parameters and sets
param N > 0 integer;
set ATOMS := 1..N;

set LAYOUTS;
set SITES;

# Sites available in each layout.
param Available {LAYOUTS, SITES} binary default 0;

# Calibrated site coordinates for each layout.
param Site_X {LAYOUTS, SITES} default 0.0;
param Site_Y {LAYOUTS, SITES} default 0.0;

# Target interaction matrix (QUBO)
param Q {ATOMS, ATOMS} default 0.0;

# Precomputed interaction strength between sites inside each layout.
# For a layout l and sites s,t:
#
#     Interaction[l,s,t] = C6 / distance(l,s,t)^6
#
# Set Interaction[l,s,s] to 0.
param Interaction {LAYOUTS, SITES, SITES} default 0.0;

# 2. Decision variables
var use_layout {LAYOUTS} binary;
var assign {LAYOUTS, ATOMS, SITES} binary;
var pair_assign {
    l in LAYOUTS,
    i in ATOMS,
    j in ATOMS,
    s in SITES,
    t in SITES:
        i < j and s != t
} binary;
var interaction_error {i in ATOMS, j in ATOMS: i < j};

# 3. Layout and assignment constraints
subject to Select_One_Layout:
    sum {l in LAYOUTS} use_layout[l] = 1;

subject to Assign_Each_Atom {i in ATOMS}:
    sum {l in LAYOUTS, s in SITES} assign[l,i,s] = 1;

subject to Use_Only_Selected_Layout_Sites {l in LAYOUTS, i in ATOMS, s in SITES}:
    assign[l,i,s] <= Available[l,s] * use_layout[l];

subject to At_Most_One_Atom_Per_Site {l in LAYOUTS, s in SITES}:
    sum {i in ATOMS} assign[l,i,s] <= Available[l,s] * use_layout[l];

subject to Pair_Assign_Upper_First {
    l in LAYOUTS,
    i in ATOMS,
    j in ATOMS,
    s in SITES,
    t in SITES:
        i < j and s != t
}:
    pair_assign[l,i,j,s,t] <= assign[l,i,s];

subject to Pair_Assign_Upper_Second {
    l in LAYOUTS,
    i in ATOMS,
    j in ATOMS,
    s in SITES,
    t in SITES:
        i < j and s != t
}:
    pair_assign[l,i,j,s,t] <= assign[l,j,t];

subject to Pair_Assign_Lower {
    l in LAYOUTS,
    i in ATOMS,
    j in ATOMS,
    s in SITES,
    t in SITES:
        i < j and s != t
}:
    pair_assign[l,i,j,s,t] >= assign[l,i,s] + assign[l,j,t] - 1;

subject to Interaction_Error_Definition {i in ATOMS, j in ATOMS: i < j}:
    interaction_error[i,j] =
        sum {l in LAYOUTS, s in SITES, t in SITES: s != t}
            Interaction[l,s,t] * pair_assign[l,i,j,s,t]
        - Q[i,j];

# 4. Objective function
#
# This objective compares the target interaction Q[i,j] with the interaction
# induced by the selected layout and atom-to-site assignment.
minimize Squared_Frobenius_Error:
    sum {i in ATOMS, j in ATOMS: i < j} interaction_error[i,j]^2;

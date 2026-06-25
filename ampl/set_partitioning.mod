# 1. Sets
set ITEMS;      # The universe of items that must be covered
set SUBSETS;    # The available subsets that can be selected

# 2. Parameters
# Incidence matrix: 1 if item 'i' belongs to subset 'j', 0 otherwise
param A {ITEMS, SUBSETS} binary; 

# Cost of each subset; in the standard problem this is often 1 to minimize the number of selected subsets
param Cost {SUBSETS} default 1;  

# 3. Decision variables
# x[j] = 1 if subset 'j' is selected, 0 otherwise
var x {SUBSETS} binary; 

# 4. Objective function
# Minimize the total cost of the selected subsets
minimize Total_Cost: 
    sum {j in SUBSETS} Cost[j] * x[j];

# 5. Classical constraints
# Exact cover constraint (Set Partitioning)
# Each item in the universe must belong to exactly one selected subset
subject to Exact_Cover {i in ITEMS}:
    sum {j in SUBSETS} A[i,j] * x[j] == 1;

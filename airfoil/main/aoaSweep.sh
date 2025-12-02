#!/bin/bash

# Define the list of X values
X_values=(0 8 10 20 28)

# Loop over each X value
for X in "${X_values[@]}"; do
    echo "Running simulation for X = $X"
    mpirun -np 16 ./airfoil "./xmlFiles/parameters_${X}.xml"
    echo "Finished simulation for X = $X"
    echo "----------------------------------------"
done

echo "All simulations completed."

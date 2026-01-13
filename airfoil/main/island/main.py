"""
This code runs shape optimization of the airfoil in order to obtain optimal lift and drag 
coefficients.

The cost function takes in as input a geometrical file describing the shape of interest, runs
palabos simulations to compute the lift and drag coefficients, and returns the lift to drag ratio.

We then run shape optimization using the scipy minimize function with the SLSQP method.
Of course, the entire question boils down to "how do we generate the geometrical file ?" 
We want to remain mostly unconstrained, benefiting from new techniques of additive manufacturing enabling us 
to print complex shapes. However, we still want to limit a few things:
1- The size of the domain that is allowed to change to contain solid parts. We do not want the optimizer to
2- The seed used for random generation (for reproducibility).
3- The minimum thickness of the solid shape (to avoid numerical instabilities in the solver and unphysical shapes).
4- The minimum thickness of the fluid channel (to avoid numerical instabilities in the solver and unphysical shapes).

A first idea is to start from a random shape, and allow it to morph using a limited number of
basis functions (e.g., Fourier modes). This way, we can control the complexity of the shape while still allowing for a wide variety of forms.
"""
import numpy as np
import random
from deap import base, creator, tools, algorithms

## Test elements
identity = np.eye(50)
identity_flat = identity.flatten()



## Define all required functions here

def simulate_ld_ratio(mask, N):
    return -1*(np.linalg.norm(mask-np.eye(N))**2) # Distance to identity as a placeholder for lift-to-drag ratio





if __name__ == "__main__":
    # --- Constants ---
    N, M = 15, 15  # Size of the mask
    POPULATION_SIZE = 1500
    GENERATIONS = 1000
    CX_PROB = 0.6  # Crossover probability
    MUT_PROB = 0.4  # Mutation probability


    
    
    # --- Fitness and Individual ---
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)

    # --- Toolbox (GA operations) ---
    toolbox = base.Toolbox()

    # Attribute generator: random boolean
    toolbox.register("attr_bool", random.randint, 0, 1)

    # Structure initializers
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, n=N*M)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    # --- Fitness function (your simulation) ---
    def evaluate(individual):
        mask = np.array(individual).reshape(N, M).astype(bool)
        ld_ratio = simulate_ld_ratio(mask, N)  # Replace with your function
        return (ld_ratio,)

    toolbox.register("evaluate", evaluate)
    
    
    toolbox.register("mate", tools.cxTwoPoint)  # Crossover
    toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)  # Mutation
    toolbox.register("select", tools.selTournament, tournsize=3)  # Selection

    # --- Main GA loop ---
    def main():
        # Initialize population
        pop = toolbox.population(n=POPULATION_SIZE)
        hof = tools.HallOfFame(1)  # Track the best individual

        # Run the GA
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", np.mean)
        stats.register("min", np.min)
        stats.register("max", np.max)

        pop, log = algorithms.eaMuCommaLambda(
            pop, toolbox, mu=POPULATION_SIZE//10, lambda_=POPULATION_SIZE, cxpb=CX_PROB, mutpb=MUT_PROB,
            ngen=GENERATIONS, stats=stats, halloffame=hof, verbose=True
        )

        # Return the best individual
        return hof[0]

    # --- Run the GA ---
    best_individual = main()
    best_mask = np.array(best_individual).reshape(N, M).astype(bool)
    print("Best mask shape:", best_mask)
    print("Best lift-to-drag ratio:", evaluate(best_individual)[0])
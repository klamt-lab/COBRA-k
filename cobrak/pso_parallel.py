from collections.abc import Callable
from copy import deepcopy
from os import cpu_count
from time import time

import numpy as np
from joblib import Parallel, delayed

from .io import json_write
from .utilities import last_n_elements_equal


# PSO class definition with genetic operators
class COBRAKPAPSO:
    def __init__(
        self,
        fitness_function: Callable[
            [list[float | int]], tuple[float, list[float | int]]
        ],
        xs_dim: int,
        gen: int,
        omega: float = 0.7298,
        eta1: float = 2.05,
        eta2: float = 2.05,
        max_vel: float = 0.5,
        extra_xs: list[list[float]] = [],
        seed: int | None = None,
        objvalue_json_path: str = "",
        max_rounds_same_objvalue: float = float("inf"),
        pop_size: int | None = None,
        diversity_threshold: float = 0.1,  # New parameter for diversity threshold
        diversity_restore_factor: float = 0.1,  # New parameter for diversity restoration
        ga_crossover_rate: float = 0.7,  # Crossover rate for genetic algorithm
        ga_mutation_rate: float = 0.01,  # Mutation rate for genetic algorithm
    ) -> None:
        # Parameters
        self.fitness_function = fitness_function
        self.xs_dim = xs_dim
        self.gen = gen
        self.omega = omega
        self.eta1 = eta1
        self.eta2 = eta2
        self.max_vel = max_vel
        self.seed = seed
        if seed is not None:
            np.random.seed(seed)

        # Initialization of random particles
        used_cpu_count = cpu_count() if pop_size is None else pop_size
        if used_cpu_count is None:
            used_cpu_count = 1
        self.xs = [
            np.random.uniform(0.0, 1.0, size=xs_dim).tolist()
            for _ in range(used_cpu_count - len(extra_xs))
        ]

        # Addition of user-defined extra particles
        if extra_xs != []:
            self.xs.extend(deepcopy(extra_xs))

        self.velocities = [
            np.random.uniform(-self.max_vel, self.max_vel, size=xs_dim).tolist()
            for _ in range(len(self.xs))
        ]
        self.personal_best_xs = self.xs.copy()
        self.personal_best_fs = [float("inf")] * len(self.xs)
        self.objvalue_json_path = objvalue_json_path
        self.objvalue_json_data: dict[float, list[float]] = {}
        self.max_rounds_same_objvalue = max_rounds_same_objvalue
        self.diversity_threshold = diversity_threshold
        self.diversity_restore_factor = diversity_restore_factor
        self.ga_crossover_rate = ga_crossover_rate
        self.ga_mutation_rate = ga_mutation_rate

    def run(self) -> tuple[float, list[float]]:
        with Parallel(n_jobs=-1) as parallel:
            # Initialization
            particle_xs = self.xs.copy()
            print("INIT START")
            particle_fs = parallel(
                delayed(self.fitness_function)(x) for x in particle_xs
            )
            if particle_fs is None:
                print("ERROR: Something went wrong during initialization!")
                raise ValueError
            particle_fs = [particle_fs[i][0] for i in range(len(particle_fs))]
            print("INIT END")
            best_swarm_x = particle_xs[particle_fs.index(min(particle_fs))].copy()
            best_f = min(particle_fs)

            # Update personal bests
            for i in range(len(particle_fs)):
                self.personal_best_fs[i] = particle_fs[i]
                self.personal_best_xs[i] = particle_xs[i].copy()

            if self.objvalue_json_path != "":
                start_time = time()
                self.objvalue_json_data[0.0] = deepcopy(self.personal_best_fs)
                json_write(self.objvalue_json_path, self.objvalue_json_data)

            # Actual algorithm
            max_objvalues = []
            for gen in range(self.gen):
                max_objvalues.append(max(self.personal_best_fs))
                if last_n_elements_equal(max_objvalues, self.max_rounds_same_objvalue):
                    break

                # Check diversity and restore if necessary
                if self.check_diversity(particle_xs) < self.diversity_threshold:
                    self.restore_diversity(particle_xs)

                # Genetic operators
                self.apply_genetic_operators(particle_xs, particle_fs)

                # Update velocities and positions in parallel
                results = parallel(
                    delayed(self.update_particle)(
                        self.personal_best_xs[current_particle],
                        particle_xs[current_particle],
                        self.velocities[current_particle],
                        best_swarm_x,
                    )
                    for current_particle in range(len(self.xs))
                )
                if results is None:
                    print("ERROR: Something went wrong during fitness calculations!")
                    raise ValueError

                # Unpack results
                for current_particle, (x, velocity, fitness) in enumerate(results):
                    particle_xs[current_particle] = x
                    self.velocities[current_particle] = velocity
                    particle_fs[current_particle] = fitness

                    # Update personal best
                    if (
                        particle_fs[current_particle]
                        < self.personal_best_fs[current_particle]
                    ):
                        self.personal_best_fs[current_particle] = particle_fs[
                            current_particle
                        ]
                        self.personal_best_xs[current_particle] = x.copy()

                    # Update the best swarm position if necessary
                    if particle_fs[current_particle] < best_f:
                        best_swarm_x = x.copy()
                        best_f = particle_fs[current_particle]

                if self.objvalue_json_path != "":
                    self.objvalue_json_data[time() - start_time] = deepcopy(
                        self.personal_best_fs
                    )
                    json_write(self.objvalue_json_path, self.objvalue_json_data)

        return best_f, best_swarm_x

    def update_particle(
        self,
        personal_best_x: list[float],
        old_x: list[float],
        velocity: list[float],
        best_swarm_x: list[float],
    ) -> tuple[list[float], list[float], float]:
        phi1 = np.random.uniform(0.0, 1.0, size=self.xs_dim)
        phi2 = np.random.uniform(0.0, 1.0, size=self.xs_dim)

        # Update velocity
        velocity = (
            self.omega * np.array(velocity)
            + self.eta1 * phi1 * (np.array(personal_best_x) - np.array(old_x))
            + self.eta2 * phi2 * (np.array(best_swarm_x) - np.array(old_x))
        ).tolist()

        # Clamp velocity
        velocity = [max(min(v, self.max_vel), -self.max_vel) for v in velocity]

        # Update position
        x = (np.array(old_x) + np.array(velocity)).tolist()

        # Boundary handling (reflection)
        x = [max(min(x_i, 1.0), 0.0) if x_i < 0.0 or x_i > 1.0 else x_i for x_i in x]

        # Update particle fitness
        fitness, _ = self.fitness_function(x)

        return x, velocity, fitness

    def check_diversity(self, particle_xs: list[list[float]]) -> float:
        # Calculate the average distance between particles
        distances = []
        for i in range(len(particle_xs)):
            for j in range(i + 1, len(particle_xs)):
                dist = np.linalg.norm(
                    np.array(particle_xs[i]) - np.array(particle_xs[j])
                )
                distances.append(dist)
        return np.mean(distances) if distances else 0.0

    def restore_diversity(self, particle_xs: list[list[float]]) -> None:
        # Perturb particles to restore diversity
        for i in range(len(particle_xs)):
            perturbation = np.random.uniform(
                -self.diversity_restore_factor,
                self.diversity_restore_factor,
                size=self.xs_dim,
            )
            particle_xs[i] = (np.array(particle_xs[i]) + perturbation).tolist()
            # Boundary handling (reflection)
            particle_xs[i] = [
                max(min(x_i, 1.0), 0.0) if x_i < 0.0 or x_i > 1.0 else x_i
                for x_i in particle_xs[i]
            ]

    def apply_genetic_operators(
        self, particle_xs: list[list[float]], particle_fs: list[float]
    ) -> None:
        # Selection: Select particles based on fitness
        sorted_indices = np.argsort(particle_fs)
        selected_indices = sorted_indices[: len(particle_xs) // 2]

        # Crossover: Combine selected particles to create new offspring
        offspring = []
        for i in range(0, len(selected_indices) - 1, 2):
            parent1 = particle_xs[selected_indices[i]]
            parent2 = particle_xs[selected_indices[i + 1]]
            if np.random.rand() < self.ga_crossover_rate:
                child1, child2 = self.crossover(parent1, parent2)
                offspring.extend((child1, child2))
            else:
                offspring.extend((parent1, parent2))

        # Mutation: Introduce random changes to some particles
        for i in range(len(offspring)):
            if np.random.rand() < self.ga_mutation_rate:
                offspring[i] = self.mutation(offspring[i])

        # Replace some particles with offspring
        for i in range(len(offspring)):
            particle_xs[i] = offspring[i]

    def crossover(
        self, parent1: list[float], parent2: list[float]
    ) -> tuple[list[float], list[float]]:
        # Single-point crossover
        crossover_point = np.random.randint(1, self.xs_dim - 1)
        child1 = parent1[:crossover_point] + parent2[crossover_point:]
        child2 = parent2[:crossover_point] + parent1[crossover_point:]
        return child1, child2

    def mutation(self, particle: list[float]) -> list[float]:
        # Random mutation
        mutated_particle = particle.copy()
        for i in range(self.xs_dim):
            if np.random.rand() < self.ga_mutation_rate:
                mutated_particle[i] += np.random.uniform(
                    -self.ga_mutation_rate, self.ga_mutation_rate
                )
                mutated_particle[i] = max(
                    min(mutated_particle[i], 1.0), 0.0
                )  # Boundary handling
        return mutated_particle

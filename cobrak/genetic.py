import operator
from collections.abc import Callable
from copy import deepcopy
from os import cpu_count
from random import choice, randint, uniform
from time import time

import numpy as np
from joblib import Parallel, delayed

from .io import json_write
from .utilities import count_last_equal_elements, last_n_elements_equal


class COBRAKGENETIC:
    def __init__(
        self,
        fitness_function: Callable[
            [list[float | int]], tuple[float, list[float | int]]
        ],
        xs_dim: int,
        gen: int,
        extra_xs: list[list[int]] = [],
        seed: int | None = None,
        objvalue_json_path: str = "",
        max_rounds_same_objvalue: float = float("inf"),
        pop_size: int | None = None,
    ) -> None:
        # Parameters
        self.fitness_function = fitness_function
        self.xs_dim = xs_dim
        self.gen = gen
        self.seed = seed
        if seed is not None:
            np.random.seed(seed)

        # Initialization of random particles
        cpu_count_value = cpu_count() if pop_size is None else pop_size
        if cpu_count_value is None:
            self.cpu_count = 1
        else:
            self.cpu_count = cpu_count_value
        self.init_xs = [
            [randint(0, 1) for _ in range(xs_dim)]
            for _ in range(self.cpu_count - len(extra_xs))
        ]

        # Addition of user-defined extra particles
        if extra_xs != []:
            self.init_xs.extend(extra_xs)

        self.tested_xs: dict[tuple[int, ...], float] = {}
        self.all_xs: dict[tuple[int, ...], float] = {}
        self.objvalue_json_path = objvalue_json_path
        self.objvalue_json_data: dict[float, list[float]] = {}
        self.max_rounds_same_objvalue = max_rounds_same_objvalue

    def _get_sorted_list_from_tested_xs(self) -> list[tuple[float, tuple[int, ...]]]:
        return sorted(
            [(fitness, x) for (x, fitness) in self.tested_xs.items()],
            key=operator.itemgetter(0),
        )

    def run(self) -> tuple[float, tuple[int, ...]]:
        init_fitnesses = Parallel(n_jobs=-1)(
            delayed(self.fitness_function)(x) for x in self.init_xs
        )
        if init_fitnesses is not None:
            self.tested_xs = {}
            for init_fitness in init_fitnesses:
                for fitness, xs in init_fitness:
                    self.tested_xs[tuple(xs)] = fitness
        else:
            print("ERROR: Something went wrong during initialization")
            raise ValueError

        if self.objvalue_json_path:
            start_time = time()
            self.objvalue_json_data[0.0] = sorted(self.tested_xs.values())
            json_write(self.objvalue_json_path, self.objvalue_json_data)

        # Actual algorithm
        max_objvalues = []
        for _ in range(self.gen):
            max_objvalues.append(max(self.tested_xs.values()))
            if last_n_elements_equal(max_objvalues, self.max_rounds_same_objvalue):  # type: ignore
                break

            xs_list = self._get_sorted_list_from_tested_xs()
            xs_list = [xs for xs in xs_list if len(xs) > 0]

            chosen_xs: list[tuple[int, ...]] = []
            # Choose some of the top 3
            for _ in range(self.cpu_count // 4):
                chosen_xs.append(choice([x[1] for x in xs_list if len(x[1]) > 0][:3]))

            # Choose some of the 25% best
            for _ in range(self.cpu_count // 2):
                chosen_xs.append(
                    choice([x[1] for x in xs_list][: round(len(xs_list) * 0.25)])
                )

            # Choose some other 75% worst
            for _ in range(self.cpu_count // 4):
                addlength = 0
                while not addlength:
                    added_xs = deepcopy(
                        choice([x[1] for x in xs_list][round(len(xs_list) * 0.25) :])
                    )
                    addlength = len(added_xs)
                chosen_xs.append(added_xs)

            # Random crossovers
            for _ in range(round(len(chosen_xs) * 0.2)):
                target = randint(0, len(chosen_xs) - 1)
                source = randint(0, len(chosen_xs) - 1)
                cut = randint(0, self.xs_dim - 1)
                chosen_xs[target] = deepcopy(chosen_xs[target][:cut]) + deepcopy(
                    chosen_xs[source][cut:]
                )

            # Test Xs in parallel
            results = Parallel(n_jobs=-1, verbose=10)(
                delayed(self.update_particle)(
                    chosen_x,
                    count_last_equal_elements(max_objvalues),
                )
                for chosen_x in chosen_xs
            )

            if results is None:
                print("ERROR: Something went wrong during fitness calculations")
                raise ValueError

            # Unpack results
            for fitnesses_and_active_xs, mutated_x in results:
                for fitness, active_x in fitnesses_and_active_xs:
                    if active_x is None or active_x == []:
                        continue
                    self.tested_xs[tuple(active_x)] = fitness
                    self.all_xs[tuple(active_x)] = fitness
                self.all_xs[tuple(mutated_x)] = max(
                    fitness for (fitness, _) in fitnesses_and_active_xs
                )

            if self.objvalue_json_path:
                self.objvalue_json_data[time() - start_time] = sorted(
                    self.tested_xs.values()
                )
                json_write(self.objvalue_json_path, self.objvalue_json_data)

        best_f_and_x = self._get_sorted_list_from_tested_xs()[0]
        return best_f_and_x[0], best_f_and_x[1]

    def update_particle(
        self,
        chosen_x: list[int],
        num_rounds_without_best_change: int,
    ) -> tuple[float, list[int], list[int]]:
        if not len(chosen_x):
            return [[1_000_000, []]], []

        min_change_p = 0.1 * 0.95**num_rounds_without_best_change
        max_change_p = 0.1 * 1.05**num_rounds_without_best_change
        change_p = uniform(min_change_p, max_change_p)
        change_p = max(0.001, change_p)
        change_p = min(0.999, change_p)

        mutation_tries = 0
        while True:
            mutated_x: list[float] = []
            match randint(0, 2):
                case 0:  # Extend
                    for x in chosen_x:
                        if x == 1:
                            mutated_x.append(1)
                            continue
                        if uniform(0.0, 1.0) < change_p:
                            mutated_x.append(1)
                        else:
                            mutated_x.append(x)
                case 1:  # Decrease
                    for x in chosen_x:
                        if x == 0:
                            mutated_x.append(0)
                            continue
                        if uniform(0.0, 1.0) < change_p:
                            mutated_x.append(0)
                        else:
                            mutated_x.append(x)
                case 2:  # Extend and decrease
                    for x in chosen_x:
                        if x == 1:
                            if uniform(0.0, 1.0) < change_p:
                                mutated_x.append(0)
                            else:
                                mutated_x.append(x)
                        else:
                            if uniform(0.0, 1.0) < change_p:
                                mutated_x.append(1)
                            else:
                                mutated_x.append(x)
                # case 3:  # Random
                #    mutated_x = [randint(0, 1) for _ in range(len(chosen_x))]

            if mutation_tries > 250:
                return [[1_000_000, []]], []
            if tuple(mutated_x) in self.all_xs:
                mutation_tries += 1
            else:
                break

        # Evaluate new position
        fitnesses_and_active_xs = self.fitness_function(mutated_x)

        return fitnesses_and_active_xs, mutated_x

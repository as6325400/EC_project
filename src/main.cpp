#include <bits/stdc++.h>

#include "ea.hpp"
#include "problem.hpp"

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    ProblemData problem = load_problem(std::cin);

    EAConfig config;
    config.population_size = 200;
    config.generations = 500;
    config.crossover_rate = 0.9;
    config.mutation_rate = 0.05;
    config.mutation_sigma = 0.02;

    AssignmentResult best = run_evolutionary_solver(problem, config);

    std::cout << best.assignments.size() << '\n';
    for (const auto& [ud_id, ris_id] : best.assignments) {
        std::cout << ud_id << ' ' << ris_id << '\n';
    }
    return 0;
}

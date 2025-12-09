#pragma once

#include "evaluator.hpp"
#include "problem.hpp"

struct EAConfig {
    int population_size = 40;
    int generations = 150;
    double crossover_rate = 0.9;
    double mutation_rate = 0.05;
    double mutation_sigma = 0.05;
    double min_cost_penalty = 0.0;
    double max_cost_penalty = 2.0;
    double penalty_mutation_sigma = 0.05;
    double usage_penalty = 0.0;
};

AssignmentResult run_evolutionary_solver(
    const ProblemData& problem,
    const EAConfig& config
);

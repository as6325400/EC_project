#pragma once

#include "evaluator.hpp"
#include "problem.hpp"

struct EAConfig {
    int population_size = 40;
    int generations = 150;
    double crossover_rate = 0.9;
    double mutation_rate = 0.05;
    double mutation_sigma = 0.05;
};

AssignmentResult run_evolutionary_solver(
    const ProblemData& problem,
    const EAConfig& config
);

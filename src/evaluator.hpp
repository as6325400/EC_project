#pragma once

#include <vector>

#include "problem.hpp"

struct AssignmentResult {
    std::vector<std::pair<int, int>> assignments;  // <UD id, RIS id/-1>
    std::vector<int> chosen_ris_index;             // per UD index, -2 unserved, -1 direct, >=0 specific RIS
    std::vector<double> used_cost_per_ud;          // consumed qubits per UD (0 if unserved)
    double total_profit{};
    double fitness{};
    double remaining_qubits{};
    double used_qubits{};
};

struct EvaluationParams {
    double cost_penalty = 0.0;
    double usage_penalty = 0.0;
};

AssignmentResult evaluate_priority(
    const ProblemData& problem,
    const std::vector<double>& priority,
    const EvaluationParams& params
);

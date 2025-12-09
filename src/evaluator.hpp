#pragma once

#include <vector>

#include "problem.hpp"

struct AssignmentResult {
    std::vector<std::pair<int, int>> assignments;  // <UD id, RIS id/-1>
    std::vector<int> chosen_ris_index;             // per UD index, -2 unserved, -1 direct, >=0 specific RIS
    double total_profit{};
    double fitness{};
    double remaining_qubits{};
};

AssignmentResult evaluate_priority(
    const ProblemData& problem,
    const std::vector<double>& priority
);

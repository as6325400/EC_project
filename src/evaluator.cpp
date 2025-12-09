#include "evaluator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace {

constexpr double EPS = 1e-9;

struct Candidate {
    int ris_index{-2};  // -1 direct, >=0 actual RIS index
    double cost{std::numeric_limits<double>::infinity()};
    double score{-std::numeric_limits<double>::infinity()};
};

bool better_candidate(const Candidate& lhs, const Candidate& rhs) {
    if (lhs.score != rhs.score) return lhs.score > rhs.score;
    if (lhs.cost != rhs.cost) return lhs.cost < rhs.cost;
    return lhs.ris_index < rhs.ris_index;
}

Candidate make_candidate(
    const ProblemData& problem,
    int ud_index,
    int ris_index,
    double cost_penalty
) {
    const int direct_idx = problem.direct_ris_index;
    double cost = std::numeric_limits<double>::infinity();
    if (ris_index == -1 && direct_idx >= 0) {
        cost = problem.cost_table[ud_index][direct_idx];
    } else if (ris_index >= 0 && ris_index < problem.actual_ris_count) {
        cost = problem.cost_table[ud_index][ris_index];
    }
    if (!std::isfinite(cost)) {
        return Candidate{};
    }
    Candidate cand;
    cand.ris_index = ris_index;
    cand.cost = cost;
    cand.score = static_cast<double>(problem.uds[ud_index].profit) - cost_penalty * cost;
    return cand;
}

Candidate best_available_candidate(
    const ProblemData& problem,
    int ud_index,
    double cost_penalty,
    const std::vector<bool>& ris_available
) {
    Candidate best = make_candidate(problem, ud_index, -1, cost_penalty);
    for (int ris_idx = 0; ris_idx < problem.actual_ris_count; ++ris_idx) {
        if (!ris_available[ris_idx]) continue;
        Candidate cand = make_candidate(problem, ud_index, ris_idx, cost_penalty);
        if (better_candidate(cand, best)) {
            best = cand;
        }
    }
    return best;
}

void assign_ud(
    AssignmentResult& result,
    int ud_index,
    const Candidate& choice,
    const ProblemData& problem,
    std::vector<bool>& ris_available,
    std::vector<int>& ris_owner,
    std::vector<double>& used_cost
) {
    result.remaining_qubits -= choice.cost;
    result.total_profit += problem.uds[ud_index].profit;
    result.chosen_ris_index[ud_index] = choice.ris_index;
    used_cost[ud_index] = choice.cost;
    if (choice.ris_index >= 0) {
        ris_available[choice.ris_index] = false;
        ris_owner[choice.ris_index] = ud_index;
    }
}

void unassign_ud(
    AssignmentResult& result,
    int ud_index,
    const ProblemData& problem,
    std::vector<bool>& ris_available,
    std::vector<int>& ris_owner,
    std::vector<double>& used_cost
) {
    const int ris_idx = result.chosen_ris_index[ud_index];
    if (ris_idx == -2) {
        return;
    }
    result.remaining_qubits += used_cost[ud_index];
    result.total_profit -= problem.uds[ud_index].profit;
    if (ris_idx >= 0) {
        ris_available[ris_idx] = true;
        ris_owner[ris_idx] = -1;
    }
    result.chosen_ris_index[ud_index] = -2;
    used_cost[ud_index] = 0.0;
}

}  // namespace

AssignmentResult evaluate_priority(
    const ProblemData& problem,
    const std::vector<double>& priority,
    const EvaluationParams& params
) {
    AssignmentResult result;
    const size_t ud_count = problem.uds.size();
    result.chosen_ris_index.assign(ud_count, -2);
    result.remaining_qubits = problem.qan.ent_gen_rate;

    std::vector<int> order(ud_count);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        if (priority[lhs] != priority[rhs]) {
            return priority[lhs] > priority[rhs];
        }
        if (problem.uds[lhs].profit != problem.uds[rhs].profit) {
            return problem.uds[lhs].profit > problem.uds[rhs].profit;
        }
        return problem.uds[lhs].id < problem.uds[rhs].id;
    });

    std::vector<bool> ris_available(problem.actual_ris_count, true);
    std::vector<int> ris_owner(problem.actual_ris_count, -1);
    std::vector<double> used_cost(ud_count, 0.0);

    for (int ud_index : order) {
        Candidate best = best_available_candidate(
            problem,
            ud_index,
            params.cost_penalty,
            ris_available
        );
        if (!std::isfinite(best.cost)) {
            continue;
        }
        if (best.cost - result.remaining_qubits > EPS) {
            continue;
        }
        assign_ud(
            result,
            ud_index,
            best,
            problem,
            ris_available,
            ris_owner,
            used_cost
        );
    }

    std::vector<int> unserved;
    for (size_t idx = 0; idx < ud_count; ++idx) {
        if (result.chosen_ris_index[idx] == -2) {
            unserved.push_back(static_cast<int>(idx));
        }
    }
    std::sort(unserved.begin(), unserved.end(), [&](int lhs, int rhs) {
        return problem.uds[lhs].profit > problem.uds[rhs].profit;
    });

    for (int ud_index : unserved) {
        Candidate best = best_available_candidate(
            problem,
            ud_index,
            params.cost_penalty,
            ris_available
        );
        if (std::isfinite(best.cost) && best.cost - result.remaining_qubits <= EPS) {
            assign_ud(
                result,
                ud_index,
                best,
                problem,
                ris_available,
                ris_owner,
                used_cost
            );
            continue;
        }

        bool improved = false;
        // Try taking an occupied RIS if current UD is more profitable.
        for (int ris_idx = 0; ris_idx < problem.actual_ris_count && !improved; ++ris_idx) {
            int occupant = ris_owner[ris_idx];
            if (occupant == -1) continue;
            Candidate cand = make_candidate(problem, ud_index, ris_idx, params.cost_penalty);
            if (!std::isfinite(cand.cost)) continue;
            const double freed_qubits = used_cost[occupant];
            if (cand.cost - (result.remaining_qubits + freed_qubits) > EPS) continue;
            if (problem.uds[ud_index].profit <= problem.uds[occupant].profit) continue;
            unassign_ud(
                result,
                occupant,
                problem,
                ris_available,
                ris_owner,
                used_cost
            );
            assign_ud(
                result,
                ud_index,
                cand,
                problem,
                ris_available,
                ris_owner,
                used_cost
            );
            improved = true;
        }
        if (improved) {
            continue;
        }

        // Try dropping a low-profit served UD to free qubits.
        Candidate feasible = best;
        if (!std::isfinite(feasible.cost)) {
            continue;
        }
        int victim = -1;
        double worst_ratio = std::numeric_limits<double>::infinity();
        for (size_t idx = 0; idx < ud_count; ++idx) {
            if (result.chosen_ris_index[idx] == -2) continue;
            if (problem.uds[ud_index].profit <= problem.uds[idx].profit) continue;
            if (used_cost[idx] <= EPS) continue;
            const double available_qubits = result.remaining_qubits + used_cost[idx];
            if (feasible.cost - available_qubits > EPS) continue;
            const double ratio = static_cast<double>(problem.uds[idx].profit)
                                 / (used_cost[idx] + EPS);
            if (ratio < worst_ratio) {
                worst_ratio = ratio;
                victim = static_cast<int>(idx);
            }
        }
        if (victim != -1) {
            unassign_ud(
                result,
                victim,
                problem,
                ris_available,
                ris_owner,
                used_cost
            );
            assign_ud(
                result,
                ud_index,
                feasible,
                problem,
                ris_available,
                ris_owner,
                used_cost
            );
        }
    }

    result.assignments.clear();
    for (size_t idx = 0; idx < ud_count; ++idx) {
        const int choice = result.chosen_ris_index[idx];
        if (choice == -2) continue;
        if (choice >= 0) {
            result.assignments.emplace_back(
                problem.uds[idx].id,
                problem.riss[choice].id
            );
        } else {
            result.assignments.emplace_back(problem.uds[idx].id, -1);
        }
    }

    result.used_qubits = problem.qan.ent_gen_rate - result.remaining_qubits;
    if (result.used_qubits < 0.0) {
        result.used_qubits = 0.0;
    }
    result.fitness = result.total_profit
                     - params.usage_penalty * result.used_qubits;
    return result;
}

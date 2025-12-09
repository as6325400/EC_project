#include "evaluator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace {

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

}  // namespace

AssignmentResult evaluate_priority(
    const ProblemData& problem,
    const std::vector<double>& priority
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

    std::vector<bool> ris_used(problem.actual_ris_count, false);
    const int direct_idx = problem.direct_ris_index;

    for (int ud_index : order) {
        const auto& ud = problem.uds[ud_index];
        Candidate best;

        // Direct connection candidate
        if (direct_idx >= 0) {
            const double cost = problem.cost_table[ud_index][direct_idx];
            if (std::isfinite(cost)) {
                const double score = (cost > 0.0)
                                     ? static_cast<double>(ud.profit) / cost
                                     : static_cast<double>(ud.profit);
                best = {-1, cost, score};
            }
        }

        // RIS candidates
        for (int ris_idx = 0; ris_idx < problem.actual_ris_count; ++ris_idx) {
            if (ris_used[ris_idx]) continue;
            const double cost = problem.cost_table[ud_index][ris_idx];
            if (!std::isfinite(cost)) continue;
            const double score = (cost > 0.0)
                                 ? static_cast<double>(ud.profit) / cost
                                 : static_cast<double>(ud.profit);
            Candidate cand{ris_idx, cost, score};
            if (better_candidate(cand, best)) {
                best = cand;
            }
        }

        if (!std::isfinite(best.cost)) {
            continue;  // no valid option -> not served
        }

        if (best.cost - result.remaining_qubits > 1e-9) {
            continue;  // insufficient qubits -> not served
        }

        result.remaining_qubits -= best.cost;
        result.total_profit += ud.profit;
        if (best.ris_index >= 0) {
            ris_used[best.ris_index] = true;
            result.assignments.emplace_back(ud.id, problem.riss[best.ris_index].id);
            result.chosen_ris_index[ud_index] = best.ris_index;
        } else {
            result.assignments.emplace_back(ud.id, -1);
            result.chosen_ris_index[ud_index] = -1;
        }
    }

    result.fitness = result.total_profit;
    return result;
}

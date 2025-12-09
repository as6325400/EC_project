#include "problem.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

constexpr double EPS = 1e-12;
constexpr int MAX_PURIFICATION_ROUNDS = 60;

double distance(const site& a, const site& b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return std::sqrt(dx * dx + dy * dy);
}

double link_success_prob(double dis, double alpha) {
    return std::exp(-alpha * dis);
}

double link_fidelity(double dis, double beta) {
    return 0.5 + 0.5 * std::exp(-beta * dis);
}

double purification_success_prob(double f1, double f2) {
    return f1 * f2 + (1.0 - f1) * (1.0 - f2);
}

double purification_fidelity(double f1, double f2) {
    const double numerator = f1 * f2;
    return numerator / (numerator + (1.0 - f1) * (1.0 - f2));
}

int required_rounds(double dis, double fidelity_threshold, double beta) {
    const double base_fidelity = link_fidelity(dis, beta);
    double current = base_fidelity;
    int rounds = 0;
    while (current + EPS < fidelity_threshold && current + EPS < 1.0) {
        current = purification_fidelity(current, base_fidelity);
        ++rounds;
        if (rounds > MAX_PURIFICATION_ROUNDS) {
            return -1;
        }
    }
    return rounds;
}

double purification_probability_product(double dis, int rounds, double beta) {
    if (rounds == 0) return 1.0;
    const double base_fidelity = link_fidelity(dis, beta);
    double probability = 1.0;
    double current = base_fidelity;
    for (int r = 1; r <= rounds; ++r) {
        const double success = purification_success_prob(current, base_fidelity);
        probability *= success;
        current = purification_fidelity(current, base_fidelity);
    }
    return probability;
}

double compute_cost(const ProblemData& problem, const UD& ud, const RIS& ris) {
    double dis{};
    if (ris.id == -1) {
        dis = distance(ud, problem.qan);
    } else {
        dis = distance(problem.qan, ris) + distance(ris, ud);
    }
    const double pe = link_success_prob(dis, problem.alpha);
    const int rounds = required_rounds(dis, ud.fidelity_th, problem.beta);
    if (rounds < 0) {
        return std::numeric_limits<double>::infinity();
    }
    const double purify_prob = purification_probability_product(dis, rounds, problem.beta);
    const double required = ud.exp_rate / pe * (rounds + 1.0) / purify_prob;
    return std::ceil(required - EPS);
}

}  // namespace

ProblemData load_problem(std::istream& in) {
    ProblemData problem;
    int UDs = 0, RISs = 0;
    if (!(in >> UDs >> RISs >> problem.alpha >> problem.beta >> problem.qan.ent_gen_rate)) {
        return problem;
    }
    in >> problem.qan.x >> problem.qan.y >> problem.qan.coverage_nums;
    for (int i = 0; i < problem.qan.coverage_nums; ++i) {
        int id;
        in >> id;
        problem.qan.coverage_UDs.insert(id);
    }

    problem.uds.resize(UDs);
    for (auto& ud : problem.uds) {
        in >> ud.id >> ud.x >> ud.y >> ud.profit >> ud.exp_rate >> ud.fidelity_th;
        ud.can_coverage = problem.qan.coverage_UDs.count(ud.id) > 0;
    }

    problem.riss.resize(RISs);
    for (auto& ris : problem.riss) {
        in >> ris.id >> ris.x >> ris.y >> ris.coverage_nums;
        for (int i = 0; i < ris.coverage_nums; ++i) {
            int id;
            in >> id;
            ris.coverage_UDs.insert(id);
        }
    }

    problem.actual_ris_count = RISs;
    RIS direct;
    direct.id = -1;
    direct.x = problem.qan.x;
    direct.y = problem.qan.y;
    direct.coverage_nums = problem.qan.coverage_nums;
    direct.coverage_UDs = problem.qan.coverage_UDs;
    problem.direct_ris_index = static_cast<int>(problem.riss.size());
    problem.riss.push_back(direct);

    const double INF = std::numeric_limits<double>::infinity();
    problem.cost_table.assign(
        problem.uds.size(),
        std::vector<double>(problem.riss.size(), INF)
    );
    problem.feasible_ris.assign(problem.uds.size(), {});

    for (size_t u = 0; u < problem.uds.size(); ++u) {
        for (size_t r = 0; r < problem.riss.size(); ++r) {
            const auto& ris = problem.riss[r];
            const bool can_use = (ris.id == -1)
                                 ? problem.uds[u].can_coverage
                                 : ris.coverage_UDs.count(problem.uds[u].id) > 0;
            if (!can_use) continue;
            problem.cost_table[u][r] = compute_cost(problem, problem.uds[u], ris);
            if (std::isfinite(problem.cost_table[u][r])) {
                problem.feasible_ris[u].push_back(static_cast<int>(r));
            }
        }
    }

    return problem;
}

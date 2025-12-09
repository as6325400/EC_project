#include "ea.hpp"

#include <algorithm>
#include <random>

namespace {

struct Individual {
    std::vector<double> priority;
    double cost_penalty{};
    AssignmentResult result;
};

std::mt19937_64& global_rng() {
    static std::mt19937_64 rng(std::random_device{}());
    return rng;
}

double gene_sigma(double base_sigma, const UD& ud) {
    return base_sigma * std::max(1.0, static_cast<double>(ud.profit));
}

Individual evaluate_individual(
    const ProblemData& problem,
    std::vector<double> priority,
    double cost_penalty,
    const EAConfig& config
) {
    Individual ind;
    ind.priority = std::move(priority);
    ind.cost_penalty = cost_penalty;
    EvaluationParams params;
    params.cost_penalty = cost_penalty;
    params.usage_penalty = config.usage_penalty;
    ind.result = evaluate_priority(problem, ind.priority, params);
    return ind;
}

const Individual& tournament_select(
    const std::vector<Individual>& population,
    int tournament_size,
    std::mt19937_64& rng
) {
    std::uniform_int_distribution<int> dist(0, static_cast<int>(population.size()) - 1);
    int best_idx = dist(rng);
    for (int i = 1; i < tournament_size; ++i) {
        int idx = dist(rng);
        if (population[idx].result.fitness > population[best_idx].result.fitness) {
            best_idx = idx;
        }
    }
    return population[best_idx];
}

void mutate_individual(
    Individual& ind,
    const ProblemData& problem,
    const EAConfig& config,
    std::mt19937_64& rng
) {
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
    std::normal_distribution<double> gauss(0.0, 1.0);
    for (size_t i = 0; i < ind.priority.size(); ++i) {
        if (prob_dist(rng) <= config.mutation_rate) {
            const double sigma = gene_sigma(config.mutation_sigma, problem.uds[i]);
            ind.priority[i] += gauss(rng) * sigma;
        }
    }
    // 交換兩個基因以探索排序空間。
    if (prob_dist(rng) <= config.mutation_rate) {
        std::uniform_int_distribution<int> idx_dist(0, static_cast<int>(ind.priority.size()) - 1);
        int a = idx_dist(rng);
        int b = idx_dist(rng);
        if (a != b) {
            std::swap(ind.priority[a], ind.priority[b]);
        }
    }
    if (prob_dist(rng) <= config.mutation_rate) {
        ind.cost_penalty += gauss(rng) * config.penalty_mutation_sigma;
        ind.cost_penalty = std::clamp(
            ind.cost_penalty,
            config.min_cost_penalty,
            config.max_cost_penalty
        );
    }
}

void local_search_swap(
    Individual& ind,
    const ProblemData& problem,
    const EAConfig& config,
    std::mt19937_64& rng
) {
    // 簡單爬山：嘗試少量交換並保留更佳解。
    std::uniform_int_distribution<int> idx_dist(0, static_cast<int>(ind.priority.size()) - 1);
    for (int step = 0; step < 2; ++step) {
        int a = idx_dist(rng);
        int b = idx_dist(rng);
        if (a == b) continue;
        std::swap(ind.priority[a], ind.priority[b]);
        EvaluationParams params;
        params.cost_penalty = ind.cost_penalty;
        params.usage_penalty = config.usage_penalty;
        AssignmentResult candidate = evaluate_priority(problem, ind.priority, params);
        if (candidate.fitness > ind.result.fitness) {
            ind.result = std::move(candidate);
        } else {
            std::swap(ind.priority[a], ind.priority[b]);
        }
    }
}

std::vector<double> random_priority(const ProblemData& problem, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> dist(-0.5, 0.5);
    std::vector<double> genes(problem.uds.size());
    for (size_t i = 0; i < problem.uds.size(); ++i) {
        genes[i] = static_cast<double>(problem.uds[i].profit) + dist(rng);
    }
    return genes;
}

double random_penalty(const EAConfig& config, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> dist(
        config.min_cost_penalty,
        config.max_cost_penalty
    );
    return dist(rng);
}

}  // namespace

AssignmentResult run_evolutionary_solver(
    const ProblemData& problem,
    const EAConfig& config
) {
    if (problem.uds.empty()) {
        return AssignmentResult{};
    }

    auto& rng = global_rng();
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);

    std::vector<Individual> population;
    population.reserve(config.population_size);
    for (int i = 0; i < config.population_size; ++i) {
        population.push_back(
            evaluate_individual(
                problem,
                random_priority(problem, rng),
                random_penalty(config, rng),
                config
            )
        );
    }

    Individual best = *std::max_element(
        population.begin(),
        population.end(),
        [](const Individual& lhs, const Individual& rhs) {
            return lhs.result.fitness < rhs.result.fitness;
        }
    );

    for (int gen = 0; gen < config.generations; ++gen) {
        std::vector<Individual> next_population;
        next_population.reserve(config.population_size);

        // Elitism
        next_population.push_back(best);

        while (static_cast<int>(next_population.size()) < config.population_size) {
            const Individual& parent1 = tournament_select(population, 3, rng);
            const Individual& parent2 = tournament_select(population, 3, rng);

            std::vector<double> child1 = parent1.priority;
            std::vector<double> child2 = parent2.priority;
            double child_penalty1 = parent1.cost_penalty;
            double child_penalty2 = parent2.cost_penalty;

            if (prob_dist(rng) <= config.crossover_rate) {
                std::uniform_real_distribution<double> alpha_dist(0.0, 1.0);
                const double alpha = alpha_dist(rng);
                for (size_t i = 0; i < child1.size(); ++i) {
                    const double a = parent1.priority[i];
                    const double b = parent2.priority[i];
                    child1[i] = alpha * a + (1.0 - alpha) * b;
                    child2[i] = (1.0 - alpha) * a + alpha * b;
                }
                child_penalty1 = alpha * parent1.cost_penalty
                                 + (1.0 - alpha) * parent2.cost_penalty;
                child_penalty2 = (1.0 - alpha) * parent1.cost_penalty
                                 + alpha * parent2.cost_penalty;
            }

            Individual offspring1;
            offspring1.priority = std::move(child1);
            offspring1.cost_penalty = child_penalty1;
            Individual offspring2;
            offspring2.priority = std::move(child2);
            offspring2.cost_penalty = child_penalty2;

            mutate_individual(offspring1, problem, config, rng);
            mutate_individual(offspring2, problem, config, rng);

            EvaluationParams params1;
            params1.cost_penalty = offspring1.cost_penalty;
            params1.usage_penalty = config.usage_penalty;
            EvaluationParams params2;
            params2.cost_penalty = offspring2.cost_penalty;
            params2.usage_penalty = config.usage_penalty;
            offspring1.result = evaluate_priority(problem, offspring1.priority, params1);
            offspring2.result = evaluate_priority(problem, offspring2.priority, params2);

            next_population.push_back(std::move(offspring1));
            if (static_cast<int>(next_population.size()) < config.population_size) {
                next_population.push_back(std::move(offspring2));
            }
        }

        // 對當代最優的少數個體做局部搜尋微調。
        const int elite_ls = std::min(3, static_cast<int>(next_population.size()));
        std::partial_sort(
            next_population.begin(),
            next_population.begin() + elite_ls,
            next_population.end(),
            [](const Individual& lhs, const Individual& rhs) {
                return lhs.result.fitness > rhs.result.fitness;
            }
        );
        for (int i = 0; i < elite_ls; ++i) {
            local_search_swap(next_population[i], problem, config, rng);
        }

        population = std::move(next_population);
        const Individual& generation_best = *std::max_element(
            population.begin(),
            population.end(),
            [](const Individual& lhs, const Individual& rhs) {
                return lhs.result.fitness < rhs.result.fitness;
            }
        );
        if (generation_best.result.fitness > best.result.fitness) {
            best = generation_best;
        }
    }

    return best.result;
}

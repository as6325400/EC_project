#include "ea.hpp"

#include <algorithm>
#include <random>

namespace {

struct Individual {
    std::vector<double> priority;
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
    std::vector<double> priority
) {
    Individual ind;
    ind.priority = std::move(priority);
    ind.result = evaluate_priority(problem, ind.priority);
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
}

std::vector<double> random_priority(const ProblemData& problem, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> dist(-0.5, 0.5);
    std::vector<double> genes(problem.uds.size());
    for (size_t i = 0; i < problem.uds.size(); ++i) {
        genes[i] = static_cast<double>(problem.uds[i].profit) + dist(rng);
    }
    return genes;
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
        population.push_back(evaluate_individual(problem, random_priority(problem, rng)));
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

            if (prob_dist(rng) <= config.crossover_rate) {
                std::uniform_real_distribution<double> alpha_dist(0.0, 1.0);
                const double alpha = alpha_dist(rng);
                for (size_t i = 0; i < child1.size(); ++i) {
                    const double a = parent1.priority[i];
                    const double b = parent2.priority[i];
                    child1[i] = alpha * a + (1.0 - alpha) * b;
                    child2[i] = (1.0 - alpha) * a + alpha * b;
                }
            }

            Individual offspring1;
            offspring1.priority = std::move(child1);
            Individual offspring2;
            offspring2.priority = std::move(child2);

            mutate_individual(offspring1, problem, config, rng);
            mutate_individual(offspring2, problem, config, rng);

            offspring1.result = evaluate_priority(problem, offspring1.priority);
            offspring2.result = evaluate_priority(problem, offspring2.priority);

            next_population.push_back(std::move(offspring1));
            if (static_cast<int>(next_population.size()) < config.population_size) {
                next_population.push_back(std::move(offspring2));
            }
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

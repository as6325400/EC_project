// ==================================================
// Headers
// ==================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdbool.h>

// ==================================================
// Abbreviations / Notations
// ==================================================

// KM     : Kuhn-Munkres Algorithm
// A, B   : The bipartition sets of the complete bipartite graph
// a, b   : Individual vertex in set A / set B

// ==================================================
// Constants
// ==================================================

#define INF 1e18
#define EPSILON 1e-9

// ==================================================
// Parameters
// ==================================================

#define QUEUE_INITIAL_MAX_SIZE 256
#define MAX_PURIFICATION_ROUND 25
#define PROFIT_SCALE 10000
#define WEIRD_ERROR_TERM_SCALE 100

// ==================================================
// Aliases
// ==================================================

typedef int64_t i64;
typedef struct general_params_t general_params_t;
typedef struct QAN_t QAN_t;
typedef struct UD_t UD_t;
typedef struct RIS_t RIS_t;
typedef struct i64_queue i64_queue;
typedef struct KM_data_t KM_data_t;

// ==================================================
// Component Parameter Structs
// ==================================================

struct general_params_t {
    // from the input
    i64 UD_num;
    i64 RIS_num;
    double alpha;
    double beta;
    i64 entangled_pairs_generation_rate;

    // Needs to be calculated
    i64 effective_RIS_num;
};

struct QAN_t {
    // from the input
    i64 x;
    i64 y;
    i64 covered_UD_num;
    i64* covered_UD_ID_array;
};

struct UD_t {
    // from the input
    i64 ID;
    i64 x;
    i64 y;
    i64 profit;
    i64 expected_rate;
    double fidelity_threshold;
};

struct RIS_t {
    // from the input
    i64 ID;
    i64 x;
    i64 y;
    i64 covered_UD_num;
    i64* covered_UD_ID_array;
};

// ==================================================
// Structs
// ==================================================

// This is a circular queue with dynamic queue max size allocation
struct i64_queue {
    i64 size;
    i64 max_size;
    i64 head;
    i64 tail;
    i64* queue;  // interval: [tail, head]
};

// Stores the data needed for Kuhn-Munkres algorithm.
struct KM_data_t {
    i64 set_size;                     // the size of each bipartition set
    i64 final_penalty_scale;          // If this equals to
                                      // penalty_scale_upper_bound + 1,
                                      // it means we won't select any RIS and UD
    i64 penalty_scale_upper_bound;    // used in binary search
    i64** least_initial_rate_matrix;
    i64** utility_matrix;
    i64* match_A;
    i64* match_B;
    i64* label_A;
    i64* label_B;
    i64* visited_A;
    i64* visited_B;
    i64_queue* queue;
    i64* even_vertex_parent_A;
    i64* slack_delta_B;
    i64* slack_parent_B;
};

// ==================================================
// Function Category Catalog
// ==================================================

// Utils
// Math Utils
// Queue
// Core Algorithm Utils
// Core Algorithm
// Input
// Initialize Input Data
// Output
// Memory Release

// ==================================================
// Utils
// ==================================================

i64* malloc_i64_1D_array(i64 size);

i64** malloc_i64_2D_array(i64 row_size, i64 col_size);

void fill_i64_1D_array(
    i64* i64_2D_array,
    i64 size,
    i64 value
);

void fill_i64_2D_array(
    i64** i64_2D_array,
    i64 row_size,
    i64 col_size,
    i64 value
);

// ==================================================
// Math Utils
// ==================================================

i64 max_i64(i64 num1, i64 num2);

i64 min_i64(i64 num1, i64 num2);

double calc_Euclidean_distance(double x1, double y1, double x2, double y2);

double calc_base_probability(double Euclidean_distance, double alpha);

double calc_base_fidelity(double Euclidean_distance, double beta);

double calc_purification_probability(double fidelity1, double fidelity2);

// Depends on the function "calc_purification_probability"
double calc_purification_fidelity(double fidelity1, double fidelity2);

// Depends on the functions
// "calc_base_probability", "calc_base_fidelity",
// "calc_purification_probability", "calc_purification_fidelity"
i64 calc_least_initial_rate(
    i64 expected_rate,
    double fidelity_threshold,
    double Euclidean_distance,
    double alpha,
    double beta
);

// Combines the "profit" and "least_initial_rate" by Lagrangian Relaxation
i64 calc_utility_metric(
    i64 profit,
    i64 least_initial_rate,
    i64 penalty_scale
);

// ==================================================
// Queue
// ==================================================

// Remember to use this function after you create a "i64_queue" immediately
void init_i64_queue(i64_queue* queue);

void _double_i64_queue_max_size(i64_queue* queue);

// Depends on the internal function "_double_i64_queue_max_size"
void push_i64_queue(i64_queue* queue, i64 value);

// Pops the front element and returns it
i64 pop_i64_queue(i64_queue* queue);

// Returns the front element without popping it
i64 get_i64_queue_front(const i64_queue* queue);

i64 get_i64_queue_size(const i64_queue* queue);

// ==================================================
// Core Algorithm Utils
// ==================================================

// KM_data->set_size = max_i64(
//     general_params->UD_num,
//     general_params->effective_RIS_num
// );
void init_KM_data_set_size(
    const general_params_t* general_params,
    KM_data_t* KM_data
);

// Depends on the functions
// "malloc_i64_2D_array", "fill_i64_2D_array",
// "calc_Euclidean_distance", "calc_least_initial_rate"
i64** gen_KM_data_least_initial_rate_matrix_and_set_penalty_scale_upper_bound(
    const general_params_t* general_params,
    const QAN_t* QAN,
    const UD_t* UDs,
    const RIS_t* RISs,
    KM_data_t* KM_data
);

// Depends on the functions
// "malloc_i64_2D_array", "fill_i64_2D_array", "calc_utility_metric"
i64** gen_KM_data_utility_matrix(
    const general_params_t* general_params,
    const QAN_t* QAN,
    const UD_t* UDs,
    const RIS_t* RISs,
    const KM_data_t* KM_data,
    i64 penalty_scale
);

// Depends on the functions
// "gen_KM_data_utility_matrix", "malloc_i64_1D_array"
KM_data_t* reset_KM_data_without_least_initial_rate_matrix_and_numeric_members(
    const general_params_t* general_params,
    const QAN_t* QAN,
    const UD_t* UDs,
    const RIS_t* RISs,
    KM_data_t* KM_data,
    i64 penalty_scale
);

i64 calc_KM_data_result_total_consumption_rate(
    const general_params_t* general_params,
    const KM_data_t* KM_data
);

// ==================================================
// Core Algorithm
// ==================================================
// This section will not have function dependency comments

void KM_relax_slack_delta_B(KM_data_t* KM_data, i64 a);

void KM_relabel(KM_data_t* KM_data);

void KM_augment_augmenting_path(KM_data_t* KM_data, i64 a, i64 b);

bool KM_traverse_alternating_tree_with_unvisited_equality_edge(
    KM_data_t* KM_data
);

bool KM_traverse_alternating_tree_with_new_equality_edge(KM_data_t* KM_data);

// Considers the bipartite graph as a complete bipartite graph.
// Nonexistent vertices are treated as dummy vertices.
// Nonexistent edges are assigned a weight of 0.
// Rows correspond to the vertices of RISs.
// Columns correspond to the vertices of UDs.
// Reference: https://web.ntnu.edu.tw/~algo/Matching.html
void Kuhn_Munkres(KM_data_t* KM_data);

// Binary searches the heuristic penalty scale.
// Applies Kuhn-Munkres Algorithm in each search round.
KM_data_t* binary_search_heuristic_penalty_scale_with_Kuhn_Munkres(
    const general_params_t* general_params,
    const QAN_t* QAN,
    const UD_t* UDs,
    const RIS_t* RISs
);

// ==================================================
// Input
// ==================================================

general_params_t* input_general_params();

// Depends on the function "malloc_i64_1D_array"
QAN_t* input_QAN();

UD_t* input_UDs(const general_params_t* general_params);

// Depends on the function "malloc_i64_1D_array"
RIS_t* input_RISs(const general_params_t* general_params, const QAN_t* QAN);

// ==================================================
// Initialize Input Data
// ==================================================

// general_params->effective_RIS_num =
//     general_params->RIS_num + QAN->covered_UD_num;
void init_general_params_effective_RIS_num(
    general_params_t* general_params,
    const QAN_t* QAN
);

// Depends on the function "malloc_i64_1D_array"
// Each UD directly reachable from the QAN is assigned a RIS dummy vertex
// located at the same position as the QAN
void init_RISs_dummy_vertices(
    const general_params_t* general_params,
    const QAN_t* QAN,
    const UD_t* UDs,
    RIS_t* RISs
);

// ==================================================
// Output
// ==================================================

void output_heuristic_solution_of_KM_data(
    const general_params_t* general_params,
    const UD_t* UDs,
    const KM_data_t* KM_data
);

// ==================================================
// Memory Release
// ==================================================

// "general_params" itself will be released
void free_general_params(general_params_t* general_params);

// "QAN" itself will be released
void free_QAN(QAN_t* QAN);

// "UDs" itself will be released
void free_UDs(UD_t* UDs);

// "RISs" itself will be released
void free_RISs(RIS_t* RISs, const general_params_t* general_params);

// "i64_2D_array" itself will be released
void free_i64_2D_array(i64** i64_2D_array, i64 row_num);

// "queue" itself will be released
void free_i64_queue(i64_queue* queue);

// "KM_data" itself won't be released
// Depends on the function "free_i64_2D_array"
void free_least_initial_rate_matrix_in_Kuhn_Munkres_data(
    KM_data_t* KM_data,
    i64 KM_data_size
);

// "KM_data" itself won't be released
// Depends on the function "free_i64_queue"
// The queue is released after each iteration of x in KM,
// so KM_data->queue is not released in the function
// "free_Kuhn_Munkres_data_without_least_initial_rate_matrix_and_queue"
void free_i64_queue_in_Kuhn_Munkres_data(KM_data_t* KM_data);

// "KM_data" itself won't be released
// Depends on the function "free_i64_2D_array"
void free_Kuhn_Munkres_data_without_least_initial_rate_matrix_and_queue(
    KM_data_t* KM_data
);

// Only "KM_data" itself will be released
void free_only_Kuhn_Munkres_data_itself(KM_data_t* KM_data);

// ==================================================
// Main
// ==================================================

int main() {
    // Input
    general_params_t* general_params = input_general_params();
    QAN_t* QAN = input_QAN();
    UD_t* UDs = input_UDs(general_params);
    RIS_t* RISs = input_RISs(general_params, QAN);

    // Initialize Input Data
    init_general_params_effective_RIS_num(general_params, QAN);
    init_RISs_dummy_vertices(general_params, QAN, UDs, RISs);

    // Core Algorithm
    KM_data_t* KM_data =
        binary_search_heuristic_penalty_scale_with_Kuhn_Munkres(
            general_params,
            QAN,
            UDs,
            RISs
        );

    // Output
    output_heuristic_solution_of_KM_data(general_params, UDs, KM_data);

    // Memory Release
    free_least_initial_rate_matrix_in_Kuhn_Munkres_data(
        KM_data,
        KM_data->set_size
    );
    free_Kuhn_Munkres_data_without_least_initial_rate_matrix_and_queue(KM_data);
    free_only_Kuhn_Munkres_data_itself(KM_data);

    free_RISs(RISs, general_params);
    free_general_params(general_params);
    free_QAN(QAN);
    free_UDs(UDs);

    return 0;
}

// ==================================================
// Utils
// ==================================================

i64* malloc_i64_1D_array(i64 size) {
    return (i64*)malloc(size * sizeof(i64));
}

i64** malloc_i64_2D_array(i64 row_size, i64 col_size) {
    i64** i64_2D_array = (i64**)malloc(row_size * sizeof(i64*));
    for (i64 i = 0; i < row_size; ++i) {
        i64_2D_array[i] = (i64*)malloc(col_size * sizeof(i64));
    }
    return i64_2D_array;
}

void fill_i64_1D_array(
    i64* i64_2D_array,
    i64 size,
    i64 value
) {
    for (i64 i = 0; i < size; ++i) {
        i64_2D_array[i] = value;
    }
}

void fill_i64_2D_array(
    i64** i64_2D_array,
    i64 row_size,
    i64 col_size,
    i64 value
) {
    for (i64 i = 0; i < row_size; ++i) {
        for (i64 j = 0; j < col_size; ++j) {
            i64_2D_array[i][j] = value;
        }
    }
}

// ==================================================
// Math Utils
// ==================================================

i64 max_i64(i64 num1, i64 num2) {
    return (num1 >= num2 ? num1 : num2);
}

i64 min_i64(i64 num1, i64 num2) {
    return (num1 <= num2 ? num1 : num2);
}

double calc_Euclidean_distance(double x1, double y1, double x2, double y2) {
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

double calc_base_probability(double Euclidean_distance, double alpha) {
    return exp((-1.0) * alpha * Euclidean_distance);
}

double calc_base_fidelity(double Euclidean_distance, double beta) {
    return 0.5 + 0.5 * exp((-1.0) * beta * Euclidean_distance);
}

double calc_purification_probability(double fidelity1, double fidelity2) {
    return fidelity1 * fidelity2 + (1.0 - fidelity1) * (1.0 - fidelity2);
}

double calc_purification_fidelity(double fidelity1, double fidelity2) {
    return fidelity1 * fidelity2
           / calc_purification_probability(fidelity1, fidelity2);
}

i64 calc_least_initial_rate(
    i64 expected_rate,
    double fidelity_threshold,
    double Euclidean_distance,
    double alpha,
    double beta
) {
    double base_probability = calc_base_probability(Euclidean_distance, alpha);
    double base_fidelity = calc_base_fidelity(Euclidean_distance, beta);

    i64 least_purification_num = 0;
    double cur_probability = 1.0;
    double cur_fidelity = base_fidelity;
    while (fabs(cur_fidelity - fidelity_threshold) >= EPSILON
           && cur_fidelity < fidelity_threshold) {
        least_purification_num += 1;
        if (least_purification_num > MAX_PURIFICATION_ROUND) return -1;

        cur_probability *= calc_purification_probability(
            cur_fidelity,
            base_fidelity
        );
        cur_fidelity = calc_purification_fidelity(
            cur_fidelity,
            base_fidelity
        );
    }

    double least_initial_rate =
        expected_rate * (least_purification_num + 1)
        / (base_probability * cur_probability);

    if (fabs(least_initial_rate - 0) < EPSILON) return 1;

    return (i64)ceil(least_initial_rate);
}

i64 calc_utility_metric(
    i64 profit,
    i64 least_initial_rate,
    i64 penalty_scale
) {
    /* return max_i64(  // 110975
        PROFIT_SCALE * profit - penalty_scale * least_initial_rate,
        0
    ); */

    // weird error term method
    return max_i64(  // 111065
        PROFIT_SCALE * profit - penalty_scale * least_initial_rate
        + WEIRD_ERROR_TERM_SCALE * (i64)round(
            sin((double)profit / (double)least_initial_rate)
        ),
        0
    );
}

// ==================================================
// Queue
// ==================================================

void init_i64_queue(i64_queue* queue) {
    queue->size = 0;
    queue->max_size = QUEUE_INITIAL_MAX_SIZE;
    queue->head = -1;
    queue->tail = 0;
    queue->queue = malloc_i64_1D_array(QUEUE_INITIAL_MAX_SIZE);
}

void _double_i64_queue_max_size(i64_queue* queue) {
    i64 old_queue_max_size = queue->max_size;

    i64* new_queue = (i64*)malloc((2 * old_queue_max_size) * sizeof(i64));

    i64 tail = queue->tail;
    for (i64 i = 0; i < queue->size; ++i) {
        new_queue[i] = queue->queue[tail];
        tail = (tail + 1) % old_queue_max_size;
    }

    queue->max_size = 2 * old_queue_max_size;

    queue->head = queue->size - 1;
    queue->tail = 0;

    free(queue->queue);
    queue->queue = new_queue;
}

void push_i64_queue(i64_queue* queue, i64 value) {
    if (queue->size == queue->max_size) {
        _double_i64_queue_max_size(queue);
    }
    queue->size += 1;
    queue->head = (queue->head + 1) % queue->max_size;
    queue->queue[queue->head] = value;
}

i64 pop_i64_queue(i64_queue* queue) {
    queue->size -= 1;
    i64 sought_value = queue->queue[queue->tail];
    queue->tail = (queue->tail + 1) % queue->max_size;
    return sought_value;
}

i64 get_i64_queue_front(const i64_queue* queue) {
    return queue->queue[queue->tail];
}

i64 get_i64_queue_size(const i64_queue* queue) {
    return queue->size;
}

// ==================================================
// Core Algorithm Utils
// ==================================================

void init_KM_data_set_size(
    const general_params_t* general_params,
    KM_data_t* KM_data
) {
    KM_data->set_size = max_i64(
        general_params->UD_num,
        general_params->effective_RIS_num
    );
}

i64** gen_KM_data_least_initial_rate_matrix_and_set_penalty_scale_upper_bound(
    const general_params_t* general_params,
    const QAN_t* QAN,
    const UD_t* UDs,
    const RIS_t* RISs,
    KM_data_t* KM_data
) {
    // Initialize the least initial rate matrix
    i64** least_initial_rate_matrix = malloc_i64_2D_array(
        KM_data->set_size,
        KM_data->set_size
    );
    fill_i64_2D_array(
        least_initial_rate_matrix,
        KM_data->set_size,
        KM_data->set_size,
        0
    );

    // Compute the least initial rate for each position (i, j)
    bool firstTime = true;
    i64 min_least_initial_rate = -1;
    i64 max_profit = -1;
    for (i64 a = 0; a < general_params->effective_RIS_num; ++a) {
        for (i64 k = 0; k < RISs[a].covered_UD_num; ++k) {
            i64 b = RISs[a].covered_UD_ID_array[k];
            
            double Euclidean_distance = calc_Euclidean_distance(
                UDs[b].x,
                UDs[b].y,
                RISs[a].x,
                RISs[a].y
            ) + calc_Euclidean_distance(
                RISs[a].x,
                RISs[a].y,
                QAN->x,
                QAN->y
            );

            least_initial_rate_matrix[a][b] = calc_least_initial_rate(
                UDs[b].expected_rate,
                UDs[b].fidelity_threshold,
                Euclidean_distance,
                general_params->alpha,
                general_params->beta
            );

            if (least_initial_rate_matrix[a][b] < 0) continue;
            if (firstTime
                || least_initial_rate_matrix[a][b] < min_least_initial_rate) {
                min_least_initial_rate = least_initial_rate_matrix[a][b];
            }
            if (firstTime
                || UDs[b].profit > max_profit) {
                max_profit = UDs[b].profit;
            }
            firstTime = false;
        }
    }

    // Calculate the penalty scale upper bound
    KM_data->penalty_scale_upper_bound = (
        min_least_initial_rate == -1 && max_profit == -1 ?
        0 :
        (i64)ceil(
            PROFIT_SCALE * max_profit / min_least_initial_rate
        ) + 1
    );

    return least_initial_rate_matrix;
}

i64** gen_KM_data_utility_matrix(
    const general_params_t* general_params,
    const QAN_t* QAN,
    const UD_t* UDs,
    const RIS_t* RISs,
    const KM_data_t* KM_data,
    i64 penalty_scale
) {
    // Initialize the utility matrix
    i64** utility_matrix = malloc_i64_2D_array(
        KM_data->set_size,
        KM_data->set_size
    );
    fill_i64_2D_array(
        utility_matrix,
        KM_data->set_size,
        KM_data->set_size,
        0
    );

    // Compute the utility metric for each position (i, j)
    for (i64 a = 0; a < general_params->effective_RIS_num; ++a) {
        for (i64 k = 0; k < RISs[a].covered_UD_num; ++k) {
            i64 b = RISs[a].covered_UD_ID_array[k];

            utility_matrix[a][b] = (
                KM_data->least_initial_rate_matrix[a][b] < 0 ?
                0 :
                calc_utility_metric(
                    UDs[b].profit,
                    KM_data->least_initial_rate_matrix[a][b],
                    penalty_scale
                )
            );
        }
    }

    return utility_matrix;
}

KM_data_t* reset_KM_data_without_least_initial_rate_matrix_and_numeric_members(
    const general_params_t* general_params,
    const QAN_t* QAN,
    const UD_t* UDs,
    const RIS_t* RISs,
    KM_data_t* KM_data,
    i64 penalty_scale
) {
    KM_data->utility_matrix = gen_KM_data_utility_matrix(
        general_params,
        QAN,
        UDs,
        RISs,
        KM_data,
        penalty_scale
    );
    KM_data->match_A = malloc_i64_1D_array(KM_data->set_size);
    KM_data->match_B = malloc_i64_1D_array(KM_data->set_size);
    KM_data->label_A = malloc_i64_1D_array(KM_data->set_size);
    KM_data->label_B = malloc_i64_1D_array(KM_data->set_size);
    KM_data->visited_A = malloc_i64_1D_array(KM_data->set_size);
    KM_data->visited_B = malloc_i64_1D_array(KM_data->set_size);
    KM_data->queue = NULL;
    KM_data->even_vertex_parent_A = malloc_i64_1D_array(KM_data->set_size);
    KM_data->slack_delta_B = malloc_i64_1D_array(KM_data->set_size);
    KM_data->slack_parent_B = malloc_i64_1D_array(KM_data->set_size);

    return KM_data;
}

i64 calc_KM_data_result_total_consumption_rate(
    const general_params_t* general_params,
    const KM_data_t* KM_data
) {
    i64 total_consumption_rate = 0;
    for (i64 effective_RIS_ID = 0;
         effective_RIS_ID < general_params->effective_RIS_num;
         ++effective_RIS_ID) {
        if (KM_data->match_A[effective_RIS_ID] == -1) continue;

        i64 cur_selected_UD_ID = KM_data->match_A[effective_RIS_ID];
        
        if (KM_data->utility_matrix
                [effective_RIS_ID]
                [cur_selected_UD_ID] == 0) continue;  // nonexistent edge

        total_consumption_rate +=
            KM_data->least_initial_rate_matrix
                [effective_RIS_ID]
                [cur_selected_UD_ID];
    }
    return total_consumption_rate;
}

// ==================================================
// Core Algorithm
// ==================================================

void KM_relax_slack_delta_B(KM_data_t* KM_data, i64 a) {
    for (i64 b = 0; b < KM_data->set_size; ++b) {
        if (KM_data->visited_B[b]) continue;
        if (KM_data->label_A[a]
            + KM_data->label_B[b]
            - KM_data->utility_matrix[a][b]
            < KM_data->slack_delta_B[b]) {
            KM_data->slack_delta_B[b] =
                KM_data->label_A[a]
                + KM_data->label_B[b]
                - KM_data->utility_matrix[a][b];
            KM_data->slack_parent_B[b] = a;
        }
    }
}

void KM_relabel(KM_data_t* KM_data) {
    i64 min_slack_delta = -100000;
    for (i64 b = 0; b < KM_data->set_size; ++b) {
        if (!KM_data->visited_B[b]) {
            min_slack_delta = (
                min_slack_delta == -100000 ?
                KM_data->slack_delta_B[b] :
                min_i64(
                    min_slack_delta,
                    KM_data->slack_delta_B[b]
                )
            );
        }
    }

    for (i64 a = 0; a < KM_data->set_size; ++a) {
        if (KM_data->visited_A[a]) {
            KM_data->label_A[a] -= min_slack_delta;
        }
    }
    for (i64 b = 0; b < KM_data->set_size; ++b) {
        if (KM_data->visited_B[b]) {
            KM_data->label_B[b] += min_slack_delta;
        }
    }
    for (i64 b = 0; b < KM_data->set_size; ++b) {
        if (!KM_data->visited_B[b]) {
            KM_data->slack_delta_B[b] -= min_slack_delta;
        }
    }
}

void KM_augment_augmenting_path(KM_data_t* KM_data, i64 a, i64 b) {
    if (a == -1 && b == -1) return;

    i64 previous_vertex_of_a = KM_data->match_A[a];
    KM_data->match_A[a] = b;
    KM_data->match_B[b] = a;

    KM_augment_augmenting_path(
        KM_data,
        KM_data->even_vertex_parent_A[a],
        previous_vertex_of_a
    );
}

bool KM_traverse_alternating_tree_with_unvisited_equality_edge(
    KM_data_t* KM_data
) {
    while (get_i64_queue_size(KM_data->queue) > 0) {
        i64 a = pop_i64_queue(KM_data->queue);
        for (i64 b = 0; b < KM_data->set_size; ++b) {
            if (!KM_data->visited_B[b]
                && KM_data->label_A[a] + KM_data->label_B[b]
                   == KM_data->utility_matrix[a][b]) {
                
                KM_data->visited_B[b] = true;

                if (KM_data->match_B[b] == -1) {  // found augmenting path
                    KM_augment_augmenting_path(
                        KM_data,
                        a,
                        b
                    );
                    return true;
                }

                i64 next_vertex = KM_data->match_B[b];
                push_i64_queue(KM_data->queue, next_vertex);
                KM_data->even_vertex_parent_A[next_vertex] = a;
                KM_data->visited_A[next_vertex] = true;
                KM_relax_slack_delta_B(KM_data, next_vertex);
            }
        }
    }
    return false;
}

bool KM_traverse_alternating_tree_with_new_equality_edge(KM_data_t* KM_data) {
    for (i64 b = 0; b < KM_data->set_size; ++b) {
        if (!KM_data->visited_B[b] && KM_data->slack_delta_B[b] == 0) {
            KM_data->visited_B[b] = true;

            if (KM_data->match_B[b] == -1) {  // found augmenting path
                KM_augment_augmenting_path(
                    KM_data,
                    KM_data->slack_parent_B[b],
                    b
                );
                return true;
            }

            i64 next_vertex = KM_data->match_B[b];
            push_i64_queue(KM_data->queue, next_vertex);
            KM_data->even_vertex_parent_A[next_vertex]
                = KM_data->slack_parent_B[b];
            KM_data->visited_A[next_vertex] = true;
            KM_relax_slack_delta_B(KM_data, next_vertex);
        }
    }
    return false;
}

void Kuhn_Munkres(KM_data_t* KM_data) {
    // Initialize label_A[a] as the max utility from a to b,
    // where there is an edge between a and b
    for (i64 a = 0; a < KM_data->set_size; ++a) {
        KM_data->label_A[a] = KM_data->utility_matrix[a][0];
        for (i64 b = 1; b < KM_data->set_size; ++b) {
            KM_data->label_A[a] = max_i64(
                KM_data->label_A[a],
                KM_data->utility_matrix[a][b]
            );
        }
    }

    // Initialize all label_A[b] as 0
    fill_i64_1D_array(KM_data->label_B, KM_data->set_size, 0);

    fill_i64_1D_array(KM_data->match_A, KM_data->set_size, -1);
    fill_i64_1D_array(KM_data->match_B, KM_data->set_size, -1);

    for (i64 a = 0; a < KM_data->set_size; ++a) {
        if (KM_data->match_A[a] != -1) continue;
        
        fill_i64_1D_array(KM_data->visited_A, KM_data->set_size, false);
        fill_i64_1D_array(KM_data->visited_B, KM_data->set_size, false);

        fill_i64_1D_array(KM_data->slack_delta_B, KM_data->set_size, INF);

        KM_data->queue = (i64_queue*)malloc(1 * sizeof(i64_queue));
        init_i64_queue(KM_data->queue);

        push_i64_queue(KM_data->queue, a);
        KM_data->even_vertex_parent_A[a] = -1;
        KM_data->visited_A[a] = true;
        KM_relax_slack_delta_B(KM_data, a);

        while (true) {
            if (KM_traverse_alternating_tree_with_unvisited_equality_edge(
                    KM_data
                )) break;
            
            KM_relabel(KM_data);

            if (KM_traverse_alternating_tree_with_new_equality_edge(
                    KM_data
                )) break;
        }

        free_i64_queue_in_Kuhn_Munkres_data(KM_data);
    }
}

KM_data_t* binary_search_heuristic_penalty_scale_with_Kuhn_Munkres(
    const general_params_t* general_params,
    const QAN_t* QAN,
    const UD_t* UDs,
    const RIS_t* RISs
) {
    KM_data_t* KM_data = (KM_data_t*)malloc(1 * sizeof(KM_data_t));
    init_KM_data_set_size(general_params, KM_data);
    KM_data->least_initial_rate_matrix =
        gen_KM_data_least_initial_rate_matrix_and_set_penalty_scale_upper_bound(
            general_params,
            QAN,
            UDs,
            RISs,
            KM_data
        );

    bool firstTime = true;
    i64 low_penalty_scale = 0;
    i64 high_penalty_scale = KM_data->penalty_scale_upper_bound;
    while (low_penalty_scale <= high_penalty_scale) {
        i64 middle_penalty_scale =
            low_penalty_scale + (high_penalty_scale - low_penalty_scale) / 2;

        if (!firstTime) {
            free_Kuhn_Munkres_data_without_least_initial_rate_matrix_and_queue(
                KM_data
            );
        }
        firstTime = false;

        reset_KM_data_without_least_initial_rate_matrix_and_numeric_members(
            general_params,
            QAN,
            UDs,
            RISs,
            KM_data,
            middle_penalty_scale
        );

        Kuhn_Munkres(KM_data);

        i64 cur_total_consumption_rate =
            calc_KM_data_result_total_consumption_rate(
                general_params,
                KM_data
            );
        
        if (cur_total_consumption_rate
            > general_params->entangled_pairs_generation_rate) {
            low_penalty_scale = middle_penalty_scale + 1;
        } else {
            high_penalty_scale = middle_penalty_scale - 1;
        }
    }

    KM_data->final_penalty_scale = low_penalty_scale;

    // To ensure that KM_data is successfully updated by the KM Algorithm
    if (KM_data->final_penalty_scale
        != KM_data->penalty_scale_upper_bound + 1) {
        if (!firstTime) {
            free_Kuhn_Munkres_data_without_least_initial_rate_matrix_and_queue(
                KM_data
            );
        }
        
        reset_KM_data_without_least_initial_rate_matrix_and_numeric_members(
            general_params,
            QAN,
            UDs,
            RISs,
            KM_data,
            KM_data->final_penalty_scale
        );

        Kuhn_Munkres(KM_data);
    }

    return KM_data;
}

// ==================================================
// Input
// ==================================================

general_params_t* input_general_params() {
    general_params_t* general_params =
        (general_params_t*)malloc(sizeof(general_params_t));
    
    scanf(
        "%" SCNd64 " %" SCNd64 " %lf" " %lf" " %" SCNd64,
        &general_params->UD_num,
        &general_params->RIS_num,
        &general_params->alpha,
        &general_params->beta,
        &general_params->entangled_pairs_generation_rate
    );

    return general_params;
}

QAN_t* input_QAN() {
    QAN_t* QAN = (QAN_t*)malloc(sizeof(QAN_t));

    scanf(
        "%" SCNd64 " %" SCNd64 " %" SCNd64,
        &QAN->x,
        &QAN->y,
        &QAN->covered_UD_num
    );

    QAN->covered_UD_ID_array = malloc_i64_1D_array(QAN->covered_UD_num);

    for (i64 i = 0; i < QAN->covered_UD_num; ++i) {
        scanf("%" SCNd64, &QAN->covered_UD_ID_array[i]);
    }

    return QAN;
}

UD_t* input_UDs(const general_params_t* general_params) {
    UD_t* UDs = (UD_t*)malloc(general_params->UD_num * sizeof(UD_t));

    for (i64 UD_ID = 0; UD_ID < general_params->UD_num; ++UD_ID) {
        scanf(
            "%" SCNd64 " %" SCNd64 " %" SCNd64 " %" SCNd64 " %" SCNd64 " %lf",
            &UDs[UD_ID].ID,
            &UDs[UD_ID].x,
            &UDs[UD_ID].y,
            &UDs[UD_ID].profit,
            &UDs[UD_ID].expected_rate,
            &UDs[UD_ID].fidelity_threshold
        );
    }

    return UDs;
}

RIS_t* input_RISs(const general_params_t* general_params, const QAN_t* QAN) {
    RIS_t* RISs = (RIS_t*)malloc(
        (general_params->RIS_num + QAN->covered_UD_num) * sizeof(RIS_t)
    );

    for (i64 RIS_ID = 0; RIS_ID < general_params->RIS_num; ++RIS_ID) {
        scanf(
            "%" SCNd64 " %" SCNd64 " %" SCNd64 " %" SCNd64,
            &RISs[RIS_ID].ID,
            &RISs[RIS_ID].x,
            &RISs[RIS_ID].y,
            &RISs[RIS_ID].covered_UD_num
        );

        RISs[RIS_ID].covered_UD_ID_array =
            malloc_i64_1D_array(RISs[RIS_ID].covered_UD_num);

        for (i64 i = 0; i < RISs[RIS_ID].covered_UD_num; ++i) {
            scanf("%" SCNd64, &RISs[RIS_ID].covered_UD_ID_array[i]);
        }
    }

    return RISs;
}

// ==================================================
// Initialize Input Data
// ==================================================

void init_general_params_effective_RIS_num(
    general_params_t* general_params,
    const QAN_t* QAN
) {
    general_params->effective_RIS_num =
        general_params->RIS_num + QAN->covered_UD_num;
}

void init_RISs_dummy_vertices(
    const general_params_t* general_params,
    const QAN_t* QAN,
    const UD_t* UDs,
    RIS_t* RISs
) {
    for (i64 RIS_dummy_vertex_ID = general_params->RIS_num;
         RIS_dummy_vertex_ID < general_params->effective_RIS_num;
         ++RIS_dummy_vertex_ID) {
        RISs[RIS_dummy_vertex_ID].ID = RIS_dummy_vertex_ID;
        RISs[RIS_dummy_vertex_ID].x = QAN->x;
        RISs[RIS_dummy_vertex_ID].y = QAN->y;
        RISs[RIS_dummy_vertex_ID].covered_UD_num = 1;
        RISs[RIS_dummy_vertex_ID].covered_UD_ID_array = malloc_i64_1D_array(1);
        RISs[RIS_dummy_vertex_ID].covered_UD_ID_array[0] =
            QAN->covered_UD_ID_array[
                RIS_dummy_vertex_ID - general_params->RIS_num
            ];
    }
}

// ==================================================
// Output
// ==================================================

void output_heuristic_solution_of_KM_data(
    const general_params_t* general_params,
    const UD_t* UDs,
    const KM_data_t* KM_data
) {
    if (KM_data->final_penalty_scale
        == KM_data->penalty_scale_upper_bound + 1) {
        printf("0\n");
        return;
    }

    i64 selected_UD_num = 0;
    for (i64 effective_RIS_ID = 0;
         effective_RIS_ID < general_params->effective_RIS_num;
         ++effective_RIS_ID) {
        if (KM_data->match_A[effective_RIS_ID] == -1) continue;

        i64 cur_selected_UD_ID = KM_data->match_A[effective_RIS_ID];

        if (KM_data->utility_matrix
                [effective_RIS_ID]
                [cur_selected_UD_ID] == 0) continue;  // nonexistent edge

        selected_UD_num += 1;
    }
    printf("%" PRId64 "\n", selected_UD_num);

    for (i64 effective_RIS_ID = 0;
         effective_RIS_ID < general_params->effective_RIS_num;
         ++effective_RIS_ID) {
        if (KM_data->match_A[effective_RIS_ID] == -1) continue;

        i64 cur_selected_UD_ID = KM_data->match_A[effective_RIS_ID];
        
        if (KM_data->utility_matrix
                [effective_RIS_ID]
                [cur_selected_UD_ID] == 0) continue;  // nonexistent edge

        printf("%" PRId64 " ", cur_selected_UD_ID);

        if (effective_RIS_ID < general_params->RIS_num) {
            printf("%" PRId64 "\n", effective_RIS_ID);
        } else {
            printf("-1\n");
        }
    }
}

// ==================================================
// Memory Release
// ==================================================

void free_general_params(general_params_t* general_params) {
    free(general_params);
}

void free_QAN(QAN_t* QAN) {
    free(QAN->covered_UD_ID_array);
    free(QAN);
}

void free_UDs(UD_t* UDs) {
    free(UDs);
}

void free_RISs(RIS_t* RISs, const general_params_t* general_params) {
    for (int effective_RIS_ID = 0;
         effective_RIS_ID < general_params->effective_RIS_num;
         ++effective_RIS_ID) {
        free(RISs[effective_RIS_ID].covered_UD_ID_array);
    }
    free(RISs);
}

void free_i64_2D_array(i64** i64_2D_array, i64 row_num) {
    for (i64 i = 0; i < row_num; ++i) {
        free(i64_2D_array[i]);
    }
    free(i64_2D_array);
}

void free_i64_queue(i64_queue* queue) {
    free(queue->queue);
    free(queue);
}

void free_least_initial_rate_matrix_in_Kuhn_Munkres_data(
    KM_data_t* KM_data,
    i64 KM_data_size
) {
    free_i64_2D_array(KM_data->least_initial_rate_matrix, KM_data_size);
}

void free_i64_queue_in_Kuhn_Munkres_data(KM_data_t* KM_data) {
    free_i64_queue(KM_data->queue);
}

void free_Kuhn_Munkres_data_without_least_initial_rate_matrix_and_queue(
    KM_data_t* KM_data
) {
    free_i64_2D_array(KM_data->utility_matrix, KM_data->set_size);
    free(KM_data->match_A);
    free(KM_data->match_B);
    free(KM_data->label_A);
    free(KM_data->label_B);
    free(KM_data->visited_A);
    free(KM_data->visited_B);
    free(KM_data->even_vertex_parent_A);
    free(KM_data->slack_delta_B);
    free(KM_data->slack_parent_B);
}

void free_only_Kuhn_Munkres_data_itself(KM_data_t* KM_data) {
    free(KM_data);
}

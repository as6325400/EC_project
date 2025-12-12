#include <bits/stdc++.h>

using namespace std;

struct Point {
    double x{};
    double y{};
};

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int ud_count = 0;
    int ris_count = 0;
    double qan_ent_rate = 0.0;
    uint64_t seed = chrono::high_resolution_clock::now().time_since_epoch().count();

    if (argc >= 4) {
        ud_count = stoi(argv[1]);
        ris_count = stoi(argv[2]);
        qan_ent_rate = stod(argv[3]);
        if (argc >= 5) seed = stoull(argv[4]);
    } else {
        if (!(cin >> ud_count >> ris_count >> qan_ent_rate)) {
            cerr << "Usage: generator <UDs> <RISs> <QAN_ent_rate> [seed]\n";
            return 1;
        }
        if (!(cin >> seed)) {
            seed = chrono::high_resolution_clock::now().time_since_epoch().count();
        }
    }

    mt19937 rng(seed);
    uniform_real_distribution<double> uni01(0.0, 1.0);

    const double alpha = 0.000400;
    const double beta = 0.001200;

    const double max_coord = (ud_count > 400 || ris_count > 250) ? 1500.0 : 1000.0;

    auto clamp_int = [](int v, int lo, int hi) {
        return max(lo, min(hi, v));
    };

    // UD placement: cluster UDs to mimic prior cases.
    const int cluster_count = max(1, static_cast<int>(sqrt(ud_count / 12.0)));
    vector<Point> cluster_centers(cluster_count);
    uniform_real_distribution<double> coord_dist(0.0, max_coord);
    for (auto& c : cluster_centers) {
        c.x = coord_dist(rng);
        c.y = coord_dist(rng);
    }

    normal_distribution<double> cluster_offset(0.0, max_coord * 0.08);
    vector<int> ud_x(ud_count), ud_y(ud_count), ud_profit(ud_count), ud_exp(ud_count);
    vector<double> ud_fid(ud_count);
    for (int i = 0; i < ud_count; ++i) {
        const int cid = uniform_int_distribution<int>(0, cluster_count - 1)(rng);
        double x = cluster_centers[cid].x + cluster_offset(rng);
        double y = cluster_centers[cid].y + cluster_offset(rng);
        x = clamp(x, 0.0, max_coord);
        y = clamp(y, 0.0, max_coord);
        ud_x[i] = clamp_int(static_cast<int>(round(x)), 0, static_cast<int>(max_coord));
        ud_y[i] = clamp_int(static_cast<int>(round(y)), 0, static_cast<int>(max_coord));

        double profit = normal_distribution<double>(140.0, 60.0)(rng);
        if (uni01(rng) < 0.08) profit += normal_distribution<double>(120.0, 40.0)(rng);
        ud_profit[i] = clamp_int(static_cast<int>(round(profit)), 8, 650);

        double exp_rate = normal_distribution<double>(70.0, 30.0)(rng);
        if (uni01(rng) < 0.06) exp_rate += normal_distribution<double>(80.0, 30.0)(rng);
        ud_exp[i] = clamp_int(static_cast<int>(round(exp_rate)), 3, 311);

        double fid = normal_distribution<double>(0.70, 0.07)(rng);
        if (uni01(rng) < 0.12) fid += normal_distribution<double>(0.08, 0.02)(rng);
        ud_fid[i] = clamp(fid, 0.55, 0.95);
    }

    // QAN near the map center with a bit of noise.
    Point qan_pos;
    normal_distribution<double> center_noise(0.0, max_coord * 0.1);
    qan_pos.x = clamp(max_coord * 0.5 + center_noise(rng), 0.0, max_coord);
    qan_pos.y = clamp(max_coord * 0.5 + center_noise(rng), 0.0, max_coord);

    // Direct QAN coverage: prefer nearby UDs.
    double direct_frac = (ud_count < 60) ? 0.25 : (ud_count < 200) ? 0.12
                                        : (ud_count < 500) ? 0.06 : 0.025;
    if (uni01(rng) < 0.15) direct_frac *= 0.4;
    if (uni01(rng) < 0.05) direct_frac = 0.0;
    int direct_cnt = min(ud_count, static_cast<int>(round(direct_frac * ud_count)));

    vector<int> direct_ids;
    if (direct_cnt > 0) {
        vector<pair<double, int>> dist;
        dist.reserve(ud_count);
        for (int i = 0; i < ud_count; ++i) {
            const double dx = ud_x[i] - qan_pos.x;
            const double dy = ud_y[i] - qan_pos.y;
            dist.emplace_back(dx * dx + dy * dy, i);
        }
        nth_element(dist.begin(), dist.begin() + direct_cnt, dist.end(),
                    [](const auto& a, const auto& b) { return a.first < b.first; });
        direct_ids.reserve(direct_cnt);
        for (int i = 0; i < direct_cnt; ++i) {
            direct_ids.push_back(dist[i].second);
        }
        sort(direct_ids.begin(), direct_ids.end());
    }

    // RIS placement.
    vector<Point> ris_pos(ris_count);
    for (int i = 0; i < ris_count; ++i) {
        const bool use_cluster = uni01(rng) < 0.7;
        if (use_cluster) {
            const int cid = i % cluster_count;
            ris_pos[i].x = clamp(cluster_centers[cid].x + cluster_offset(rng), 0.0, max_coord);
            ris_pos[i].y = clamp(cluster_centers[cid].y + cluster_offset(rng), 0.0, max_coord);
        } else {
            ris_pos[i].x = coord_dist(rng);
            ris_pos[i].y = coord_dist(rng);
        }
    }

    // Coverage planning.
    const bool sparse_mode = (ud_count > 600 && uni01(rng) < 0.22);
    double edges_per_ud = sparse_mode ? uniform_real_distribution<double>(0.6, 1.3)(rng)
                                      : uniform_real_distribution<double>(1.8, 5.8)(rng);
    if (!sparse_mode && ud_count < 120 && uni01(rng) < 0.25) {
        edges_per_ud += uniform_real_distribution<double>(0.0, 2.0)(rng);
    }
    if (!sparse_mode && ud_count > 800 && uni01(rng) < 0.15) {
        edges_per_ud *= 0.7;
    }

    const double total_edges = edges_per_ud * ud_count;
    const int avg_cover = max(1, static_cast<int>(round(total_edges / max(1, ris_count))));
    const int max_cover = max(6, min(22, avg_cover + 8));

    vector<unordered_set<int>> ris_cov(ris_count);
    vector<int> ris_cap(ris_count, max_cover);

    auto nearest_candidates = [&](const Point& p, int need) {
        const int k = min(ud_count, need * 3 + 15);
        vector<pair<double, int>> dist;
        dist.reserve(ud_count);
        for (int i = 0; i < ud_count; ++i) {
            const double dx = ud_x[i] - p.x;
            const double dy = ud_y[i] - p.y;
            dist.emplace_back(dx * dx + dy * dy, i);
        }
        if (k < ud_count) {
            nth_element(dist.begin(), dist.begin() + k, dist.end(),
                        [](const auto& a, const auto& b) { return a.first < b.first; });
            dist.resize(k);
        }
        shuffle(dist.begin(), dist.end(), rng);
        vector<int> res;
        res.reserve(min(need, static_cast<int>(dist.size())));
        for (int i = 0; i < need && i < static_cast<int>(dist.size()); ++i) {
            res.push_back(dist[i].second);
        }
        return res;
    };

    normal_distribution<double> cover_noise(static_cast<double>(avg_cover),
                                            max(1.0, avg_cover * 0.25));
    for (int i = 0; i < ris_count; ++i) {
        int cov_sz = clamp_int(static_cast<int>(round(cover_noise(rng))), 1, max_cover);
        vector<int> pick = nearest_candidates(ris_pos[i], cov_sz);
        ris_cov[i].insert(pick.begin(), pick.end());
        ris_cap[i] = max_cover - static_cast<int>(ris_cov[i].size());
    }

    // Ensure a reasonable portion of UDs are coverable unless sparse mode is intended.
    auto recompute_cover = [&](vector<int>& covered) {
        covered.assign(ud_count, 0);
        for (int id : direct_ids) {
            ++covered[id];
        }
        for (const auto& cov : ris_cov) {
            for (int id : cov) {
                ++covered[id];
            }
        }
    };

    vector<int> covered;
    recompute_cover(covered);
    const double desired_ratio = sparse_mode ? 0.0 : 0.6;
    double covered_ratio = count_if(covered.begin(), covered.end(),
                                    [](int v) { return v > 0; })
                           / static_cast<double>(max(1, ud_count));

    if (covered_ratio + 1e-9 < desired_ratio) {
        vector<int> uncovered;
        uncovered.reserve(ud_count);
        for (int i = 0; i < ud_count; ++i) {
            if (covered[i] == 0) uncovered.push_back(i);
        }
        // Greedy: attach uncovered UDs to their nearest RIS with room.
        for (int uid : uncovered) {
            int best_ris = -1;
            double best_dist = numeric_limits<double>::infinity();
            for (int r = 0; r < ris_count; ++r) {
                if (ris_cap[r] <= 0) continue;
                const double dx = ud_x[uid] - ris_pos[r].x;
                const double dy = ud_y[uid] - ris_pos[r].y;
                const double d2 = dx * dx + dy * dy;
                if (d2 < best_dist) {
                    best_dist = d2;
                    best_ris = r;
                }
            }
            if (best_ris != -1) {
                if (ris_cov[best_ris].insert(uid).second) {
                    --ris_cap[best_ris];
                }
            }
        }
        recompute_cover(covered);
    }

    cout << ud_count << " " << ris_count << " " << fixed << setprecision(6) << alpha << " "
         << beta << " " << setprecision(0) << qan_ent_rate << "\n";
    cout << defaultfloat;
    cout << static_cast<int>(round(qan_pos.x)) << " " << static_cast<int>(round(qan_pos.y))
         << " " << direct_ids.size() << "\n";
    if (!direct_ids.empty()) {
        for (size_t i = 0; i < direct_ids.size(); ++i) {
            if (i) cout << " ";
            cout << direct_ids[i];
        }
        cout << "\n";
    }

    cout << fixed << setprecision(2);
    for (int i = 0; i < ud_count; ++i) {
        cout << i << " " << ud_x[i] << " " << ud_y[i] << " " << ud_profit[i] << " "
             << ud_exp[i] << " " << ud_fid[i] << "\n";
    }

    for (int i = 0; i < ris_count; ++i) {
        cout << i << " " << static_cast<int>(round(ris_pos[i].x)) << " "
             << static_cast<int>(round(ris_pos[i].y)) << " " << ris_cov[i].size() << "\n";
        if (!ris_cov[i].empty()) {
            vector<int> sorted_ids(ris_cov[i].begin(), ris_cov[i].end());
            sort(sorted_ids.begin(), sorted_ids.end());
            size_t idx = 0;
            for (int id : sorted_ids) {
                if (idx++) cout << " ";
                cout << id;
            }
            cout << "\n";
        }
    }

    return 0;
}

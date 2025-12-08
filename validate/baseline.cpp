#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdbool.h>

#define MAX_COV 1000
#define MAX_PUMP_ROUNDS 50  // 安全上限
#define EPS_FID 1e-9
#define PE_MIN 1e-15
#define PUN_MIN 1e-18
#define VAL_SAFE_LIMIT 9e18  // double -> long long 安全界限

typedef struct {
    long long id;
    double x, y;
    double profit;
    double exp_rate;      // d_u
    double fidelity_th;   // F_th^u
    int accepted;
    long long used_ris;   // -1 = direct, >=0 = RIS id
} UD;

typedef struct {
    long long id;
    double x, y;
    int covered_cnt;
    long long *covered;   // 改成動態配置
    int used;
} RIS;

typedef struct {
    double x, y;
    int covered_cnt;
    long long *covered;   // 改成動態配置
} QAN;

// ========== 工具函式區 ==========

// 計算兩點距離
double dist(double x1, double y1, double x2, double y2) {
    return hypot(x1 - x2, y1 - y2);
}

// Pe = exp(-alpha * l)
double calc_Pe(double alpha, double l) {
    if (!isfinite(alpha) || !isfinite(l)) return 0.0;
    double arg = -alpha * l;
    if (arg < -700) arg = -700;
    double v = exp(arg);
    if (!isfinite(v)) return 0.0;
    if (v < PE_MIN) v = 0.0;
    return v;
}

// Fe = 0.5 + 0.5 * exp(-beta * l)
double calc_Fe(double beta, double l) {
    if (!isfinite(beta) || !isfinite(l)) return 0.5;
    double arg = -beta * l;
    if (arg < -700) arg = -700;
    double v = 0.5 + 0.5 * exp(arg);
    if (!isfinite(v)) return 0.5;
    return v;
}

// P_p(F1,F2) = F1F2 + (1-F1)(1-F2)
double P_p(double F1, double F2) {
    if (!isfinite(F1) || !isfinite(F2)) return 0.0;
    double v = F1 * F2 + (1.0 - F1) * (1.0 - F2);
    if (!isfinite(v)) return 0.0;
    return v;
}

// F_p(F1,F2) = (F1F2) / [F1F2 + (1-F1)(1-F2)]
double F_p(double F1, double F2) {
    double a = F1 * F2;
    double b = (1.0 - F1) * (1.0 - F2);
    double denom = a + b;
    if (!isfinite(denom) || fabs(denom) < 1e-18) return 0.5;
    double v = a / denom;
    if (!isfinite(v)) return 0.5;
    return v;
}

// 計算最小 n 與累積成功率
int compute_min_n_and_Pn(double Fe, double Pe, double F_th, int *out_n, double *out_Pn) {
    if (!isfinite(Fe) || !isfinite(Pe) || !isfinite(F_th)) return 0;
    if (Fe + EPS_FID >= F_th) {
        *out_n = 0;
        *out_Pn = 1.0;
        return 1;
    }

    double F_prev = Fe, P_prev = 1.0;
    for (int n = 1; n <= MAX_PUMP_ROUNDS; ++n) {
        double Pp = P_p(F_prev, Fe);
        double Fnew = F_p(F_prev, Fe);
        double Pnew = Pp * P_prev;
        if (!isfinite(Fnew) || !isfinite(Pnew)) return 0;
        if (Fnew + EPS_FID >= F_th) {
            *out_n = n;
            *out_Pn = Pnew;
            return 1;
        }
        F_prev = Fnew;
        P_prev = Pnew;
        if (P_prev < PUN_MIN) break;
    }
    return 0;
}

// 計算所需 s 值
long long calc_s_required(double d_u, double Pe, int n, double P_un) {
    if (!isfinite(d_u) || !isfinite(Pe) || !isfinite(P_un)) return LLONG_MAX;
    if (Pe <= PE_MIN || P_un <= PUN_MIN) return LLONG_MAX;

    double denom = Pe * P_un;
    if (!isfinite(denom) || denom <= 0.0) return LLONG_MAX;

    double val = (d_u / denom) * (double)(n + 1);
    if (!isfinite(val)) return LLONG_MAX;
    if (val < 0.0) return 0;
    if (val > VAL_SAFE_LIMIT) return LLONG_MAX;

    long long ret = (long long)ceil(val - 1e-12);
    if (ret < 0) ret = 0;
    return ret;
}

// ========== 主程式 ==========

int main() {
    int n_ud, n_ris;
    double alpha, beta, Ent_Gen_Rate;
    QAN qan;

    // 讀 header
    if (scanf("%d %d %lf %lf %lf", &n_ud, &n_ris, &alpha, &beta, &Ent_Gen_Rate) != 5) {
        fprintf(stderr, "Header read error\n");
        return 0;
    }

    // 動態配置記憶體
    UD *uds = (UD *)calloc(n_ud, sizeof(UD));
    RIS *riss = (RIS *)calloc(n_ris, sizeof(RIS));
    if (!uds || !riss) {
        fprintf(stderr, "Memory allocation failed\n");
        return 0;
    }

    // QAN
    if (scanf("%lf %lf %d", &qan.x, &qan.y, &qan.covered_cnt) != 3) {
        fprintf(stderr, "QAN read error\n");
        return 0;
    }
    if (qan.covered_cnt < 0) qan.covered_cnt = 0;
    qan.covered = (long long *)calloc(qan.covered_cnt, sizeof(long long));
    for (int i = 0; i < qan.covered_cnt; ++i)
        scanf("%lld", &qan.covered[i]);

    // 讀 UDs
    for (int i = 0; i < n_ud; ++i) {
        long long id; double x, y, pf, er, ft;
        if (scanf("%lld %lf %lf %lf %lf %lf", &id, &x, &y, &pf, &er, &ft) != 6) {
            fprintf(stderr, "UD input error at %d\n", i);
            return 0;
        }
        uds[i].id = id;
        uds[i].x = x; uds[i].y = y;
        uds[i].profit = pf;
        uds[i].exp_rate = er;
        uds[i].fidelity_th = ft;
        uds[i].accepted = 0;
        uds[i].used_ris = -1;
    }

    // 讀 RISs
    for (int i = 0; i < n_ris; ++i) {
        long long id; double x, y; int cov;
        if (scanf("%lld %lf %lf %d", &id, &x, &y, &cov) != 4) {
            fprintf(stderr, "RIS input error at %d\n", i);
            return 0;
        }
        riss[i].id = id;
        riss[i].x = x; riss[i].y = y;
        riss[i].covered_cnt = cov;
        riss[i].used = 0;
        riss[i].covered = (long long *)calloc(cov, sizeof(long long));
        for (int j = 0; j < cov; ++j)
            scanf("%lld", &riss[i].covered[j]);
    }

    double residual_rate = Ent_Gen_Rate;
    if (!isfinite(residual_rate) || residual_rate < 0.0) residual_rate = 0.0;

    // 主迴圈
    for (int ui = 0; ui < n_ud; ++ui) {
        long long best_s = LLONG_MAX;
        long long best_ris_id = -2;
        int found_path = 0;

        // 檢查 direct path
        int in_qan_cov = 0;
        for (int k = 0; k < qan.covered_cnt; ++k)
            if (qan.covered[k] == uds[ui].id) in_qan_cov = 1;
        if (in_qan_cov) {
            double l = dist(qan.x, qan.y, uds[ui].x, uds[ui].y);
            double Pe = calc_Pe(alpha, l);
            double Fe = calc_Fe(beta, l);
            int n_min; double P_un;
            if (compute_min_n_and_Pn(Fe, Pe, uds[ui].fidelity_th, &n_min, &P_un)) {
                long long s_req = calc_s_required(uds[ui].exp_rate, Pe, n_min, P_un);
                if (s_req != LLONG_MAX && s_req <= (long long)(residual_rate + 1e-9)) {
                    best_s = s_req;
                    best_ris_id = -1;
                    found_path = 1;
                }
            }
        }

        // 檢查 RIS path
        for (int ri = 0; ri < n_ris; ++ri) {
            if (riss[ri].used) continue;
            int covers = 0;
            for (int c = 0; c < riss[ri].covered_cnt; ++c)
                if (riss[ri].covered[c] == uds[ui].id) covers = 1;
            if (!covers) continue;

            double l = dist(qan.x, qan.y, riss[ri].x, riss[ri].y)
                     + dist(riss[ri].x, riss[ri].y, uds[ui].x, uds[ui].y);
            double Pe = calc_Pe(alpha, l);
            double Fe = calc_Fe(beta, l);
            int n_min; double P_un;
            if (compute_min_n_and_Pn(Fe, Pe, uds[ui].fidelity_th, &n_min, &P_un)) {
                long long s_req = calc_s_required(uds[ui].exp_rate, Pe, n_min, P_un);
                if (s_req != LLONG_MAX && s_req <= (long long)(residual_rate + 1e-9)) {
                    if (!found_path || s_req < best_s ||
                        (s_req == best_s && (best_ris_id == -1 || riss[ri].id < best_ris_id))) {
                        best_s = s_req;
                        best_ris_id = riss[ri].id;
                        found_path = 1;
                    }
                }
            }
        }

        if (found_path && best_s != LLONG_MAX && best_s <= (long long)(residual_rate + 1e-9)) {
            uds[ui].accepted = 1;
            uds[ui].used_ris = best_ris_id;
            residual_rate -= (double)best_s;
            if (residual_rate < 0.0) residual_rate = 0.0;
            if (best_ris_id >= 0) {
                for (int ri = 0; ri < n_ris; ++ri)
                    if (riss[ri].id == best_ris_id) riss[ri].used = 1;
            }
        }
    }

    // 輸出結果
    int acc_cnt = 0;
    for (int i = 0; i < n_ud; ++i)
        if (uds[i].accepted) acc_cnt++;
    printf("%d\n", acc_cnt);
    for (int i = 0; i < n_ud; ++i)
        if (uds[i].accepted)
            printf("%lld %lld\n", uds[i].id, uds[i].used_ris);

    // 釋放記憶體
    free(qan.covered);
    for (int i = 0; i < n_ris; ++i) free(riss[i].covered);
    free(riss);
    free(uds);
    return 0;
}

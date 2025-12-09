#pragma once

#include <istream>
#include <set>
#include <vector>

class site {
  public:
    double x{}, y{};
    int id{};
  protected:
    site() = default;
};

class cov {
  public:
    int coverage_nums{};
    std::set<int> coverage_UDs;
  protected:
    cov() = default;
};

class UD : public site {
  public:
    bool can_coverage{};
    int profit{};
    double exp_rate{};
    double fidelity_th{};
};

class RIS : public site, public cov {
};

class QAN : public site, public cov {
  public:
    double ent_gen_rate{};
};

struct ProblemData {
    std::vector<UD> uds;
    std::vector<RIS> riss;
    QAN qan;
    double alpha{};
    double beta{};
    int actual_ris_count{};
    int direct_ris_index{};
    std::vector<std::vector<double>> cost_table;
};

ProblemData load_problem(std::istream& in);

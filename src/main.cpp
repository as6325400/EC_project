#include <bits/stdc++.h>

using namespace std;

class site{
  public:
    double x, y;
    int id;
  protected:
    site() = default;
};

class cov
{
  public:
    int coverage_nums;
    set<int> coverage_UDs;
  protected:
    cov() = default;
};

class UD : public site
{ 
  public:
    bool can_coverage;
    int profit;
    double exp_rate, fidelity_th;
};

class RIS : public site, public cov
{ 

};

class QAN : public site, public cov
{ 
  public:
    double ent_gen_rate;
};

vector<UD> UD_Array;
vector<RIS> RIS_Array;
double alpha, beta;
QAN qan;

double dis(const site &a, const site &b){
  return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

double Pe(const UD &ud, const RIS *ris){
  double d;
  if(ris == nullptr) d = dis(ud, qan);
  else d = dis(qan, *ris) + dis(*ris, ud);
  return exp(-1 * alpha * d);
}

double Fe(const UD &ud, const RIS *ris){
  double d;
  if(ris == nullptr) d = dis(ud, qan);
  else d = dis(qan, *ris) + dis(*ris, ud);
  return 0.5 + 0.5 * exp(-1 * ::beta * d);
}

double Pp(const double &f1, const double &f2){
  return f1 * f2 + (1 - f1) * (1 - f2);
}

double Fp(const double &f1, const double &f2){
  return (f1 * f2) / (f1 * f2 + (1 - f1) * (1 - f2));
}

double F(const UD &ud, const RIS *ris, int n){
  if(n == 0) return Fe(ud, ris);
  return Fp(F(ud, ris, n - 1), Fe(ud, ris));
}

double P(const UD &ud, const RIS *ris, int n){
  if(n == 0) return 1;
  return Pp(F(ud, ris, n - 1), Fe(ud, ris)) * P(ud, ris, n - 1);
}

int N(const UD &ud, const RIS *ris){
  int ans = 0;
  for(int i = 0; ;i++){
    if(F(ud, ris, i) >= ud.fidelity_th){
      ans = i;
      break;
    }
  }
  return ans;
}

double S(const UD &ud, const RIS *ris){
  int rounds = N(ud, ris);
  return (ud.exp_rate / Pe(ud, ris)) * ((rounds + 1) / P(ud, ris, rounds)); 
}

int main(){
  
  int UDs, RISs;

  cin >> UDs >> RISs >> ::alpha >> ::beta >> qan.ent_gen_rate;
  cin >> qan.x >>qan.y >> qan.coverage_nums;
  
  for(int i = 0; i < qan.coverage_nums; i++){
    int id;
    cin >> id;
    qan.coverage_UDs.insert(id);
  }

  UD_Array.resize(UDs);
  RIS_Array.resize(RISs);

  for(auto &ud:UD_Array){
    cin >> ud.id >> ud.x >> ud.y >> ud.profit >> ud.exp_rate >> ud.fidelity_th;
    if(qan.coverage_UDs.count(ud.id)) ud.can_coverage = true;
    else ud.can_coverage = false;
  }

  for(auto &ris:RIS_Array){
    cin >> ris.id >> ris.x >> ris.y >> ris.coverage_nums;
    for(int i = 0; i < ris.coverage_nums; i++){
      int id;
      cin >> id;
      ris.coverage_UDs.insert(id);
    }
  }

  int accept_UD = 0;
  vector<pair<int, int>> accept_UD_list;
  vector<bool> used(RISs, false);
  sort(UD_Array.begin(), UD_Array.end(),
     [](const UD &a, const UD &b){
        return a.profit > b.profit;
     });
  for(int i = 0; i < UDs; i++){
    auto &ud = UD_Array[i];
    double cost = 1e9;
    int min_cost_id = -1;
    if(ud.can_coverage) cost = S(ud, nullptr);
    for(int j = 0; j < RISs; j++){
      auto &ris = RIS_Array[j];
      if(j >= 0 and ris.coverage_UDs.count(ud.id) == 0) continue;
      if(j >= 0 and used[ris.id]) continue;
      double candidate_cost = S(ud, &ris);
      if(candidate_cost < cost){
        cost = candidate_cost;
        min_cost_id = ris.id;
      }
    }
    if(qan.ent_gen_rate >= ceil(cost)){
      accept_UD++;
      qan.ent_gen_rate -= ceil(cost);
      accept_UD_list.push_back({ud.id, min_cost_id});
      if(min_cost_id != -1) used[min_cost_id] = true;
    }
  }
  cout << accept_UD << '\n';
  for(auto [idx, rid]:accept_UD_list){
    cout << idx << ' ' << rid << '\n';
  }

  return 0;
}

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

double dis(site &a, site &b){
  return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

double Pe(int UID, int RID){
  double d;
  if(RID == -1) d = dis(UD_Array[UID], qan);
  else d = dis(qan, RIS_Array[RID]) + dis(RIS_Array[RID], UD_Array[UID]);
  return exp(-1 * alpha * d);
}

double Fe(int UID, int RID){
  double d;
  if(RID == -1) d = dis(UD_Array[UID], qan);
  else d = dis(qan, RIS_Array[RID]) + dis(RIS_Array[RID], UD_Array[UID]);
  return 0.5 + 0.5 * exp(-1 * ::beta * d);
}

double Pp(double f1, double f2){
  return f1 * f2 + (1 - f1) * (1 - f2);
}

double Fp(double f1, double f2){
  return (f1 * f2) / (f1 * f2 + (1 - f1) * (1 - f2));
}

double F(int UID, int RID, int n){
  if(n == 0) return Fe(UID, RID);
  return Fp(F(UID, RID, n - 1), Fe(UID, RID));
}

double P(int UID, int RID, int n){
  if(n == 0) return 1;
  return Pp(F(UID, RID, n - 1), Fe(UID, RID)) * P(UID, RID, n - 1);
}

int N(int UID, int RID){
  int ans = 0;
  for(int i = 0; ;i++){
    if(F(UID, RID, i) >= UD_Array[UID].fidelity_th){
      ans = i;
      break;
    }
  }
  return ans;
}

double S(int UID, int RID){
  return (UD_Array[UID].exp_rate / Pe(UID, RID)) * ((N(UID, RID) + 1) / P(UID, RID, N(UID, RID))); 
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
  for(int i = 0; i < UDs; i++){
    auto ud = UD_Array[i];
    double cost = 1e18;
    int min_cost_id = -1;
    for(int j = -1; j < RISs; j++){
      if(j == -1 and !ud.can_coverage) continue;
      if(j >= 0 and RIS_Array[j].coverage_UDs.count(i) == 0) continue;
      if(j >= 0 and used[j]) continue;
      if(S(i, j) < cost){
        cost = S(i, j);
        min_cost_id = j;
      }
    }
    if(qan.ent_gen_rate >= ceil(cost)){
      accept_UD++;
      qan.ent_gen_rate -= ceil(cost);
      accept_UD_list.push_back({i, min_cost_id});
      if(min_cost_id != -1) used[min_cost_id] = true;
    }
  }

  cout << accept_UD << '\n';
  for(auto [idx, rid]:accept_UD_list){
    cout << idx << ' ' << rid << '\n';
  }

  return 0;
}
#include <bits/stdc++.h>

using namespace std;

struct UD
{
  int x, y, id;
  bool can_coverage;
  int profit;
  double exp_rate, fth;
};

struct RIS
{
  int x, int y;
  set<int> coverage_UDs;
}

double fidelity(int dis){
  return 0.5 + 0.5 * exp(-dis * 12e-4);
}

double Probability(int dis){
  return exp(-dis * 4e-4);
}

int main(){
  int UDs, RISs, X, Y, coverage_nums;
  double alpha, beta, genrate;
  cin >> UDs >> RISs >> alpha >> beta >> genrate;
  cin >> X >> Y >> coverage_nums;
  set<int> st;
  for(int i = 0; i < coverage_nums; i++){
    int id;
    cin >> id;
    st.insert(id);
  }
  vector<UD> UD_Array(UDs);
  for(auto &ud:UD_Array){
    cin >> ud.id >> ud.x >> ud.y >> ud.profit >> ud.exp_rate >> ud.fth;
    if(st.count(ud.id)) ud.can_coverage = true;
    else ud.can_coverage = false;
  }
  vector<RIS> RIS_Array(RISs);
  for(auto &ris:RIS_Array){
    int num;
    cin >> ris.id >> ris.x >> ris.y >> num;
    for(int i = 0; i < num; i++){
      int id;
      cin >> id;
      ris.coverage_UDs.insert(id);
    }
  }
  return 0;
}
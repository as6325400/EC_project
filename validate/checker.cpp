#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
using namespace std;
#define int long long 
#define x first 
#define y second
void wa() { exit(43); }
void accept() { exit(42); }
const double EPS=1e-12;
class device{
    public:
    int id,coverage;
    pair<int,int> pos;
    int qubit_num;
    vector<int> coverUD;
};
class UD{
    public:
    int id,Exp_Rate,profit;
    double F_th;
    pair<int,int> pos;
};
int numTestcase = 20;
int UDs, RISs,Ent_Gen_rate;
int total_profit,BaseLine_profit,cur_profit;
//string msg;// baseline data
double alpha, beta_p;
void readInput(ifstream &inFile);
void readAns(ifstream &ansFile);
void readUser();
bool isInt(string &str);
int s2i(string str);
vector<int> ss2vector(string &s);
vector<string> ss2string(string &s);
vector<device> RIS;
vector<UD> UD_list;
device QAN;
void readInput(ifstream &inFile) {
  inFile>>UDs>>RISs>>alpha>>beta_p>>QAN.qubit_num;
  QAN.id=0;
  inFile>>QAN.pos.x>>QAN.pos.y>>QAN.coverage;
  for(int i=0;i<QAN.coverage;i++){
      int ud_id;
      inFile>>ud_id;
      QAN.coverUD.push_back(ud_id);
  }
  for(int i=0;i<UDs;i++){
    UD tmp;
    inFile>>tmp.id>>tmp.pos.x>>tmp.pos.y>>tmp.profit>>tmp.Exp_Rate>>tmp.F_th;
    UD_list.push_back(tmp);
  }
  for(int i=0;i<RISs;i++){
    device tmp;
    tmp.qubit_num=1;
    inFile>>tmp.id>>tmp.pos.x>>tmp.pos.y>>tmp.coverage;
    for(int j=0;j<tmp.coverage;j++){
        int ud_id;
        inFile>>ud_id;
        tmp.coverUD.push_back(ud_id);
    }
    RIS.push_back(tmp);
  }
  return;
}
double link_success_prob(double dis) {
  return exp(-alpha * dis);
}

double link_fidelity(double dis) {
  return 0.5 + (0.5 * exp(-beta_p * (dis)));
}
double purification_fidelity(double f1, double f2) {
  return (f1 * f2) / ((f1 * f2) + (1 - f1) * (1 - f2));
}

double purification_success_prob(double f1, double f2) {
  return (f1 * f2) + ((1 - f1) * (1 - f2));
}
double purify_fidelity(double dis,int round){
    if(round==0) return link_fidelity(dis);
    double link_f1=purify_fidelity(dis,round-1);
    double link_f2=link_fidelity(dis);
    double purify_fid=purification_fidelity(link_f1,link_f2);
    return purify_fid;
}
double purify_prob(double dis,int round){
    if(round==0) return 1.0;
    double link_f1=purify_fidelity(dis,round-1);
    double link_f2=link_fidelity(dis);
    double purify_p=purification_success_prob(link_f1,link_f2)*purify_prob(dis,round-1);
    return purify_p;
}
int purify_num(double dis,double F_th){  //d(u) dudududu max verstappen
    for(int i=0;;i++){
      double cur_fid=purify_fidelity(dis,i);
      if(cur_fid+EPS>=1.0) return i;
      if(cur_fid>=F_th) return i;
    }
}
int satisfy(double dis,double num,double F_th){
    int round=purify_num(dis,F_th);
    double sat=num/link_success_prob(dis)*(round+1)/purify_prob(dis,round);
    return ceil(sat);
}
double dist(pair<int,int> a,pair<int,int> b){
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}
double RIS_pass_dist(int UDid,int RISid){
    if(RISid==-1) return dist(QAN.pos,UD_list[UDid].pos);
    return dist(QAN.pos,RIS[RISid].pos)+dist(RIS[RISid].pos,UD_list[UDid].pos);
}
void readAns(ifstream &ansFile) {
  // check this is first case or last case;
  //ansFile>>msg;
  ansFile>>cur_profit;
  BaseLine_profit+=cur_profit;
  return;
}
void readUser(){
  try{
    string stringLine;
    set<int> usedUD;
    set<int> usedRIS;
    vector<pair<int,int>> user_output; //UDid,RISid
    int res_qubit=QAN.qubit_num;
    int flag=0;
    while(getline(cin,stringLine)){
        vector<int> in_data=ss2vector(stringLine);
        if(flag==0){
            if((int)in_data.size()!=1){
                cerr<<"user output format of AccUDs error!"<<endl;
                wa();
            }
            flag=1;
            int acc_num=in_data[0];
            if(acc_num<0||acc_num>UDs){
                cerr<<"user output AccUDs out of range!"<<endl;
                wa(); 
            }
            if(acc_num==0){
                cerr<<"each testcase must accept at least one UD!"<<endl;
                wa();
            }
        }
        else{
            if((int)in_data.size()!=2){
                cerr<<"user output of <Acc_UD_ID,Used_RIS_ID_> error!"<<endl;
                wa();
            }
            int UDid=in_data[0],RISid=in_data[1];
            if(UDid<0||UDid>=UDs){
                cerr<<"user output UDid out of range!"<<endl;
                wa(); 
            }
            if(RISid<-1||RISid>=RISs){
                cerr<<"user output RISid out of range!"<<endl;
                wa();
            }
            if(UDid>=0&&usedUD.count(UDid)){
                cerr<<"user output repeat UDid!"<<endl;
                wa();
            }
            if(RISid>=0&&usedRIS.count(RISid)){
                cerr<<"user use same RIS again!"<<endl;
                wa(); 
            }
            usedUD.insert(UDid);
            if(RISid!=-1) usedRIS.insert(RISid);
            if(RISid!=-1&&find(RIS[RISid].coverUD.begin(),RIS[RISid].coverUD.end(),UDid)==RIS[RISid].coverUD.end()){
                cerr<<"user output RISid can't cover UDid!"<<endl;
                wa();
            }
            if(RISid==-1&&find(QAN.coverUD.begin(),QAN.coverUD.end(),UDid)==QAN.coverUD.end()){
                cerr<<"user output QAN can't cover UDid!"<<endl;
                wa();
            }
            user_output.push_back({UDid,RISid});
            int need_qubit=satisfy(RIS_pass_dist(UDid,RISid),UD_list[UDid].Exp_Rate,UD_list[UDid].F_th);
            //cerr<<"UDid: "<<UDid<<" RISid: "<<RISid<<" need_qubit: "<<need_qubit<<endl;
            if(need_qubit>res_qubit){
                cerr<<"user output qubit_num exceed!"<<endl;
                wa();
            }
            
            res_qubit-=need_qubit;
            total_profit+=UD_list[UDid].profit;
        }
    }
  }catch (const runtime_error &e) {
    cerr << "run error: " << e.what() << endl;
    exit(1);
  } catch (const exception &e) {
    cerr << "sys error: " << e.what() << endl;
    exit(1);
  }
}
#undef int
int main(int argc, char *argv[]){
  ifstream in_file, ans_file;
  in_file.open(argv[1]);
  ans_file.open(argv[2]);

  readInput(in_file);
  readAns(ans_file);
  readUser();
  //cerr<<fixed<<setprecision(6)<<"the total profit: "<<total_profit<<endl;
  /* if(total_profit<cur_profit){
    cerr<<fixed<<setprecision(6)<<"user total profit: "<<total_profit<<endl;
    cerr<<fixed<<setprecision(6)<<"baseline total profit: "<<cur_profit<<endl;
    cerr<<"Profit is lower than baseline"<<endl;
    wa();
  } */
  if(total_profit<=0){
    cerr<<"the student score is zero!\n";
    wa();
  }
  cerr<<"OPT_SCORE="<<total_profit<<endl;
  accept();
  return 0;
}
#define int long long

bool isInt(string &str) {
    if (str.empty()) return false;  // 空字串不是整數
    if (str.size() > 10) return false; // 可再依需求調整上限

    int start = 0;
    if (str[0] == '-') {
        if (str.size() == 1) return false; // 只有 "-" 不是數字
        start = 1; // 從下一個字元開始檢查
    }

    for (int i = start; i < (int)str.size(); i++) {
        if (!isdigit(str[i]))
            return false;
    }
    return true;
}


int s2i(string str) {
  if (!isInt(str)) {
    cerr << "node's id or RIS's id isn't number " << str << endl;
    wa();
  }
  return stoi(str);
}

int s2i_AccUDs(string str) {
  if (!isInt(str)) {
    cerr << "AccUD isn't number " << str << endl;
    wa();
  }
  return stoi(str);
}

vector<int> ss2vector(string &s) {
  vector<int> res;
  stringstream ss;
  ss << s;
  try {
    string tmp;
    while (ss >> tmp) {
      res.push_back(s2i(tmp));
    }
    return res;
  } catch (const runtime_error &e) {
    cerr << "(str2vector) run error: " << e.what() << endl;
  } catch (const exception &e) {
    cerr << "(str2vector) sys error: " << e.what() << endl;
  }
  return res;
}

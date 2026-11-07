#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
struct Witness {vector<int> rows; int root=-1;};
map<int,Witness> cycles,endpoint_values; map<int,long long> hist;
int dp[512][9],rootdp[512][9];long long seen=0;
void visit(const vector<int>&a){
 int n=a.size(),N=1<<n;
 for(auto &row:dp) fill(begin(row),end(row),0);
 for(auto &row:rootdp) fill(begin(row),end(row),0);
 for(int i=0;i<n;i++)dp[1<<i][i]=1;
 rootdp[1][0]=1;
 for(int S=1;S<N;S++)for(int i=0;i<n;i++)if(dp[S][i]||rootdp[S][i]){
  int todo=a[i]&~S;
  while(todo){int bit=todo&-todo,j=__builtin_ctz((unsigned)bit);todo-=bit;
   dp[S|bit][j]+=dp[S][i]; rootdp[S|bit][j]+=rootdp[S][i];
  }
 }
 int H=0,C=0;
 for(int i=0;i<n;i++){H+=dp[N-1][i];if(n>=3 && (a[i]&1))C+=rootdp[N-1][i];if(!endpoint_values.count(dp[N-1][i]))endpoint_values[dp[N-1][i]]={a,i};}
 if(!(H%2) || n*C>H)throw runtime_error("parity/rotation control failed");
 if(!cycles.count(C)){cycles[C]={a,-1};} hist[C]++;seen++;
}
void write(const string&name,const map<int,Witness>&data){ofstream o(name);o<<"value\troot\trows\n";for(auto &[v,w]:data){o<<v<<'\t'<<w.root<<'\t';for(int i=0;i<(int)w.rows.size();i++)o<<(i?",":"")<<w.rows[i];o<<'\n';}}
int main(int argc,char**argv){
 if(argc!=2){throw runtime_error("usage census DIRECTORY");} string dir=argv[1];
 for(int n=1;n<=9;n++){
  cycles.clear();endpoint_values.clear();hist.clear();seen=0;
  int m=min(n,8);ifstream f(dir+"/classes"+to_string(m)+".rows");if(!f)throw runtime_error("missing class bank");
  vector<int>a(m);while(f>>a[0]){for(int j=1;j<m;j++)if(!(f>>a[j]))throw runtime_error("short row");
   if(n<9)visit(a);else for(int cut=0;cut<256;cut++){vector<int>b=a;for(int i=0;i<8;i++)b[i]|=((cut>>i)&1)<<8;b.push_back(255^cut);visit(b);}
  }
  string base=dir+"/n"+to_string(n);write(base+"_cycles.tsv",cycles);write(base+"_endpoints.tsv",endpoint_values);
  ofstream out(base+"_hist.tsv");out<<"cycles\t"<<(n<=8?"isomorphism_classes":"extension_presentations")<<'\n';for(auto [c,mult]:hist)out<<c<<'\t'<<mult<<'\n';
  cout<<"n="<<n<<" rows="<<seen<<" cycle_values="<<cycles.size()<<" cycle_max="<<cycles.rbegin()->first<<" endpoint_values="<<endpoint_values.size()<<" endpoint_max="<<endpoint_values.rbegin()->first<<"\n"<<flush;
 }
}

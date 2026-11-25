#include <iostream>
#include <vector>
#include <tuple>
using namespace std;
int n;vector<pair<int,int>> E;
int run(int m,bool *sc=nullptr,int a0=-1,int b0=-1){
 int R[8]={};for(int i=0;i<n-1;i++)R[i]|=1<<(i+1);
 for(int e=0;e<(int)E.size();e++){auto[a,b]=E[e];if(m>>e&1)R[b]|=1<<a;else R[a]|=1<<b;}
 if(a0>=0){R[a0]^=1<<b0;R[b0]^=1<<a0;}
 if(sc){int S=1;for(int k=0;k<n;k++)for(int x=0;x<n;x++)if(S>>x&1)S|=R[x];
 *sc=S==(1<<n)-1;S=1;for(int k=0;k<n;k++)for(int x=0;x<n;x++)if(S&R[x])S|=1<<x;
 *sc&=S==(1<<n)-1;}
 int dp[256][8]={}; for(int i=0;i<n;i++)dp[1<<i][i]=1;
 for(int M=1;M<(1<<n);M++)for(int i=0;i<n;i++)if(dp[M][i])
   for(int a=R[i]&~M;a;a&=a-1){int j=__builtin_ctz(a);dp[M|1<<j][j]+=dp[M][i];}
 int H=0;for(int i=0;i<n;i++)H+=dp[(1<<n)-1][i];return H;
}


int main(){
 vector<vector<int>> cores;
 for(n=4;n<=7;n++){
  E.clear();for(int i=0;i<n;i++)for(int j=i+2;j<n;j++)E.emplace_back(i,j);
  int counts[128]={};
  for(int m=0;m<1<<(int)E.size();m++){
   int h=run(m);if(h<128)counts[h]++;
   if(n==7&&h==55){
    for(int a=0;a<n;a++)for(int b=a+1;b<n;b++)if(run(m,nullptr,a,b)==53||run(m,nullptr,a,b)==57)return 1;
    vector<int> cut(n);for(int i=0;i<n-1;i++)cut[i]|=1<<(i+1);
    for(int e=0;e<(int)E.size();e++){auto[a,b]=E[e];if(m>>e&1)cut[b]|=1<<a;else cut[a]|=1<<b;}
    cores.push_back(move(cut));
   }
  }
  cout<<"order "<<n<<" / H13,H33,H55 "<<counts[13]<<" "<<counts[33]<<" "<<counts[55]<<endl;
 }
 if(cores.size()!=220)return 2;
 n=7; if(run(81)!=55||run(593)!=65||run(597)!=53||run(595)!=57)return 5;
 cout<<"Two-edit fork: 55 at mask 81, 65 at 593, then 53 at 597 or 57 at 595."<<endl;

 n=8;E.clear();for(int i=0;i<n;i++)for(int j=i+2;j<n;j++)E.emplace_back(i,j);
 for(const auto &T:cores)for(int side=0;side<2;side++){
  // The old forward Hamilton path stays a Hamilton path after adding a source or a sink.
  int m=0;
  for(int i=0;i<(int)E.size();i++){
    auto[a,b]=E[i];
    if(side){if(b!=7&&T[b]>>a&1)m|=1<<i;}
    else if(a>0&&T[b-1]>>(a-1)&1)m|=1<<i;
  }
  if(run(m)!=55)return 3;
  for(int a=0;a<n;a++)for(int b=a+1;b<n;b++)if(run(m,nullptr,a,b)==53||run(m,nullptr,a,b)==57)return 4;
 }
 int all55=0,strong55=0;
 for(int m=0;m<(1<<(int)E.size());m++){bool sc;if(run(m,&sc)==55){all55++;strong55+=sc;}}
 if(all55!=440||strong55!=0)return 6;
 cout<<"Order 8 independent sweep: 440 spine masks at 55, none strong."<<endl;
 cout<<"No 55→53 or 55→57 at order seven or among singleton pads."<<endl;
}

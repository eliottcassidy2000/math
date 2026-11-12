// Optional bounded discovery. Counts h=15 deletion-core attachments, not all tournaments.
#include <bits/stdc++.h>
using namespace std;using ll=long long;
vector<tuple<ll,ll,int,int>> records;
vector<pair<string,vector<int>>> sources;
vector<int> join(vector<vector<int>> cs){
 vector<int>D;int n=0;for(auto&c:cs)n+=c.size();int off=0;
 for(auto&c:cs){int nn=c.size();for(int a:c)D.push_back((a<<off)|(((1<<n)-1)^((1<<(off+nn))-1)));off+=nn;}return D;
}
vector<array<ll,20>> dp(vector<int>D){
 int n=D.size(),N=1<<n;vector<array<ll,20>>a(N);
 for(int i=0;i<n;i++)a[1<<i][i]=1;
 for(int s=1;s<N;s++)for(int i=0;i<n;i++)if(a[s][i]){
  int v=D[i]&~s;while(v){int j=__builtin_ctz(v);v&=v-1;a[s|(1<<j)][j]+=a[s][i];}
 }return a;
}
void scan(vector<int>D,string label){
 int sid=sources.size();sources.push_back({label,D});int n=D.size(),N=1<<n,F=N-1;
 auto e=dp(D);vector<int>R(n);for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(D[i]>>j&1)R[j]|=1<<i;
 auto st=dp(R);vector<ll>H(N);H[0]=1;
 for(int s=1;s<N;s++)for(int i=0;i<n;i++)H[s]+=e[s][i];ll h=H[F];
 if(h!=15)throw runtime_error("deletion h");
 ll E[20]={},W[20][20]={},B[20][20]={},L[20]={};
 for(int s=1;s<N;s++)for(int i=0;i<n;i++){
  E[i]+=e[s][i]*H[F^s];if(e[s][i])for(int j=0;j<n;j++)W[i][j]+=e[s][i]*st[F^s][j];
 }
 for(int i=0;i<n;i++){L[i]=e[F][i]-st[F][i];for(int j=0;j<n;j++){L[i]+=W[i][j];B[i][j]=W[i][j]+W[j][i];}}
 ll hc=h,v=h;int old=0;
 for(int t=0;t<N;t++){
  int x=t^(t>>1);if(t){int k=__builtin_ctz(x^old);int sign=(x>>k&1)?1:-1;ll delta=L[k];
   for(int j=0;j<n;j++)if(old>>j&1)delta-=B[k][j];hc+=sign*delta;v+=sign*E[k];}
  old=x;records.emplace_back(v,hc,sid,x);
 }
}
int main(int argc,char**argv){
 int total=argc>1?atoi(argv[1]):14;if(total<7||total>18){cerr<<"Choose deletion order7..18\n";return 2;}
 vector<int>A={6,12,24,17,3},B={2,4,1},C={2,4,9,3};
 for(int a=0;a<=total-5;a++){
  vector<vector<int>>cs;for(int j=0;j<a;j++)cs.push_back({0});cs.push_back(A);
  for(int j=a;j<total-5;j++)cs.push_back({0});scan(join(cs),"C5:"+to_string(a));
 }
 for(int sw=0;sw<2;sw++)for(int a=0;a<=total-7;a++)for(int b=0;b<=total-7-a;b++){
  vector<vector<int>>cs;for(int j=0;j<a;j++)cs.push_back({0});cs.push_back(sw?C:B);
  for(int j=0;j<b;j++)cs.push_back({0});cs.push_back(sw?B:C);
  for(int j=a+b;j<total-7;j++)cs.push_back({0});scan(join(cs),"C3C4:"+to_string(sw)+":"+to_string(a)+":"+to_string(b));
 }
 for(auto core:vector<vector<int>>{{36,1,2,7,15,30},{44,1,2,6,15,30}})for(int a=0;a<=total-6;a++){
  vector<vector<int>>cs;for(int j=0;j<a;j++)cs.push_back({0});cs.push_back(core);
  for(int j=a;j<total-6;j++)cs.push_back({0});scan(join(cs),"C6:"+to_string(core[0])+":"+to_string(a));
 }
 sort(records.begin(),records.end());ll lastV=-1,lastH=-10;int run=0,best=0;vector<tuple<ll,ll,int,int>>chain;
 for(auto rec:records){auto[v,h,sid,x]=rec;if(v==lastV&&h==lastH)continue;
  if(v==lastV&&h==lastH+2){run++;chain.push_back(rec);}else{run=1;chain={rec};}lastV=v;lastH=h;
  if(run>best){best=run;cout<<"BEST "<<best<<" V "<<v<<"\n";
   for(auto[vv,hh,ss,xx]:chain){cout<<sources[ss].first<<" H "<<hh<<" X "<<xx<<" D";for(auto a:sources[ss].second)cout<<" "<<a;cout<<"\n";}cout<<flush;}
 }
 cerr<<"skeletons "<<sources.size()<<" cuts "<<records.size()<<" max_run "<<best<<"\n";
}

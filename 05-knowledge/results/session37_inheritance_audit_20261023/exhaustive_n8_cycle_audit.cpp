// Independent n=8 census via rooted cycles and odd-cycle packing identity.
// No row-assignment determinant/permanent DP and no all-starts Hamiltonian DP.
// Edge convention: upper-pair-lex bit 1 means i->j for i<j.
// Universe: all 2^21 presentations containing path 0->1->...->7.
// Every tournament has a Hamiltonian path, so this covers all isomorphism types.
#include <array>
#include <vector>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <fstream>
#include <algorithm>
using namespace std;
constexpr int n=8,N=1<<n,full=N-1;
struct Key { array<int,9> a{}; bool operator==(const Key&o)const{return a==o.a;} };
struct Hash { size_t operator()(const Key&x) const {size_t h=1937;for(int v:x.a) h=h*1000003+unsigned(v);return h;} };
struct Value {uint32_t mask;int H;int D35;};
struct Pair {int s,t,type;};
int main(int argc,char**argv) {
 int limit=argc>1?stoi(argv[1]):1<<21;
 vector<pair<int,int>> freeedges;uint32_t fixedmask=0;int b=0;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,b++) {if(j==i+1)fixedmask|=1u<<b;else freeedges.push_back({i,j});}
 int size[N],root[N],afterroot[N];
 for(int s=1;s<N;s++){size[s]=__builtin_popcount((unsigned)s);root[s]=__builtin_ctz((unsigned)s);afterroot[s]=full^((1<<(root[s]+1))-1);}
 vector<Pair> pairs;
 for(int s=1;s<N;s++)if(size[s]>=3)for(int t=s+1;t<N;t++)if(!(s&t)&&size[t]>=3) {
  int a=size[s],b=size[t];if(a>b)swap(a,b);
  int type=(a==3&&b==3)?0:(a==3&&b==4)?1:(a==3&&b==5)?2:(a==4&&b==4)?3:-1;
  if(type>=0)pairs.push_back({s,t,type});
 }
 unordered_map<Key,Value,Hash> seen;seen.reserve(4000);
 long long presentation_H_sum=0;int maxH=0,minH=100000;
 for(int it=0;it<limit;it++) {
  int out[n]={};uint32_t mask=fixedmask;int freebit=0,edgebit=0;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,edgebit++){
   bool forward=j==i+1 || ((it>>freebit++)&1);
   // Short-circuit evaluation increments freebit only for free edges.
   if(forward){out[i]|=1<<j;mask|=1u<<edgebit;}else out[j]|=1<<i;
  }
  int dp[N][n]={},cycles[N]={},totals[n+1]={},pairtotals[4]={};
  for(int i=0;i<n;i++)dp[1<<i][i]=1;
  for(int s=1;s<N;s++) {
   int r=root[s];
   for(int v=0;v<n;v++) if(dp[s][v]) {
    int ways=dp[s][v];
    if(size[s]>=3 && ((out[v]>>r)&1))cycles[s]+=ways;
    int todo=out[v]&afterroot[s]&~s;
    while(todo){int w=__builtin_ctz((unsigned)todo);todo&=todo-1;dp[s|(1<<w)][w]+=ways;}
   }
   totals[size[s]]+=cycles[s];
  }
  for(auto p:pairs)pairtotals[p.type]+=cycles[p.s]*cycles[p.t];
  Key key;for(int m=3;m<=8;m++)key.a[m-3]=totals[m];
  key.a[6]=pairtotals[0];key.a[7]=pairtotals[1];key.a[8]=pairtotals[2]+pairtotals[3];
  int H=1+2*(totals[3]+totals[5]+totals[7])+4*(pairtotals[0]+pairtotals[2]);
  auto [pos,added]=seen.emplace(key,Value{mask,H,pairtotals[2]});
  if(!added&&pos->second.H!=H){cout<<"COLLISION masks "<<pos->second.mask<<' '<<mask<<" H "<<pos->second.H<<' '<<H<<'\n';return 1;}
  if(!added&&pos->second.D35!=pairtotals[2]){cerr<<"INTERNAL ERROR: D35 changes but H unchanged\n";return 2;}
  presentation_H_sum+=H;if(H>maxH)maxH=H;if(H<minH)minH=H;
  if((it+1)%262144==0)cerr<<it+1<<" presentations; "<<seen.size()<<" keys\n";
 }
 if(argc>2){
  vector<pair<Key,Value>> rows(seen.begin(),seen.end());
  sort(rows.begin(),rows.end(),[](const auto&a,const auto&b){return a.first.a<b.first.a;});
  ofstream f(argv[2]);f<<"mask\tH\tc3\tc4\tc5\tc6\tc7\tc8\tD33\tD34\tE8\tD35\n";
  for(const auto &p:rows){f<<p.second.mask<<'\t'<<p.second.H;for(int v:p.first.a)f<<'\t'<<v;f<<'\t'<<p.second.D35<<'\n';}
 }
 cout<<"NO COLLISION count "<<limit<<" keys "<<seen.size()<<"\n";
 cout<<"Presentation H sum "<<presentation_H_sum<<" min "<<minH<<" max "<<maxH<<"\n";
 cout<<"Cycle-pair placements "<<pairs.size()<<"\n";
 cout<<"All realized fibers also fix D35 and D44 separately.\n";
}

// Full simple-graph endpoint certificate. Two independent exact kernels.
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <string>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <stdexcept>
using namespace std;
constexpr int N=10, FULL=(1<<N)-1;
int adjm[N],ei[N][N],stamp[1<<N]={},val[1<<N],haf[1<<N],epoch=0;
vector<pair<int,int>> edges; vector<uint64_t> patterns[5];
void require(bool ok,const string&s){if(!ok)throw runtime_error(s);}
void make_patterns(int mask,int k,uint64_t em,int total){
 if(!k){patterns[total].push_back(em);return;}
 if(__builtin_popcount((unsigned)mask)<2*k)return;
 int v=__builtin_ctz((unsigned)mask),t=mask^(1<<v);
 make_patterns(t,k,em,total);
 for(int b=t;b;b&=b-1){int w=__builtin_ctz((unsigned)b);make_patterns(t^(1<<w),k-1,em|(1ull<<ei[v][w]),total);}
}
void recur(int s){
 if(stamp[s]==epoch)return;stamp[s]=epoch;
 if(!s){val[s]=haf[s]=1;return;}
 int v=__builtin_ctz((unsigned)s),t=s^(1<<v),k=__builtin_popcount((unsigned)s);
 recur(t);int z=(k%2?1:4)*val[t],h=0;
 for(int b=adjm[v]&t;b;b&=b-1){int u=t^(b&-b);recur(u);z-=val[u];h+=haf[u];}
 val[s]=z;haf[s]=(k%2?0:h);
}
int p21_direct(){
 int ans=0;
 for(auto [i,j]:edges)if(adjm[i]>>j&1){
  int rem=FULL^(1<<i)^(1<<j);
  for(int b=rem;b;b&=b-1){int v=__builtin_ctz((unsigned)b);int d=__builtin_popcount((unsigned)(adjm[v]&rem));ans+=d*(d-1)/2;}
 }
 return ans;
}
string quoted(const string&s){string r="\"";for(char c:s){if(c=='\\'||c=='\"')r+='\\';r+=c;}return r+'\"';}
int main(int argc,char**argv){try{
 string mode=argc>1?argv[1]:"dual";require(mode=="recursive"||mode=="patterns"||mode=="dual","kernel");
 long long expected=argc>2?stoll(argv[2]):12005168;
 for(int j=1;j<N;j++)for(int i=0;i<j;i++){ei[i][j]=ei[j][i]=edges.size();edges.push_back({i,j});}
 for(int k=2;k<=4;k++)make_patterns(FULL,k,0,k);
 require(patterns[2].size()==630&&patterns[3].size()==3150&&patterns[4].size()==4725,"pattern universe");
 long long count=0,minimum=numeric_limits<long long>::max();array<long long,46>hist{};vector<string>zeros;
 string line;
 while(getline(cin,line)){
  require(line.size()==9 && line[0]=='I',"strict order-ten graph6 row");
  fill(adjm,adjm+N,0);uint64_t graph=0;
  for(int e=0;e<45;e++){int c=(unsigned char)line[1+e/6]-63;require(c>=0&&c<=63,"graph6 alphabet");if(c>>(5-e%6)&1){auto[i,j]=edges[e];adjm[i]|=1<<j;adjm[j]|=1<<i;graph|=1ull<<e;}}
  require(((line.back()-63)&7)==0,"graph6 padding");
  int d[N],t[N]={},m=__builtin_popcountll(graph),w=0;
  for(int i=0;i<N;i++){d[i]=__builtin_popcount((unsigned)adjm[i]);w+=d[i]*(d[i]-1)/2;}
  for(int i=0;i<N;i++)for(int b=adjm[i];b;b&=b-1)t[i]+=d[__builtin_ctz((unsigned)b)]-1;
  int p21=0;
  for(auto[i,j]:edges)if(adjm[i]>>j&1)p21+=w-d[i]*(d[i]-1)/2-d[j]*(d[j]-1)/2-t[i]-t[j]+d[i]+d[j]-2+__builtin_popcount((unsigned)(adjm[i]&adjm[j]));
  long long r=-1,r2=-1;
  if(mode!="patterns"){
   ++epoch;recur(FULL);
   long long scaled=1024-val[FULL]-haf[FULL]-156*m+16*p21;
   require(scaled%4==0,"scaled divisibility");r=scaled/4;
  }
  if(mode!="recursive"){
   int mk[5]={1,m,0,0,0};
   for(int k=2;k<=4;k++)for(uint64_t p:patterns[k])mk[k]+=(graph&p)==p;
   int p2=p21_direct();require(p2==p21,"complement wedge count");
   r2=25*m-16*mk[2]+4*mk[3]-mk[4]+4*p2;
   if(mode=="dual")require(r2==r,"independent kernels disagree");else r=r2;
  }
  require(r>=0,"negative residual at "+line);
  minimum=min(minimum,r);hist[m]++;count++;
  if(r==0){zeros.push_back(line);require(m==0||(m==5&&all_of(d,d+N,[](int x){return x==1;})),"unexpected zero type");}
 }
 require(count==expected,"wrong universe count");
 cout<<"{\"status\":\"passed\",\"n\":10,\"kernel\":"<<quoted(mode)<<",\"graphs\":"<<count<<",\"minimum\":"<<minimum<<",\"zero_graph6\":[";
 for(size_t i=0;i<zeros.size();i++){if(i)cout<<',';cout<<quoted(zeros[i]);}cout<<"],\"edge_histogram\":[";
 for(int i=0;i<=45;i++){if(i)cout<<',';cout<<hist[i];}cout<<"],\"matching_pattern_counts\":[630,3150,4725]}\n";
}catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <tuple>
#include <vector>
using I=long long; using Row=std::array<I,6>;
using P=std::array<I,5>;
std::map<P,std::array<I,3>> lib;
void visit(int n,I code){
 const int V=1<<n,F=V-1;
 int out[9]={}; int k=0;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
  if(code>>k++&1)out[i]|=1<<j;else out[j]|=1<<i;
 }
 I en[512][9]={}, st[512][9]={}, all[512]={};
 for(int i=0;i<n;i++){en[1<<i][i]=1;st[1<<i][i]=1;}
 all[0]=1;
 for(int s=1;s<V;s++)for(int i=0;i<n;i++)if(s>>i&1){
  all[s]+=en[s][i];
  int To=F^s;
  for(int j=0;j<n;j++)if(To>>j&1){
   if(out[i]>>j&1)en[s|(1<<j)][j]+=en[s][i];
   else st[s|(1<<j)][j]+=st[s][i];
  }
 }
 I q=0,E[9]={},L[9]={},J[9][9]={};
 for(int s=0;s<V;s++){
  q+=all[s]*all[F^s];
  for(int u=0;u<n;u++){
   E[u]+=en[s][u]*all[F^s];L[u]+=st[s][u]*all[F^s];
   for(int v=0;v<n;v++)J[u][v]+=en[s][u]*st[F^s][v];
  }
 }
 // All masks, including the empty one, are legal for guarded tapes.
 I e[512]={},l[512]={},j[512]={};
 for(int x=0;x<V;x++){
  if(x){ int u=__builtin_ctz(x), y=x&(x-1); e[x]=e[y]+E[u];l[x]=l[y]+L[u]; j[x]=j[y];
   for(int v=0;v<n;v++) if(y>>v&1) j[x]+=J[u][v]+J[v][u];
  }
  if(e[x]>=l[x]) lib.emplace(P{all[F],q,e[x],l[x],j[x]},std::array<I,3>{n,code,x});
  else { // Keep a normalized copy under reversal of every edge, not just the signature.
   I mirrored=code^((1LL<<(n*(n-1)/2))-1);
   lib.emplace(P{all[F],q,l[x],e[x],j[x]},std::array<I,3>{n,mirrored,x});
  }
 }
}
int main(int argc,char**argv){
 int maxN=argc>1?std::atoi(argv[1]):6;
 std::mt19937_64 rng(20261115);
 for(int n=1;n<=maxN;n++){
  I lim=1LL<<(n*(n-1)/2);
  if(n<=7){
   std::vector<bool> used(lim);
   for(I c=0;c<lim;c++)if(!used[c]){
    visit(n,c);
    bool a[9][9]={};int ix=0;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){a[i][j]=(c>>ix++)&1;a[j][i]=!a[i][j];}
    std::vector<int> pi(n);std::iota(pi.begin(),pi.end(),0);
    do{I y=0;int ix=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)y|=(I)a[pi[i]][pi[j]]<<ix++;used[y]=true;}while(std::next_permutation(pi.begin(),pi.end()));
   }
  }else{
   for(int t=0;t<6000;t++)visit(n,rng()%lim);
  }
  fprintf(stderr,"n=%d distinct partial profiles=%zu\n",n,lib.size());
 }
 // A common t factor grows both guards: a and b are relative shifts.
 using Coef=std::pair<I,I>; using V=std::tuple<I,I,I,I,I,I,I,I,I,I>;
 std::map<Coef,std::map<I,V>> sh;
 for(const auto& [p,r]:lib){
  auto [h,q,e,l,j]=p; auto [n,c,v]=r;
  for(int flip=0;flip<(e!=l?2:1);flip++){
   for(int a=0;a<=4;a++)for(int b=0;b<=4;b++){
    I X=1LL<<a,Y=1LL<<b;
    sh[{q*X*Y,e*X-l*Y}].emplace(h-j,V{h,q,e,l,j,n,c,v,a,b});
   }
   std::swap(e,l);c^=(1LL<<(n*(n-1)/2))-1;
  }
 }
 auto dump=[](const V&v){auto [h,q,e,l,j,n,c,m,a,b]=v;
   printf("%lld %lld %lld %lld %lld ; %lld %lld %lld ; %lld %lld\n",h,q,e,l,j,n,c,m,a,b);};
 int longest=0,moving=0, printed=0; 
 for(const auto& [coef,table]:sh){
  for(auto it=table.begin();it!=table.end();it++){
   int len=1; for(auto nx=std::next(it);nx!=table.end()&&nx->first==it->first+2*len;nx++,len++){}
   if(len>longest){
    longest=len; printf("RUN %d high=(%lld,%lld) start=%lld\n",len,coef.first,coef.second,it->first);
    auto nx=it; for(int k=0;k<len;k++,nx++)dump(nx->second);
   }
   if(auto nx=table.find(it->first+2);nx!=table.end()&&std::get<0>(nx->second)!=std::get<0>(it->second)&&std::get<1>(nx->second)!=std::get<1>(it->second)&&std::get<2>(nx->second)!=std::get<2>(it->second)&&std::get<3>(nx->second)!=std::get<3>(it->second)){
    moving++;if(printed<12 && coef.second!=0){printed++;printf("MOVING high=(%lld,%lld) offsets=%lld,%lld\n",coef.first,coef.second,it->first,nx->first);dump(it->second);dump(nx->second);}
   }
  }
 }
 printf("total %d moving pairs preferred by table (not an exhaustive pair multiplicity)\n",moving);
}

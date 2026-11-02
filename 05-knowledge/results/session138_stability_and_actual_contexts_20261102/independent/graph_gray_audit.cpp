// Independent full-universe graph certificate using Gray-code local updates.
// No producer counting loop is used. Each toggle adds/removes patterns that
// contain that edge, partitioned by its role in the pattern.
#include <array>
#include <vector>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <cstdint>
#include <algorithm>
using namespace std;
void need(bool b,const char*s){if(!b)throw runtime_error(s);}
struct Audit {
 int n; vector<pair<int,int>> edges;array<unsigned,7>a{};
 int m=0,w=0,m2=0,m3=0,p21=0;
 explicit Audit(int nn):n(nn){for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)edges.push_back({i,j});}
 int pop(unsigned x){return __builtin_popcount(x);}
 pair<int,int> induced(unsigned vs){int es=0,ws=0;for(int i=0;i<n;i++)if(vs>>i&1){int d=pop(a[i]&vs);es+=d;ws+=d*(d-1)/2;}return {es/2,ws};}
 void change(int edge,bool add){
  auto [u,v]=edges[edge];unsigned ub=1u<<u,vb=1u<<v,all=(1u<<n)-1;
  if(!add){a[u]^=vb;a[v]^=ub;}
  unsigned rem=all^ub^vb;auto [em,wm]=induced(rem);
  int dw=pop(a[u])+pop(a[v]);int dm2=em,dm3=em*(em-1)/2-wm;
  int dp21=wm;
  for(int x=0;x<n;x++)if(rem>>x&1){
   int mult=((a[u]>>x)&1)+((a[v]>>x)&1);
   if(mult)dp21+=mult*induced(rem^(1u<<x)).first;
  }
  int s=add?1:-1;m+=s;w+=s*dw;m2+=s*dm2;m3+=s*dm3;p21+=s*dp21;
  if(add){a[u]^=vb;a[v]^=ub;}
 }
 void brute_check(){
  vector<int> present;for(int i=0;i<(int)edges.size();i++){auto[u,v]=edges[i];if(a[u]>>v&1)present.push_back(i);}
  int bm2=0,bm3=0,bw=0,bp=0;
  auto verts=[&](int i){auto[u,v]=edges[i];return (1u<<u)|(1u<<v);};
  for(int i=0;i<(int)present.size();i++)for(int j=i+1;j<(int)present.size();j++){
   unsigned vi=verts(present[i]),vj=verts(present[j]);bool disjoint=!(vi&vj);
   if(disjoint)bm2++;else bw++;
   for(int k=j+1;k<(int)present.size();k++){
    unsigned vk=verts(present[k]);
    if(disjoint&&!(vi&vk)&&!(vj&vk))bm3++;
    if((vi&vj)&&!((vi|vj)&vk))bp++;
    if((vi&vk)&&!((vi|vk)&vj))bp++;
    if((vj&vk)&&!((vj|vk)&vi))bp++;
   }
  }
  need(m==(int)present.size()&&m2==bm2&&m3==bm3&&w==bw&&p21==bp,"independent brute snapshot mismatch");
 }
 void run(){
  int total=1<<edges.size(),minimum=1000000,zeros=0;unsigned witness=0;
  int maxnum=0,maxden=1;unsigned alphawitness=0;uint64_t digest=1469598103934665603ULL;int snapshots=0;
  for(int t=0;t<total;t++){
   unsigned g=(unsigned)t^((unsigned)t>>1);
   if(t){int e=__builtin_ctz((unsigned)t);change(e,(g>>e)&1);}
   if(t<16 || t%4093==0 || t==total-1){brute_check();snapshots++;}
   need(w==m*(m-1)/2-m2,"wedge partition mismatch");
   int base=22*m-24*m2+6*m3+6*p21;
   int residual=n==7?base+8*w:11*m-12*m2+3*m3+3*p21;
   if(residual<minimum){minimum=residual;witness=g;}if(!residual)zeros++;
   if(n==7&&w&&(-base)*maxden>maxnum*(24*w)){maxnum=-base;maxden=24*w;alphawitness=g;int d=gcd(maxnum,maxden);maxnum/=d;maxden/=d;}
   digest^=(uint64_t)(residual+100000);digest*=1099511628211ULL;
  }
  need(minimum==0,"graph certificate failed");
  if(n==7)need(maxnum==1&&maxden==3,"optimal graph alpha mismatch");
  cout<<"{\"n\":"<<n<<",\"graphs\":"<<total<<",\"minimum\":"<<minimum<<",\"zero_graphs\":"<<zeros<<",\"witness\":"<<witness<<",\"independent_brute_snapshots\":"<<snapshots<<",\"gray_residual_digest\":"<<digest;
  if(n==7)cout<<",\"sharp_alpha\":\""<<maxnum<<"/"<<maxden<<"\",\"alpha_witness\":"<<alphawitness;
  cout<<"}\n";
 }
};
int main(){Audit(6).run();Audit(7).run();}

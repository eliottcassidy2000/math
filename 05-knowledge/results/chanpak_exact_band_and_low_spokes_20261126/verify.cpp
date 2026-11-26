#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using I=int64_t;
void need(bool ok,const char* s){if(!ok) throw runtime_error(s);}
int bits(I u){int n=0;for(;u;u>>=1)++n;return n;}
I independent(I u) {
 // Weighted inversion sum, scanning most significant to least significant.
 I pref=0,H=1; int prev=0;
 for (int j=bits(u)-1;j>=0;--j) {
    int b=(u>>j)&1; if(!b) H+=pref+prev;
    pref=2*pref+b; prev=b;
 }
 return H;
}
struct G {I h,e,R;int n;};
int main(int argc,char**argv){
 need(argc==5,"usage: verify first_odd last_odd witnesses rooted_graphs");
 I lo=stoll(argv[1]),hi=stoll(argv[2]);
 ifstream gs(argv[4]); need(bool(gs),"no graph file");
 map<I,G> graph;
 for(string line;getline(gs,line);){
   istringstream ss(line);I id,h,e,R;need(bool(ss>>id>>h>>e>>R),"bad profile");vector<int> r;
   for(int a;ss>>a;)r.push_back(a);
   int n=r.size(),M=(1<<n)-1;need(n>=1&&n<=7,"order");
   for(int i=0;i<n;i++){
     need(r[i]>=0&&r[i]<=M&&!(r[i]>>i&1),"row");
     for(int j=0;j<i;j++)need((r[i]>>j&1)+(r[j]>>i&1)==1,"not tournament");
   }
   vector<I> H(M+1),S(M+1);I end=0;H[0]=1;
   // Literal path-word enumeration on every induced subset of the small bank.
   for(int s=1;s<=M;s++){
     vector<int> p;for(int i=0;i<n;i++)if(s>>i&1)p.push_back(i);
     do {
       bool ok=true;for(size_t i=1;i<p.size();i++)if(!(r[p[i-1]]>>p[i]&1)){ok=false;break;}
       if(ok){++H[s];if(p.front()==0)++S[s];if(s==M&&p.back()==0)++end;}
     }while(next_permutation(p.begin(),p.end()));
   }
   I two=0;for(int s=0;s<=M;s+=2)two+=H[s]*S[M^s];
   need(H[M]==h&&end==e&&two==R,"bad literal bank profile");
   need(graph.emplace(id,G{h,e,R,n}).second,"duplicate profile index");
 }
 ifstream in(argv[3]);need(bool(in),"no witness file");I which=lo,n,a,v,f1=0,f2=0;int cap1=0,cap2=0;
 setbuf(stdout,nullptr);
 while(in>>n>>a>>v){
   need(n==which&&n<=hi&&a&&v>=0,"missing or unordered target");
   I f=independent(v);
   if(a>0){
     need(a*(v+1)+independent(a-1)*f-(a&1)*!(v&1)==n,"F1 height");
     ++f1;cap1=max(cap1,bits(a-1)+bits(v)+2);
   }else{
     const auto& g=graph.at(-a);need(g.h*f-g.e*!(v&1)+g.R*(v+1)==n,"F2 height");
     ++f2;cap2=max(cap2,bits(v)+1+g.n);
   }
   which+=2;
 }
 need(which==hi+2,"incomplete band");
 cout<<"verified "<<lo<<" .. "<<hi<<" F1 "<<f1<<" F2 "<<f2<<" graph profiles "<<graph.size()
     <<" largest order F1 "<<cap1<<" F2 "<<cap2<<"\n";
}

// Exact longest-path certificate on defect residue x quadratic phase.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using I=std::int64_t;
void need(bool b,const char*m){if(!b)throw std::runtime_error(m);}
struct Node {int parent,a;};
int main(int argc,char**argv){
 int E=argc>1?std::stoi(argv[1]):4;need(E>=1&&E<=9,"levels one through nine");
 int q=1;for(int j=0;j<E;++j)q*=5;
 int cap=4*E+8,S=2*q;I scale=8LL*q,cost=5LL*q,NEG=-(I(1)<<60);
 std::vector<int> height(q,E),power(E+1,1),inv(cap+1,1),two(cap+1,1),kap(cap+1,0);
 for(int j=1;j<=E;++j)power[j]=5*power[j-1];
 for(int z=1;z<q;++z){int a=z,h=0;while(a%5==0){++h;a/=5;}height[z]=h;}
 for(int a=1;a<=cap;++a){inv[a]=I(inv[a-1])*((q+1)/2)%q;two[a]=2*two[a-1]%q;kap[a]=height[(two[a]-2+q)%q];}
 std::vector<I>V(S,NEG);std::vector<int>path(S,-1);std::vector<unsigned char>in(S,0);std::deque<int>queue;
 std::vector<Node>nodes;
 auto transition=[&](int s,int a){int z=s/2,active=s%2;return int(((3LL*z+two[a]-2+q)*inv[a])%q)*2+(active^(a%2==0));};
 auto delta=[&](int s,int a,int t)->I{
  if(a==1)return s%2? -cost/power[height[s/2]]:-cost;
  return (t%2?scale*std::max(height[t/2],kap[a]):0)-cost*a;
 };
 auto update=[&](int t,I score,int parent,int a){
  if(score<=V[t])return;
  V[t]=score;nodes.push_back({parent,a});path[t]=int(nodes.size())-1;
  if(score>0){std::vector<int>w;int j=path[t];while(j>=0){w.push_back(nodes[j].a);j=nodes[j].parent;}std::reverse(w.begin(),w.end());
   std::cout<<"{\"counterexample\":[";for(size_t k=0;k<w.size();++k){if(k)std::cout<<",";std::cout<<w[k];}std::cout<<"],\"scaled_debt\":"<<score<<",\"scale\":"<<scale<<"}\n";std::exit(2);}
  if(!in[t]){in[t]=1;queue.push_back(t);}
 };
 // Initial ones have zero defect and do not alter quadratic phase.
 for(int a=2;a<=cap;++a){int t=transition(1,a);update(t,delta(1,a,t),-1,a);}
 I popped=0;
 while(!queue.empty()){
  int s=queue.front();queue.pop_front();in[s]=0;++popped;need(popped<200000000,"search guard");
  I score=V[s];int parent=path[s];
  for(int a=1;a<=cap;++a){int t=transition(s,a);update(t,score+delta(s,a,t),parent,a);}
 }
 I lo=*std::min_element(V.begin(),V.end()),hi=*std::max_element(V.begin(),V.end());
 need(lo>NEG,"every state reached");need(hi<=0,"nonpositive budget debt");
 need(cost*(cap+1)>=scale*E-lo,"unbounded letter and seed tail");
 I checked=0;
 for(int s=0;s<S;++s)for(int a=1;a<=cap;++a){int t=transition(s,a);need(V[t]>=V[s]+delta(s,a,t),"edge certificate");++checked;}
 for(int a=2;a<=cap;++a){int t=transition(1,a);need(V[t]>=delta(1,a,t),"seed certificate");}
 if(argc>2){std::ofstream f(argv[2],std::ios::binary);need(bool(f),"output");for(I x:V)for(int j=0;j<8;++j)f.put(char((std::uint64_t(x)>>(8*j))&255));need(bool(f),"complete write");}
 int best=int(std::max_element(V.begin(),V.end())-V.begin());std::vector<int>w;int j=path[best];while(j>=0){w.push_back(nodes[j].a);j=nodes[j].parent;}std::reverse(w.begin(),w.end());
 std::cout<<"{\"E\":"<<E<<",\"q\":"<<q<<",\"states\":"<<S<<",\"scale\":"<<scale<<",\"cap\":"<<cap<<",\"min\":"<<lo<<",\"max\":"<<hi<<",\"popped\":"<<popped<<",\"updates\":"<<nodes.size()<<",\"inequalities\":"<<checked<<",\"best_word\":[";for(size_t i=0;i<w.size();++i){if(i)std::cout<<",";std::cout<<w[i];}std::cout<<"]}\n";
}

// Independent compact graph reconstruction and Kahn pair propagation.
// No producer includes, scheduling, reverse graph, or certificate receipt imported.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>
using namespace std;using I=int64_t;using U=uint64_t;using P=pair<I,I>;
I modpow(I a,I n,I q){I r=1;for(;n;n>>=1,a=(__int128)a*a%q)if(n&1)r=(__int128)r*a%q;return r;}
struct Node{U key;uint32_t first,count;};
struct Arc{uint32_t source,target,next;uint8_t letter;};
int main(int argc,char**argv){try{
 if(argc!=3)throw runtime_error("usage E binary");int E=stoi(argv[1]);if(E<2||E>10)throw runtime_error("E");const int R=2,B=2;
 vector<int>p(E+1,1);for(int i=1;i<=E;i++)p[i]=5*p[i-1];int q=p[E];
 auto height=[&](int z){if(!z)return E;int e=0;for(;z%5==0;z/=5)++e;return e;};
 vector<uint32_t>base(E+2),units(E+1);for(int e=0;e<=E;e++){units[e]=4*p[max(R,e)-1];U n=e==E?1:4*q/p[e+1];U end=base[e]+n*units[e];if(end>numeric_limits<uint32_t>::max())throw runtime_error("index overflow");base[e+1]=end;}
 uint32_t states=base[E+1];auto rank=[](int u){return u-u/5-1;};auto unit=[](int u){return u+1+u/4;};
 auto sid=[&](int z,int m){int e=height(z);return base[e]+(z?rank(z/p[e])*units[e]:0)+rank(m);};
 auto decode=[&](uint32_t s){int e=upper_bound(base.begin(),base.end(),s)-base.begin()-1;int i=s-base[e];return pair<int,int>{e==E?0:p[e]*unit(i/units[e]),unit(i%units[e])};};
 ifstream f(argv[2],ios::binary|ios::ate);if(!f||f.tellg()!=8ULL*states)throw runtime_error("binary length");f.seekg(0);vector<int16_t>pot(states);array<uint8_t,65536>buf;
 for(U at=0;at<states;){int n=min<U>(8192,states-at);f.read((char*)buf.data(),8*n);if(!f)throw runtime_error("binary read");for(int i=0;i<n;i++){U u=0;for(int j=0;j<8;j++)u|=U(buf[8*i+j])<<(8*j);I x=(I)u;if(x<=-30000||x>0)throw runtime_error("encoded potential");pot[at+i]=x;}at+=n;}
 int low=*min_element(pot.begin(),pot.end()),high=*max_element(pot.begin(),pot.end());
 if(low-high-8*E+325<=B||low-8*E+325<=B)throw runtime_error("large letter gap");
 I inv3=modpow(3,4LL*q/5-1,q);vector<int>twopow(65,1);for(int a=1;a<=64;a++)twopow[a]=2LL*twopow[a-1]%q;
 vector<Node>nodes;vector<Arc>arcs;unordered_map<U,uint32_t>index;
 auto add=[&](U k){auto it=index.find(k);if(it!=index.end())return it->second;if(nodes.size()>=6000000)throw runtime_error("node resource bound");uint32_t i=nodes.size();index.emplace(k,i);nodes.push_back({k,0,0});return i;};
 array<uint32_t,3>roots;uint32_t terminal=sid(0,1);for(int b=0;b<=2;b++)roots[b]=add(U(terminal)*3+b);
 U candidates=0;for(uint32_t i=0;i<nodes.size();i++){
  U key=nodes[i].key;uint32_t state=key/3;int budget=key%3;auto [z,m]=decode(state);int h=height(z),mod=p[max(R,h)],charge=min(h,height(m-1));uint32_t start=arcs.size();
  for(int a=1;a<=64;a++){
   int required=pot[state]-8*charge+5*a;if(required<low||required-budget>high)continue;
   I v=(I(twopow[a])*(z-1)+2)%q;if(v<0)v+=q;int oldz=v*inv3%q,oldmod=p[max(R,height(oldz))],common=min(mod,oldmod);
   int oldresidue=(I(m)*twopow[a]%common)*inv3%common;
   for(int oldm=oldresidue;oldm<oldmod;oldm+=common){++candidates;uint32_t old=sid(oldz,oldm);int slack=required-pot[old];if(slack<0||slack>budget)continue;uint32_t j=add(U(old)*3+budget-slack);arcs.push_back({j,i,0,uint8_t(a)});}
  }
  nodes[i].first=start;nodes[i].count=arcs.size()-start;
 }
 // Independently enumerate every canonical and one-digit extra entry.
 unordered_map<U,vector<P>>entry;U canonical=0,longentries=0;
 auto ones=[&](I k){I sum=0;for(int e=1;e<=E;e++)sum+=k/(2*p[e-1]);return sum;};
 for(int a=2;a<=64;a++){
  I inv2=modpow(twopow[a],4LL*q/5-1,q);int z=(q+1-2*inv2%q)%q,h=height(z),ell=max(R,h),mod=p[ell],T=2*p[ell-1],rho=3LL*((mod+1)/2)%mod;
  for(int s=0;s<T;s++){
   int m=3LL*modpow(rho,s,mod)%mod*inv2%mod;uint32_t state=sid(z,m);int reward=min(h,height(m-1));I score=8*(ones(s)+reward)-5LL*(s+a),gap=pot[state]-score;
   if(gap<0)throw runtime_error("invalid entry");if(gap<=2)entry[U(state)*3+gap].push_back({s+a,s+1});
   if(gap==0){for(int e=ell;e<=E;e++){I k=s+2*p[e-1];I longscore=8*(ones(k)+reward)-5*(k+a);if(longscore!=score-2||3LL*modpow(rho,k,mod)%mod*inv2%mod!=m)throw runtime_error("long entry exact check");entry[U(state)*3+2].push_back({k+a,k+1});++longentries;}}
   ++canonical;
  }
 }
 // Convert to outgoing adjacency and process source-before-target with Kahn.
 const uint32_t NONE=numeric_limits<uint32_t>::max();vector<uint32_t>head(nodes.size(),NONE),indegree(nodes.size());
 for(uint32_t e=0;e<arcs.size();e++){auto &a=arcs[e];a.next=head[a.source];head[a.source]=e;++indegree[a.target];}
 vector<vector<P>>pairs(nodes.size());vector<U>paths(nodes.size(),0);U retained=0;vector<uint32_t>queue;queue.reserve(nodes.size());
 for(uint32_t i=0;i<nodes.size();i++){auto it=entry.find(nodes[i].key);if(it!=entry.end()){pairs[i]=it->second;paths[i]=it->second.size();retained+=paths[i];}if(!indegree[i])queue.push_back(i);}
 for(size_t at=0;at<queue.size();at++){
  uint32_t s=queue[at];auto &ps=pairs[s];sort(ps.begin(),ps.end());ps.erase(unique(ps.begin(),ps.end()),ps.end());
  for(uint32_t e=head[s];e!=NONE;e=arcs[e].next){auto a=arcs[e];if(numeric_limits<U>::max()-paths[a.target]<paths[s])throw runtime_error("path count overflow");paths[a.target]+=paths[s];for(auto [A,L]:ps)pairs[a.target].push_back({A+a.letter,L+1});if(pairs[a.target].size()>1000000)throw runtime_error("pair resource bound");if(!--indegree[a.target])queue.push_back(a.target);}
 }
 if(queue.size()!=nodes.size())throw runtime_error("graph cycle");
 cout<<"{\"status\":\"PASS\",\"E\":"<<E<<",\"R\":2,\"budget\":2,\"nodes\":"<<nodes.size()<<",\"edges\":"<<arcs.size()<<",\"scanned_candidates\":"<<candidates<<",\"entries_checked\":"<<canonical<<",\"long_entry_controls\":"<<longentries<<",\"retained_entries\":"<<retained<<",\"acyclic\":true,\"records\":[";
 for(int b=0;b<=2;b++){if(b)cout<<",";auto i=roots[b];cout<<"{\"slack\":"<<b<<",\"score\":"<<-pot[terminal]+b<<",\"abstract_paths\":"<<paths[i]<<",\"pairs\":[";bool first=true;U survivors=0;for(auto [A,L]:pairs[i]){if(!(0<L&&L<A))throw runtime_error("invalid word pair");if(!first)cout<<",";first=false;cout<<"["<<A<<","<<L<<"]";if(modpow(3,L,q)==modpow(2,A,q))++survivors;}cout<<"],\"multiplier_survivors_mod_5E\":"<<survivors<<"}";}
 cout<<"]}\n";
 }catch(exception&e){cerr<<e.what()<<"\n";return 1;}}

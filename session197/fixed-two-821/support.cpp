#include <bits/stdc++.h>
using namespace std;
struct State{array<uint16_t,27>x;bool operator==(State const&b)const{return x==b.x;}};
struct Hasher{size_t operator()(State const&s)const{uint64_t h=0x9e3779b97f4a7c15ULL;for(auto x:s.x){h^=x+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2);}return h;}};
#include "automaton.inc"
int main(int argc,char**argv){
 int X=argc>1?atoi(argv[1]):173;if(X<1||X>65534)return 1;
 int cycle_limit=argc>3?atoi(argv[3]):3; if(cycle_limit<1||cycle_limit>X)return 1;
 int cap=X+1;size_t stop=argc>2?atoll(argv[2]):5000000; if(stop<2||stop>UINT32_MAX)return 1;
 unordered_set<State,Hasher>seen;vector<State>Q;vector<uint32_t>parent,depth;vector<uint8_t>letter;
 seen.insert(seed);Q.push_back(seed);parent.push_back(0);depth.push_back(0);letter.push_back(255);
 struct Witness{uint32_t p;int g,c;};map<int,Witness>support;
 size_t processed=0,edges=0,cycle_rejected=0;bool complete=true;
 for(size_t pos=0;pos<Q.size()&&complete;++pos){
  State st=Q[pos];processed=pos+1;
  for(int g=2;g<12;g+=3){
   int S=g/3;int C=S==0?0:S==1?st.x[1]+st.x[2]:S==2?st.x[3]+st.x[5]:st.x[2]+st.x[5];
   if(C>cycle_limit){cycle_rejected++;continue;}
   State ns{};
   for(int i=0;i<27;i++)if(st.x[i])for(auto [j,m]:tr[g][i])ns.x[j]=min(cap,(int)ns.x[j]+st.x[i]*m);
   int h=0;for(int i:hc)h+=ns.x[i];if(h>X)continue;edges++;
   if(C&& !support.count(h))support[h]={(uint32_t)pos,g,C};
   if(seen.insert(ns).second){Q.push_back(ns);parent.push_back(pos);depth.push_back(depth[pos]+1);letter.push_back(g);if(Q.size()>=stop){complete=false;break;}}
  }
 }
 cout<<"{\"bound\":"<<X<<",\"cycle_limit\":"<<cycle_limit<<",\"complete\":"<<(complete?"true":"false")<<",\"states\":"<<Q.size()<<",\"processed\":"<<processed<<",\"admissible_edges\":"<<edges<<",\"cycle_rejected_edges\":"<<cycle_rejected<<",\"max_shortest_word\":"<<*max_element(depth.begin(),depth.end())<<",\"missing_odds\":[";
 bool first=true;for(int h=3;h<=X;h+=2)if(!support.count(h)){if(!first)cout<<',';cout<<h;first=false;}
 cout<<"],\"witnesses\":{";first=true;
 for(auto [h,w]:support){if(!first)cout<<',';first=false;vector<int>v={w.g};for(uint32_t p=w.p;p;p=parent[p])v.push_back(letter[p]);reverse(v.begin(),v.end());cout<<'"'<<h<<"\":{\"C\":"<<w.c<<",\"word\":[";for(size_t i=0;i<v.size();i++){if(i)cout<<',';cout<<'['<<v[i]/3<<','<<v[i]%3<<']';}cout<<"]}";}
 cout<<"}}\n";return complete?0:2;
}

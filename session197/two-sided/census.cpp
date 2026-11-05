#include <array>
#include <cstdint>
#include <vector>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <set>
using namespace std;
struct State { array<uint16_t,27>x;bool operator==(const State&o)const{return x==o.x;} };
struct Hash {size_t operator()(State const&s)const{uint64_t h=1469598103934665603ULL;for(auto v:s.x){h^=v;h*=1099511628211ULL;}return h;}};
#include "automaton.inc"
int main(int argc,char**argv){if(argc!=2)return 1;ifstream f(argv[1]);if(!f)return 2;int p,q;bool firstcase=true;cout<<"[";while(f>>p>>q){State init{};int v;for(int i=0;i<27;i++){if(!(f>>v)||v<0||v>152)return 3;init.x[i]=v;}unordered_set<State,Hash>seen;seen.insert(init);vector<State>queue={init};set<int>support;
for(size_t pos=0;pos<queue.size();pos++){State st=queue[pos];int h0=0;for(int i:hc)h0+=st.x[i];support.insert(h0);for(int g=0;g<12;g++){State next{};for(int i=0;i<27;i++)for(auto [j,m]:tr[g][i])next.x[j]=min(152,(int)next.x[j]+m*st.x[i]);int h=0;for(int i:hc)h+=next.x[i];if(h>151)continue;if(seen.insert(next).second)queue.push_back(next);}}
sort(queue.begin(),queue.end(),[](State const&a,State const&b){return a.x<b.x;});if(!firstcase)cout<<',';firstcase=false;cout<<"{\"ports\":["<<p<<","<<q<<"],\"support\":[";bool sep=false;for(int h:support){if(sep)cout<<',';sep=true;cout<<h;}cout<<"],\"states\":[";sep=false;for(auto const&st:queue){if(sep)cout<<',';sep=true;cout<<'[';for(int i=0;i<27;i++){if(i)cout<<',';cout<<st.x[i];}cout<<']';}cout<<"]}";}cout<<"]\n";return 0;}

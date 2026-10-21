#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
constexpr int N=12, FULL=(1<<N)-1;
struct Orbit { vector<pair<int,int>> edges; int forced=-1; };
vector<Orbit> orbits;
array<int,N> adj,deg;
array<int,N+1> groups;
array<array<unsigned, (1<<N)*N>, N+1> seen{};
array<unsigned,N+1> stamp{};
struct State { uint16_t mask; uint8_t last; uint64_t path; };
vector<uint64_t> witnesses;
array<int,N> scores;
uint64_t prefixes=0;
void insist(bool x,const string &s){if(!x)throw runtime_error(s);}
void init(){
    bool done[N][N]{};
    for(int a=0;a<N;a++)for(int b=a+1;b<N;b++)if(!done[a][b]){
        Orbit o;
        for(int k=0;k<3;k++){
            int x=(a+4*k)%N,y=(b+4*k)%N;
            done[min(x,y)][max(x,y)]=true;
            o.edges.push_back({x,y});
            int v=((x+1)%N==y)?1:((y+1)%N==x?0:-1);
            if(v>=0){insist(o.forced<0||o.forced==v,"inconsistent cycle");o.forced=v;}
        }
        orbits.push_back(o);
    }
    insist(orbits.size()==22,"orbit count");
    insist(count_if(orbits.begin(),orbits.end(),[](const Orbit&o){return o.forced<0;})==18,"free count");
}
void frame(unsigned mask){
    adj.fill(0);int k=0;
    for(auto&o:orbits){int v=o.forced<0?((mask>>k++)&1):o.forced;
        for(auto[x,y]:o.edges){if(!v)swap(x,y);adj[x]|=1<<y;}}
    groups.fill(0);
    for(int i=0;i<N;i++){deg[i]=__builtin_popcount((unsigned)adj[i]);groups[deg[i]]|=1<<i;}
}
bool seek(int depth,const vector<State>&states){
    prefixes++;
    if(depth==N){
        bool palindrome=true;
        for(int i=0;i<N;i++)if(scores[i]+scores[N-1-i]!=N-1)palindrome=false;
        if(!palindrome)witnesses.push_back(states.front().path);
        return witnesses.size()>=8;
    }
    for(int d=0;d<N;d++)if(groups[d]){
        vector<State> next;auto tick=++stamp[depth];
        for(auto s:states){
            unsigned rem=groups[d]&~s.mask;
            if(depth)rem&=adj[s.last];
            while(rem){unsigned bit=rem&-rem;rem-=bit;int v=__builtin_ctz(bit),mask=s.mask|bit;
                int ix=mask*N+v;
                if(seen[depth][ix]==tick)continue;
                seen[depth][ix]=tick;next.push_back({(uint16_t)mask,(uint8_t)v,s.path|((uint64_t)v<<(4*depth))});
            }
        }
        if(!next.empty()){scores[depth]=d;if(seek(depth+1,next))return true;}
    }
    return false;
}
int main(int argc,char**argv){
    init();unsigned limit=(argc>1?stoul(argv[1]):1<<18);
    insist(limit>0&&limit<=(1u<<18),"frame limit");
    ofstream bank;
    if(argc>2){bank.open(argv[2],ios::binary);insist(bool(bank),"bank open");bank.write("D12SCR01",8);}
    int nonsc=0,passed=0;vector<unsigned> residual;auto start=chrono::steady_clock::now();
    for(unsigned mask=0;mask<limit;mask++){
        frame(mask);bool symmetric=true;
        for(int d=0;d<N;d++)if(__builtin_popcount((unsigned)groups[d])!=__builtin_popcount((unsigned)groups[N-1-d]))symmetric=false;
        if(!symmetric){nonsc++;continue;}
        witnesses.clear();vector<State> first{{0,0,0}};
        if(seek(0,first)){
            passed++;
            if(bank.is_open())for(auto p:witnesses){
                for(int k=0;k<6;k++)bank.put(char((p>>(8*k))&255));
                bool kept[N][N]{};
                for(int k=0;k<N-1;k++){
                    int a=(p>>(4*k))&15,b=(p>>(4*(k+1)))&15;
                    kept[a][b]=kept[b][a]=true;
                }
                array<uint8_t,9> target{};int k=0;
                for(int a=0;a<N;a++)for(int b=a+1;b<N;b++,k++)
                    if(bool((adj[a]>>b)&1)^!kept[a][b])target[k/8]|=1u<<(k%8);
                bank.write(reinterpret_cast<const char*>(target.data()),9);
            }
        }else{residual.push_back(mask);cerr<<"residual "<<mask<<" witnesses "<<witnesses.size()<<'\n';}
    }
    double sec=chrono::duration<double>(chrono::steady_clock::now()-start).count();
    cout<<"frames "<<limit<<" asymmetric-score "<<nonsc<<" eight-black "<<passed<<" residual "<<residual.size()<<" prefixes "<<prefixes<<" seconds "<<sec<<'\n';
    if(bank.is_open()){bank.close();insist(bool(bank),"bank write");}
    return residual.empty()?0:2;
}

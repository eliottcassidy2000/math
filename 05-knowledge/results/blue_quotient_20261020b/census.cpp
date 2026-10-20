// Exact n=9 tournament universe, complete blue graph, certified black tree.
// Inherits the fixed-path/score conventions and n=8 class input of
// 04-computation/experiments/tiling_line_metagraph_census_20261020.cpp.
// Canonicalization is independently implemented ordered equitable refinement
// plus individualization; every inherited n=8 Aut value is checked.
// A black tree is a SUBGRAPH CERTIFICATE, not a full black-line census.
// Build: g++ -O3 -std=c++17 census.cpp -o census.exe
// Run: census.exe inherited_classes_n8.tsv output_prefix
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
using Code = uint64_t;
using Rows = array<uint16_t, 10>;
using Partition = vector<vector<int>>;

static void need(bool x, const string& s) { if (!x) throw runtime_error(s); }
static int pop(unsigned x) { return __builtin_popcount(x); }
static Code encode(const Rows& a, const vector<int>& p) {
    Code z=0;
    for (int i=0;i<(int)p.size();i++)
        for (int j=i+1;j<(int)p.size();j++) z=(z<<1)|((a[p[i]]>>p[j])&1);
    return z;
}
static Rows decode(Code z, int n) {
    Rows a{}; int bit=n*(n-1)/2-1;
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++,bit--)
        if((z>>bit)&1) a[i]|=1u<<j; else a[j]|=1u<<i;
    return a;
}
static Rows converse(const Rows& a,int n) {
    Rows b{}; for(int i=0;i<n;i++) b[i]=a[i]^((1u<<n)-1)^(1u<<i);
    return b;
}
static Code packed(const vector<int>& p) {
    Code z=0; for(int i=0;i<(int)p.size();i++) z|=Code(p[i])<<(4*i); return z;
}
static vector<int> unpack(Code z,int n) {
    vector<int> p(n); for(int i=0;i<n;i++) p[i]=(z>>(4*i))&15; return p;
}
static Partition refine(const Rows& a,Partition p) {
    while(true) {
        Partition q;
        for(auto cell:p) {
            map<vector<int>,vector<int>> groups;
            for(int v:cell) {
                vector<int> signature;
                for(auto other:p) {
                    int count=0; for(int w:other) count+=(a[v]>>w)&1;
                    signature.push_back(count);
                }
                groups[signature].push_back(v);
            }
            for(auto& entry:groups) q.push_back(entry.second);
        }
        if(q.size()==p.size()) return q;
        p=move(q);
    }
}
struct Canon { Code code=~Code(0), order=0; int aut=0; };
static Canon canonical(const Rows& a,int n) {
    map<int,vector<int>> scores;
    for(int v=0;v<n;v++) scores[pop(a[v])].push_back(v);
    Partition initial; for(auto& e:scores) initial.push_back(e.second);
    Canon answer;
    function<void(Partition)> visit=[&](Partition p) {
        p=refine(a,move(p));
        int branch=-1;
        for(int i=0;i<(int)p.size();i++) if(p[i].size()>1) {branch=i;break;}
        if(branch<0) {
            vector<int> order; for(auto cell:p) order.push_back(cell[0]);
            Code code=encode(a,order);
            if(code<answer.code) {answer.code=code;answer.order=packed(order);answer.aut=1;}
            else if(code==answer.code) answer.aut++;
            return;
        }
        for(int v:p[branch]) {
            Partition q=p;
            vector<int> rest; for(int w:p[branch]) if(w!=v) rest.push_back(w);
            q[branch]={v}; q.insert(q.begin()+branch+1,rest);
            visit(move(q));
        }
    };
    visit(initial); return answer;
}
static uint64_t ham(const Rows& a,int n) {
    vector<array<uint64_t,10>> dp(1u<<n);
    for(int v=0;v<n;v++) dp[1u<<v][v]=1;
    for(unsigned mask=1;mask<(1u<<n);mask++) for(int v=0;v<n;v++) if(dp[mask][v]) {
        unsigned avail=a[v]&~mask;
        while(avail) {int w=__builtin_ctz(avail);avail&=avail-1;dp[mask|(1u<<w)][w]+=dp[mask][v];}
    }
    uint64_t total=0; for(int v=0;v<n;v++) total+=dp.back()[v];return total;
}
struct DSU {
    vector<int> parent,size;
    DSU(int n):parent(n),size(n,1) {iota(parent.begin(),parent.end(),0);}
    int find(int x) {while(parent[x]!=x) {parent[x]=parent[parent[x]];x=parent[x];}return x;}
    bool join(int a,int b) {a=find(a);b=find(b);if(a==b)return false;if(size[a]<size[b])swap(a,b);parent[b]=a;size[a]+=size[b];return true;}
};
struct Class {
    Code code,dualOrder=0; int aut,dual=-1,bm=0; uint64_t H=0;
};
static vector<string> split(const string& line) {
    vector<string> fields;string x;istringstream s(line);while(getline(s,x,'\t'))fields.push_back(x);return fields;
}
int main(int argc,char**argv) {
    try {
        need(argc==3,"usage: census classes_n8.tsv output_prefix");
        const int n=9; string prefix=argv[2];
        ifstream input(argv[1]);need(bool(input),"missing inherited input");
        vector<Rows> parents;string line;getline(input,line);set<Code> parentCodes;
        uint64_t parentMass=0;
        while(getline(input,line)) {
            auto f=split(line);Rows a=decode(stoull(f[1]),8);
            Canon c=canonical(a,8);
            need(c.aut==stoi(f[3]),"inherited n8 automorphism mismatch");
            need(parentCodes.insert(c.code).second,"inherited classes merged by new canonicalizer");
            parents.push_back(a);parentMass+=40320/c.aut;
        }
        need(parents.size()==6880 && parentMass==(1ull<<28),"inherited n8 complete orbit mass");
        cerr<<"n8 input validated: 6880 classes, all Aut values and full labelled mass\n";
        unordered_map<Code,int> catalog;catalog.reserve(220000);
        for(int j=0;j<(int)parents.size();j++) {
            for(unsigned pat=0;pat<256;pat++) {
                Rows a=parents[j];a[8]=0;
                for(int v=0;v<8;v++) if((pat>>v)&1) a[v]|=1u<<8;else a[8]|=1u<<v;
                Canon c=canonical(a,n);
                auto result=catalog.emplace(c.code,c.aut);
                need(result.second||result.first->second==c.aut,"Aut not class invariant");
            }
            if(j%512==511) cerr<<"augmentation "<<j+1<<"/6880 classes "<<catalog.size()<<"\n";
        }
        vector<Class> classes;classes.reserve(catalog.size());
        for(auto& e:catalog) classes.push_back({e.first,0,e.second,-1,0,0});
        sort(classes.begin(),classes.end(),[](auto&a,auto&b){return a.code<b.code;});
        unordered_map<Code,int> index;index.reserve(classes.size()*2);
        for(int i=0;i<(int)classes.size();i++) index[classes[i].code]=i;
        uint64_t mass=0;int sc=0;
        for(int i=0;i<(int)classes.size();i++) {
            Class& c=classes[i];mass+=362880/c.aut;
            Rows a=decode(c.code,n);Canon rev=canonical(converse(a,n),n);
            need(index.count(rev.code),"converse missing from complete universe");
            c.dual=index.at(rev.code);c.dualOrder=rev.order;
            if(c.dual==i) sc++;
        }
        need(classes.size()==191536 && mass==(1ull<<36) && sc==2752,"n9 class/mass/SC control");
        for(int i=0;i<(int)classes.size();i++) need(classes[classes[i].dual].dual==i,"converse involution");
        cerr<<"n9 complete classes "<<classes.size()<<" SC "<<sc<<" labelled mass "<<mass<<"\n";
        auto mid=[&](int i){return min(i,classes[i].dual);};
        vector<pair<int,int>> tiles;
        for(int i=0;i<n;i++)for(int j=i+2;j<n;j++)tiles.push_back({i,j});
        vector<int> mirror(tiles.size());vector<uint32_t> orbit;
        for(int k=0;k<(int)tiles.size();k++) {
            auto e=make_pair(n-1-tiles[k].second,n-1-tiles[k].first);
            mirror[k]=find(tiles.begin(),tiles.end(),e)-tiles.begin();
            if(mirror[k]>=k)orbit.push_back((1u<<k)|(1u<<mirror[k]));
        }
        need(orbit.size()==16,"blue cube dimension");
        const uint32_t all=(1u<<tiles.size())-1;
        auto fromTiling=[&](uint32_t t) {
            Rows a{};for(int i=0;i+1<n;i++)a[i]|=1u<<(i+1);
            for(int k=0;k<(int)tiles.size();k++) {
                int i=tiles[k].first,j=tiles[k].second;
                if((t>>k)&1)a[i]|=1u<<j;else a[j]|=1u<<i;
            }
            return a;
        };
        map<pair<int,int>,int> blue;
        for(unsigned q=0;q<(1u<<orbit.size());q++) {
            uint32_t t=0;for(int k=0;k<(int)orbit.size();k++)if((q>>k)&1)t|=orbit[k];
            Canon c=canonical(fromTiling(t),n);
            int i=index.at(c.code);need(classes[i].dual==i,"blue not self-converse");classes[i].bm++;
            if(t<(t^all)) {
                int j=index.at(canonical(fromTiling(t^all),n).code);
                blue[{min(i,j),max(i,j)}]++;
            }
        }
        int pb=0,leaf=0,branch=0;vector<bool> active(classes.size(),false);
        for(int i=0;i<(int)classes.size();i++)if(mid(i)==i) {
            Class& c=classes[i];
            if(c.dual==i) {
                c.H=ham(decode(c.code,n),n);
                need(c.H%c.aut==0 && c.bm%2==1,"SC fibre parity");
                active[i]=(c.H/c.aut!=uint64_t(c.bm));
                need(c.H/c.aut>=uint64_t(c.bm),"blue exceeds all fibre");
                if(!active[i]) {pb++;leaf+=c.bm==1;branch+=c.bm==3;}
            } else active[i]=true;
        }
        need(pb==6 && leaf==5 && branch==1,"pure-blue all-order control");
        int activeCount=count(active.begin(),active.end(),true);
        DSU black(classes.size()); int components=activeCount;
        ofstream witness(prefix+"_black_tree.tsv");
        witness<<"source_code\ttarget_code\tpath_order\tpartner_canonical_order\n";
        vector<int> processed(classes.size(),0);
        uint64_t testedPaths=0,blackPaths=0;vector<int> remainingByRound;
        for(int budget: {1,4,16,64,256,1024,1000000}) {
            int largest=-1;for(int i=0;i<(int)classes.size();i++)if(active[i]&&(largest<0||black.size[black.find(i)]>black.size[black.find(largest)]))largest=i;
            for(int i=0;i<(int)classes.size() && components>1;i++) if(active[i]) {
                if(budget>1 && black.find(i)==black.find(largest)) continue;
                Rows a=decode(classes[i].code,n);vector<int> path;int seen=0;bool reached=false;
                function<void(unsigned,int)> dfs=[&](unsigned used,int last) {
                    if(seen>=budget||reached)return;
                    if(used==((1u<<n)-1)) {
                        seen++;if(seen<=processed[i])return;testedPaths++;
                        uint32_t t=0;
                        for(int k=0;k<(int)tiles.size();k++)
                            if((a[path[tiles[k].first]]>>path[tiles[k].second])&1)t|=1u<<k;
                        uint32_t reflected=0;
                        for(int k=0;k<(int)tiles.size();k++)if((t>>k)&1)reflected|=1u<<mirror[k];
                        if(t==reflected)return;
                        blackPaths++;
                        Rows partner=converse(a,n);
                        for(int k=0;k+1<n;k++) {int u=path[k],v=path[k+1];partner[u]^=1u<<v;partner[v]^=1u<<u;}
                        Canon c=canonical(partner,n);int j=index.at(c.code),merged=mid(j);
                        need(active[merged],"black line hits pure-blue node");
                        if(black.join(i,merged)) {
                            components--;
                            witness<<classes[i].code<<'\t'<<c.code<<'\t'<<packed(path)<<'\t'<<c.order<<'\n';
                        }
                        if(budget>1 && black.find(i)==black.find(largest)) reached=true;
                        return;
                    }
                    unsigned avail=last<0?((1u<<n)-1):(a[last]&~used);
                    while(avail&&!reached&&seen<budget) {
                        int v=__builtin_ctz(avail);avail&=avail-1;path.push_back(v);dfs(used|(1u<<v),v);path.pop_back();
                    }
                };
                dfs(0,-1);processed[i]=max(processed[i],seen);
            }
            remainingByRound.push_back(components);
            cerr<<"black witness budget "<<budget<<" components "<<components<<" tested paths "<<testedPaths<<"\n";
            if(components==1)break;
        }
        witness.close();
        ofstream cls(prefix+"_classes.tsv");cls<<"code\tAut\tconverse_code\tconverse_order\tblue_mass\tSC_H\n";
        for(auto& c:classes)cls<<c.code<<'\t'<<c.aut<<'\t'<<classes[c.dual].code<<'\t'<<c.dualOrder<<'\t'<<c.bm<<'\t'<<c.H<<'\n';
        ofstream bl(prefix+"_blue.tsv");bl<<"source_code\ttarget_code\tmultiplicity\n";
        for(auto&e:blue)bl<<classes[e.first.first].code<<'\t'<<classes[e.first.second].code<<'\t'<<e.second<<'\n';
        ofstream summary(prefix+"_summary.txt");
        summary<<"ordinary_classes "<<classes.size()<<"\nSC_classes "<<sc<<"\nmerged_nodes "<<(classes.size()+sc)/2
               <<"\nlabelled_orbit_mass "<<mass<<"\nblue_tilings "<<(1u<<orbit.size())<<"\nblue_lines "<<(1u<<(orbit.size()-1))
               <<"\nblue_node_pairs "<<blue.size()<<"\npure_blue "<<pb<<"\npure_blue_leaves "<<leaf<<"\npure_blue_trivalent "<<branch
               <<"\nblack_active_nodes "<<activeCount<<"\nblack_active_components_certified "<<components
               <<"\nblack_spanning_forest_edges "<<activeCount-components<<"\ntested_Hamiltonian_paths "<<testedPaths
               <<"\nblack_paths_tested "<<blackPaths<<"\nremaining_components_by_round";
        for(int x:remainingByRound)summary<<' '<<x;summary<<"\n";
        need(components==1,"black connectivity not certified: retain forest and inspect components");
        cerr<<"COMPLETE n9 universe + full blue graph + authenticated black spanning tree\n";
    }catch(const exception&e){cerr<<"ERROR "<<e.what()<<"\n";return 1;}
}

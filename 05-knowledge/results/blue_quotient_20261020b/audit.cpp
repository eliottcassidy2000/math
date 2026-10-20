// Independent replay using the inherited score-block exhaustive canonicalizer.
// This deliberately does not use producer equitable-refinement code.
// Build g++ -O3 -std=c++17 audit.cpp -o audit.exe
// Run audit.exe prefix
#define main inherited_census_main
#include "../../../04-computation/experiments/tiling_line_metagraph_census_20261020.cpp"
#undef main
#include <fstream>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <set>
#include <iostream>
static void require(bool b,const string&s){if(!b)throw runtime_error(s);}
static vector<string> fields(const string&s){vector<string>r;string x;istringstream in(s);while(getline(in,x,'\t'))r.push_back(x);return r;}
static void matrix(uint64_t code,int a[16][16]) {
    memset(a,0,16*16*sizeof(int));int bit=n*(n-1)/2-1;
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,bit--) {a[i][j]=(code>>bit)&1;a[j][i]=1-a[i][j];}
}
static vector<int> permutation(uint64_t x) {
    vector<int>p;for(int i=0;i<n;i++)p.push_back((x>>(4*i))&15);
    auto q=p;sort(q.begin(),q.end());for(int i=0;i<n;i++)require(q[i]==i,"invalid witness permutation");
    return p;
}
static uint64_t codeOf(const int a[16][16],const vector<int>&p) {
    uint64_t code=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)code=(code<<1)|a[p[i]][p[j]];
    return code;
}
struct AuditClass {uint64_t code,dual,order,H;int aut,bm;};
int main(int argc,char**argv) {
    try {
        require(argc==2,"usage audit prefix");n=9;string prefix=argv[1],line;
        ifstream input(prefix+"_classes.tsv");require(bool(input),"missing classes");getline(input,line);
        map<uint64_t,AuditClass> classes;
        while(getline(input,line)) {
            auto f=fields(line);AuditClass c{stoull(f[0]),stoull(f[2]),stoull(f[3]),stoull(f[5]),stoi(f[1]),stoi(f[4])};
            require(classes.emplace(c.code,c).second,"duplicate declared class");
        }
        unordered_map<uint64_t,uint64_t> oldToNew;oldToNew.reserve(classes.size()*2);
        uint64_t orbitMass=0;int sc=0,a[16][16],b[16][16],score[16],aut,done=0;
        for(auto&entry:classes) {
            auto c=entry.second;matrix(c.code,a);
            uint64_t old=canonForm(a,aut,score);
            require(aut==c.aut,"independent Aut mismatch");
            require(oldToNew.emplace(old,c.code).second,"independent canonical collision");
            require(362880%aut==0,"Aut does not divide9!");
            orbitMass+=362880/aut;
            for(int i=0;i<n;i++)for(int j=0;j<n;j++)b[i][j]=a[j][i];
            require(codeOf(b,permutation(c.order))==c.dual,"converse explicit isomorphism");
            require(classes.count(c.dual)&&classes.at(c.dual).dual==c.code,"converse closure");
            if(c.code==c.dual) {sc++;require(uint64_t(hamPaths(a))==c.H,"independent SC Hamiltonian count");}
            if(++done%20000==0)cerr<<"legacy canonical replay "<<done<<"/"<<classes.size()<<"\n";
        }
        require(classes.size()==191536&&orbitMass==(1ull<<36)&&sc==2752,"complete universe control");
        cerr<<"ALL191536 classes independently distinct; orbit mass2^36\n";
        vector<uint64_t> nodes;map<uint64_t,int> nodeIndex;
        int pureBlue=0,leaf=0,branch=0;
        for(auto&entry:classes)if(entry.first<=entry.second.dual) {
            nodeIndex[entry.first]=nodes.size();nodes.push_back(entry.first);
            auto c=entry.second;if(c.code==c.dual && c.H/c.aut==uint64_t(c.bm)) {pureBlue++;leaf+=c.bm==1;branch+=c.bm==3;}
        }
        auto merged=[&](uint64_t code){return min(code,classes.at(code).dual);};
        vector<int> parent(nodes.size());iota(parent.begin(),parent.end(),0);
        function<int(int)> findroot=[&](int x){if(parent[x]!=x)parent[x]=findroot(parent[x]);return parent[x];};
        ifstream tree(prefix+"_black_tree.tsv");require(bool(tree),"missing tree");getline(tree,line);
        size_t edges=0;
        while(getline(tree,line)) {
            auto f=fields(line);uint64_t source=stoull(f[0]),target=stoull(f[1]);
            require(classes.count(source)&&classes.count(target),"tree endpoint outside universe");
            matrix(source,a);vector<int>path=permutation(stoull(f[2]));
            for(int i=0;i+1<n;i++)require(a[path[i]][path[i+1]],"not a Hamiltonian path");
            bool symmetric=true;
            for(int i=0;i<n;i++)for(int j=i+2;j<n;j++)
                if(a[path[i]][path[j]]!=a[path[n-1-j]][path[n-1-i]])symmetric=false;
            require(!symmetric,"tree edge is blue");
            for(int i=0;i<n;i++)for(int j=0;j<n;j++)b[i][j]=a[j][i];
            for(int i=0;i+1<n;i++){int u=path[i],v=path[i+1];b[u][v]^=1;b[v][u]^=1;}
            require(codeOf(b,permutation(stoull(f[3])))==target,"partner explicit isomorphism");
            int u=findroot(nodeIndex.at(merged(source))),v=findroot(nodeIndex.at(merged(target)));
            require(u!=v,"witness is not a forest");parent[v]=u;edges++;
        }
        set<int> activeComponents;
        for(uint64_t code:nodes) {
            auto c=classes.at(code);
            if(c.code==c.dual && c.H/c.aut==uint64_t(c.bm))continue;
            activeComponents.insert(findroot(nodeIndex.at(code)));
        }
        require(pureBlue==6&&leaf==5&&branch==1,"PB exact check");
        require(edges==97137 && activeComponents.size()==1,"complete black tree replay");
        cerr<<"ALL97137 black edges verified by HP and explicit isomorphism\n";
        vector<pair<int,int>> cells;
        for(int i=0;i<n;i++)for(int j=i+2;j<n;j++)cells.push_back({i,j});
        vector<uint32_t> orbits;
        for(int k=0;k<(int)cells.size();k++) {
            auto e=make_pair(n-1-cells[k].second,n-1-cells[k].first);
            int mirror=find(cells.begin(),cells.end(),e)-cells.begin();
            if(k<=mirror)orbits.push_back((1u<<k)|(1u<<mirror));
        }
        vector<uint64_t> blueClass(1u<<orbits.size());map<uint64_t,int>blueMass;
        for(unsigned q=0;q<blueClass.size();q++) {
            uint32_t t=0;for(int j=0;j<(int)orbits.size();j++)if((q>>j)&1)t|=orbits[j];
            memset(a,0,sizeof(a));for(int i=0;i+1<n;i++)a[i][i+1]=1;
            for(int j=0;j<(int)cells.size();j++){auto e=cells[j];a[e.first][e.second]=(t>>j)&1;a[e.second][e.first]=1-a[e.first][e.second];}
            for(int i=0;i+1<n;i++)a[i+1][i]=0;
            uint64_t old=canonForm(a,aut,score);
            require(oldToNew.count(old),"blue class absent");
            blueClass[q]=oldToNew.at(old);blueMass[blueClass[q]]++;
            if((q+1)%8192==0)cerr<<"legacy blue replay "<<q+1<<"/"<<blueClass.size()<<"\n";
        }
        map<pair<uint64_t,uint64_t>,int> blue,expected;
        for(unsigned q=0;q<blueClass.size()/2;q++) {
            uint64_t u=blueClass[q],v=blueClass[q^(blueClass.size()-1)];
            blue[{min(u,v),max(u,v)}]++;
        }
        ifstream blueInput(prefix+"_blue.tsv");getline(blueInput,line);
        while(getline(blueInput,line)){auto f=fields(line);expected[{stoull(f[0]),stoull(f[1])}]=stoi(f[2]);}
        require(blue==expected,"independent complete blue matrix mismatch");
        for(auto&entry:classes)require(entry.second.bm==blueMass[entry.first],"independent blue fibre mismatch");
        ofstream report(prefix+"_independent_audit.txt");
        report<<"PASS\nindependent_canonicalizer inherited_score_block_full_permutations\n"
              <<"ordinary_classes "<<classes.size()<<"\nlabelled_orbit_mass "<<orbitMass
              <<"\nSC_classes "<<sc<<"\nmerged_nodes "<<nodes.size()
              <<"\nblack_tree_witnesses_replayed "<<edges<<"\nblack_active_components "<<activeComponents.size()
              <<"\nblue_tilings_recomputed "<<blueClass.size()<<"\nblue_node_pairs "<<blue.size()
              <<"\npure_blue_leaves "<<leaf<<"\npure_blue_trivalent "<<branch<<"\n";
        cerr<<"INDEPENDENT AUDIT PASSED\n";
    }catch(const exception&e){cerr<<"ERROR "<<e.what()<<"\n";return 1;}
}

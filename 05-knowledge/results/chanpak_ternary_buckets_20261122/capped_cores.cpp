// Height-capped augmentation of all strong tournaments, with fork classification.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
#include <nauty.h>

using Rows = std::vector<int>;
using Count = std::int64_t;
using Edge = std::pair<int,int>;

int low(int mask) { return __builtin_ctz((unsigned)mask); }

Rows canonical(const Rows& rows) {
    int n=rows.size();
    graph g[16]{}, canon[16]{};
    int labels[16], cells[16], orbits[16];
    for (int i=0; i<n; ++i)
        for (int bits=rows[i]; bits; bits&=bits-1) ADDONEARC(g,i,low(bits),1);
    static DEFAULTOPTIONS_GRAPH(options);
    options.digraph=TRUE;
    options.getcanon=TRUE;
    statsblk stats;
    densenauty(g,labels,cells,orbits,&options,&stats,1,n,canon);
    Rows result(n);
    for (int i=0; i<n; ++i) for (int j=0; j<n; ++j)
        if (ISELEMENT(GRAPHROW(canon,i,1),j)) result[i]|=1<<j;
    return result;
}

using Vec = std::array<Count,16>;
struct Reader {
    int n, full;
    std::vector<Vec> starts, ends;
    std::array<Vec,16> links{};
    Count height=0;

    explicit Reader(const Rows& rows): n(rows.size()),full((1<<n)-1),
                                        starts(full+1),ends(full+1) {
        for (int i=0; i<n; ++i) starts[1<<i][i]=ends[1<<i][i]=1;
        for (int set=1; set<=full; ++set) for (int bits=set; bits; bits&=bits-1) {
            int i=low(bits), rest=set^(1<<i);
            for (int t=rest&rows[i]; t; t&=t-1) starts[set][i]+=starts[rest][low(t)];
            for (int t=rest&~rows[i]; t; t&=t-1) ends[set][i]+=ends[rest][low(t)];
        }
        for (int i=0; i<n; ++i) height+=ends[full][i];
        for (int s=1; s<full; ++s) {
            int t=full^s;
            for (int bits=s; bits; bits&=bits-1) {
                int i=low(bits);
                if (!ends[s][i]) continue;
                for (int right=t; right; right&=right-1) {
                    int j=low(right);
                    links[i][j]+=ends[s][i]*starts[t][j];
                }
            }
        }
    }

    Edge forks(const Rows& rows,std::map<int,int>& counts,std::array<int,4>& fibre,
               std::ostream& report) const {
        ++fibre[0];
        std::array<Edge,3> arcs{Edge{-1,-1},Edge{-1,-1},Edge{-1,-1}};
        for (int i=0; i<n; ++i) for (int t=rows[i]; t; t&=t-1) {
            int j=low(t);
            Count delta=links[j][i]-links[i][j];
            if (delta==-2*height/3) arcs[0]={i,j};
            if (delta==-2) arcs[1]={i,j};
            if (delta== 2) arcs[2]={i,j};
        }
        if (arcs[1].first>=0) ++fibre[1];
        if (arcs[2].first>=0) ++fibre[2];
        if (arcs[1].first>=0 && arcs[2].first>=0) ++fibre[3];
        if (height%3 || arcs[0].first<0 || arcs[1].first<0 || arcs[2].first<0) return arcs[2];
        ++counts[height/3];
        report<<"center "<<n<<" parent "<<height/3<<" rows";
        for (int r:rows) report<<' '<<r;
        Rows converse=rows;
        for (int i=0;i<n;++i) converse[i]^=full^(1<<i);
        report<<" self-converse "<<(canonical(converse)==canonical(rows))<<" arcs";
        for (auto [i,j]:arcs) report<<' '<<i<<','<<j;
        report<<'\n';
        return arcs[2];
    }

    void grow(const Rows& rows,Count cap,std::map<Rows,Count>& next) const {
        // mask = destinations of arcs from the new vertex. Flip its least bit
        // on in the prior mask, updating first, last, and two-slot contacts.
        std::vector<Count> values(full+1);
        values[0]=height;
        for (int mask=1; mask<full; ++mask) {
            int v=low(mask), old=mask^(1<<v);
            Count value=values[old]+starts[full][v]-ends[full][v];
            for (int t=full^mask; t; t&=t-1) value+=links[low(t)][v];
            for (int t=old; t; t&=t-1) value-=links[v][low(t)];
            values[mask]=value;
            if (value>cap) continue;
            Rows child=rows;
            for (int t=full^mask; t; t&=t-1) child[low(t)]|=1<<n;
            child.push_back(mask);
            auto [it,fresh]=next.emplace(canonical(child),value);
            if (!fresh && it->second!=value) std::abort();
        }
    }
};

int main(int argc,char**argv) {
    const bool show_fibres=argc==4 && std::string(argv[3])=="--fibres";
    const bool save_edges=argc==6 && std::string(argv[3])=="--edges";
    if (argc!=3 && !show_fibres && !save_edges) return 2;
    int cap=std::stoi(argv[1]),max_order=std::stoi(argv[2]);
    if (max_order<3 || max_order>15 || cap<3) return 2;
    int first=save_edges?std::stoi(argv[4]):0, last=save_edges?std::stoi(argv[5]):0;
    if (save_edges && (first<4 || last<first || last>=cap || first%2 || last%2)) return 2;
    std::map<Count,std::pair<Rows,Edge>> certificates;
    std::ostream& report=save_edges?std::cerr:std::cout;
    std::map<Rows,Count> level{{{2,4,1},3}};
    for (int n=3; n<=max_order && !level.empty(); ++n) {
        std::map<Rows,Count> next;
        std::map<int,int> forks;
        std::map<Count,std::array<int,4>> fibres;
        Count minimum=cap;
        for (const auto& [rows,recorded]:level) {
            Reader dp(rows);
            if (recorded!=dp.height) std::abort();
            minimum=std::min(minimum,dp.height);
            Edge upper=dp.forks(rows,forks,fibres[recorded],report);
            Count target=recorded+1;
            if (save_edges && upper.first>=0 && target>=first && target<=last)
                certificates.try_emplace(target,rows,upper);
            if (n<max_order) dp.grow(rows,cap,next);
        }
        if (show_fibres) for (const auto& [h,types]:fibres) {
            std::cout<<"fibre "<<n<<' '<<h;
            for (auto count:types) std::cout<<' '<<count;
            std::cout<<'\n';
        }
        report<<"order "<<n<<" classes "<<level.size()<<" minimum "<<minimum<<" forks";
        for (auto [h,k]:forks) report<<' '<<h<<':'<<k;
        report<<std::endl;
        if (n<max_order && next.empty()) report<<"order "<<n+1<<" classes 0; closed"<<std::endl;
        level.swap(next);
    }
    for (const auto& [target,certificate]:certificates) {
        const auto& [rows,arc]=certificate;
        std::cout<<target<<' '<<arc.first<<' '<<arc.second;
        for (int row:rows) std::cout<<' '<<row;
        std::cout<<'\n';
    }
}

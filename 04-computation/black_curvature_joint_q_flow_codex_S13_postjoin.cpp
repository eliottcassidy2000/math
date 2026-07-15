// Exact n=8 joint (source q,target q) black-boundary flow disintegration.
//
// This is an atlas-free continuation of THM-814.  It classifies every endpoint
// as self-converse/non-self-converse by exact tournament backtracking, orients
// each non-tied mixed--pure-black complement edge toward increasing C3, and
// records the two matrices
//
//   O[i][j] = mixed(q=i) -> pure(q=j),
//   R[i][j] = pure(q=i)  -> mixed(q=j).
//
// It reports raw counts, source-conditioned rates O[i][j]/M_i and R[i][j]/P_i,
// and target-conditioned rates O[i][j]/P_j and R[i][j]/M_j.  Here M_i,P_i are
// the full black mixed/pure endpoint masses from THM-814.
//
// Tournament Analysis uses q strata as vertices.  The pairwise observable is
// the larger conditional outward-minus-reverse marginal; the switch changes
// source conditioning to target conditioning, and q=0<...<6 is the tie path.
// The tiling vertices retain marked Hamiltonian-path positions; q vertices are
// only a flow diagnostic.  No LRC gap/owner/loneliness predicate is preserved.

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

struct Tile { int a,b; };

static std::vector<Tile> tiles(int n) {
    std::vector<Tile> out;
    for(int b=1;b<=n-2;++b) for(int a=n;a>=b+2;--a) out.push_back({a,b});
    return out;
}

static int tile_pos(const std::vector<Tile>& ts,int a,int b) {
    for(int i=0;i<(int)ts.size();++i) if(ts[i].a==a&&ts[i].b==b) return i;
    assert(false); return -1;
}

static std::vector<int> reflection_map(int n,const std::vector<Tile>& ts) {
    std::vector<int> sigma(ts.size());
    for(int i=0;i<(int)ts.size();++i)
        sigma[i]=tile_pos(ts,n-ts[i].b+1,n-ts[i].a+1);
    return sigma;
}

static uint32_t reflect_mask(uint32_t mask,const std::vector<int>& sigma) {
    uint32_t out=0;
    for(int i=0;i<(int)sigma.size();++i) out|=((mask>>i)&1u)<<sigma[i];
    return out;
}

struct Tournament {
    std::array<uint8_t,8> out{};
    std::array<uint8_t,8> deg{};
};

static Tournament tournament_of(uint32_t mask,const std::vector<Tile>& ts) {
    Tournament T;
    for(int i=0;i<7;++i) T.out[i]|=uint8_t(1u<<(i+1));
    for(int bit=0;bit<(int)ts.size();++bit) {
        int i=8-ts[bit].a,j=8-ts[bit].b;
        if((mask>>bit)&1u) T.out[j]|=uint8_t(1u<<i);
        else T.out[i]|=uint8_t(1u<<j);
    }
    for(int i=0;i<8;++i) T.deg[i]=uint8_t(__builtin_popcount(T.out[i]));
    return T;
}

static Tournament converse(const Tournament& T) {
    Tournament C;
    for(int i=0;i<8;++i) {
        C.out[i]=uint8_t((~T.out[i])&(0xffu^(1u<<i)));
        C.deg[i]=uint8_t(7-T.deg[i]);
    }
    return C;
}

static bool isomorphic(const Tournament& A,const Tournament& B,uint64_t& states) {
    std::array<int,8> ha{},hb{};
    for(int i=0;i<8;++i){++ha[A.deg[i]];++hb[B.deg[i]];}
    if(ha!=hb) return false;
    std::array<int,8> order{}; std::iota(order.begin(),order.end(),0);
    std::sort(order.begin(),order.end(),[&](int u,int v){
        if(ha[A.deg[u]]!=ha[A.deg[v]]) return ha[A.deg[u]]<ha[A.deg[v]];
        if(A.deg[u]!=A.deg[v]) return A.deg[u]<A.deg[v];
        return u<v;
    });
    std::array<int,8> image; image.fill(-1);
    auto dfs=[&](auto&& self,int pos,uint16_t used)->bool{
        ++states;
        if(pos==8) return true;
        int u=order[pos];
        for(int v=0;v<8;++v) {
            if(((used>>v)&1u)||A.deg[u]!=B.deg[v]) continue;
            bool ok=true;
            for(int k=0;k<pos&&ok;++k) {
                int x=order[k],y=image[x];
                ok=(((A.out[u]>>x)&1u)==((B.out[v]>>y)&1u));
            }
            if(!ok) continue;
            image[u]=v;
            if(self(self,pos+1,uint16_t(used|(1u<<v)))) return true;
            image[u]=-1;
        }
        return false;
    };
    return dfs(dfs,0,0);
}

static bool self_converse(const Tournament& T,uint64_t& states,bool& score_symmetric) {
    std::array<int,8> h{};
    for(int i=0;i<8;++i) ++h[T.deg[i]];
    score_symmetric=true;
    for(int d=0;d<8;++d) score_symmetric&=h[d]==h[7-d];
    return score_symmetric&&isomorphic(T,converse(T),states);
}

static int tile_bit(uint32_t mask,const std::vector<Tile>& ts,int a,int b) {
    return (mask>>tile_pos(ts,a,b))&1u;
}

// Returns (q(mask),q(complement),C3(complement)-C3(mask)).
static std::array<int,3> curvature(uint32_t line,uint32_t full,const std::vector<Tile>& ts) {
    uint32_t mask=(line&1u)?(line^full):line;
    int B=0,T=0,q0=0;
    for(int x=3;x<8;++x) B+=tile_bit(mask,ts,x,1);
    for(int y=2;y<7;++y) T+=tile_bit(mask,ts,8,y);
    for(int x=3;x<7;++x) q0+=tile_bit(mask,ts,x,1)*tile_bit(mask,ts,8,x);
    int q1=q0+tile_bit(mask,ts,7,1)+tile_bit(mask,ts,8,2);
    return {q0,q1,6-B-T};
}

static std::string fraction(uint64_t a,uint64_t b) {
    assert(b>0);
    uint64_t g=std::gcd(a,b);
    std::ostringstream out; out<<(a/g)<<"/"<<(b/g); return out.str();
}

struct Fingerprint {
    std::map<int,int> scores;
    int c3=0;
    std::vector<int> scc;
    uint64_t hpaths=0;
    std::array<std::array<int,7>,7> adj{};
};

template<class Better> static Fingerprint fingerprint(Better better) {
    Fingerprint f;
    for(int i=0;i<7;++i) for(int j=i+1;j<7;++j)
        if(better(i,j)) f.adj[i][j]=1; else f.adj[j][i]=1;
    for(int i=0;i<7;++i) {
        int score=0; for(int j=0;j<7;++j) score+=f.adj[i][j]; ++f.scores[score];
    }
    for(int i=0;i<7;++i)for(int j=i+1;j<7;++j)for(int k=j+1;k<7;++k)
        f.c3+=f.adj[i][j]*f.adj[j][k]*f.adj[k][i]
             +f.adj[j][i]*f.adj[k][j]*f.adj[i][k];
    std::array<int,7> seen{}; std::vector<int> order;
    auto dfs1=[&](auto&& self,int u)->void{
        seen[u]=1; for(int v=0;v<7;++v) if(f.adj[u][v]&&!seen[v]) self(self,v);
        order.push_back(u);
    };
    for(int u=0;u<7;++u) if(!seen[u]) dfs1(dfs1,u);
    seen.fill(0);
    auto dfs2=[&](auto&& self,int u)->int{
        seen[u]=1; int z=1;
        for(int v=0;v<7;++v) if(f.adj[v][u]&&!seen[v]) z+=self(self,v);
        return z;
    };
    for(auto it=order.rbegin();it!=order.rend();++it) if(!seen[*it]) f.scc.push_back(dfs2(dfs2,*it));
    std::sort(f.scc.begin(),f.scc.end(),std::greater<int>());
    std::array<std::array<uint64_t,7>,128> dp{};
    for(int v=0;v<7;++v) dp[1<<v][v]=1;
    for(int mask=1;mask<128;++mask) for(int v=0;v<7;++v) if(dp[mask][v])
        for(int w=0;w<7;++w) if(!(mask&(1<<w))&&f.adj[v][w]) dp[mask|(1<<w)][w]+=dp[mask][v];
    for(int v=0;v<7;++v) f.hpaths+=dp[127][v];
    return f;
}

static void print_fp(const char* name,const Fingerprint& f) {
    std::cout<<"  "<<name<<" scores={"; bool first=true;
    for(auto [s,c]:f.scores){if(!first)std::cout<<",";first=false;std::cout<<s<<":"<<c;}
    std::cout<<"} C3="<<f.c3<<" SCC=(";
    for(size_t i=0;i<f.scc.size();++i) std::cout<<(i?",":"")<<f.scc[i];
    std::cout<<") Hpaths="<<f.hpaths<<"\n";
}

int main() {
    const int n=8;
    auto ts=tiles(n); auto sigma=reflection_map(n,ts);
    const uint32_t masks=1u<<21,full=masks-1,lines=1u<<20;
    std::vector<uint8_t> sc(masks);
    uint64_t iso_states=0,score_candidates=0,sc_count=0;
    for(uint32_t mask=0;mask<masks;++mask) {
        bool score=false; sc[mask]=self_converse(tournament_of(mask,ts),iso_states,score);
        score_candidates+=score; sc_count+=sc[mask];
    }
    assert(score_candidates==744678&&sc_count==58712);

    std::array<uint64_t,7> mixed{},pure{},out_source{},rev_source{},out_target{},rev_target{};
    std::array<std::array<uint64_t,7>,7> outward{},reverse{};
    uint64_t black_lines=0,tied=0;
    for(uint32_t line=0;line<lines;++line) {
        uint32_t mask=(line&1u)?(line^full):line,mate=mask^full;
        if(reflect_mask(mask,sigma)==mask) continue;
        ++black_lines;
        auto c=curvature(line,full,ts); int q0=c[0],q1=c[1],flux=c[2];
        bool m0=sc[mask],m1=sc[mate];
        ++(m0?mixed[q0]:pure[q0]); ++(m1?mixed[q1]:pure[q1]);
        if(m0==m1) continue;
        if(flux==0){++tied;continue;}
        bool source_mixed=flux>0?m0:m1;
        int qs=flux>0?q0:q1,qt=flux>0?q1:q0;
        if(source_mixed){++outward[qs][qt];++out_source[qs];++out_target[qt];}
        else {++reverse[qs][qt];++rev_source[qs];++rev_target[qt];}
    }
    assert(black_lines==1046528&&tied==12584);
    const std::array<uint64_t,7> expected_out{5296,5166,1326,1820,2158,180,172};
    const std::array<uint64_t,7> expected_rev{8620,9086,1356,2496,2422,368,22};
    const std::array<uint64_t,7> expected_mixed{8282,16240,19160,7844,2734,180,176};
    const std::array<uint64_t,7> expected_pure{404710,702608,559720,274780,81682,14156,784};
    assert(out_source==expected_out&&rev_source==expected_rev&&mixed==expected_mixed&&pure==expected_pure);
    assert(std::accumulate(out_target.begin(),out_target.end(),uint64_t(0))==16118);
    assert(std::accumulate(rev_target.begin(),rev_target.end(),uint64_t(0))==24370);
    int support_cells=0,source_strict_cells=0,target_reverse_cells=0;
    for(int i=0;i<7;++i) for(int j=0;j<7;++j) {
        assert((outward[i][j]>0)==(reverse[i][j]>0));
        assert((outward[i][j]&1u)==0&&(reverse[i][j]&1u)==0);
        if(!outward[i][j]) continue;
        ++support_cells;
        // Equivalently, the raw edge-count ratio lies strictly between the
        // source mixed/pure mass odds and the target pure/mixed mass odds.
        source_strict_cells+=(__int128)outward[i][j]*pure[i]>(__int128)reverse[i][j]*mixed[i];
        target_reverse_cells+=(__int128)outward[i][j]*mixed[j]<(__int128)reverse[i][j]*pure[j];
        assert(std::abs(i-j)<=2);
    }
    assert(support_cells==17&&source_strict_cells==17&&target_reverse_cells==17);

    auto source_better=[&](int i,int j){
        __int128 ni=(__int128)out_source[i]*pure[i]-(__int128)rev_source[i]*mixed[i];
        __int128 di=(__int128)mixed[i]*pure[i];
        __int128 nj=(__int128)out_source[j]*pure[j]-(__int128)rev_source[j]*mixed[j];
        __int128 dj=(__int128)mixed[j]*pure[j];
        return ni*dj!=nj*di?ni*dj>nj*di:i<j;
    };
    auto target_better=[&](int i,int j){
        __int128 ni=(__int128)out_target[i]*mixed[i]-(__int128)rev_target[i]*pure[i];
        __int128 di=(__int128)mixed[i]*pure[i];
        __int128 nj=(__int128)out_target[j]*mixed[j]-(__int128)rev_target[j]*pure[j];
        __int128 dj=(__int128)mixed[j]*pure[j];
        return ni*dj!=nj*di?ni*dj>nj*di:i<j;
    };
    auto sfp=fingerprint(source_better),tfp=fingerprint(target_better);
    int flips=0;for(int i=0;i<7;++i)for(int j=i+1;j<7;++j)flips+=sfp.adj[i][j]!=tfp.adj[i][j];

    std::cout<<"N=8 JOINT SOURCE-q / TARGET-q BLACK BOUNDARY FLOW\n"<<std::string(78,'=')<<"\n";
    std::cout<<"black lines="<<black_lines<<" non-tied outward/reverse=16118/24370 tied="<<tied<<"\n";
    std::cout<<"self-converse classifier score-candidates/SC/states="<<score_candidates<<"/"<<sc_count<<"/"<<iso_states<<"\n\n";
    std::cout<<"RAW MATRICES (rows=source q, columns=target q)\n";
    for(const auto* label:{"O mixed->pure","R pure->mixed"}) {
        bool is_out=label[0]=='O'; std::cout<<label<<"\n";
        for(int i=0;i<7;++i){std::cout<<"  q"<<i<<":";for(int j=0;j<7;++j)std::cout<<" "<<(is_out?outward[i][j]:reverse[i][j]);std::cout<<"\n";}
    }
    std::cout<<"\nSOURCE MARGINAL (THM-814 replay)\n";
    for(int q=0;q<7;++q)
        std::cout<<"  q="<<q<<" O/M="<<out_source[q]<<"/"<<mixed[q]<<"="<<fraction(out_source[q],mixed[q])
                 <<" R/P="<<rev_source[q]<<"/"<<pure[q]<<"="<<fraction(rev_source[q],pure[q])<<"\n";
    std::cout<<"\nTARGET MARGINAL\n";
    int target_outward_strict=0;
    for(int q=0;q<7;++q) {
        __int128 lhs=(__int128)out_target[q]*mixed[q],rhs=(__int128)rev_target[q]*pure[q];
        target_outward_strict+=lhs>rhs;
        std::cout<<"  q="<<q<<" O/P="<<out_target[q]<<"/"<<pure[q]<<"="<<fraction(out_target[q],pure[q])
                 <<" R/M="<<rev_target[q]<<"/"<<mixed[q]<<"="<<fraction(rev_target[q],mixed[q])
                 <<" comparison="<<(lhs>rhs?"outward":lhs<rhs?"reverse":"tie")<<"\n";
    }
    std::cout<<"  strict outward target strata="<<target_outward_strict<<"/7\n\n";
    std::cout<<"NONZERO JOINT CELLS WITH BOTH CONDITIONINGS\n";
    for(int i=0;i<7;++i)for(int j=0;j<7;++j)if(outward[i][j]||reverse[i][j])
        std::cout<<"  "<<i<<"->"<<j<<" O="<<outward[i][j]<<" R="<<reverse[i][j]
                 <<" source(O/M,R/P)=("<<fraction(outward[i][j],mixed[i])<<","<<fraction(reverse[i][j],pure[i])<<")"
                 <<" target(O/P,R/M)=("<<fraction(outward[i][j],pure[j])<<","<<fraction(reverse[i][j],mixed[j])<<")\n";
    std::cout<<"  common support cells="<<support_cells
             <<" source-strict="<<source_strict_cells
             <<" target-reverse="<<target_reverse_cells<<"; every |target q-source q|<=2\n";
    std::cout<<"\nTOURNAMENT ANALYSIS (q strata as vertices)\n";
    std::cout<<"  observable=conditional outward-minus-reverse marginal\n";
    std::cout<<"  switches=source conditioning / target conditioning; tie path=q0<...<q6; flips="<<flips<<"\n";
    print_fp("source",sfp); print_fp("target",tfp);
    std::cout<<"\nPreservation: marked endpoint q, category, C3 direction, literal multiplicity.\n";
    std::cout<<"Destroyed: node identity, layer position, owner/metric/LRC predicate, continuation.\n";
}

// Atlas-free exact n=8 black-edge boundary for THM-814.
//
// Two audits share the literal complement-line stream.  First, the
// reflection-orbit (B2,B3) key is sorted and tournament isomorphism is run
// only inside its collision cells.  Second, endpoint self-converseness gives
// mixed versus pure-black directly and disintegrates boundary flow by q.
// Tournament Analysis uses q strata as vertices, raw versus source-normalized
// outward bias as switches, and increasing q as the tie Hamiltonian path.

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

struct Tile { int a, b; };

static std::vector<Tile> tiles(int n) {
    std::vector<Tile> out;
    for (int b = 1; b <= n - 2; ++b)
        for (int a = n; a >= b + 2; --a) out.push_back({a,b});
    return out;
}

static int tile_pos(const std::vector<Tile>& ts, int a, int b) {
    for (int i = 0; i < (int)ts.size(); ++i)
        if (ts[i].a == a && ts[i].b == b) return i;
    throw std::runtime_error("tile not found");
}

static std::vector<int> reflection_map(int n, const std::vector<Tile>& ts) {
    std::vector<int> sigma(ts.size());
    for (int i = 0; i < (int)ts.size(); ++i)
        sigma[i] = tile_pos(ts, n-ts[i].b+1, n-ts[i].a+1);
    return sigma;
}

static uint32_t reflect_mask(uint32_t mask, const std::vector<int>& sigma) {
    uint32_t out = 0;
    for (int i = 0; i < (int)sigma.size(); ++i)
        out |= ((mask >> i) & 1u) << sigma[i];
    return out;
}

static bool symmetric(uint32_t mask, const std::vector<int>& sigma) {
    return reflect_mask(mask, sigma) == mask;
}

static int composition_rank(const std::array<int,4>& c, int total) {
    int rank = 0;
    for (int a = 0; a <= total; ++a)
        for (int b = 0; b <= total-a; ++b)
            for (int d = 0; d <= total-a-b; ++d) {
                std::array<int,4> x{a,b,d,total-a-b-d};
                if (x == c) return rank;
                ++rank;
            }
    throw std::runtime_error("invalid weak composition");
}

static uint32_t s2_rank(uint32_t mask, const std::vector<Tile>& ts,
                        const std::vector<int>& sigma) {
    const std::array<int,5> layer_size{1,1,2,2,3};
    const std::array<int,6> radix{4,4,10,10,20,4};
    std::array<std::array<int,4>,5> comp{};
    int fixed_ones = 0;
    for (int bit = 0; bit < (int)ts.size(); ++bit) {
        int tau = ts[bit].a + ts[bit].b - 1;
        if (tau < 8) {
            int x = (mask >> bit) & 1u;
            int y = (mask >> sigma[bit]) & 1u;
            ++comp[tau-3][2*x+y];
        } else if (tau == 8) {
            assert(sigma[bit] == bit);
            fixed_ones += (mask >> bit) & 1u;
        }
    }
    uint32_t rank = 0;
    for (int i = 0; i < 5; ++i)
        rank = rank * radix[i] + composition_rank(comp[i], layer_size[i]);
    rank = rank * radix[5] + fixed_ones;
    assert(rank < 128000);
    return rank;
}

static int s3_stratum(const Tile& tile, int n) {
    bool A = tile.a < n;
    bool B = tile.a - tile.b >= 3;
    bool C = tile.b >= 2;
    if (A && !B && !C) return 0;
    if (!A && B && !C) return 1;
    if (!A && !B && C) return 2;
    if (A && B && !C) return 3;
    if (A && !B && C) return 4;
    if (!A && B && C) return 5;
    if (A && B && C) return 6;
    throw std::runtime_error("empty B3 membership");
}

static uint32_t s3_rank(uint32_t mask, const std::vector<Tile>& ts, int n) {
    const std::array<int,7> radix{2,2,2,5,5,5,7};
    std::array<int,7> counts{};
    for (int bit = 0; bit < (int)ts.size(); ++bit)
        counts[s3_stratum(ts[bit],n)] += (mask >> bit) & 1u;
    uint32_t rank = 0;
    for (int i = 0; i < 7; ++i) rank = rank * radix[i] + counts[i];
    assert(rank < 7000);
    return rank;
}

static uint32_t literal_b23(uint32_t mask, const std::vector<Tile>& ts,
                            const std::vector<int>& sigma) {
    return s2_rank(mask,ts,sigma) * 7000u + s3_rank(mask,ts,8);
}

struct Tournament {
    std::array<uint8_t,8> out{};
    std::array<uint8_t,8> deg{};
};

static Tournament tournament_of(uint32_t mask, const std::vector<Tile>& ts) {
    Tournament T;
    for (int i = 0; i < 7; ++i) T.out[i] |= uint8_t(1u << (i+1));
    for (int bit = 0; bit < (int)ts.size(); ++bit) {
        int i = 8 - ts[bit].a, j = 8 - ts[bit].b;
        assert(i < j && j >= i+2);
        if ((mask >> bit) & 1u) T.out[j] |= uint8_t(1u << i);
        else T.out[i] |= uint8_t(1u << j);
    }
    for (int i = 0; i < 8; ++i) T.deg[i] = uint8_t(__builtin_popcount(T.out[i]));
    return T;
}

static Tournament converse(const Tournament& T) {
    Tournament C;
    for (int i = 0; i < 8; ++i) {
        C.out[i] = uint8_t((~T.out[i]) & (0xffu ^ (1u << i)));
        C.deg[i] = uint8_t(7 - T.deg[i]);
    }
    return C;
}

struct IsoCounter { uint64_t calls=0, states=0; };

static bool isomorphic(const Tournament& A, const Tournament& B, IsoCounter& counter) {
    ++counter.calls;
    std::array<int,8> ha{}, hb{};
    for (int i = 0; i < 8; ++i) { ++ha[A.deg[i]]; ++hb[B.deg[i]]; }
    if (ha != hb) return false;
    std::array<int,8> order{};
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int u, int v) {
        if (ha[A.deg[u]] != ha[A.deg[v]]) return ha[A.deg[u]] < ha[A.deg[v]];
        if (A.deg[u] != A.deg[v]) return A.deg[u] < A.deg[v];
        return u < v;
    });
    std::array<int,8> image;
    image.fill(-1);
    uint16_t used = 0;
    auto dfs = [&](auto&& self, int pos, uint16_t used_now) -> bool {
        ++counter.states;
        if (pos == 8) return true;
        int u = order[pos];
        for (int v = 0; v < 8; ++v) {
            if ((used_now >> v) & 1u || A.deg[u] != B.deg[v]) continue;
            bool ok = true;
            for (int k = 0; k < pos && ok; ++k) {
                int x = order[k], y = image[x];
                bool au = (A.out[u] >> x) & 1u;
                bool bv = (B.out[v] >> y) & 1u;
                ok = au == bv;
            }
            if (!ok) continue;
            image[u] = v;
            if (self(self, pos+1, uint16_t(used_now | (1u << v)))) return true;
            image[u] = -1;
        }
        return false;
    };
    return dfs(dfs,0,used);
}

static bool self_converse(const Tournament& T, IsoCounter& counter,
                          bool& score_symmetric) {
    std::array<int,8> h{};
    for (int i = 0; i < 8; ++i) ++h[T.deg[i]];
    score_symmetric = true;
    for (int d = 0; d < 8; ++d) score_symmetric &= h[d] == h[7-d];
    if (!score_symmetric) return false;
    Tournament C = converse(T);
    return isomorphic(T,C,counter);
}

static bool same_unordered_node_pair(uint32_t line_a, uint32_t line_b,
                                     uint32_t full, const std::vector<Tile>& ts,
                                     IsoCounter& counter) {
    uint32_t a = (line_a & 1u) ? (line_a ^ full) : line_a;
    uint32_t b = (line_b & 1u) ? (line_b ^ full) : line_b;
    Tournament A = tournament_of(a,ts), Ac = tournament_of(a^full,ts);
    Tournament B = tournament_of(b,ts), Bc = tournament_of(b^full,ts);
    // Each row already represents the full staircase-reflection orbit, so
    // converse choices are present geometrically.  Compare the two endpoint
    // tournaments up to ordinary isomorphism, allowing the endpoints to swap.
    return (isomorphic(A,B,counter) && isomorphic(Ac,Bc,counter)) ||
           (isomorphic(A,Bc,counter) && isomorphic(Ac,B,counter));
}

static int tile_bit(uint32_t mask, const std::vector<Tile>& ts, int a, int b) {
    return (mask >> tile_pos(ts,a,b)) & 1u;
}

static std::array<int,4> curvature(uint32_t line, uint32_t full,
                                   const std::vector<Tile>& ts) {
    uint32_t mask = (line & 1u) ? (line ^ full) : line;
    int B=0,T=0,q0=0;
    for (int x = 3; x < 8; ++x) B += tile_bit(mask,ts,x,1);
    for (int y = 2; y < 7; ++y) T += tile_bit(mask,ts,8,y);
    for (int x = 3; x < 7; ++x)
        q0 += tile_bit(mask,ts,x,1) * tile_bit(mask,ts,8,x);
    int q1 = q0 + tile_bit(mask,ts,7,1) + tile_bit(mask,ts,8,2);
    return {q0,q1,6-B-T,B-T};
}

static int fixed_moment(uint32_t line, uint32_t full, const std::vector<Tile>& ts) {
    uint32_t mask = (line & 1u) ? (line ^ full) : line;
    int moment=0, pos=0;
    for (int bit = 0; bit < (int)ts.size(); ++bit)
        if (ts[bit].a + ts[bit].b - 1 == 8)
            moment += pos++ * int((mask >> bit) & 1u);
    return moment;
}

struct Row { uint64_t key; uint32_t line; };

struct FlowRow {
    uint64_t out=0, rev=0, mixed=0, black=0;
};

struct TournamentFingerprint {
    std::map<int,int> score_hist;
    int c3=0;
    std::vector<int> scc;
    uint64_t hpaths=0;
    std::array<std::array<int,7>,7> adj{};
};

template<class Better>
static TournamentFingerprint fingerprint(Better better) {
    TournamentFingerprint f;
    for (int i=0;i<7;++i) for (int j=i+1;j<7;++j)
        if (better(i,j)) f.adj[i][j]=1; else f.adj[j][i]=1;
    for (int i=0;i<7;++i) {
        int s=0; for(int j=0;j<7;++j)s+=f.adj[i][j]; ++f.score_hist[s];
    }
    for(int i=0;i<7;++i)for(int j=i+1;j<7;++j)for(int k=j+1;k<7;++k) {
        int cyc=f.adj[i][j]*f.adj[j][k]*f.adj[k][i]+f.adj[j][i]*f.adj[k][j]*f.adj[i][k];
        f.c3+=cyc;
    }
    // Kosaraju SCC sizes.
    std::array<int,7> seen{}; std::vector<int> order;
    auto dfs1=[&](auto&& self,int u)->void{seen[u]=1;for(int v=0;v<7;++v)if(f.adj[u][v]&&!seen[v])self(self,v);order.push_back(u);};
    for(int u=0;u<7;++u)if(!seen[u])dfs1(dfs1,u);
    seen.fill(0);
    auto dfs2=[&](auto&& self,int u)->int{seen[u]=1;int z=1;for(int v=0;v<7;++v)if(f.adj[v][u]&&!seen[v])z+=self(self,v);return z;};
    for(auto it=order.rbegin();it!=order.rend();++it)if(!seen[*it])f.scc.push_back(dfs2(dfs2,*it));
    std::sort(f.scc.begin(),f.scc.end(),std::greater<int>());
    // Hamiltonian path subset DP.
    std::array<std::array<uint64_t,7>,128> dp{};
    for(int v=0;v<7;++v)dp[1<<v][v]=1;
    for(int mask=1;mask<128;++mask)for(int v=0;v<7;++v)if(dp[mask][v])
        for(int w=0;w<7;++w)if(!(mask&(1<<w))&&f.adj[v][w])dp[mask|(1<<w)][w]+=dp[mask][v];
    for(int v=0;v<7;++v)f.hpaths+=dp[127][v];
    return f;
}

static void print_fingerprint(const char* name, const TournamentFingerprint& f) {
    std::cout << "  " << name << " score_hist={";
    bool first=true; for(auto [s,c]:f.score_hist){if(!first)std::cout<<",";first=false;std::cout<<s<<":"<<c;}
    std::cout << "} C3=" << f.c3 << " SCC=(";
    for(size_t i=0;i<f.scc.size();++i)std::cout<<(i?",":"")<<f.scc[i];
    std::cout << ") Hpaths=" << f.hpaths << "\n";
}

int main() {
    const int n=8;
    auto ts=tiles(n); assert(ts.size()==21);
    auto sigma=reflection_map(n,ts);
    const uint32_t full=(1u<<21)-1, line_count=1u<<20, mask_count=1u<<21;

    // Direct self-converse classifier for every endpoint; no class atlas.
    std::vector<uint8_t> sc(mask_count);
    IsoCounter sc_counter;
    uint64_t score_symmetric_candidates=0, sc_tilings=0, blue_tilings=0, black_mixed=0;
    for(uint32_t mask=0;mask<mask_count;++mask) {
        bool score_symmetric=false;
        bool value=self_converse(tournament_of(mask,ts),sc_counter,score_symmetric);
        score_symmetric_candidates += score_symmetric;
        sc[mask]=value;
        sc_tilings += value;
        bool blue=symmetric(mask,sigma);
        blue_tilings += blue;
        black_mixed += value && !blue;
    }
    assert(score_symmetric_candidates==744678);
    assert(sc_tilings==58712 && blue_tilings==4096 && black_mixed==54616);

    // Stream black reflection-orbit rows and source-q flow.
    std::vector<Row> rows; rows.reserve(523264);
    std::array<FlowRow,7> flow{};
    std::array<uint64_t,7> black_q{};
    uint64_t black_lines=0, cross_phase=0, tied_cross_phase=0;
    for(uint32_t line=0;line<line_count;++line) {
        uint32_t mask=(line&1u)?(line^full):line;
        if(symmetric(mask,sigma)) continue;
        ++black_lines;
        uint32_t mate=mask^full;
        auto cur=curvature(line,full,ts);
        int q0=cur[0],q1=cur[1],flux=cur[2];
        ++black_q[q0]; ++black_q[q1];
        bool mixed0=sc[mask], mixed1=sc[mate];
        if(mixed0) ++flow[q0].mixed; else ++flow[q0].black;
        if(mixed1) ++flow[q1].mixed; else ++flow[q1].black;
        if(mixed0 != mixed1) {
            ++cross_phase;
            if(flux==0) ++tied_cross_phase;
            else {
                // flux=C3(kappa mask)-C3(mask), so positive flux makes the
                // apex-zero endpoint the increasing-C3 source.
                bool source_mixed = flux>0 ? mixed0 : mixed1;
                int source_q = flux>0 ? q0 : q1;
                if(source_mixed) ++flow[source_q].out; else ++flow[source_q].rev;
            }
        }

        uint32_t reflected=reflect_mask(mask,sigma);
        uint32_t reflected_line=std::min(reflected,reflected^full);
        if(line>reflected_line) continue;
        uint32_t lit=literal_b23(mask,ts,sigma);
        uint32_t rlit=literal_b23(reflected,ts,sigma);
        uint32_t lo=std::min(lit,rlit), hi=std::max(lit,rlit);
        rows.push_back({(uint64_t(lo)<<32)|hi,line});
    }
    assert(black_lines==1046528 && rows.size()==523264);
    assert(cross_phase==53072 && tied_cross_phase==12584);
    const std::array<uint64_t,7> predicted_black_q{412992,718848,578880,282624,84416,14336,960};
    assert(black_q==predicted_black_q);
    const std::array<FlowRow,7> expected{{
        {5296,8620,8282,404710}, {5166,9086,16240,702608},
        {1326,1356,19160,559720}, {1820,2496,7844,274780},
        {2158,2422,2734,81682}, {180,368,180,14156}, {172,22,176,784}
    }};
    for(int q=0;q<7;++q) {
        assert(flow[q].out==expected[q].out && flow[q].rev==expected[q].rev);
        assert(flow[q].mixed==expected[q].mixed && flow[q].black==expected[q].black);
        __int128 lhs=(__int128)flow[q].out*flow[q].black;
        __int128 rhs=(__int128)flow[q].rev*flow[q].mixed;
        assert(lhs>rhs);
    }
    assert(std::accumulate(flow.begin(),flow.end(),uint64_t(0),[](uint64_t s,const FlowRow&r){return s+r.out;})==16118);
    assert(std::accumulate(flow.begin(),flow.end(),uint64_t(0),[](uint64_t s,const FlowRow&r){return s+r.rev;})==24370);

    std::sort(rows.begin(),rows.end(),[](const Row&a,const Row&b){return a.key!=b.key?a.key<b.key:a.line<b.line;});
    uint64_t key_cells=0,key_collision_cells=0,key_excess=0,candidate_pairs=0,max_mult=0;
    uint64_t augmented_cells=0,augmented_collision_cells=0,augmented_excess=0;
    std::vector<std::tuple<uint32_t,uint32_t,uint64_t>> witnesses;
    IsoCounter pair_counter;
    for(size_t i=0;i<rows.size();) {
        size_t j=i+1; while(j<rows.size()&&rows[j].key==rows[i].key)++j;
        uint64_t mult=j-i; ++key_cells; key_collision_cells+=mult>1; key_excess+=mult-1;
        candidate_pairs+=mult*(mult-1)/2; max_mult=std::max(max_mult,mult);
        std::vector<std::vector<size_t>> groups;
        for(size_t k=i;k<j;++k) {
            bool placed=false;
            for(auto& group:groups) if(same_unordered_node_pair(rows[k].line,rows[group[0]].line,full,ts,pair_counter)) {
                group.push_back(k); placed=true; break;
            }
            if(!placed) groups.push_back({k});
        }
        augmented_cells+=groups.size();
        for(auto& group:groups) if(group.size()>1) {
            ++augmented_collision_cells; augmented_excess+=group.size()-1;
            assert(group.size()==2);
            witnesses.push_back({rows[group[0]].line,rows[group[1]].line,rows[i].key});
        }
        i=j;
    }
    assert(key_cells==331500 && key_collision_cells==144500 && key_excess==191764);
    assert(candidate_pairs==275368 && max_mult==10);
    assert(augmented_cells==523248 && augmented_collision_cells==16 && augmented_excess==16);
    assert(witnesses.size()==16);
    std::sort(witnesses.begin(),witnesses.end());
    std::array<int,4> epsilon_groups{};
    IsoCounter witness_loop_counter;
    for(auto [a,b,key]:witnesses) {
        assert((a^b)==0x02080u);
        auto ca=curvature(a,full,ts), cb=curvature(b,full,ts);
        assert(ca[0]==0&&ca[1]==2&&ca[2]==0&&cb[0]==0&&cb[1]==2&&cb[2]==0);
        assert(fixed_moment(a,full,ts)!=fixed_moment(b,full,ts));
        for(uint32_t line:{a,b}) {
            uint32_t mask=(line&1u)?(line^full):line;
            assert(!sc[mask] && !sc[mask^full]);
            Tournament A=tournament_of(mask,ts), B=tournament_of(mask^full,ts);
            Tournament Bc=converse(B);
            assert(!isomorphic(A,B,witness_loop_counter));
            assert(!isomorphic(A,Bc,witness_loop_counter));
        }
        int slot=(ca[3]+4)/2; assert(slot>=0&&slot<4); ++epsilon_groups[slot];
    }
    const std::array<int,4> expected_epsilon_groups{4,4,4,4};
    assert(epsilon_groups==expected_epsilon_groups);

    auto raw_better=[&](int i,int j){int64_t a=int64_t(flow[i].out)-int64_t(flow[i].rev),b=int64_t(flow[j].out)-int64_t(flow[j].rev);return a!=b?a>b:i<j;};
    auto normalized_better=[&](int i,int j){
        __int128 ni=(__int128)flow[i].out*flow[i].black-(__int128)flow[i].rev*flow[i].mixed;
        __int128 di=(__int128)flow[i].mixed*flow[i].black;
        __int128 nj=(__int128)flow[j].out*flow[j].black-(__int128)flow[j].rev*flow[j].mixed;
        __int128 dj=(__int128)flow[j].mixed*flow[j].black;
        return ni*dj!=nj*di?ni*dj>nj*di:i<j;
    };
    auto raw_fp=fingerprint(raw_better), norm_fp=fingerprint(normalized_better);
    int gauge_flips=0;for(int i=0;i<7;++i)for(int j=i+1;j<7;++j)gauge_flips+=raw_fp.adj[i][j]!=norm_fp.adj[i][j];

    std::cout << "THM-814 N=8 BLACK CURVATURE/ORBIT-CODEC BOUNDARY\n" << std::string(82,'=') << "\n";
    std::cout << "lines="<<line_count<<" black="<<black_lines<<" black reflection orbits="<<rows.size()<<"\n";
    std::cout << "atlas-free SC: score-symmetric candidates="<<score_symmetric_candidates
              <<" SC tilings="<<sc_tilings<<" blue="<<blue_tilings
              <<" black-mixed endpoints="<<black_mixed<<"\n\n";
    std::cout << "ORBIT-SYMMETRIZED B2+B3 CODEC\n";
    std::cout << "  B23 cells="<<key_cells<<" collision_cells="<<key_collision_cells
              <<" excess="<<key_excess<<" candidate_pairs="<<candidate_pairs<<" max="<<max_mult<<"\n";
    std::cout << "  + unordered merged node pair: cells="<<augmented_cells
              <<" collision_cells="<<augmented_collision_cells<<" excess="<<augmented_excess<<"\n";
    std::cout << "  iso calls/states="<<pair_counter.calls<<"/"<<pair_counter.states<<"\n";
    std::cout << "  witnesses (line,line,key_lo,key_hi,epsilon):\n";
    for(auto [a,b,key]:witnesses) {
        auto c=curvature(a,full,ts);
        std::cout << "    "<<a<<" "<<b<<" 0x"<<std::hex<<uint32_t(key>>32)
                  <<" 0x"<<uint32_t(key)<<std::dec<<" eps="<<c[3]<<"\n";
    }
    std::cout << "  every xor=0x02080; qpair=(0,2); flux=0; epsilon groups=(-4,-2,0,2)x4\n";
    std::cout << "  fixed-layer first moment separates every collision.\n\n";
    std::cout << "SOURCE-q BLACK BOUNDARY FLOW\n";
    std::cout << "  black endpoint q polynomial coefficients=(";
    for(int q=0;q<7;++q) std::cout<<(q?",":"")<<black_q[q];
    std::cout << ")\n";
    for(int q=0;q<7;++q) {
        std::cout << "  q="<<q<<": outward="<<flow[q].out<<"/"<<flow[q].mixed
                  <<" reverse="<<flow[q].rev<<"/"<<flow[q].black<<" strict=True\n";
    }
    std::cout << "  raw non-tied outward/reverse=16118/24370; tied cross-phase="<<tied_cross_phase<<"\n";
    std::cout << "  all-size denominator: 2^(M-2n+5)(4+(1+x)^2)(3+x)^(n-4)-2^(R-n+2)G_n(x)\n\n";
    std::cout << "TOURNAMENT ANALYSIS (q strata as vertices)\n";
    std::cout << "  observable=larger outward-minus-reverse signal\n";
    std::cout << "  switches=raw count bias / source-normalized rate bias; tie path=q0<...<q6; flips="<<gauge_flips<<"\n";
    print_fingerprint("raw",raw_fp); print_fingerprint("normalized",norm_fp);
    std::cout << "\nThe count-only fixed ABC layer is the exact codec kernel. Curvature controls\n"
                 "flow but cannot repair these balanced positional collisions.\n";
    return 0;
}

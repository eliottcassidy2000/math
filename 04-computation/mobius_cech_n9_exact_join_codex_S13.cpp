// Exact n=9 lower-codec collision join for THM-828 (using THM-818's R8).
//
// Input is the certified raw little-endian uint16 merged-node rank atlas for
// all 2^21 n=8 masks.  The engine reconstructs H8, its off-diagonal kernel
// lists L[d]={x:H8(x)=H8(x xor d)}, enumerates the 2^14 apex-zero S2 parity
// syndromes, and uses the exact A/C two-face factorization.  For every
// nonzero syndrome D, dA and dC are nonzero; A and C cover every n=9 cell
// except the apex, so an overlap-compatible pair of faces glues uniquely once
// the apex is fixed to zero.  The B condition, upper colour, and literal S2
// word are then checked directly.

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

struct Tile { int a, b; };

static std::vector<Tile> tiles(int n) {
    std::vector<Tile> out;
    for (int b=1;b<=n-2;++b) for (int a=n;a>=b+2;--a) out.push_back({a,b});
    return out;
}

static int pos(const std::vector<Tile>& t,int a,int b) {
    for(int i=0;i<(int)t.size();++i) if(t[i].a==a&&t[i].b==b) return i;
    throw std::runtime_error("tile not found");
}

static std::vector<int> reflection(int n,const std::vector<Tile>& t) {
    std::vector<int> s(t.size());
    for(int i=0;i<(int)t.size();++i) s[i]=pos(t,n-t[i].b+1,n-t[i].a+1);
    return s;
}

static uint32_t reflect_mask(uint32_t x,const std::vector<int>& s) {
    uint32_t y=0; for(int i=0;i<(int)s.size();++i) y|=((x>>i)&1u)<<s[i]; return y;
}

struct FaceMaps {
    // upper bit -> local lower bit; -1 means absent
    std::array<std::vector<int>,3> down;
    std::array<uint32_t,3> cover{};
};

static FaceMaps face_maps(int n,const std::vector<Tile>& up,const std::vector<Tile>& lo) {
    FaceMaps f; for(auto& v:f.down) v.assign(up.size(),-1);
    for(int i=0;i<(int)up.size();++i) {
        int a=up[i].a,b=up[i].b;
        if(a<n) f.down[0][i]=pos(lo,a,b);
        if(a-b>=3) f.down[1][i]=pos(lo,a-1,b);
        if(b>=2) f.down[2][i]=pos(lo,a-1,b-1);
        for(int r=0;r<3;++r) if(f.down[r][i]>=0) f.cover[r]|=1u<<i;
    }
    return f;
}

static uint32_t restrict_face(uint32_t up,int role,const FaceMaps& f) {
    uint32_t x=0;
    for(int i=0;i<(int)f.down[role].size();++i) if(f.down[role][i]>=0)
        x|=((up>>i)&1u)<<f.down[role][i];
    return x;
}

static uint32_t lift_face(uint32_t x,int role,const FaceMaps& f) {
    uint32_t up=0;
    for(int i=0;i<(int)f.down[role].size();++i) if(f.down[role][i]>=0)
        up|=((x>>f.down[role][i])&1u)<<i;
    return up;
}

static uint16_t overlap_key(uint32_t x,int role_a,int role_c,const FaceMaps& f) {
    uint16_t key=0; int k=0;
    for(int i=0;i<(int)f.down[0].size();++i)
        if(f.down[role_a][i]>=0 && f.down[role_c][i]>=0)
            key|=((x>>f.down[role_a][i])&1u)<<k++;
    assert(k==15); return key;
}

static int comp_rank(const std::array<int,4>& want,int total) {
    int r=0;
    for(int a=0;a<=total;++a) for(int b=0;b<=total-a;++b)
      for(int c=0;c<=total-a-b;++c,++r)
        if(std::array<int,4>{a,b,c,total-a-b-c}==want) return r;
    throw std::runtime_error("bad composition");
}

static std::array<int,7> s2_word(uint32_t x,const std::vector<Tile>& t,
                                 const std::vector<int>& sigma) {
    std::array<std::array<int,4>,6> c{}; int fixed=0;
    for(int i=0;i<(int)t.size();++i) {
        int tau=t[i].a+t[i].b-1;
        if(tau<9) {
            int u=(x>>i)&1u,v=(x>>sigma[i])&1u;
            ++c[tau-3][2*u+v];
        } else if(tau==9) { assert(sigma[i]==i); fixed+=(x>>i)&1u; }
    }
    const std::array<int,6> sizes{1,1,2,2,3,3};
    std::array<int,7> w{};
    for(int i=0;i<6;++i) w[i]=comp_rank(c[i],sizes[i]);
    w[6]=fixed; return w;
}

static std::array<std::array<int,4>,7> s2_first_moments(
        uint32_t x,const std::vector<Tile>& t,const std::vector<int>& sigma) {
    std::array<std::array<int,4>,7> m{};
    std::array<int,7> position{};
    for(int i=0;i<(int)t.size();++i) {
        int tau=t[i].a+t[i].b-1;
        if(tau<9) {
            int layer=tau-3, state=2*int((x>>i)&1u)+int((x>>sigma[i])&1u);
            m[layer][state]+=position[layer]++;
        } else if(tau==9) {
            if((x>>i)&1u) m[6][1]+=position[6];
            ++position[6];
        }
    }
    return m;
}

static std::array<int,6> antisymmetric_moment(
        uint32_t x,const std::vector<Tile>& t,const std::vector<int>& sigma) {
    std::array<int,6> j{},position{};
    for(int i=0;i<(int)t.size();++i) {
        int tau=t[i].a+t[i].b-1;
        if(tau<9) {
            int layer=tau-3;
            j[layer]+=position[layer]++*(int((x>>i)&1u)-int((x>>sigma[i])&1u));
        }
    }
    return j;
}

static int lexicographic_orientation(const std::array<int,6>& j) {
    for(int z:j) if(z) return z>0?1:-1;
    return 0;
}

static std::array<int,4> curvature(uint32_t x,int n,const std::vector<Tile>& t) {
    int B=0,T=0,q0=0;
    auto bit=[&](int a,int b){return int((x>>pos(t,a,b))&1u);};
    for(int a=3;a<n;++a) B+=bit(a,1);
    for(int b=2;b<n-1;++b) T+=bit(n,b);
    for(int a=3;a<n-1;++a) q0+=bit(a,1)*bit(n,a);
    int q1=q0+bit(n-1,1)+bit(n,2);
    return {q0,q1,n-2-B-T,B-T};
}

static std::vector<uint32_t> admissible_differences(const std::vector<Tile>& t,
                                                     const std::vector<int>& sigma) {
    std::array<std::vector<std::pair<int,int>>,6> pairs;
    std::vector<int> fixed;
    for(int i=0;i<(int)t.size();++i) {
        int tau=t[i].a+t[i].b-1;
        if(sigma[i]==i) fixed.push_back(i);
        else if(i<sigma[i]) {
            int j=sigma[i], tau_j=t[j].a+t[j].b-1;
            int low=(tau<9)?i:j, high=(tau<9)?j:i;
            int low_tau=std::min(tau,tau_j);
            assert(low_tau>=3 && low_tau<=8);
            pairs[low_tau-3].push_back({low,high});
        }
    }
    assert(fixed.size()==4 && std::find(fixed.begin(),fixed.end(),0)!=fixed.end());
    std::vector<std::vector<uint32_t>> local;
    for(const auto& ps:pairs) {
        std::vector<uint32_t> choices; int s=ps.size();
        for(int lo=0;lo<(1<<s);++lo) if(__builtin_parity((unsigned)lo)==0)
          for(int hi=0;hi<(1<<s);++hi) if(__builtin_parity((unsigned)hi)==0) {
            uint32_t d=0; for(int j=0;j<s;++j) {
                if((lo>>j)&1) d|=1u<<ps[j].first;
                if((hi>>j)&1) d|=1u<<ps[j].second;
            } choices.push_back(d);
          }
        local.push_back(std::move(choices));
    }
    std::vector<uint32_t> fc;
    for(int z=0;z<16;++z) if(__builtin_parity((unsigned)z)==0) {
        // fixed is in tile order, and apex is not assumed first in this loop.
        uint32_t d=0; bool apex_ok=true;
        for(int j=0;j<4;++j) { if(fixed[j]==0 && ((z>>j)&1)) apex_ok=false; if((z>>j)&1)d|=1u<<fixed[j]; }
        if(apex_ok) fc.push_back(d);
    }
    std::vector<uint32_t> ds{0};
    for(const auto& choices:local) {
        std::vector<uint32_t> next; next.reserve(ds.size()*choices.size());
        for(uint32_t a:ds) for(uint32_t b:choices) next.push_back(a|b);
        ds.swap(next);
    }
    std::vector<uint32_t> all; all.reserve(ds.size()*fc.size());
    for(uint32_t a:ds) for(uint32_t b:fc) all.push_back(a|b);
    std::sort(all.begin(),all.end()); all.erase(std::unique(all.begin(),all.end()),all.end());
    assert(all.size()==16384 && all.front()==0); return all;
}

struct Entry { uint32_t h, x; };
struct DX { uint32_t d,x; };
static bool operator<(const DX&a,const DX&b){return std::tie(a.d,a.x)<std::tie(b.d,b.x);}

struct Stats {
    uint64_t syndromes=0, ac_matches=0, b_kernel=0, upper_colour=0;
    std::array<uint64_t,7> s2_prefix{};
    uint64_t d_face_support=0,d_ac_support=0,d_b_support=0,d_full_support=0;
    size_t max_la=0,max_lc=0,max_c_projection_bucket=0;
    uint64_t max_ac_matches_per_d=0;
};

struct Witness {
    uint32_t D,u,v,xA,xC,xB,dA,dC,dB;
    int first_m1_layer;
};

int main(int argc,char**argv) {
    std::string atlas="/tmp/n8_merged_node_rank_u16.bin", witnesses="/tmp/n9_join_witnesses.tsv",summary,json;
    uint32_t shard_count=1,shard_index=0;
    for(int i=1;i<argc;++i) {
        std::string a=argv[i];
        if(a=="--node-atlas"&&i+1<argc) atlas=argv[++i];
        else if(a=="--witnesses"&&i+1<argc) witnesses=argv[++i];
        else if(a=="--summary"&&i+1<argc) summary=argv[++i];
        else if(a=="--json"&&i+1<argc) json=argv[++i];
        else if(a=="--shard-count"&&i+1<argc) shard_count=std::stoul(argv[++i]);
        else if(a=="--shard-index"&&i+1<argc) shard_index=std::stoul(argv[++i]);
        else throw std::runtime_error("usage: [--node-atlas PATH] [--witnesses PATH] [--summary PATH] [--json PATH] [--shard-count N --shard-index I]");
    }
    if(shard_index>=shard_count) throw std::runtime_error("bad shard");
    const uint32_t N=1u<<21,FULL8=N-1;
    std::vector<uint16_t> node(N);
    std::ifstream in(atlas,std::ios::binary); in.read((char*)node.data(),node.size()*sizeof(uint16_t));
    if(!in || in.peek()!=EOF) throw std::runtime_error("node atlas must be exactly 2^21 uint16 values");
    if(*std::max_element(node.begin(),node.end())!=3527) throw std::runtime_error("unexpected node ranks");

    auto t8=tiles(8),t9=tiles(9); auto s8=reflection(8,t8),s9=reflection(9,t9);
    auto fm=face_maps(9,t9,t8);
    std::vector<uint32_t> H(N); std::vector<Entry> e; e.reserve(N);
    for(uint32_t x=0;x<N;++x) {
        uint32_t blue=(reflect_mask(x,s8)==x);
        H[x]=((uint32_t(node[x])*3528u+node[x^FULL8])<<1)|blue;
        e.push_back({H[x],x});
    }
    std::sort(e.begin(),e.end(),[](const Entry&a,const Entry&b){return std::tie(a.h,a.x)<std::tie(b.h,b.x);});
    std::vector<DX> dx; dx.reserve(3900264); uint64_t rrows=0; size_t support=0,maxf=0;
    for(size_t i=0;i<e.size();) {
        size_t j=i+1; while(j<e.size()&&e[j].h==e[i].h)++j;
        size_t f=j-i; ++support; maxf=std::max(maxf,f); rrows+=f*f;
        for(size_t p=i;p<j;++p) for(size_t q=i;q<j;++q) if(p!=q)
            dx.push_back({e[p].x^e[q].x,e[p].x});
        i=j;
    }
    if(support!=876512||maxf!=26||rrows!=5997416||dx.size()!=3900264)
        throw std::runtime_error("R8 census mismatch");
    std::sort(dx.begin(),dx.end());
    size_t d_support=0,d_max=0; for(size_t i=0;i<dx.size();) {size_t j=i+1;while(j<dx.size()&&dx[j].d==dx[i].d)++j;++d_support;d_max=std::max(d_max,j-i);i=j;}
    if(d_support!=249149||d_max!=5360) throw std::runtime_error("L[d] census mismatch");

    auto range_for=[&](uint32_t d){
        DX lo{d,0},hi{d+1,0};
        auto a=std::lower_bound(dx.begin(),dx.end(),lo), b=std::lower_bound(dx.begin(),dx.end(),hi);
        return std::pair(a,b);
    };
    auto syndromes=admissible_differences(t9,s9);
    const std::set<uint32_t> expected_absent{0x1026286,0x4dd3c9e,0x54a5692,0x5537214};
    std::map<uint32_t,std::array<uint64_t,5>> absent_gate_counts;
    Stats st;
    std::vector<Witness> found;
    std::array<uint64_t,7> first_m1{};
    std::map<uint32_t,uint64_t> full_by_d;
    for(size_t ordinal=1;ordinal<syndromes.size();++ordinal) {
        if((ordinal-1)%shard_count!=shard_index) continue;
        uint32_t D=syndromes[ordinal]; ++st.syndromes;
        uint32_t dA=restrict_face(D,0,fm),dB=restrict_face(D,1,fm),dC=restrict_face(D,2,fm);
        if(!dA||!dC) throw std::runtime_error("nonzero syndrome has zero A/C face difference");
        auto ra=range_for(dA),rc=range_for(dC);
        if(expected_absent.count(D)) absent_gate_counts[D]={uint64_t(ra.second-ra.first),uint64_t(rc.second-rc.first),0,0,0};
        if(ra.first==ra.second||rc.first==rc.second) continue;
        ++st.d_face_support;
        st.max_la=std::max(st.max_la,size_t(ra.second-ra.first));
        st.max_lc=std::max(st.max_lc,size_t(rc.second-rc.first));
        std::vector<std::pair<uint16_t,uint32_t>> ci; ci.reserve(rc.second-rc.first);
        for(auto it=rc.first;it!=rc.second;++it) ci.push_back({overlap_key(it->x,2,0,fm),it->x});
        std::sort(ci.begin(),ci.end());
        for(size_t i=0;i<ci.size();) {size_t j=i+1;while(j<ci.size()&&ci[j].first==ci[i].first)++j;st.max_c_projection_bucket=std::max(st.max_c_projection_bucket,j-i);i=j;}
        uint64_t local_ac=0,local_b=0,local_full=0;
        for(auto ia=ra.first;ia!=ra.second;++ia) {
            uint32_t xA=ia->x; if(xA>(xA^dA)) continue; // unique endpoint orientation: A is always first nonzero face
            uint16_t key=overlap_key(xA,0,2,fm);
            auto lo=std::lower_bound(ci.begin(),ci.end(),std::pair<uint16_t,uint32_t>{key,0});
            auto hi=std::upper_bound(ci.begin(),ci.end(),std::pair<uint16_t,uint32_t>{key,UINT32_MAX});
            for(auto ic=lo;ic!=hi;++ic) {
                ++st.ac_matches; ++local_ac; uint32_t xC=ic->second;
                uint32_t u=lift_face(xA,0,fm)|lift_face(xC,2,fm), v=u^D;
                assert((u&1u)==0&&(v&1u)==0&&u!=v&&(u^v)==D);
                assert(restrict_face(u,0,fm)==xA&&restrict_face(u,2,fm)==xC);
                assert(H[xA]==H[xA^dA]&&H[xC]==H[xC^dC]);
                uint32_t xB=restrict_face(u,1,fm);
                if(H[xB]!=H[xB^dB]) continue;
                ++st.b_kernel; ++local_b;
                if((reflect_mask(u,s9)==u)!=(reflect_mask(v,s9)==v)) continue;
                ++st.upper_colour;
                auto wu=s2_word(u,t9,s9),wv=s2_word(v,t9,s9); bool ok=true;
                for(int k=0;k<7;++k) { ok=ok&&(wu[k]==wv[k]); if(ok)++st.s2_prefix[k]; }
                if(!ok) continue;
                ++local_full;
                auto mu=s2_first_moments(u,t9,s9),mv=s2_first_moments(v,t9,s9);
                int sep=-1; for(int k=0;k<7&&sep<0;++k) if(mu[k]!=mv[k]) sep=k;
                if(sep<0) throw std::runtime_error("S2+M1 failed to separate a collision");
                ++first_m1[sep];
                if(u>v) {std::swap(u,v);xA^=dA;xB^=dB;xC^=dC;}
                found.push_back({D,u,v,xA,xC,xB,dA,dC,dB,sep});
            }
        }
        if(local_ac)++st.d_ac_support;
        if(local_b)++st.d_b_support;
        if(local_full){++st.d_full_support;full_by_d[D]=local_full;}
        if(expected_absent.count(D)) {
            absent_gate_counts[D][2]=local_ac;
            absent_gate_counts[D][3]=local_b;
            absent_gate_counts[D][4]=local_full;
        }
        st.max_ac_matches_per_d=std::max(st.max_ac_matches_per_d,local_ac);
    }

    std::sort(found.begin(),found.end(),[](const Witness&a,const Witness&b){return std::tie(a.u,a.v)<std::tie(b.u,b.v);});
    for(size_t i=1;i<found.size();++i) if(found[i-1].u==found[i].u&&found[i-1].v==found[i].v)
        throw std::runtime_error("duplicate canonical endpoint pair");
    std::set<uint32_t> endpoints; std::set<uint64_t> pairs;
    uint64_t digest=14695981039346656037ull;
    auto mix=[&](uint32_t z){for(int j=0;j<4;++j){digest^=(z>>(8*j))&255u;digest*=1099511628211ull;}};
    uint64_t reflection_fixed=0,b_diagonal=0,pure_black=0,three_cross=0;
    std::map<int,uint64_t> weight_hist;
    uint64_t zero_total_j=0,zero_clock_j=0;
    std::set<std::array<int,6>> j_vectors;
    std::map<uint32_t,std::array<int,6>> j_by_d;
    std::set<std::array<uint16_t,6>> face_node_triples;
    std::map<std::pair<int,int>,uint64_t> curvature_hist;
    uint64_t curvature_tuple_preserved=0,smith_balance=0;
    for(const auto&w:found) {
        endpoints.insert(w.u);endpoints.insert(w.v);pairs.insert((uint64_t(w.u)<<28)|w.v);
        mix(w.u);mix(w.v);++weight_hist[__builtin_popcount(w.D)];
        b_diagonal+=(w.dB==0);
        if(reflect_mask(w.u,s9)!=w.v||reflect_mask(w.v,s9)!=w.u)
            throw std::runtime_error("collision is not its grid-reflection mate");
        ++reflection_fixed;
        assert(s2_word(w.u,t9,s9)==s2_word(w.v,t9,s9));
        bool u_black=reflect_mask(w.u,s9)!=w.u;
        bool a_black=reflect_mask(w.xA,s8)!=w.xA,b_black=reflect_mask(w.xB,s8)!=w.xB,
             c_black=reflect_mask(w.xC,s8)!=w.xC;
        pure_black+=(u_black&&a_black&&b_black&&c_black);
        bool a_cross=node[w.xA]!=node[w.xA^FULL8],b_cross=node[w.xB]!=node[w.xB^FULL8],
             c_cross=node[w.xC]!=node[w.xC^FULL8];
        three_cross+=(a_cross&&b_cross&&c_cross);
        face_node_triples.insert({node[w.xA],node[w.xA^FULL8],node[w.xB],node[w.xB^FULL8],
                                  node[w.xC],node[w.xC^FULL8]});
        auto cu=curvature(w.u,9,t9),cv=curvature(w.v,9,t9);
        curvature_tuple_preserved+=(cu==cv);
        smith_balance+=(cu[3]==0&&cv[3]==0);
        ++curvature_hist[{cu[0],cu[1]}];
        auto ju=antisymmetric_moment(w.u,t9,s9),jv=antisymmetric_moment(w.v,t9,s9);
        for(int k=0;k<6;++k) if(jv[k]!=-ju[k]) throw std::runtime_error("J does not reverse under reflection");
        if(lexicographic_orientation(ju)==0||lexicographic_orientation(ju)!=-lexicographic_orientation(jv))
            throw std::runtime_error("lexicographic J sign failed to orient collision");
        j_vectors.insert(ju);
        auto [it,inserted]=j_by_d.emplace(w.D,ju);
        if(!inserted&&it->second!=ju) throw std::runtime_error("one D sector has multiple oriented J vectors");
        int total=0,clock=0;for(int k=0;k<6;++k){total+=ju[k];clock+=(k+3)*ju[k];}
        zero_total_j+=(total==0);zero_clock_j+=(clock==0);
    }
    std::vector<uint32_t> sectors; for(auto [d,c]:full_by_d) sectors.push_back(d);
    uint64_t cube_edges=0; int cube_diameter=0;
    std::array<uint64_t,4> cube_coordinate_edges{};
    std::map<int,uint64_t> cube_degree_hist;
    if(shard_count==1) {
        const std::array<uint64_t,7> expected{636,636,88,82,68,58,58};
        if(st.ac_matches!=9540||st.b_kernel!=636||st.upper_colour!=636||st.s2_prefix!=expected||found.size()!=58)
            throw std::runtime_error("full exact join census mismatch");
        const std::array<uint64_t,7> expected_m1{0,0,38,12,8,0,0};
        if(endpoints.size()!=116||reflection_fixed!=58||first_m1!=expected_m1||b_diagonal!=0||
           pure_black!=58||three_cross!=58||face_node_triples.size()!=58)
            throw std::runtime_error("collision classification mismatch");
        if(weight_hist!=std::map<int,uint64_t>{{8,44},{12,8},{16,6}}||full_by_d.size()!=11)
            throw std::runtime_error("difference-mask classification mismatch");
        if(j_vectors.size()!=11||j_by_d.size()!=11||zero_total_j!=42||zero_clock_j!=8)
            throw std::runtime_error("antisymmetric-moment classification mismatch");
        if(st.d_face_support!=319||st.d_ac_support!=45||st.d_b_support!=16||st.d_full_support!=11||
           st.max_la!=4104||st.max_lc!=4104||st.max_c_projection_bucket!=16||st.max_ac_matches_per_d!=7944||
           digest!=0x53b4b074be8ae851ull)
            throw std::runtime_error("join sparsity or digest mismatch");
        const std::map<std::pair<int,int>,uint64_t> expected_curvature{
            {{0,0},13},{{0,1},5},{{0,2},8},{{1,1},5},{{1,2},2},{{1,3},2},
            {{2,3},6},{{3,3},3},{{3,4},2},{{3,5},10},{{4,5},1},{{5,7},1}};
        if(curvature_tuple_preserved!=58||smith_balance!=58||curvature_hist!=expected_curvature)
            throw std::runtime_error("black-curvature disintegration mismatch");
        const std::array<uint32_t,4> defect_basis{0x192486,0x8c2c0c,0x11b4600,0x4483414};
        std::set<uint32_t> defect_span;
        for(int z=0;z<16;++z){uint32_t d=0;for(int i=0;i<4;++i)if((z>>i)&1)d^=defect_basis[i];defect_span.insert(d);}
        std::set<uint32_t> absent=defect_span;
        absent.erase(0); for(auto [d,c]:full_by_d) absent.erase(d);
        if(defect_span.size()!=16||absent!=expected_absent)
            throw std::runtime_error("rank-four reflection-defect cube mismatch");
        std::map<uint32_t,int> sector_index;
        for(int i=0;i<(int)sectors.size();++i) sector_index[sectors[i]]=i;
        std::vector<std::vector<int>> dist(sectors.size(),std::vector<int>(sectors.size(),99));
        std::vector<int> degree(sectors.size());
        for(int i=0;i<(int)sectors.size();++i){
            dist[i][i]=0;
            for(int k=0;k<4;++k){
                auto it=sector_index.find(sectors[i]^defect_basis[k]);
                if(it!=sector_index.end()&&i<it->second){
                    int j=it->second;dist[i][j]=dist[j][i]=1;++degree[i];++degree[j];
                    ++cube_edges;++cube_coordinate_edges[k];
                }
            }
        }
        for(int k=0;k<(int)sectors.size();++k)for(int i=0;i<(int)sectors.size();++i)
          for(int j=0;j<(int)sectors.size();++j)dist[i][j]=std::min(dist[i][j],dist[i][k]+dist[k][j]);
        for(int i=0;i<(int)sectors.size();++i){++cube_degree_hist[degree[i]];
          for(int j=0;j<(int)sectors.size();++j)cube_diameter=std::max(cube_diameter,dist[i][j]);}
        if(cube_edges!=14||cube_diameter!=5||cube_degree_hist!=std::map<int,uint64_t>{{1,1},{2,4},{3,5},{4,1}}||
           cube_coordinate_edges!=std::array<uint64_t,4>{4,3,3,4})
            throw std::runtime_error("occupied reflection-defect cube graph mismatch");
        const std::map<uint32_t,std::array<uint64_t,5>> expected_gates{
            {0x1026286,{0,0,0,0,0}}, {0x4dd3c9e,{8,8,4,4,0}},
            {0x54a5692,{1028,1028,7944,504,0}}, {0x5537214,{12,12,4,4,0}}};
        if(absent_gate_counts!=expected_gates)
            throw std::runtime_error("missing-sector gate genealogy mismatch");
    }
    uint64_t sector_gauge_flips=0;
    auto first_j=[&](uint32_t d){int k=0;while(k<6&&j_by_d.at(d)[k]==0)++k;return k+3;};
    for(size_t i=0;i<sectors.size();++i)for(size_t j=i+1;j<sectors.size();++j){
        uint32_t a=sectors[i],b=sectors[j];
        auto structural_a=std::tuple{first_j(a),__builtin_popcount(a),full_by_d[a],a};
        auto structural_b=std::tuple{first_j(b),__builtin_popcount(b),full_by_d[b],b};
        auto empirical_a=std::tuple{-int64_t(full_by_d[a]),first_j(a),__builtin_popcount(a),a};
        auto empirical_b=std::tuple{-int64_t(full_by_d[b]),first_j(b),__builtin_popcount(b),b};
        sector_gauge_flips+=((structural_a<structural_b)!=(empirical_a<empirical_b));
    }
    if(shard_count==1&&sector_gauge_flips!=22)
        throw std::runtime_error("reflection-defect Tournament Analysis mismatch");
    std::ofstream wit(witnesses); if(!wit) throw std::runtime_error("cannot open witnesses");
    wit<<"D\tu\tv\txA\txC\txB\tdA\tdC\tdB\tfirst_M1_layer\tJ3,J4,J5,J6,J7,J8\teta\n";
    for(const auto&w:found) {
        auto j=antisymmetric_moment(w.u,t9,s9);
        wit<<std::hex<<w.D<<'\t'<<w.u<<'\t'<<w.v<<'\t'<<w.xA<<'\t'<<w.xC<<'\t'<<w.xB<<'\t'<<w.dA<<'\t'<<w.dC<<'\t'<<w.dB<<std::dec<<'\t'<<w.first_m1_layer+3<<'\t';
        for(int k=0;k<6;++k)wit<<(k?",":"")<<j[k];
        wit<<'\t'<<lexicographic_orientation(j)<<'\n';
    }
    std::ostringstream report;
    report<<"THM-828 EXACT N=9 TWO-FACE JOIN\n"
             <<"R8 rows/offdiag="<<rrows<<'/'<<dx.size()<<" H support/max="<<support<<'/'<<maxf<<"\n"
             <<"L[d] support/max="<<d_support<<'/'<<d_max<<" syndromes(total/nonzero)="<<syndromes.size()<<"/"<<syndromes.size()-1<<"\n"
             <<"shard="<<shard_index<<'/'<<shard_count<<" processed="<<st.syndromes<<"\n"
             <<"AC matches="<<st.ac_matches<<" B-kernel="<<st.b_kernel<<" upper-colour="<<st.upper_colour<<"\n"
             <<"D support face-pair/AC/B/full="<<st.d_face_support<<'/'<<st.d_ac_support<<'/'<<st.d_b_support<<'/'<<st.d_full_support
             <<" max L_A/L_C/C-bucket/AC-per-D="<<st.max_la<<'/'<<st.max_lc<<'/'<<st.max_c_projection_bucket<<'/'<<st.max_ac_matches_per_d<<"\n"
             <<"S2 prefix tau3,tau4,tau5,tau6,tau7,tau8,fixed9=";
    for(int k=0;k<7;++k) report<<(k?",":"")<<st.s2_prefix[k];
    report<<"\nfull witnesses="<<found.size()<<" endpoints="<<endpoints.size()
             <<" reflection-fixed-pairs="<<reflection_fixed<<" pair-FNV1a64=0x"<<std::hex<<digest<<std::dec<<"\n"
             <<"face-difference-nonzero ABC=111:"<<(found.size()-b_diagonal)<<",101:"<<b_diagonal<<"\n"
             <<"pure-black UABC=KKKK="<<pure_black<<" three-cross="<<three_cross
             <<" distinct-six-node-triples="<<face_node_triples.size()<<"\n"
             <<"D multiplicities=";for(auto [d,c]:full_by_d)report<<std::hex<<d<<std::dec<<':'<<c<<',';
    report<<" weights=";for(auto [w,c]:weight_hist)report<<w<<':'<<c<<',';
    report<<" first-M1 tau5/tau6/tau7="<<first_m1[2]<<'/'<<first_m1[3]<<'/'<<first_m1[4]
             <<" J-vectors="<<j_vectors.size()<<" zero-sum/clock="<<zero_total_j<<'/'<<zero_clock_j
             <<" file="<<witnesses<<"\n"
             <<"reflection-defect span rank=4 occupied/missing=11/4 missing=1026286,4dd3c9e,54a5692,5537214\n"
             <<"occupied-cube graph edges/diameter=14/5 degrees=1:1,2:4,3:5,4:1 coordinate-edges=4,3,3,4\n"
             <<"missing-sector gates d:LA/LC/AC/B/S2=";
    for(auto [d,a]:absent_gate_counts)report<<std::hex<<d<<std::dec<<':'<<a[0]<<'/'<<a[1]<<'/'<<a[2]<<'/'<<a[3]<<'/'<<a[4]<<',';
    report<<"\n"
             <<"black curvature tuple-preserved/Smith-wall="<<curvature_tuple_preserved<<'/'<<smith_balance<<" q-pairs=";
    for(auto [q,c]:curvature_hist)report<<'('<<q.first<<','<<q.second<<"):"<<c<<',';
    report<<"\nraw-codec fibres singletons/doubletons/cells/excess/max="
             <<((1u<<27)-116)<<"/58/"<<((1u<<27)-58)<<"/58/2\n"
             <<"TOURNAMENT ANALYSIS vertices=11 occupied reflection-defect sectors\n"
             <<"  observable=(first skew layer,Hamming weight,fibre multiplicity); switches=structural/empirical; tie-path=hex D\n"
             <<"  both gauges: score histogram 0..10 once, C3=0, SCC=1^11, Hamiltonian paths=1; edge flips="<<sector_gauge_flips<<"\n";
    std::cout<<report.str();
    if(!summary.empty()){std::ofstream o(summary);if(!o)throw std::runtime_error("cannot open summary");o<<report.str();}
    if(!json.empty()) {
        std::ofstream o(json);if(!o)throw std::runtime_error("cannot open json");
        o<<"{\n  \"schema_version\": 1,\n  \"theorem\": \"THM-828\",\n"
         <<"  \"R8_rows\": "<<rrows<<",\n  \"R8_off_diagonal\": "<<dx.size()<<",\n"
         <<"  \"syndrome_dimension\": 14,\n  \"nonzero_syndromes\": "<<st.syndromes<<",\n"
         <<"  \"face_pair_supported_syndromes\": "<<st.d_face_support<<",\n"
         <<"  \"AC_matches\": "<<st.ac_matches<<",\n  \"B_kernel\": "<<st.b_kernel<<",\n"
         <<"  \"upper_colour\": "<<st.upper_colour<<",\n  \"S2_prefix\": [";
        for(int k=0;k<7;++k)o<<(k?", ":"")<<st.s2_prefix[k];
        o<<"],\n  \"collision_pairs\": "<<found.size()<<",\n  \"distinct_endpoints\": "<<endpoints.size()<<",\n"
         <<"  \"reflection_mate_pairs\": "<<reflection_fixed<<",\n  \"difference_sectors\": "<<full_by_d.size()<<",\n"
         <<"  \"pure_black_KKKK\": "<<pure_black<<",\n  \"three_cross_pairs\": "<<three_cross<<",\n"
         <<"  \"distinct_six_node_face_triples\": "<<face_node_triples.size()<<",\n"
         <<"  \"first_J_layer_tau5_tau6_tau7\": ["<<first_m1[2]<<", "<<first_m1[3]<<", "<<first_m1[4]<<"],\n"
         <<"  \"reflection_defect_rank\": 4,\n  \"occupied_reflection_defect_sectors\": 11,\n"
         <<"  \"missing_reflection_defect_sectors\": [\"0x1026286\", \"0x4dd3c9e\", \"0x54a5692\", \"0x5537214\"],\n"
         <<"  \"occupied_cube_graph\": {\"connected\": true, \"edges\": "<<cube_edges<<", \"diameter\": "<<cube_diameter<<", \"degree_histogram\": {\"1\": 1, \"2\": 4, \"3\": 5, \"4\": 1}, \"coordinate_edge_counts\": [4, 3, 3, 4]},\n"
         <<"  \"black_curvature_tuple_preserved\": "<<curvature_tuple_preserved<<",\n  \"smith_balance_wall\": "<<smith_balance<<",\n"
         <<"  \"raw_codec_fibres\": {\"singletons\": "<<((1u<<27)-116)<<", \"doubletons\": 58, \"cells\": "<<((1u<<27)-58)<<", \"collision_excess\": 58, \"maximum\": 2},\n"
         <<"  \"tournament_analysis\": {\"vertices\": \"11 occupied reflection-defect sectors\", \"observable\": [\"first skew layer\", \"Hamming weight\", \"fibre multiplicity\"], \"switches\": [\"structural\", \"empirical\"], \"tie_hamiltonian_path\": \"hex D\", \"directed_3cycles\": 0, \"scc_sizes\": [1,1,1,1,1,1,1,1,1,1,1], \"hamiltonian_paths\": 1, \"edge_flips\": "<<sector_gauge_flips<<"},\n"
         <<"  \"pair_fnv1a64\": \"0x"<<std::hex<<digest<<std::dec<<"\"\n}\n";
    }
}

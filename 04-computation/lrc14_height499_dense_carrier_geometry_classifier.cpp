#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

using V = std::array<int, 3>;

static int igcd3(int x, int y, int z) {
    return std::gcd(std::gcd(std::abs(x), std::abs(y)), std::abs(z));
}

static V canonical_direction(V v) {
    int g = igcd3(v[0], v[1], v[2]);
    for (int &x : v) x /= g;
    for (int x : v) {
        if (x != 0) {
            if (x < 0) for (int &y : v) y = -y;
            break;
        }
    }
    return v;
}

static bool parallel(const V &u, const V &v) {
    return int64_t(u[0]) * v[1] == int64_t(u[1]) * v[0]
        && int64_t(u[0]) * v[2] == int64_t(u[2]) * v[0]
        && int64_t(u[1]) * v[2] == int64_t(u[2]) * v[1];
}

static V add(const V &u, const V &v, int sign) {
    return {u[0] + sign*v[0], u[1] + sign*v[1], u[2] + sign*v[2]};
}

static bool equal_up_to_sign(const V &u, const V &v) {
    return u == v || (u[0] == -v[0] && u[1] == -v[1] && u[2] == -v[2]);
}

static int affine_rank(const std::vector<V> &p) {
    if (p.size() <= 1) return 0;
    V first{p[1][0]-p[0][0],p[1][1]-p[0][1],p[1][2]-p[0][2]};
    for (size_t i=2;i<p.size();++i) {
        V d{p[i][0]-p[0][0],p[i][1]-p[0][1],p[i][2]-p[0][2]};
        V cr{first[1]*d[2]-first[2]*d[1],first[2]*d[0]-first[0]*d[2],first[0]*d[1]-first[1]*d[0]};
        if (cr != V{0,0,0}) return 2;
    }
    return 1;
}

struct CarrierData {
    int count = 0;
    bool rank_one = true;
    bool norm_four_ray = false;
    bool a2_six = false;
    V first{0,0,0};
    std::set<V> directions;
    std::set<V> points;
};

static CarrierData carriers(const V &w) {
    CarrierData out;
    V n{0,0,0};
    const std::set<V> norm4{{1,-2,1},{1,2,-1},{2,1,-1}};
    while (n[0] <= w[0] && n[1] <= w[1] && n[2] <= w[2]) {
        V lo, hi;
        for (int i=0; i<3; ++i) {
            lo[i] = std::max(0, 14*n[i]-3);
            hi[i] = std::min(14*w[i], 14*n[i]+3);
        }
        int il = 0, iu = 0;
        for (int i=1; i<3; ++i) {
            if (int64_t(lo[i])*w[il] > int64_t(lo[il])*w[i]) il = i;
            if (int64_t(hi[i])*w[iu] < int64_t(hi[iu])*w[i]) iu = i;
        }
        if (int64_t(lo[il])*w[iu] < int64_t(hi[iu])*w[il]) {
            int labels[3];
            for (int i=0; i<3; ++i) labels[i] = (3 - (w[i]*n[i])%3)%3;
            if (labels[0] != labels[1] && labels[0] != labels[2]
                    && labels[1] != labels[2]) {
                V c{
                    w[1]*n[2] - w[2]*n[1],
                    w[2]*n[0] - w[0]*n[2],
                    w[0]*n[1] - w[1]*n[0]
                };
                ++out.count;
                out.points.insert(c);
                V d = canonical_direction(c);
                out.directions.insert(d);
                if (out.count == 1) out.first = c;
                else if (!parallel(out.first, c)) out.rank_one = false;
            }
        }
        // Advance every interval attaining the earliest right endpoint.
        for (int i=0; i<3; ++i) {
            if (int64_t(hi[i])*w[iu] == int64_t(hi[iu])*w[i]) ++n[i];
        }
    }
    if (out.count == 0) out.rank_one = false;
    if (out.rank_one) out.norm_four_ray = norm4.count(canonical_direction(out.first));
    if (out.count == 6 && out.directions.size() == 3) {
        std::vector<V> d(out.directions.begin(), out.directions.end());
        for (int i=0; i<3; ++i) for (int j=i+1; j<3; ++j) {
            int k = 3-i-j;
            for (int s : {-1,1}) {
                if (equal_up_to_sign(add(d[i],d[j],s),d[k])) out.a2_six = true;
            }
        }
    }
    return out;
}

static std::set<V> carriers_lattice_box(const V &w) {
    V B{
        (3*(w[1]+w[2])-1)/14,
        (3*(w[0]+w[2])-1)/14,
        (3*(w[0]+w[1])-1)/14
    };
    std::set<V> out;
    for (int x=-B[0]; x<=B[0]; ++x) for (int y=-B[1]; y<=B[1]; ++y) {
        int64_t nz=-int64_t(w[0])*x-int64_t(w[1])*y;
        if (nz%w[2]) continue;
        int z=int(nz/w[2]);
        if (std::abs(z)>B[2] || x%3==0 || y%3==0 || z%3==0) continue;
        out.insert({x,y,z});
    }
    return out;
}

static std::pair<int,V> minimum_l1_relation(const V &w) {
    int g01=std::gcd(w[0],w[1]), g02=std::gcd(w[0],w[2]), g12=std::gcd(w[1],w[2]);
    std::array<std::pair<int,V>,3> seed{{
        {(w[0]+w[1])/g01,{w[1]/g01,-w[0]/g01,0}},
        {(w[0]+w[2])/g02,{w[2]/g02,0,-w[0]/g02}},
        {(w[1]+w[2])/g12,{0,w[2]/g12,-w[1]/g12}}
    }};
    auto initial=*std::min_element(seed.begin(),seed.end(),[](const auto &x,const auto &y){return x.first<y.first;});
    int upper=initial.first;
    V best=initial.second;
    for (int x=-upper; x<=upper; ++x) {
        int rem=upper-std::abs(x);
        for (int y=-rem; y<=rem; ++y) {
            int64_t nz=-int64_t(w[0])*x-int64_t(w[1])*y;
            if (nz%w[2]) continue;
            int z=int(nz/w[2]);
            int norm=std::abs(x)+std::abs(y)+std::abs(z);
            if (norm && norm<upper) { upper=norm; best={x,y,z}; }
        }
    }
    return {upper,canonical_direction(best)};
}

static std::string vs(const V &v) {
    return "("+std::to_string(v[0])+","+std::to_string(v[1])+","+std::to_string(v[2])+")";
}

int main(int argc, char **argv) {
#ifdef _WIN32
    _setmode(_fileno(stdout), _O_BINARY);
#endif
    int H = argc > 1 ? std::stoi(argv[1]) : 499;
    bool audit_dense = argc > 2 && std::string(argv[2]) == "--audit-dense";
    std::vector<int> values;
    for (int x=1; x<=H; x+=2) if (x%3) values.push_back(x);
    int64_t triples=0, dense=0, ray=0, rank2=0, other_rank1=0, a2=0;
    std::map<int,int64_t> dense_by_height, rank2_by_height, other_ray_by_height;
    std::vector<std::tuple<V,CarrierData>> exceptional;
    std::vector<std::tuple<int,int,V,bool,int>> nonnorm4_rows;
    int64_t rank2_half_cores=0;
    int64_t tight_rank2_half_cores=0, zero_sum_half_triangles=0;
    std::vector<std::tuple<V,int,int>> rank2_half_examples;
    std::tuple<int,int,V,int> max_half{0,1,V{0,0,0},0};
    for (size_t k=2; k<values.size(); ++k) {
        int c=values[k];
        for (size_t j=1; j<k; ++j) for (size_t i=0; i<j; ++i) {
            V w{values[i],values[j],c};
            if (igcd3(w[0],w[1],w[2]) != 1) continue;
            ++triples;
            CarrierData cd=carriers(w);
            if (cd.count) {
                std::vector<V> half;
                for (const V &C:cd.points) {
                    bool positive=true, inner=true;
                    for (int q=0;q<3;++q) {
                        int cm=(C[q]%3+3)%3, wm=w[q]%3;
                        positive &= cm==wm;
                        inner &= 28*std::abs(C[q]) < 3*(w[0]+w[1]+w[2]-w[q]);
                    }
                    if (positive && inner) half.push_back(C);
                }
                int m=cd.count/2, h=int(half.size()), ar=affine_rank(half);
                auto [mh,mm,mw,mr]=max_half;
                if (int64_t(h)*mm > int64_t(mh)*m) max_half={h,m,w,ar};
                if (ar==2) {
                    ++rank2_half_cores;
                    if (m==3*h-3) ++tight_rank2_half_cores;
                    if (h==3) {
                        V total{0,0,0};
                        for (const V &q:half) for (int j=0;j<3;++j) total[j]+=q[j];
                        if (total==V{0,0,0}) ++zero_sum_half_triangles;
                    }
                    if (rank2_half_examples.size()<20) rank2_half_examples.push_back({w,h,m});
                }
            }
            if (cd.count && !(cd.rank_one && cd.norm_four_ray))
                nonnorm4_rows.push_back({cd.count,c,w,cd.rank_one,int(cd.directions.size())});
            if (int64_t(11)*cd.count <= int64_t(2)*c) continue;
            if (audit_dense && carriers_lattice_box(w) != cd.points) {
                std::cerr << "dense sweep/lattice mismatch " << vs(w) << "\n";
                return 2;
            }
            ++dense; ++dense_by_height[c];
            if (cd.rank_one && cd.norm_four_ray) { ++ray; continue; }
            if (cd.rank_one) { ++other_rank1; ++other_ray_by_height[c]; }
            else { ++rank2; ++rank2_by_height[c]; }
            if (cd.a2_six) ++a2;
            exceptional.push_back({w,cd});
        }
    }
    std::cout << "LRC DENSE CARRIER GEOMETRY CLASSIFIER\n";
    std::cout << "universe=sorted primitive distinct positive odd ternary-unit triples max_height=" << H << "\n";
    std::cout << "dense_definition=11*N>2*c exact; triples=" << triples << " dense=" << dense << "\n";
    std::cout << "dense_norm4_rank1=" << ray << " dense_other_rank1=" << other_rank1
              << " dense_rank2=" << rank2 << " dense_A2_six=" << a2 << "\n";
    std::cout << "exception_count=" << exceptional.size() << "\n";
    size_t shown=0;
    for (const auto &[w,cd] : exceptional) {
        if (shown++ >= 200) break;
        std::cout << "exception w=" << vs(w) << " N=" << cd.count
                  << " rank=" << (cd.rank_one?1:2) << " A2=" << cd.a2_six << " dirs=";
        bool first=true;
        for (const auto &d:cd.directions) { if(!first) std::cout << ";"; first=false; std::cout << vs(d); }
        std::cout << "\n";
    }
    std::cout << "dense_by_height=";
    for (auto [h,n]:dense_by_height) std::cout << h << ":" << n << ",";
    std::cout << "\nrank2_by_height=";
    for (auto [h,n]:rank2_by_height) std::cout << h << ":" << n << ",";
    std::cout << "\nother_rank1_by_height=";
    for (auto [h,n]:other_ray_by_height) std::cout << h << ":" << n << ",";
    auto [mh,mm,mw,mr]=max_half;
    std::cout << "\nhalf_body rank2_core_count=" << rank2_half_cores
              << " tight_rank2_count=" << tight_rank2_half_cores
              << " zero_sum_triangle_count=" << zero_sum_half_triangles
              << " max_ratio=" << mh << "/" << mm << " max_row=" << vs(mw)
              << " max_core_rank=" << mr << " examples=";
    for (const auto &[w,h,m]:rank2_half_examples) std::cout << vs(w) << ":" << h << "/" << m << ",";
    std::sort(nonnorm4_rows.begin(),nonnorm4_rows.end(),[](const auto &x,const auto &y){
        int nx=std::get<0>(x), cx=std::get<1>(x), ny=std::get<0>(y), cy=std::get<1>(y);
        if (int64_t(nx)*cy != int64_t(ny)*cx) return int64_t(nx)*cy > int64_t(ny)*cx;
        return std::get<2>(x) < std::get<2>(y);
    });
    std::cout << "\ntop_nonnorm4_density=";
    for (size_t i=0; i<std::min<size_t>(30,nonnorm4_rows.size()); ++i) {
        auto [n,c,w,r,d]=nonnorm4_rows[i];
        auto [mu,rel]=minimum_l1_relation(w);
        std::cout << "\n  w=" << vs(w) << " N=" << n << " N/c=" << n << "/" << c
                  << " rank=" << (r?1:2) << " dirs=" << d << " mu1=" << mu
                  << " minrel=" << vs(rel);
    }
    std::cout << "\nverdict=PASS audit_dense=" << audit_dense << "\n";
}

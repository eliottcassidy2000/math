// Independent interval-overlay audit for THM-4207.
//
// Unlike the primary literal joint-wall path, this implementation first
// intersects the newcomer safe-comb interval lists and then overlays that
// intersection on the fixed pool's 7,133 labelled cells.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using i64 = std::int64_t;
using i128 = __int128_t;
using u32 = std::uint32_t;

namespace {

constexpr std::array<int, 30> P = {
    8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
    120,126,132,143,145,168,170,176,190,193,240,252,264,286,290};
constexpr i64 D = INT64_C(18241159416480);

void check(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

i64 lcm_checked(i64 a, i64 b) {
    const i128 value = static_cast<i128>(a / std::gcd(a,b)) * b;
    check(value <= std::numeric_limits<i64>::max(), "lcm overflow");
    return static_cast<i64>(value);
}

std::string dec(i128 value) {
    if (value == 0) return "0";
    const bool neg = value < 0;
    if (neg) value = -value;
    std::string s;
    while (value) {
        s.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

bool safe_mid(int speed, i64 a, i64 b, i64 denom) {
    const i128 raw = static_cast<i128>(speed) * (a+b);
    const i128 residue = raw % (static_cast<i128>(2)*denom);
    return static_cast<i128>(7)*residue >= denom &&
           static_cast<i128>(7)*residue <= static_cast<i128>(13)*denom;
}

struct Interval { i64 a; i64 b; };

std::vector<Interval> comb(int speed, i64 denom) {
    check(denom % (14LL*speed) == 0, "comb denominator mismatch");
    const i64 unit = denom/(14LL*speed);
    std::vector<Interval> out;
    out.reserve(speed);
    for (int k=0;k<speed;++k) {
        out.push_back({(14LL*k+1)*unit,(14LL*k+13)*unit});
    }
    return out;
}

std::vector<Interval> intersect_lists(const std::vector<Interval>& x,
                                      const std::vector<Interval>& y) {
    std::vector<Interval> out;
    std::size_t i=0,j=0;
    while (i<x.size() && j<y.size()) {
        const i64 a=std::max(x[i].a,y[j].a);
        const i64 b=std::min(x[i].b,y[j].b);
        if (a<b) out.push_back({a,b});
        if (x[i].b<y[j].b) ++i;
        else if (y[j].b<x[i].b) ++j;
        else { ++i; ++j; }
    }
    return out;
}

struct FixedCell { i64 a; i64 b; u32 failed; };

std::vector<FixedCell> fixed_cells() {
    std::vector<i64> walls={0,D};
    for (int speed:P) {
        const i64 unit=D/(14LL*speed);
        for (int k=0;k<speed;++k) {
            walls.push_back((14LL*k+1)*unit);
            walls.push_back((14LL*k+13)*unit);
        }
    }
    std::sort(walls.begin(),walls.end());
    walls.erase(std::unique(walls.begin(),walls.end()),walls.end());
    check(walls.size()==7134,"fixed wall count");
    std::vector<FixedCell> cells;
    for (std::size_t i=0;i+1<walls.size();++i) {
        u32 failed=0;
        for (int v=0;v<30;++v) {
            if (!safe_mid(P[v],walls[i],walls[i+1],D)) failed|=u32{1}<<v;
        }
        cells.push_back({walls[i],walls[i+1],failed});
    }
    return cells;
}

struct AtomGeometry {
    i64 denom;
    std::unordered_map<u32,i64> mass;
};

AtomGeometry atoms(const std::vector<FixedCell>& fixed,
                   const std::vector<int>& newcomers) {
    i64 denom=D;
    for (int q:newcomers) denom=lcm_checked(denom,14LL*q);
    std::vector<Interval> allowed={{0,denom}};
    for (int q:newcomers) allowed=intersect_lists(allowed,comb(q,denom));
    const i64 scale=denom/D;
    AtomGeometry result{denom,{}};
    std::size_t j=0;
    i128 allowed_total=0;
    for (const FixedCell& c:fixed) {
        const i64 ca=c.a*scale, cb=c.b*scale;
        while (j<allowed.size() && allowed[j].b<=ca) ++j;
        std::size_t k=j;
        i64 overlap=0;
        while (k<allowed.size() && allowed[k].a<cb) {
            const i64 a=std::max(ca,allowed[k].a);
            const i64 b=std::min(cb,allowed[k].b);
            if (a<b) overlap+=b-a;
            if (allowed[k].b>=cb) break;
            ++k;
        }
        if (overlap) result.mass[c.failed]+=overlap;
        allowed_total+=overlap;
    }
    i64 direct=0;
    for (const Interval& x:allowed) direct+=x.b-x.a;
    check(allowed_total==direct,"allowed sweep lost mass");
    return result;
}

i64 repaired_mass(u32 deletion,const AtomGeometry& g) {
    i64 result=0;
    u32 subset=deletion;
    while (true) {
        const auto it=g.mass.find(subset);
        if (it!=g.mass.end()) result+=it->second;
        if (subset==0) break;
        subset=(subset-1)&deletion;
    }
    return result;
}

i128 delta(u32 deletion,const AtomGeometry& g) {
    return static_cast<i128>(63)*repaired_mass(deletion,g)-
           static_cast<i128>(4)*g.denom;
}

i64 body_mass(u32 body,const AtomGeometry& g) {
    i64 result=0;
    for (const auto& [failed,length]:g.mass) {
        if ((failed&body)==0) result+=length;
    }
    return result;
}

i128 body_delta(u32 body,const AtomGeometry& g) {
    return static_cast<i128>(63)*body_mass(body,g)-
           static_cast<i128>(4)*g.denom;
}

u32 mask(std::initializer_list<int> labels) {
    u32 out=0;
    for (int label:labels) {
        const auto it=std::find(P.begin(),P.end(),label);
        check(it!=P.end(),"mask label");
        out|=u32{1}<<static_cast<int>(it-P.begin());
    }
    return out;
}

std::string labels(u32 x) {
    std::ostringstream out;
    bool first=true;
    for (int v=0;v<30;++v) if (x&(u32{1}<<v)) {
        if (!first) out << ',';
        first=false;
        out << P[v];
    }
    return out.str();
}

u32 next(u32 x) {
    const u32 low=x&(~x+1u), ripple=x+low;
    return ripple|(((x^ripple)>>2)/low);
}

struct FullLayer {
    std::vector<u32> edges;
    std::uint64_t equalities=0;
    i128 best_delta=0;
    u32 best_edge=0;
};

FullLayer full_layer(const AtomGeometry& g,int d) {
    FullLayer out;
    out.best_delta=-static_cast<i128>(4)*g.denom;
    u32 x=(u32{1}<<d)-1, limit=u32{1}<<30;
    while (x && x<limit) {
        const i128 value=delta(x,g);
        if (value==0) ++out.equalities;
        if (value>=0) out.edges.push_back(x);
        if (value>out.best_delta) {
            out.best_delta=value;
            out.best_edge=x;
        }
        x=next(x);
    }
    return out;
}

bool covers(u32 body,const std::vector<u32>& edges,std::uint64_t& checks) {
    for (u32 edge:edges) {
        ++checks;
        if ((body&edge)==0) return false;
    }
    return true;
}

void audit_tau(const std::vector<u32>& edges,int expected) {
    std::uint64_t candidates=0,checks=0;
    for (int k=0;k<=expected;++k) {
        u32 x=k==0?0:(u32{1}<<k)-1, limit=u32{1}<<30;
        while (true) {
            ++candidates;
            if (covers(x,edges,checks)) {
                check(k==expected,"transversal number changed");
                std::cout << "TAU " << k << " WITNESS " << labels(x)
                          << " CANDIDATES " << candidates
                          << " CHECKS " << checks << '\n';
                return;
            }
            if (k==0) break;
            x=next(x);
            if (!x || x>=limit) break;
        }
    }
    throw std::runtime_error("tau upper bound false");
}

std::uint64_t independent_key(std::uint64_t value) {
    value^=UINT64_C(0xd6e8feb86659fd93);
    value^=value>>32;
    value*=UINT64_C(0xd6e8feb86659fd93);
    value^=value>>32;
    value*=UINT64_C(0xd6e8feb86659fd93);
    return value^(value>>32);
}

struct SevenAudit {
    std::uint64_t candidates=0;
    std::uint64_t checks=0;
    std::uint64_t max_checks=0;
    u32 closest=0;
    u32 missed=0;
    u32 cover=0;
};

void recurse_seven(int need,int next_vertex,u32 chosen,
                   const std::vector<u32>& edges,SevenAudit& audit) {
    if (need==0) {
        ++audit.candidates;
        std::uint64_t row=0;
        u32 missed=0;
        for (u32 edge:edges) {
            ++row;
            if ((chosen&edge)==0) {
                missed=edge;
                break;
            }
        }
        audit.checks+=row;
        if (row>audit.max_checks) {
            audit.max_checks=row;
            audit.closest=chosen;
            audit.missed=missed;
        }
        if (missed==0) audit.cover=chosen;
        return;
    }
    for (int vertex=next_vertex;vertex<=30-need;++vertex) {
        recurse_seven(need-1,vertex+1,chosen|(u32{1}<<vertex),edges,audit);
    }
}

void independent_depth8_audit(std::vector<u32> edges,u32 explicit_cover) {
    std::uint64_t control_checks=0;
    check(covers(explicit_cover,edges,control_checks),"explicit d8 cover failed");
    std::sort(edges.begin(),edges.end(),[](u32 a,u32 b) {
        const auto ka=independent_key(a),kb=independent_key(b);
        if (ka!=kb) return ka<kb;
        return a<b;
    });
    SevenAudit audit;
    recurse_seven(7,0,0,edges,audit);
    check(audit.candidates==UINT64_C(2035800),"seven universe changed");
    check(audit.cover==0,"unexpected seven-cover");
    std::cout << "INDEPENDENT_D8_TAU 8 COVER " << labels(explicit_cover)
              << " SEVEN_CANDIDATES " << audit.candidates
              << " CHECKS " << audit.checks
              << " MAX_CHECKS " << audit.max_checks
              << " CLOSEST " << labels(audit.closest)
              << " MISSED " << labels(audit.missed) << '\n';
}

} // namespace

int main() {
    try {
        const auto fixed=fixed_cells();
        const auto q50=atoms(fixed,{50});
        const auto q51=atoms(fixed,{51});
        const auto pair=atoms(fixed,{50,51});
        check(q50.denom==INT64_C(91205797082400) &&
                  q51.denom==INT64_C(18241159416480) &&
                  pair.denom==INT64_C(91205797082400),
              "natural denominators changed");
        const u32 r4=mask({88,95,176,193});
        std::cout << "R4_Q50 DEN " << q50.denom << " NUM " << repaired_mass(r4,q50)
                  << " DELTA63 " << dec(delta(r4,q50)) << '\n';
        std::cout << "R4_Q51 DEN " << q51.denom << " NUM " << repaired_mass(r4,q51)
                  << " DELTA63 " << dec(delta(r4,q51)) << '\n';
        std::cout << "R4_PAIR DEN " << pair.denom << " NUM " << repaired_mass(r4,pair)
                  << " DELTA63 " << dec(delta(r4,pair)) << '\n';
        const std::array<u32,4> all_e4={
            mask({88,95,176,193}), mask({88,145,176,193}),
            mask({88,145,193,290}), mask({145,168,193,290})};
        const FullLayer q50e4=full_layer(q50,4);
        const FullLayer paire4=full_layer(pair,4);
        const FullLayer paire5=full_layer(pair,5);
        check(q50e4.edges.size()==4 && q50e4.equalities==0,
              "q50 depth-four layer changed");
        check(paire4.edges.empty() && paire4.equalities==0 &&
                  paire5.edges.empty() && paire5.equalities==0,
              "pair shallow empty staircase changed");
        for (std::size_t i=0;i<all_e4.size();++i) {
            check(q50e4.edges[i]==all_e4[i],"q50 depth-four edge list changed");
        }
        for (u32 edge:all_e4) {
            check(delta(edge,q51)>0 && delta(edge,pair)<0,
                  "depth-four intersection hostile changed");
            std::cout << "E4_EDGE " << labels(edge)
                      << " Q50 " << dec(delta(edge,q50))
                      << " Q51 " << dec(delta(edge,q51))
                      << " PAIR " << dec(delta(edge,pair)) << '\n';
        }
        const FullLayer e6=full_layer(pair,6);
        const FullLayer e7=full_layer(pair,7);
        check(e6.edges.size()==39 && e6.equalities==0,
              "pair depth-six layer changed");
        check(e7.edges.size()==10114 && e7.equalities==0,
              "pair depth-seven layer changed");
        std::cout << "PAIR_E6 " << e6.edges.size() << '\n';
        audit_tau(e6.edges,2);
        std::cout << "PAIR_E7 " << e7.edges.size() << '\n';
        audit_tau(e7.edges,5);
        FullLayer e8=full_layer(pair,8);
        const u32 d8_cover=mask({16,88,95,143,168,193,240,290});
        const u32 d8_body9=d8_cover|mask({8});
        const u32 d8_primary_miss=mask({8,42,88,132,145,170,176,264});
        check(e8.edges.size()==311544 && e8.equalities==0 &&
                  e8.best_delta==static_cast<i128>(97750631715282LL) &&
                  e8.best_edge==mask({80,88,95,145,168,190,193,290}),
              "pair depth-eight layer changed");
        check(delta(d8_primary_miss,pair)==
                  static_cast<i128>(1334721427452LL),
              "depth-eight missed-edge control changed");
        check(body_delta(d8_body9,pair)==
                  static_cast<i128>(614449941537852LL),
              "depth-eight safe cover-body control changed");
        std::cout << "PAIR_E8 " << e8.edges.size()
                  << " EQUALITIES " << e8.equalities
                  << " BEST_DELTA63 " << dec(e8.best_delta)
                  << " BEST " << labels(e8.best_edge) << '\n';
        independent_depth8_audit(std::move(e8.edges),d8_cover);
        std::cout << "D8_CONTROLS MISSED_EDGE_DELTA63 "
                  << dec(delta(d8_primary_miss,pair))
                  << " COVER_BODY9_DEN " << pair.denom
                  << " COVER_BODY9_DELTA63 " << dec(body_delta(d8_body9,pair))
                  << " COVER_BODY9 " << labels(d8_body9) << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL " << e.what() << '\n';
        return 1;
    }
}

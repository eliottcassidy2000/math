#define main inherited_chart_main
#include "/tmp/math-wt-lrc14-breakthrough-20260826/04-computation/lrc14_fixed_fifty_four_petal_haar_chart_thm4240.cpp"
#undef main

#include <fstream>
#include <tuple>

namespace bridge {

struct HighBody {
    std::uint16_t petals;
    u32 core;
    int cutoff;
    i64 mass;
    i64 components;
    i128 delta;
};

struct Literal {
    i64 numerator;
    i64 denominator;
    i128 scaled_margin;
};

std::string all_labels(const HighBody& body) {
    std::vector<int> labels;
    for (int j=0;j<18;++j) if (body.core&(u32{1}<<j)) labels.push_back(C[j]);
    for (int j=0;j<12;++j) if (body.petals&(std::uint16_t{1}<<j)) labels.push_back(O[j]);
    std::sort(labels.begin(),labels.end());
    std::string out;
    for (int label:labels) {
        if (!out.empty()) out+=',';
        out+=std::to_string(label);
    }
    return out;
}

Literal integrate(const HighBody& body,int outsider) {
    std::vector<int> speeds;
    for (int j=0;j<18;++j) if (body.core&(u32{1}<<j)) speeds.push_back(C[j]);
    for (int j=0;j<12;++j) if (body.petals&(std::uint16_t{1}<<j)) speeds.push_back(O[j]);
    speeds.push_back(Q);
    speeds.push_back(outsider);
    std::sort(speeds.begin(),speeds.end());
    speeds.erase(std::unique(speeds.begin(),speeds.end()),speeds.end());
    require(speeds.size()==11,"literal speed collision");
    i64 denominator=1;
    for (int speed:speeds) denominator=lcm_exact(denominator,14LL*speed);
    std::vector<i64> walls{0,denominator};
    for (int speed:speeds) {
        const i64 unit=denominator/(14LL*speed);
        for (int tooth=0;tooth<speed;++tooth) {
            walls.push_back((14LL*tooth+1)*unit);
            walls.push_back((14LL*tooth+13)*unit);
        }
    }
    std::sort(walls.begin(),walls.end());
    walls.erase(std::unique(walls.begin(),walls.end()),walls.end());
    i64 numerator=0;
    for (std::size_t j=0;j+1<walls.size();++j) {
        bool safe=true;
        for (int speed:speeds) {
            if (!safe_mid(speed,walls[j],walls[j+1],denominator)) {
                safe=false;
                break;
            }
        }
        if (safe) numerator+=walls[j+1]-walls[j];
    }
    return {numerator,denominator,static_cast<i128>(63)*numerator-
                                  static_cast<i128>(4)*denominator};
}

} // namespace bridge

int main() {
    using namespace bridge;
    const BaseArrangement arrangement=make_base_arrangement();
    require(arrangement.denominator==91205797082400LL,"base denominator changed");
    require(arrangement.cells.size()==7213,"base cells changed");
    std::vector<HighBody> high;
    std::uint64_t profiles=0;
    std::array<std::uint64_t,10> profiles_by_k{};
    std::array<int,10> high_by_k{};
    std::array<int,10> max_by_k{};
    for (int k=4;k<=6;++k) {
        const int core_size=9-k;
        for (std::uint16_t petals=0;petals<(std::uint16_t{1}<<12);++petals) {
            if (std::popcount(petals)!=k) continue;
            std::vector<i64> mass(N,0),starts(N,0),continuations(N,0);
            for (std::size_t j=0;j<arrangement.cells.size();++j) {
                const Cell& current=arrangement.cells[j];
                const Cell& previous=arrangement.cells[(j+arrangement.cells.size()-1)%arrangement.cells.size()];
                if (!current.q_safe || (current.o_safe&petals)!=petals) continue;
                mass[current.c_failed]+=current.length;
                starts[current.c_failed]+=1;
                if (previous.q_safe && (previous.o_safe&petals)==petals)
                    continuations[current.c_failed|previous.c_failed]+=1;
            }
            zeta(mass); zeta(starts); zeta(continuations);
            for (u32 core=ALL;;core=(core-1)&ALL) {
                if (std::popcount(core)==core_size) {
                    ++profiles;
                    ++profiles_by_k[k];
                    const u32 allowed=ALL^core;
                    const i64 components=starts[allowed]-continuations[allowed];
                    const i128 delta=static_cast<i128>(54)*mass[allowed]-
                                     static_cast<i128>(4)*arrangement.denominator;
                    require(components>0,"nonpositive components");
                    require(delta>0,"nonstrict base");
                    const i128 numerator=static_cast<i128>(54)*components*arrangement.denominator;
                    const i128 divisor=static_cast<i128>(7)*delta;
                    const int cutoff=static_cast<int>((numerator+divisor-1)/divisor);
                    max_by_k[k]=std::max(max_by_k[k],cutoff);
                    if (cutoff>590) {
                        high.push_back({petals,core,cutoff,mass[allowed],components,delta});
                        ++high_by_k[k];
                    }
                }
                if (core==0) break;
            }
        }
    }
    std::sort(high.begin(),high.end(),[](const HighBody& a,const HighBody& b){
        if (a.cutoff!=b.cutoff) return a.cutoff>b.cutoff;
        return all_labels(a)<all_labels(b);
    });
    std::cout << "HEADER denominator=" << arrangement.denominator
              << " cells=" << arrangement.cells.size()
              << " k4_profiles=" << profiles_by_k[4]
              << " k5_profiles=" << profiles_by_k[5]
              << " k6_profiles=" << profiles_by_k[6]
              << " total_profiles=" << profiles
              << " k4_max=" << max_by_k[4]
              << " k5_max=" << max_by_k[5]
              << " k6_max=" << max_by_k[6]
              << " k4_high=" << high_by_k[4]
              << " k5_high=" << high_by_k[5]
              << " k6_high=" << high_by_k[6]
              << " total_high=" << high.size() << '\n';
    for (std::size_t j=0;j<high.size();++j) {
        const HighBody& body=high[j];
        std::cout << "HIGH index=" << j
                  << " k=" << std::popcount(body.petals)
                  << " B={" << all_labels(body) << "}"
                  << " petals={" << o_labels(body.petals) << "}"
                  << " core={" << c_labels(body.core) << "}"
                  << " mass=" << body.mass
                  << " components=" << body.components
                  << " delta=" << decimal(body.delta)
                  << " cutoff=" << body.cutoff << '\n';
    }
    std::uint64_t literal_rows=0;
    int failures=0,equalities=0;
    bool closest_set=false;
    Literal closest{};
    std::size_t closest_index=0;
    int closest_r=0;
    for (std::size_t j=0;j<high.size();++j) {
        const HighBody& body=high[j];
        for (int r=590;r<body.cutoff;++r) {
            const Literal result=integrate(body,r);
            ++literal_rows;
            if (result.scaled_margin<0) ++failures;
            if (result.scaled_margin==0) ++equalities;
            if (!closest_set || result.scaled_margin*closest.denominator <
                                closest.scaled_margin*result.denominator) {
                closest_set=true;
                closest=result;
                closest_index=j;
                closest_r=r;
            }
            std::cout << "LITERAL index=" << j
                      << " r=" << r
                      << " numerator=" << result.numerator
                      << " denominator=" << result.denominator
                      << " scaled_margin=" << decimal(result.scaled_margin)
                      << '\n';
        }
    }
    std::cout << "SUMMARY high_bodies=" << high.size()
              << " literal_rows=" << literal_rows
              << " failures=" << failures
              << " equalities=" << equalities
              << " closest_index=" << closest_index
              << " closest_B={" << all_labels(high[closest_index]) << "}"
              << " closest_r=" << closest_r
              << " closest_numerator=" << closest.numerator
              << " closest_denominator=" << closest.denominator
              << " closest_scaled_margin=" << decimal(closest.scaled_margin)
              << '\n';
    return failures==0 ? 0 : 3;
}

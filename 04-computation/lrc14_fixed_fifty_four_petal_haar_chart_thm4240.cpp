// Exact fixed-50 fourth-petal Haar chart certificate.
//
// Starting from the THM-4234 universal triple {132,176,264}, this program
// tests each possible fourth petal.  A candidate has three newly exposed
// triple layers, with six labels chosen from C, and one quadruple layer, with
// five labels chosen from C.  The limiting pass records exact Haar surplus,
// circular component count, and the resulting THM-4170 tail cutoff.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using i64 = std::int64_t;
using u32 = std::uint32_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int,18> C = {
    8,16,40,42,80,84,85,88,95,120,126,143,145,168,193,240,252,286};
constexpr std::array<int,12> O = {
    10,15,20,30,60,63,132,170,176,190,264,290};
constexpr int Q=50;
constexpr std::array<int,3> STAR={6,8,10}; // 132,176,264 in O.
constexpr u32 ALL=(u32{1}<<18)-1;
constexpr std::size_t N=(std::size_t{1}<<18);

void require(bool ok,const char* message) {
    if (!ok) throw std::runtime_error(message);
}

i64 lcm_exact(i64 a,i64 b) {
    const i128 value=static_cast<i128>(a/std::gcd(a,b))*b;
    require(value<=std::numeric_limits<i64>::max(),"lcm overflow");
    return static_cast<i64>(value);
}

bool safe_mid(int speed,i64 left,i64 right,i64 denominator) {
    i128 residue=static_cast<i128>(speed)*(left+right)%
                 (static_cast<i128>(2)*denominator);
    if (residue<0) residue+=static_cast<i128>(2)*denominator;
    return static_cast<i128>(7)*residue>=denominator &&
           static_cast<i128>(7)*residue<=static_cast<i128>(13)*denominator;
}

std::string decimal(i128 value) {
    if (value==0) return "0";
    const bool negative=value<0;
    if (negative) value=-value;
    std::string result;
    while (value!=0) {
        result.push_back(static_cast<char>('0'+value%10));
        value/=10;
    }
    if (negative) result.push_back('-');
    std::reverse(result.begin(),result.end());
    return result;
}

std::string c_labels(u32 mask) {
    std::string result;
    for (int bit=0;bit<18;++bit) if (mask&(u32{1}<<bit)) {
        if (!result.empty()) result+=',';
        result+=std::to_string(C[bit]);
    }
    return result;
}

std::string o_labels(std::uint16_t mask) {
    std::string result;
    for (int bit=0;bit<12;++bit) if (mask&(std::uint16_t{1}<<bit)) {
        if (!result.empty()) result+=',';
        result+=std::to_string(O[bit]);
    }
    return result;
}

void zeta(std::vector<i64>& values) {
    for (int bit=0;bit<18;++bit)
        for (u32 mask=0;mask<(u32{1}<<18);++mask)
            if (mask&(u32{1}<<bit)) values[mask]+=values[mask^(u32{1}<<bit)];
}

struct Cell {
    i64 length;
    u32 c_failed;
    std::uint16_t o_safe;
    bool q_safe;
};

struct BaseArrangement {
    i64 denominator;
    std::vector<Cell> cells;
};

BaseArrangement make_base_arrangement() {
    std::vector<int> speeds(C.begin(),C.end());
    speeds.insert(speeds.end(),O.begin(),O.end());
    speeds.push_back(Q);
    std::sort(speeds.begin(),speeds.end());
    speeds.erase(std::unique(speeds.begin(),speeds.end()),speeds.end());
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
    std::vector<Cell> cells;
    cells.reserve(walls.size()-1);
    for (std::size_t i=0;i+1<walls.size();++i) {
        u32 failed=0;
        for (int bit=0;bit<18;++bit)
            if (!safe_mid(C[bit],walls[i],walls[i+1],denominator))
                failed|=u32{1}<<bit;
        std::uint16_t osafe=0;
        for (int bit=0;bit<12;++bit)
            if (safe_mid(O[bit],walls[i],walls[i+1],denominator))
                osafe|=static_cast<std::uint16_t>(1u<<bit);
        cells.push_back({walls[i+1]-walls[i],failed,osafe,
                         safe_mid(Q,walls[i],walls[i+1],denominator)});
    }
    return {denominator,std::move(cells)};
}

struct BaseStats {
    std::uint64_t profiles=0;
    int nonstrict=0;
    i128 min_delta=0;
    u32 min_body=0;
    int max_cutoff=0;
    u32 cutoff_body=0;
};

BaseStats base_stats(const BaseArrangement& arrangement,
                     std::uint16_t required,int core_size) {
    std::vector<i64> mass(N,0),starts(N,0),continuations(N,0);
    for (std::size_t i=0;i<arrangement.cells.size();++i) {
        const Cell& current=arrangement.cells[i];
        const Cell& previous=arrangement.cells[
            (i+arrangement.cells.size()-1)%arrangement.cells.size()];
        if (!current.q_safe || (current.o_safe&required)!=required) continue;
        mass[current.c_failed]+=current.length;
        starts[current.c_failed]+=1;
        if (previous.q_safe && (previous.o_safe&required)==required)
            continuations[current.c_failed|previous.c_failed]+=1;
    }
    zeta(mass); zeta(starts); zeta(continuations);
    BaseStats result;
    bool initialized=false;
    for (u32 body=ALL;;body=(body-1)&ALL) {
        if (std::popcount(body)==core_size) {
            ++result.profiles;
            const u32 allowed=ALL^body;
            const i64 components=starts[allowed]-continuations[allowed];
            require(components>0,"component count nonpositive");
            const i128 delta=static_cast<i128>(54)*mass[allowed]-
                             static_cast<i128>(4)*arrangement.denominator;
            if (!initialized || delta<result.min_delta) {
                initialized=true;
                result.min_delta=delta;
                result.min_body=body;
            }
            if (delta<=0) {
                ++result.nonstrict;
            } else {
                const i128 numerator=static_cast<i128>(54)*components*
                                     arrangement.denominator;
                const i128 divisor=static_cast<i128>(7)*delta;
                const int cutoff=static_cast<int>((numerator+divisor-1)/divisor);
                if (cutoff>result.max_cutoff) {
                    result.max_cutoff=cutoff;
                    result.cutoff_body=body;
                }
            }
        }
        if (body==0) break;
    }
    return result;
}

std::uint16_t bit(int index) {
    return static_cast<std::uint16_t>(1u<<index);
}

bool in_pool(int value) {
    return std::find(C.begin(),C.end(),value)!=C.end() ||
           std::find(O.begin(),O.end(),value)!=O.end();
}

} // namespace

int main(int argc,char** argv) {
    const BaseArrangement arrangement=make_base_arrangement();
    require(arrangement.denominator==91205797082400LL,"base denominator changed");
    require(arrangement.cells.size()==7213,"base cell count changed");

    if (argc==3 && std::string(argv[1])=="--all-petal-base") {
        const int petal_count=std::stoi(argv[2]);
        require(petal_count>=0 && petal_count<=9,"petal count outside 0..9");
        const int core_size=9-petal_count;
        std::uint64_t profiles=0;
        int lanes=0,strict=0,nonstrict=0,global_cutoff=0;
        bool best_set=false;
        i128 best_min_delta=0;
        std::uint16_t best_required=0;
        int best_cutoff=0;
        for (std::uint16_t required=0;required<(std::uint16_t{1}<<12);++required) {
            if (std::popcount(required)!=petal_count) continue;
            const BaseStats stats=base_stats(arrangement,required,core_size);
            ++lanes;
            profiles+=stats.profiles;
            global_cutoff=std::max(global_cutoff,stats.max_cutoff);
            if (stats.nonstrict==0) ++strict; else ++nonstrict;
            if (stats.nonstrict==0 &&
                (!best_set || stats.min_delta>best_min_delta)) {
                best_set=true;
                best_min_delta=stats.min_delta;
                best_required=required;
                best_cutoff=stats.max_cutoff;
            }
            std::cout << "PETAL_BASE k=" << petal_count
                      << " labels=" << o_labels(required)
                      << " profiles=" << stats.profiles
                      << " nonstrict=" << stats.nonstrict
                      << " min_delta=" << decimal(stats.min_delta)
                      << " min_body=" << c_labels(stats.min_body)
                      << " max_cutoff=" << stats.max_cutoff
                      << " cutoff_body=" << c_labels(stats.cutoff_body)
                      << '\n';
        }
        std::cout << "PETAL_BASE_SUMMARY k=" << petal_count
                  << " core_size=" << core_size
                  << " denominator=" << arrangement.denominator
                  << " cells=" << arrangement.cells.size()
                  << " lanes=" << lanes
                  << " profiles=" << profiles
                  << " strict=" << strict
                  << " nonstrict=" << nonstrict
                  << " global_max_cutoff=" << global_cutoff
                  << " best=" << o_labels(best_required)
                  << " best_min_delta=" << decimal(best_min_delta)
                  << " best_cutoff=" << best_cutoff << '\n';
        return nonstrict==0 ? 0 : 3;
    }

    if (argc==2 && std::string(argv[1])=="--all-quadruple-base") {
        std::uint64_t profiles=0;
        int quadruples=0,strict=0,nonstrict=0,global_cutoff=0;
        bool best_set=false;
        i128 best_min_delta=0;
        std::uint16_t best_required=0;
        int best_cutoff=0;
        for (int a=0;a<12;++a) for (int b=a+1;b<12;++b)
            for (int c=b+1;c<12;++c) for (int d=c+1;d<12;++d) {
                const std::uint16_t required=static_cast<std::uint16_t>(
                    bit(a)|bit(b)|bit(c)|bit(d));
                const BaseStats stats=base_stats(arrangement,required,5);
                require(stats.profiles==8568,"quadruple profile count changed");
                ++quadruples;
                profiles+=stats.profiles;
                global_cutoff=std::max(global_cutoff,stats.max_cutoff);
                if (stats.nonstrict==0) ++strict; else ++nonstrict;
                if (stats.nonstrict==0 &&
                    (!best_set || stats.min_delta>best_min_delta)) {
                    best_set=true;
                    best_min_delta=stats.min_delta;
                    best_required=required;
                    best_cutoff=stats.max_cutoff;
                }
                std::cout << "QUADRUPLE labels=" << o_labels(required)
                          << " profiles=" << stats.profiles
                          << " nonstrict=" << stats.nonstrict
                          << " min_delta=" << decimal(stats.min_delta)
                          << " min_body=" << c_labels(stats.min_body)
                          << " max_cutoff=" << stats.max_cutoff
                          << " cutoff_body=" << c_labels(stats.cutoff_body)
                          << '\n';
            }
        require(quadruples==495,"quadruple count changed");
        require(profiles==UINT64_C(4241160),"quadruple universe changed");
        std::cout << "QUADRUPLE_SUMMARY denominator=" << arrangement.denominator
                  << " cells=" << arrangement.cells.size()
                  << " quadruples=" << quadruples
                  << " profiles=" << profiles
                  << " strict=" << strict
                  << " nonstrict=" << nonstrict
                  << " global_max_cutoff=" << global_cutoff
                  << " best=" << o_labels(best_required)
                  << " best_min_delta=" << decimal(best_min_delta)
                  << " best_cutoff=" << best_cutoff << '\n';
        return nonstrict==0 ? 0 : 3;
    }

    const std::uint16_t star=static_cast<std::uint16_t>(
        bit(STAR[0])|bit(STAR[1])|bit(STAR[2]));
    std::uint64_t total_profiles=0;
    int candidates=0;
    bool best_set=false;
    i128 best_floor=0;
    int best_extra=-1,best_cutoff=0;
    for (int extra=0;extra<12;++extra) {
        if (star&bit(extra)) continue;
        ++candidates;
        i128 candidate_floor=0;
        bool floor_set=false;
        int candidate_cutoff=0;
        int candidate_nonstrict=0;
        for (int omitted=0;omitted<3;++omitted) {
            std::uint16_t required=bit(extra);
            for (int j=0;j<3;++j) if (j!=omitted) required|=bit(STAR[j]);
            const BaseStats stats=base_stats(arrangement,required,6);
            require(stats.profiles==18564,"triple profile count changed");
            total_profiles+=stats.profiles;
            candidate_cutoff=std::max(candidate_cutoff,stats.max_cutoff);
            candidate_nonstrict+=stats.nonstrict;
            if (!floor_set || stats.min_delta<candidate_floor) {
                floor_set=true;
                candidate_floor=stats.min_delta;
            }
            std::cout << "LAYER extra=" << O[extra]
                      << " kind=triple labels=" << o_labels(required)
                      << " profiles=" << stats.profiles
                      << " nonstrict=" << stats.nonstrict
                      << " min_delta=" << decimal(stats.min_delta)
                      << " min_body=" << c_labels(stats.min_body)
                      << " max_cutoff=" << stats.max_cutoff
                      << " cutoff_body=" << c_labels(stats.cutoff_body) << '\n';
        }
        const std::uint16_t required=static_cast<std::uint16_t>(star|bit(extra));
        const BaseStats stats=base_stats(arrangement,required,5);
        require(stats.profiles==8568,"quadruple profile count changed");
        total_profiles+=stats.profiles;
        candidate_cutoff=std::max(candidate_cutoff,stats.max_cutoff);
        candidate_nonstrict+=stats.nonstrict;
        if (!floor_set || stats.min_delta<candidate_floor)
            candidate_floor=stats.min_delta;
        std::cout << "LAYER extra=" << O[extra]
                  << " kind=quadruple labels=" << o_labels(required)
                  << " profiles=" << stats.profiles
                  << " nonstrict=" << stats.nonstrict
                  << " min_delta=" << decimal(stats.min_delta)
                  << " min_body=" << c_labels(stats.min_body)
                  << " max_cutoff=" << stats.max_cutoff
                  << " cutoff_body=" << c_labels(stats.cutoff_body) << '\n';
        std::cout << "CANDIDATE extra=" << O[extra]
                  << " profiles=64260"
                  << " nonstrict=" << candidate_nonstrict
                  << " min_delta=" << decimal(candidate_floor)
                  << " max_cutoff=" << candidate_cutoff << '\n';
        if (candidate_nonstrict==0 &&
            (!best_set || candidate_floor>best_floor ||
             (candidate_floor==best_floor && candidate_cutoff<best_cutoff))) {
            best_set=true;
            best_floor=candidate_floor;
            best_extra=extra;
            best_cutoff=candidate_cutoff;
        }
    }
    require(candidates==9,"candidate count changed");
    require(total_profiles==UINT64_C(578340),"total profile count changed");
    std::cout << "SUMMARY denominator=" << arrangement.denominator
              << " cells=" << arrangement.cells.size()
              << " candidates=" << candidates
              << " profiles=" << total_profiles
              << " best_extra=" << O[best_extra]
              << " best_min_delta=" << decimal(best_floor)
              << " best_max_cutoff=" << best_cutoff << '\n';

    if (!(argc==2 && std::string(argv[1])=="--selected-four-petal")) return 0;

    constexpr int EXTRA=5; // 63 in O.
    struct Lane {
        std::uint16_t required;
        int core_size;
        int cutoff;
        std::uint64_t rows=0;
        std::uint64_t checks=0;
        std::uint64_t failures=0;
        std::uint64_t equalities=0;
        bool min_set=false;
        i128 min_delta=0;
        i64 min_denominator=0;
        int min_r=0;
        u32 min_body=0;
    };
    std::array<Lane,4> lanes{};
    for (int omitted=0;omitted<3;++omitted) {
        std::uint16_t required=bit(EXTRA);
        for (int j=0;j<3;++j) if (j!=omitted) required|=bit(STAR[j]);
        const BaseStats stats=base_stats(arrangement,required,6);
        lanes[omitted].required=required;
        lanes[omitted].core_size=6;
        lanes[omitted].cutoff=stats.max_cutoff;
    }
    lanes[3].required=static_cast<std::uint16_t>(star|bit(EXTRA));
    lanes[3].core_size=5;
    lanes[3].cutoff=base_stats(arrangement,lanes[3].required,5).max_cutoff;
    require(lanes[0].cutoff==475 && lanes[1].cutoff==427 &&
            lanes[2].cutoff==477 && lanes[3].cutoff==422,
            "selected layer cutoffs changed");

    const int literal_stop=std::max({lanes[0].cutoff,lanes[1].cutoff,
                                     lanes[2].cutoff,lanes[3].cutoff});
    std::vector<int> base_speeds(C.begin(),C.end());
    base_speeds.insert(base_speeds.end(),O.begin(),O.end());
    base_speeds.push_back(Q);
    std::sort(base_speeds.begin(),base_speeds.end());
    base_speeds.erase(std::unique(base_speeds.begin(),base_speeds.end()),
                      base_speeds.end());

    for (int outsider=1;outsider<literal_stop;++outsider) {
        if (outsider==Q || in_pool(outsider)) continue;
        std::vector<int> row_speeds=base_speeds;
        row_speeds.push_back(outsider);
        std::sort(row_speeds.begin(),row_speeds.end());
        row_speeds.erase(std::unique(row_speeds.begin(),row_speeds.end()),
                         row_speeds.end());
        i64 row_denominator=1;
        for (int speed:row_speeds)
            row_denominator=lcm_exact(row_denominator,14LL*speed);
        std::vector<i64> walls{0,row_denominator};
        for (int speed:row_speeds) {
            const i64 unit=row_denominator/(14LL*speed);
            for (int tooth=0;tooth<speed;++tooth) {
                walls.push_back((14LL*tooth+1)*unit);
                walls.push_back((14LL*tooth+13)*unit);
            }
        }
        std::sort(walls.begin(),walls.end());
        walls.erase(std::unique(walls.begin(),walls.end()),walls.end());
        std::vector<Cell> row_cells;
        row_cells.reserve(walls.size()-1);
        for (std::size_t i=0;i+1<walls.size();++i) {
            u32 failed=0;
            for (int bit_index=0;bit_index<18;++bit_index)
                if (!safe_mid(C[bit_index],walls[i],walls[i+1],row_denominator))
                    failed|=u32{1}<<bit_index;
            std::uint16_t osafe=0;
            for (int bit_index=0;bit_index<12;++bit_index)
                if (safe_mid(O[bit_index],walls[i],walls[i+1],row_denominator))
                    osafe|=static_cast<std::uint16_t>(1u<<bit_index);
            const bool fixed_safe=
                safe_mid(Q,walls[i],walls[i+1],row_denominator) &&
                safe_mid(outsider,walls[i],walls[i+1],row_denominator);
            row_cells.push_back(
                {walls[i+1]-walls[i],failed,osafe,fixed_safe});
        }

        std::array<i128,4> row_min_delta{};
        std::array<u32,4> row_min_body{};
        std::array<std::uint64_t,4> row_failures{},row_equalities{};
#pragma omp parallel for schedule(static)
        for (int lane_index=0;lane_index<4;++lane_index) {
            const Lane& lane=lanes[lane_index];
            if (outsider>=lane.cutoff) continue;
            std::vector<i64> mass(N,0);
            for (const Cell& cell:row_cells)
                if (cell.q_safe &&
                    (cell.o_safe&lane.required)==lane.required)
                    mass[cell.c_failed]+=cell.length;
            zeta(mass);
            bool initialized=false;
            i128 least=0;
            u32 least_body=0;
            std::uint64_t bad=0,equal=0;
            for (u32 body=ALL;;body=(body-1)&ALL) {
                if (std::popcount(body)==lane.core_size) {
                    const i128 delta=static_cast<i128>(63)*mass[ALL^body]-
                                     static_cast<i128>(4)*row_denominator;
                    if (!initialized || delta<least) {
                        initialized=true;
                        least=delta;
                        least_body=body;
                    }
                    bad+=delta<0;
                    equal+=delta==0;
                }
                if (body==0) break;
            }
            row_min_delta[lane_index]=least;
            row_min_body[lane_index]=least_body;
            row_failures[lane_index]=bad;
            row_equalities[lane_index]=equal;
        }

        for (int lane_index=0;lane_index<4;++lane_index) {
            Lane& lane=lanes[lane_index];
            if (outsider>=lane.cutoff) continue;
            ++lane.rows;
            lane.checks+=lane.core_size==6 ? UINT64_C(18564) : UINT64_C(8568);
            lane.failures+=row_failures[lane_index];
            lane.equalities+=row_equalities[lane_index];
            if (!lane.min_set ||
                row_min_delta[lane_index]*lane.min_denominator <
                    lane.min_delta*row_denominator) {
                lane.min_set=true;
                lane.min_delta=row_min_delta[lane_index];
                lane.min_denominator=row_denominator;
                lane.min_r=outsider;
                lane.min_body=row_min_body[lane_index];
            }
        }
    }

    std::uint64_t total_rows=0,total_checks=0,total_failures=0,total_equalities=0;
    for (const Lane& lane:lanes) {
        total_rows+=lane.rows;
        total_checks+=lane.checks;
        total_failures+=lane.failures;
        total_equalities+=lane.equalities;
        std::cout << "LITERAL kind="
                  << (lane.core_size==6 ? "triple" : "quadruple")
                  << " labels=" << o_labels(lane.required)
                  << " core_size=" << lane.core_size
                  << " cutoff=" << lane.cutoff
                  << " rows=" << lane.rows
                  << " checks=" << lane.checks
                  << " failures=" << lane.failures
                  << " equalities=" << lane.equalities
                  << " min_delta=" << decimal(lane.min_delta)
                  << " min_denominator=" << lane.min_denominator
                  << " min_r=" << lane.min_r
                  << " min_body=" << c_labels(lane.min_body) << '\n';
    }
    std::cout << "LITERAL_SUMMARY labels=63,132,176,264"
              << " rows=" << total_rows
              << " checks=" << total_checks
              << " failures=" << total_failures
              << " equalities=" << total_equalities << '\n';
    return total_failures==0 ? 0 : 2;
}

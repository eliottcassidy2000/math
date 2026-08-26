// Exact fixed-50 pair and triple Haar chart certificate.
//
// For every unordered pair {p,q} of the twelve labels in P\C, and every
// L in binom(C,7), compute the exact mass and component count of
// V=G_(L union {p,q,50}).  THM-4170 can close the extra outsider r only when
// (6/7)mu(V)>4/63; the corresponding sufficient cutoff is also recorded.
// It then scans every admissible outsider below the sufficient cutoff.  The
// finite scan plus THM-4170 closes an edge exactly when no failure occurs.

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

constexpr std::array<int, 18> C = {
    8,16,40,42,80,84,85,88,95,120,126,143,145,168,193,240,252,286};
constexpr std::array<int, 12> O = {
    10,15,20,30,60,63,132,170,176,190,264,290};
constexpr int Q = 50;

void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

i64 lcm_exact(i64 a, i64 b) {
    const i128 value = static_cast<i128>(a / std::gcd(a,b)) * b;
    require(value <= std::numeric_limits<i64>::max(), "lcm overflow");
    return static_cast<i64>(value);
}

bool safe_mid(int speed, i64 left, i64 right, i64 denominator) {
    i128 residue = static_cast<i128>(speed) * (left + right) %
                   (static_cast<i128>(2) * denominator);
    if (residue < 0) residue += static_cast<i128>(2) * denominator;
    return static_cast<i128>(7) * residue >= denominator &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * denominator;
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    bool negative = value < 0;
    if (negative) value = -value;
    std::string result;
    while (value != 0) {
        result.push_back(static_cast<char>('0' + value % 10));
        value /= 10;
    }
    if (negative) result.push_back('-');
    std::reverse(result.begin(), result.end());
    return result;
}

std::string c_labels(u32 mask) {
    std::string result;
    for (int bit=0; bit<18; ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        if (!result.empty()) result += ',';
        result += std::to_string(C[bit]);
    }
    return result;
}

bool in_pool(int value) {
    return std::find(C.begin(),C.end(),value)!=C.end() ||
           std::find(O.begin(),O.end(),value)!=O.end();
}

struct Cell {
    i64 length;
    u32 c_failed;
    std::uint16_t o_safe;
    bool q_safe;
};

void zeta(std::vector<i64>& values) {
    for (int bit=0; bit<18; ++bit) {
        for (u32 mask=0; mask<(u32{1}<<18); ++mask) {
            if (mask & (u32{1}<<bit)) values[mask] += values[mask^(u32{1}<<bit)];
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    std::vector<int> speeds(C.begin(), C.end());
    speeds.insert(speeds.end(), O.begin(), O.end());
    speeds.push_back(Q);
    std::sort(speeds.begin(), speeds.end());
    speeds.erase(std::unique(speeds.begin(), speeds.end()), speeds.end());

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
        for (int bit=0;bit<18;++bit) {
            if (!safe_mid(C[bit],walls[i],walls[i+1],denominator))
                failed|=u32{1}<<bit;
        }
        std::uint16_t osafe=0;
        for (int bit=0;bit<12;++bit) {
            if (safe_mid(O[bit],walls[i],walls[i+1],denominator))
                osafe|=static_cast<std::uint16_t>(1u<<bit);
        }
        cells.push_back({walls[i+1]-walls[i],failed,osafe,
                         safe_mid(Q,walls[i],walls[i+1],denominator)});
    }

    constexpr u32 ALL=(u32{1}<<18)-1;
    constexpr std::size_t N=(std::size_t{1}<<18);
    int strict_pairs=0;
    int nonstrict_pairs=0;
    std::uint64_t profiles=0;
    int global_max_cutoff=0;
    std::array<int,66> pair_a{},pair_b{},cutoffs{};
    int pair_index=0;
    for (int a=0;a<12;++a) for (int b=a+1;b<12;++b,++pair_index) {
        pair_a[pair_index]=a;
        pair_b[pair_index]=b;
        const std::uint16_t pair=static_cast<std::uint16_t>((1u<<a)|(1u<<b));
        std::vector<i64> mass(N,0),starts(N,0),continuations(N,0);
        for (std::size_t i=0;i<cells.size();++i) {
            const Cell& current=cells[i];
            const Cell& previous=cells[(i+cells.size()-1)%cells.size()];
            if (!current.q_safe || (current.o_safe&pair)!=pair) continue;
            mass[current.c_failed]+=current.length;
            starts[current.c_failed]+=1;
            if (previous.q_safe && (previous.o_safe&pair)==pair)
                continuations[current.c_failed|previous.c_failed]+=1;
        }
        zeta(mass); zeta(starts); zeta(continuations);

        i128 min_delta=0;
        u32 min_body=0;
        int nonstrict=0;
        int max_cutoff=0;
        u32 cutoff_body=0;
        for (u32 body=ALL;;body=(body-1)&ALL) {
            if (std::popcount(body)==7) {
                ++profiles;
                const u32 allowed=ALL^body;
                const i64 components=starts[allowed]-continuations[allowed];
                require(components>0,"component count nonpositive");
                const i128 delta=static_cast<i128>(54)*mass[allowed]-
                                 static_cast<i128>(4)*denominator;
                if (min_body==0 || delta<min_delta) {
                    min_delta=delta; min_body=body;
                }
                if (delta<=0) {
                    ++nonstrict;
                } else {
                    const i128 numerator=static_cast<i128>(54)*components*denominator;
                    const i128 divisor=static_cast<i128>(7)*delta;
                    const int cutoff=static_cast<int>((numerator+divisor-1)/divisor);
                    if (cutoff>max_cutoff) {max_cutoff=cutoff;cutoff_body=body;}
                }
            }
            if (body==0) break;
        }
        if (nonstrict==0) ++strict_pairs; else ++nonstrict_pairs;
        cutoffs[pair_index]=max_cutoff;
        global_max_cutoff=std::max(global_max_cutoff,max_cutoff);
        std::cout << "PAIR p=" << O[a] << " q=" << O[b]
                  << " nonstrict=" << nonstrict
                  << " min_delta=" << decimal(min_delta)
                  << " min_body=" << c_labels(min_body)
                  << " max_cutoff=" << max_cutoff
                  << " cutoff_body=" << c_labels(cutoff_body) << '\n';
    }
    require(denominator==91205797082400LL,"base denominator changed");
    require(profiles==UINT64_C(2100384),"profile count changed");
    std::cout << "SUMMARY denominator=" << denominator
              << " cells=" << cells.size()
              << " pairs=66 profiles=" << profiles
              << " strict_pairs=" << strict_pairs
              << " nonstrict_pairs=" << nonstrict_pairs
              << " global_max_cutoff=" << global_max_cutoff << '\n';

    if (argc==2 && std::string(argv[1])=="--triple-base") {
        std::uint64_t triple_profiles=0;
        int strict_triples=0,nonstrict_triples=0;
        i128 best_min_delta=0;
        int best_a=0,best_b=0,best_c=0,best_cutoff=0;
        bool best_set=false;
        for (int a=0;a<12;++a) for (int b=a+1;b<12;++b)
            for (int c=b+1;c<12;++c) {
                const std::uint16_t required=static_cast<std::uint16_t>(
                    (1u<<a)|(1u<<b)|(1u<<c));
                std::vector<i64> mass(N,0),starts(N,0),continuations(N,0);
                for (std::size_t i=0;i<cells.size();++i) {
                    const Cell& current=cells[i];
                    const Cell& previous=cells[(i+cells.size()-1)%cells.size()];
                    if (!current.q_safe || (current.o_safe&required)!=required) continue;
                    mass[current.c_failed]+=current.length;
                    starts[current.c_failed]+=1;
                    if (previous.q_safe && (previous.o_safe&required)==required)
                        continuations[current.c_failed|previous.c_failed]+=1;
                }
                zeta(mass);zeta(starts);zeta(continuations);
                bool initialized=false;
                i128 least=0;
                u32 least_body=0,cutoff_body=0;
                int nonstrict=0,max_cutoff=0;
                for (u32 body=ALL;;body=(body-1)&ALL) {
                    if (std::popcount(body)==6) {
                        ++triple_profiles;
                        const u32 allowed=ALL^body;
                        const i64 components=starts[allowed]-continuations[allowed];
                        require(components>0,"triple component count nonpositive");
                        const i128 delta=static_cast<i128>(54)*mass[allowed]-
                                         static_cast<i128>(4)*denominator;
                        if (!initialized || delta<least) {
                            initialized=true;least=delta;least_body=body;
                        }
                        if (delta<=0) ++nonstrict;
                        else {
                            const i128 numerator=static_cast<i128>(54)*components*denominator;
                            const i128 divisor=static_cast<i128>(7)*delta;
                            const int cutoff=static_cast<int>((numerator+divisor-1)/divisor);
                            if (cutoff>max_cutoff) {max_cutoff=cutoff;cutoff_body=body;}
                        }
                    }
                    if (body==0) break;
                }
                if (nonstrict==0) ++strict_triples;else ++nonstrict_triples;
                if (nonstrict==0 && (!best_set || least>best_min_delta)) {
                    best_set=true;best_min_delta=least;
                    best_a=a;best_b=b;best_c=c;best_cutoff=max_cutoff;
                }
                std::cout << "TRIPLE p=" << O[a] << " q=" << O[b] << " u=" << O[c]
                          << " nonstrict=" << nonstrict
                          << " min_delta=" << decimal(least)
                          << " min_body=" << c_labels(least_body)
                          << " max_cutoff=" << max_cutoff
                          << " cutoff_body=" << c_labels(cutoff_body) << '\n';
            }
        require(triple_profiles==UINT64_C(4084080),"triple profile count changed");
        std::cout << "TRIPLE_SUMMARY profiles=" << triple_profiles
                  << " strict=" << strict_triples
                  << " nonstrict=" << nonstrict_triples
                  << " best=" << O[best_a] << ',' << O[best_b] << ',' << O[best_c]
                  << " best_min_delta=" << decimal(best_min_delta)
                  << " best_cutoff=" << best_cutoff << '\n';
        return 0;
    }

    if (argc==2 && std::string(argv[1])=="--selected-triple") {
        constexpr int A=6;   // 132
        constexpr int B=8;   // 176
        constexpr int CC=10; // 264
        const std::uint16_t required=static_cast<std::uint16_t>(
            (1u<<A)|(1u<<B)|(1u<<CC));
        std::vector<i64> base_mass(N,0),base_starts(N,0),base_continuations(N,0);
        for (std::size_t i=0;i<cells.size();++i) {
            const Cell& current=cells[i];
            const Cell& previous=cells[(i+cells.size()-1)%cells.size()];
            if (!current.q_safe || (current.o_safe&required)!=required) continue;
            base_mass[current.c_failed]+=current.length;
            base_starts[current.c_failed]+=1;
            if (previous.q_safe && (previous.o_safe&required)==required)
                base_continuations[current.c_failed|previous.c_failed]+=1;
        }
        zeta(base_mass);zeta(base_starts);zeta(base_continuations);
        int cutoff=0,nonstrict=0;
        i128 limiting_min=0;
        bool limiting_set=false;
        for (u32 body=ALL;;body=(body-1)&ALL) {
            if (std::popcount(body)==6) {
                const u32 allowed=ALL^body;
                const i64 components=base_starts[allowed]-base_continuations[allowed];
                const i128 delta=static_cast<i128>(54)*base_mass[allowed]-
                                 static_cast<i128>(4)*denominator;
                if (!limiting_set || delta<limiting_min) {
                    limiting_set=true;limiting_min=delta;
                }
                if (delta<=0) ++nonstrict;
                else {
                    const i128 numerator=static_cast<i128>(54)*components*denominator;
                    const i128 divisor=static_cast<i128>(7)*delta;
                    cutoff=std::max(cutoff,static_cast<int>((numerator+divisor-1)/divisor));
                }
            }
            if (body==0) break;
        }
        require(nonstrict==0 && cutoff==470,"selected triple base changed");

        std::uint64_t rows=0,checks=0,equalities=0,failures=0;
        bool min_set=false;
        i128 literal_min=0;
        i64 literal_min_denominator=0;
        int literal_min_r=0;
        u32 literal_min_body=0;
        for (int outsider=1;outsider<cutoff;++outsider) {
            if (outsider==Q || in_pool(outsider)) continue;
            std::vector<int> row_speeds=speeds;
            row_speeds.push_back(outsider);
            std::sort(row_speeds.begin(),row_speeds.end());
            row_speeds.erase(std::unique(row_speeds.begin(),row_speeds.end()),row_speeds.end());
            i64 row_denominator=1;
            for (int speed:row_speeds) row_denominator=lcm_exact(row_denominator,14LL*speed);
            std::vector<i64> row_walls{0,row_denominator};
            for (int speed:row_speeds) {
                const i64 unit=row_denominator/(14LL*speed);
                for (int tooth=0;tooth<speed;++tooth) {
                    row_walls.push_back((14LL*tooth+1)*unit);
                    row_walls.push_back((14LL*tooth+13)*unit);
                }
            }
            std::sort(row_walls.begin(),row_walls.end());
            row_walls.erase(std::unique(row_walls.begin(),row_walls.end()),row_walls.end());
            std::vector<i64> mass(N,0);
            for (std::size_t i=0;i+1<row_walls.size();++i) {
                if (!safe_mid(Q,row_walls[i],row_walls[i+1],row_denominator) ||
                    !safe_mid(outsider,row_walls[i],row_walls[i+1],row_denominator) ||
                    !safe_mid(O[A],row_walls[i],row_walls[i+1],row_denominator) ||
                    !safe_mid(O[B],row_walls[i],row_walls[i+1],row_denominator) ||
                    !safe_mid(O[CC],row_walls[i],row_walls[i+1],row_denominator)) continue;
                u32 failed=0;
                for (int bit=0;bit<18;++bit) {
                    if (!safe_mid(C[bit],row_walls[i],row_walls[i+1],row_denominator))
                        failed|=u32{1}<<bit;
                }
                mass[failed]+=row_walls[i+1]-row_walls[i];
            }
            zeta(mass);
            bool row_set=false;
            i128 row_min=0;
            u32 row_body=0;
            for (u32 body=ALL;;body=(body-1)&ALL) {
                if (std::popcount(body)==6) {
                    ++checks;
                    const i128 delta=static_cast<i128>(63)*mass[ALL^body]-
                                     static_cast<i128>(4)*row_denominator;
                    if (!row_set || delta<row_min) {row_set=true;row_min=delta;row_body=body;}
                    equalities+=delta==0;
                    failures+=delta<0;
                }
                if (body==0) break;
            }
            ++rows;
            if (!min_set || row_min*literal_min_denominator<literal_min*row_denominator) {
                min_set=true;literal_min=row_min;literal_min_denominator=row_denominator;
                literal_min_r=outsider;literal_min_body=row_body;
            }
        }
        std::cout << "SELECTED_TRIPLE labels=132,176,264"
                  << " limiting_profiles=18564"
                  << " limiting_min_delta=" << decimal(limiting_min)
                  << " cutoff=" << cutoff
                  << " literal_rows=" << rows
                  << " literal_checks=" << checks
                  << " failures=" << failures
                  << " equalities=" << equalities
                  << " min_delta=" << decimal(literal_min)
                  << " min_denominator=" << literal_min_denominator
                  << " min_r=" << literal_min_r
                  << " min_body=" << c_labels(literal_min_body) << '\n';
        return failures==0 ? 0 : 3;
    }

    std::array<std::uint64_t,66> literal_rows{},literal_checks{},equalities{},failures{};
    std::array<bool,66> min_set{};
    std::array<i128,66> min_delta{};
    std::array<i64,66> min_denominator{};
    std::array<int,66> min_r{};
    std::array<u32,66> min_body{};

    for (int outsider=1;outsider<global_max_cutoff;++outsider) {
        if (outsider==Q || in_pool(outsider)) continue;
        bool active=false;
        for (int lane=0;lane<66;++lane) active=active || outsider<cutoffs[lane];
        if (!active) continue;

        std::vector<int> row_speeds=speeds;
        row_speeds.push_back(outsider);
        std::sort(row_speeds.begin(),row_speeds.end());
        row_speeds.erase(std::unique(row_speeds.begin(),row_speeds.end()),row_speeds.end());
        i64 row_denominator=1;
        for (int speed:row_speeds) row_denominator=lcm_exact(row_denominator,14LL*speed);
        std::vector<i64> row_walls{0,row_denominator};
        for (int speed:row_speeds) {
            const i64 unit=row_denominator/(14LL*speed);
            for (int tooth=0;tooth<speed;++tooth) {
                row_walls.push_back((14LL*tooth+1)*unit);
                row_walls.push_back((14LL*tooth+13)*unit);
            }
        }
        std::sort(row_walls.begin(),row_walls.end());
        row_walls.erase(std::unique(row_walls.begin(),row_walls.end()),row_walls.end());
        std::vector<Cell> row_cells;
        row_cells.reserve(row_walls.size()-1);
        for (std::size_t i=0;i+1<row_walls.size();++i) {
            u32 failed=0;
            for (int bit=0;bit<18;++bit) {
                if (!safe_mid(C[bit],row_walls[i],row_walls[i+1],row_denominator))
                    failed|=u32{1}<<bit;
            }
            std::uint16_t osafe=0;
            for (int bit=0;bit<12;++bit) {
                if (safe_mid(O[bit],row_walls[i],row_walls[i+1],row_denominator))
                    osafe|=static_cast<std::uint16_t>(1u<<bit);
            }
            const bool both=safe_mid(Q,row_walls[i],row_walls[i+1],row_denominator) &&
                            safe_mid(outsider,row_walls[i],row_walls[i+1],row_denominator);
            row_cells.push_back({row_walls[i+1]-row_walls[i],failed,osafe,both});
        }

        std::array<i128,66> row_min_delta{};
        std::array<u32,66> row_min_body{};
        std::array<std::uint64_t,66> row_equalities{};
        std::array<std::uint64_t,66> row_failures{};

#pragma omp parallel for schedule(dynamic)
        for (int lane=0;lane<66;++lane) {
            if (outsider>=cutoffs[lane]) continue;
            const std::uint16_t pair=static_cast<std::uint16_t>(
                (1u<<pair_a[lane])|(1u<<pair_b[lane]));
            std::vector<i64> mass(N,0);
            for (const Cell& cell:row_cells) {
                if (cell.q_safe && (cell.o_safe&pair)==pair)
                    mass[cell.c_failed]+=cell.length;
            }
            zeta(mass);
            bool initialized=false;
            i128 least=0;
            u32 least_body=0;
            std::uint64_t eq=0,bad=0;
            for (u32 body=ALL;;body=(body-1)&ALL) {
                if (std::popcount(body)==7) {
                    const i128 delta=static_cast<i128>(63)*mass[ALL^body]-
                                     static_cast<i128>(4)*row_denominator;
                    if (!initialized || delta<least) {
                        initialized=true; least=delta; least_body=body;
                    }
                    eq+=delta==0;
                    bad+=delta<0;
                }
                if (body==0) break;
            }
            row_min_delta[lane]=least;
            row_min_body[lane]=least_body;
            row_equalities[lane]=eq;
            row_failures[lane]=bad;
        }

        for (int lane=0;lane<66;++lane) {
            if (outsider>=cutoffs[lane]) continue;
            ++literal_rows[lane];
            literal_checks[lane]+=UINT64_C(31824);
            equalities[lane]+=row_equalities[lane];
            failures[lane]+=row_failures[lane];
            if (!min_set[lane] ||
                row_min_delta[lane]*min_denominator[lane] <
                    min_delta[lane]*row_denominator) {
                min_set[lane]=true;
                min_delta[lane]=row_min_delta[lane];
                min_denominator[lane]=row_denominator;
                min_r[lane]=outsider;
                min_body[lane]=row_min_body[lane];
            }
        }
        if (outsider%50==0)
            std::cout << "PROGRESS outsider=" << outsider << '\n' << std::flush;
    }

    std::uint64_t total_rows=0,total_checks=0,total_equalities=0,total_failures=0;
    for (int lane=0;lane<66;++lane) {
        total_rows+=literal_rows[lane]; total_checks+=literal_checks[lane];
        total_equalities+=equalities[lane]; total_failures+=failures[lane];
        std::cout << "LITERAL p=" << O[pair_a[lane]]
                  << " q=" << O[pair_b[lane]]
                  << " rows=" << literal_rows[lane]
                  << " checks=" << literal_checks[lane]
                  << " failures=" << failures[lane]
                  << " equalities=" << equalities[lane]
                  << " min_delta=" << decimal(min_delta[lane])
                  << " min_denominator=" << min_denominator[lane]
                  << " min_r=" << min_r[lane]
                  << " min_body=" << c_labels(min_body[lane]) << '\n';
    }
    std::cout << "LITERAL_SUMMARY rows=" << total_rows
              << " checks=" << total_checks
              << " failures=" << total_failures
              << " equalities=" << total_equalities << '\n';
    return total_failures==0 ? 0 : 2;
}

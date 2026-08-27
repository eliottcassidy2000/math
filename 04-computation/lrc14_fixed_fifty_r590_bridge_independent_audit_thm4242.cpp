// Independent exact audit for the fixed-50 r>=590 direct tail.
//
// Scope: B is a nine-subset of P=C union O with |B intersect O| in {4,5,6}.
// The analytic pass is deliberately dual to the petal-first primary engine:
// it fixes the C-core first, scatters cells by the complete O-safe mask, and
// uses an O-superset zeta transform.  It independently recovers every body
// whose sufficient component-discrepancy cutoff exceeds 590.  Each recovered
// body is then replayed directly, and every 590 <= r < cutoff is integrated by
// a grouped endpoint-event sweep and by an independent literal midpoint pass.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using i64 = std::int64_t;
using u32 = std::uint32_t;
using u16 = std::uint16_t;
using i128 = __int128_t;

namespace {

constexpr std::array<int,18> C = {
    8,16,40,42,80,84,85,88,95,120,126,143,145,168,193,240,252,286};
constexpr std::array<int,12> O = {
    10,15,20,30,60,63,132,170,176,190,264,290};
constexpr int Q = 50;
constexpr int TAIL = 590;
constexpr u32 C_ALL = (u32{1} << 18) - 1;
constexpr u16 O_ALL = (u16{1} << 12) - 1;
constexpr std::size_t O_N = std::size_t{1} << 12;

void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}

i64 lcm_exact(i64 a, i64 b) {
    const i128 value = static_cast<i128>(a / std::gcd(a,b)) * b;
    require(value <= std::numeric_limits<i64>::max(), "lcm overflow");
    return static_cast<i64>(value);
}

std::string decimal(i128 value) {
    if (value == 0) return "0";
    const bool negative = value < 0;
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

std::string mask_labels(u32 core, u16 petals) {
    std::vector<int> labels;
    for (int bit=0; bit<18; ++bit)
        if (core & (u32{1} << bit)) labels.push_back(C[bit]);
    for (int bit=0; bit<12; ++bit)
        if (petals & (u16{1} << bit)) labels.push_back(O[bit]);
    std::sort(labels.begin(),labels.end());
    std::string result;
    for (int label: labels) {
        if (!result.empty()) result += ',';
        result += std::to_string(label);
    }
    return result;
}

std::vector<int> body_speeds(u32 core, u16 petals, bool include_q=true) {
    std::vector<int> speeds;
    for (int bit=0; bit<18; ++bit)
        if (core & (u32{1} << bit)) speeds.push_back(C[bit]);
    for (int bit=0; bit<12; ++bit)
        if (petals & (u16{1} << bit)) speeds.push_back(O[bit]);
    if (include_q) speeds.push_back(Q);
    std::sort(speeds.begin(),speeds.end());
    require(std::adjacent_find(speeds.begin(),speeds.end()) == speeds.end(),
            "duplicate speed in body");
    return speeds;
}

bool safe_mid(int speed, i64 left, i64 right, i64 denominator) {
    i128 residue = static_cast<i128>(speed) * (left + right) %
                   (static_cast<i128>(2) * denominator);
    if (residue < 0) residue += static_cast<i128>(2) * denominator;
    return static_cast<i128>(7) * residue >= denominator &&
           static_cast<i128>(7) * residue <= static_cast<i128>(13) * denominator;
}

struct Cell {
    i64 left;
    i64 right;
    u32 c_safe;
    u16 o_safe;
    bool q_safe;
};

struct MasterArrangement {
    i64 denominator;
    std::vector<Cell> cells;
};

MasterArrangement make_master() {
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
        for (int tooth=0; tooth<speed; ++tooth) {
            walls.push_back((14LL*tooth+1)*unit);
            walls.push_back((14LL*tooth+13)*unit);
        }
    }
    std::sort(walls.begin(),walls.end());
    walls.erase(std::unique(walls.begin(),walls.end()),walls.end());
    std::vector<Cell> cells;
    cells.reserve(walls.size()-1);
    for (std::size_t i=0;i+1<walls.size();++i) {
        u32 c_safe=0;
        u16 o_safe=0;
        for (int bit=0; bit<18; ++bit)
            if (safe_mid(C[bit],walls[i],walls[i+1],denominator))
                c_safe |= u32{1} << bit;
        for (int bit=0; bit<12; ++bit)
            if (safe_mid(O[bit],walls[i],walls[i+1],denominator))
                o_safe |= u16{1} << bit;
        cells.push_back({walls[i],walls[i+1],c_safe,o_safe,
                         safe_mid(Q,walls[i],walls[i+1],denominator)});
    }
    return {denominator,std::move(cells)};
}

void superset_zeta(std::array<i64,O_N>& values) {
    for (int bit=0; bit<12; ++bit)
        for (u16 mask=0; mask<=O_ALL; ++mask)
            if ((mask & (u16{1}<<bit)) == 0)
                values[mask] += values[mask | (u16{1}<<bit)];
}

struct Exceptional {
    int k;
    u32 core;
    u16 petals;
    i64 mass;
    i64 components;
    i128 delta;
    int cutoff;
};

struct AnalyticSummary {
    std::uint64_t cores=0;
    std::uint64_t profiles=0;
    std::uint64_t nonstrict=0;
    int max_cutoff=0;
    u32 max_core=0;
    u16 max_petals=0;
    i64 max_mass=0;
    i64 max_components=0;
    i128 max_delta=0;
    std::vector<Exceptional> exceptional;
};

AnalyticSummary enumerate_core_first(const MasterArrangement& master, int k) {
    const int core_size=9-k;
    AnalyticSummary summary;
    for (u32 core=0; core<=C_ALL; ++core) {
        if (std::popcount(core)!=core_size) continue;
        ++summary.cores;
        std::array<i64,O_N> mass{};
        std::array<i64,O_N> starts{};
        std::array<i64,O_N> continuations{};
        for (std::size_t i=0;i<master.cells.size();++i) {
            const Cell& current=master.cells[i];
            const Cell& previous=master.cells[
                (i+master.cells.size()-1)%master.cells.size()];
            const bool current_core=current.q_safe &&
                                    (current.c_safe&core)==core;
            if (!current_core) continue;
            mass[current.o_safe] += current.right-current.left;
            starts[current.o_safe] += 1;
            const bool previous_core=previous.q_safe &&
                                     (previous.c_safe&core)==core;
            if (previous_core)
                continuations[current.o_safe & previous.o_safe] += 1;
        }
        superset_zeta(mass);
        superset_zeta(starts);
        superset_zeta(continuations);
        for (u16 petals=0; petals<=O_ALL; ++petals) {
            if (std::popcount(petals)!=k) continue;
            ++summary.profiles;
            const i64 components=starts[petals]-continuations[petals];
            require(components>0,"nonpositive component count");
            const i128 delta=static_cast<i128>(54)*mass[petals]-
                             static_cast<i128>(4)*master.denominator;
            if (delta<=0) {
                ++summary.nonstrict;
                continue;
            }
            const i128 numerator=static_cast<i128>(54)*components*
                                 master.denominator;
            const i128 divisor=static_cast<i128>(7)*delta;
            const int cutoff=static_cast<int>((numerator+divisor-1)/divisor);
            if (cutoff>summary.max_cutoff) {
                summary.max_cutoff=cutoff;
                summary.max_core=core;
                summary.max_petals=petals;
                summary.max_mass=mass[petals];
                summary.max_components=components;
                summary.max_delta=delta;
            }
            if (cutoff>TAIL)
                summary.exceptional.push_back(
                    {k,core,petals,mass[petals],components,delta,cutoff});
        }
    }
    std::sort(summary.exceptional.begin(),summary.exceptional.end(),
              [](const Exceptional& a,const Exceptional& b) {
                  return std::tie(a.cutoff,a.core,a.petals) >
                         std::tie(b.cutoff,b.core,b.petals);
              });
    return summary;
}

struct Event {
    i64 x;
    int delta;
};

struct Integral {
    i64 mass;
    i64 components;
    std::size_t grouped_walls;
    std::size_t raw_events;
};

std::vector<Event> make_events(const std::vector<int>& speeds,i64 denominator) {
    std::vector<Event> events;
    std::size_t reserve=0;
    for (int speed:speeds) reserve+=static_cast<std::size_t>(2*speed);
    events.reserve(reserve);
    for (int speed:speeds) {
        require(denominator%(14LL*speed)==0,"denominator misses speed wall");
        const i64 unit=denominator/(14LL*speed);
        for (int tooth=0; tooth<speed; ++tooth) {
            events.push_back({(14LL*tooth+1)*unit,+1});
            events.push_back({(14LL*tooth+13)*unit,-1});
        }
    }
    std::sort(events.begin(),events.end(),[](const Event& a,const Event& b) {
        return std::tie(a.x,a.delta)<std::tie(b.x,b.delta);
    });
    return events;
}

Integral grouped_event_integral(const std::vector<int>& speeds,i64 denominator) {
    const std::vector<Event> events=make_events(speeds,denominator);
    i64 mass=0;
    i64 position=0;
    int active_count=0;
    bool left_active=false;
    i64 components=0;
    std::size_t groups=0;
    std::size_t index=0;
    while (index<events.size()) {
        const i64 x=events[index].x;
        require(x>position,"event coordinates not strictly advancing by group");
        if (left_active) mass+=x-position;
        int net=0;
        while (index<events.size() && events[index].x==x) {
            net+=events[index].delta;
            ++index;
        }
        active_count+=net;
        require(active_count>=0 && active_count<=static_cast<int>(speeds.size()),
                "invalid event active count");
        const bool right_active=active_count==static_cast<int>(speeds.size());
        if (right_active && !left_active) ++components;
        left_active=right_active;
        position=x;
        ++groups;
    }
    if (left_active) mass+=denominator-position;
    require(active_count==0,"endpoint sweep did not return to empty state");
    return {mass,components,groups,events.size()};
}

Integral midpoint_integral(const std::vector<int>& speeds,i64 denominator) {
    const std::vector<Event> events=make_events(speeds,denominator);
    std::vector<i64> walls{0,denominator};
    walls.reserve(events.size()+2);
    for (const Event& event:events) walls.push_back(event.x);
    std::sort(walls.begin(),walls.end());
    walls.erase(std::unique(walls.begin(),walls.end()),walls.end());
    i64 mass=0;
    i64 components=0;
    bool previous=false;
    for (std::size_t i=0;i+1<walls.size();++i) {
        bool active=true;
        for (int speed:speeds)
            if (!safe_mid(speed,walls[i],walls[i+1],denominator)) {
                active=false;
                break;
            }
        if (active) mass+=walls[i+1]-walls[i];
        if (active && !previous) ++components;
        previous=active;
    }
    // Every positive speed is unsafe on the first and last open cells, so the
    // linear component count is already the circular positive-length count.
    require(!safe_mid(speeds.front(),0,walls[1],denominator),
            "unexpected active first cell");
    return {mass,components,walls.size()-2,events.size()};
}

struct LiteralSummary {
    std::uint64_t bodies=0;
    std::uint64_t rows=0;
    std::uint64_t grouped_checks=0;
    std::uint64_t midpoint_checks=0;
    std::uint64_t replay_failures=0;
    std::uint64_t integrator_mismatches=0;
    std::uint64_t failures=0;
    std::uint64_t equalities=0;
    bool min_set=false;
    i128 min_delta=0;
    i64 min_denominator=1;
    i64 min_mass=0;
    int min_r=0;
    Exceptional min_body{};
    std::size_t max_grouped_walls=0;
    std::size_t max_raw_events=0;
};

LiteralSummary literal_audit(const MasterArrangement& master,
                             const std::vector<Exceptional>& exceptional) {
    LiteralSummary summary;
    summary.bodies=exceptional.size();
    for (const Exceptional& body:exceptional) {
        std::vector<int> base=body_speeds(body.core,body.petals,true);
        const Integral replay_event=grouped_event_integral(base,master.denominator);
        const Integral replay_mid=midpoint_integral(base,master.denominator);
        if (replay_event.mass!=body.mass ||
            replay_event.components!=body.components ||
            replay_mid.mass!=body.mass ||
            replay_mid.components!=body.components)
            ++summary.replay_failures;
        for (int r=TAIL;r<body.cutoff;++r) {
            std::vector<int> speeds=base;
            speeds.push_back(r);
            std::sort(speeds.begin(),speeds.end());
            require(std::adjacent_find(speeds.begin(),speeds.end())==speeds.end(),
                    "literal outsider duplicates a fixed speed");
            const i64 denominator=lcm_exact(master.denominator,14LL*r);
            const Integral event_value=grouped_event_integral(speeds,denominator);
            const Integral midpoint_value=midpoint_integral(speeds,denominator);
            ++summary.rows;
            ++summary.grouped_checks;
            ++summary.midpoint_checks;
            summary.max_grouped_walls=std::max(summary.max_grouped_walls,
                                                event_value.grouped_walls);
            summary.max_raw_events=std::max(summary.max_raw_events,
                                             event_value.raw_events);
            if (event_value.mass!=midpoint_value.mass ||
                event_value.components!=midpoint_value.components ||
                event_value.grouped_walls!=midpoint_value.grouped_walls)
                ++summary.integrator_mismatches;
            const i128 delta=static_cast<i128>(63)*event_value.mass-
                             static_cast<i128>(4)*denominator;
            summary.failures+=delta<0;
            summary.equalities+=delta==0;
            if (!summary.min_set ||
                delta*summary.min_denominator <
                    summary.min_delta*denominator) {
                summary.min_set=true;
                summary.min_delta=delta;
                summary.min_denominator=denominator;
                summary.min_mass=event_value.mass;
                summary.min_r=r;
                summary.min_body=body;
            }
        }
    }
    return summary;
}

} // namespace

int main() {
    const MasterArrangement master=make_master();
    require(master.denominator==INT64_C(91205797082400),
            "master denominator changed");
    require(master.cells.size()==7213,"master cell count changed");
    std::cout << "AUDIT architecture=dual-core-first-O-superset"
              << " literal=grouped-endpoint-plus-midpoint"
              << " threshold=" << TAIL
              << " denominator=" << master.denominator
              << " cells=" << master.cells.size() << '\n';

    std::array<AnalyticSummary,3> analytics;
    std::vector<Exceptional> all_exceptional;
    for (int k=4;k<=6;++k) {
        AnalyticSummary result=enumerate_core_first(master,k);
        std::cout << "ANALYTIC k=" << k
                  << " core_size=" << 9-k
                  << " cores=" << result.cores
                  << " profiles=" << result.profiles
                  << " nonstrict=" << result.nonstrict
                  << " max_cutoff=" << result.max_cutoff
                  << " max_body=" << mask_labels(result.max_core,result.max_petals)
                  << " max_mass=" << result.max_mass
                  << " max_components=" << result.max_components
                  << " max_delta=" << decimal(result.max_delta)
                  << " exceptional_bodies=" << result.exceptional.size()
                  << '\n';
        for (const Exceptional& body:result.exceptional)
            std::cout << "EXCEPTION k=" << body.k
                      << " body=" << mask_labels(body.core,body.petals)
                      << " mass=" << body.mass
                      << " components=" << body.components
                      << " delta=" << decimal(body.delta)
                      << " cutoff=" << body.cutoff
                      << " literal_rows=" << body.cutoff-TAIL << '\n';
        all_exceptional.insert(all_exceptional.end(),
                               result.exceptional.begin(),result.exceptional.end());
        analytics[k-4]=std::move(result);
    }

    std::sort(all_exceptional.begin(),all_exceptional.end(),
              [](const Exceptional& a,const Exceptional& b) {
                  return std::tie(a.k,a.core,a.petals)<
                         std::tie(b.k,b.core,b.petals);
              });
    const LiteralSummary literal=literal_audit(master,all_exceptional);
    for (int k=4;k<=6;++k) {
        std::vector<Exceptional> lane;
        for (const Exceptional& body:all_exceptional)
            if (body.k==k) lane.push_back(body);
        const LiteralSummary row=literal_audit(master,lane);
        std::cout << "LITERAL k=" << k
                  << " bodies=" << row.bodies
                  << " rows=" << row.rows
                  << " grouped_checks=" << row.grouped_checks
                  << " midpoint_checks=" << row.midpoint_checks
                  << " replay_failures=" << row.replay_failures
                  << " integrator_mismatches=" << row.integrator_mismatches
                  << " failures=" << row.failures
                  << " equalities=" << row.equalities;
        if (row.min_set)
            std::cout << " min_delta=" << decimal(row.min_delta)
                      << " min_denominator=" << row.min_denominator
                      << " min_mass=" << row.min_mass
                      << " min_r=" << row.min_r
                      << " min_body="
                      << mask_labels(row.min_body.core,row.min_body.petals);
        else
            std::cout << " min_delta=NA min_denominator=NA min_mass=NA"
                      << " min_r=NA min_body=NA";
        std::cout << " max_grouped_walls=" << row.max_grouped_walls
                  << " max_raw_events=" << row.max_raw_events << '\n';
    }

    std::uint64_t analytic_profiles=0;
    std::uint64_t analytic_nonstrict=0;
    int analytic_max=0;
    for (const AnalyticSummary& row:analytics) {
        analytic_profiles+=row.profiles;
        analytic_nonstrict+=row.nonstrict;
        analytic_max=std::max(analytic_max,row.max_cutoff);
    }
    std::cout << "SUMMARY analytic_profiles=" << analytic_profiles
              << " analytic_nonstrict=" << analytic_nonstrict
              << " analytic_max_cutoff=" << analytic_max
              << " exceptional_bodies=" << all_exceptional.size()
              << " literal_rows=" << literal.rows
              << " grouped_checks=" << literal.grouped_checks
              << " midpoint_checks=" << literal.midpoint_checks
              << " replay_failures=" << literal.replay_failures
              << " integrator_mismatches=" << literal.integrator_mismatches
              << " failures=" << literal.failures
              << " equalities=" << literal.equalities;
    if (literal.min_set)
        std::cout << " min_delta=" << decimal(literal.min_delta)
                  << " min_denominator=" << literal.min_denominator
                  << " min_mass=" << literal.min_mass
                  << " min_r=" << literal.min_r
                  << " min_body="
                  << mask_labels(literal.min_body.core,literal.min_body.petals);
    std::cout << " max_grouped_walls=" << literal.max_grouped_walls
              << " max_raw_events=" << literal.max_raw_events << '\n';

    require(analytic_nonstrict==0,"analytic base has nonstrict profile");
    require(literal.replay_failures==0,"direct base replay mismatch");
    require(literal.integrator_mismatches==0,"literal integrator mismatch");
    require(literal.failures==0,"literal tail failure");
    return 0;
}

// Independent exact audit of the fixed-50 high-petal limiting atlas.
//
// This implementation deliberately reverses the primary computation.  It
// fixes a C-core body, scatters grouped wall intervals by their 12-bit
// O-safe mask, and applies a superset zeta transform on O.  Thus all petal
// lanes for a core body are recovered simultaneously.  It imports neither
// the primary midpoint ledger nor its petal-first 18-bit subset transform.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using u32 = std::uint32_t;
using u64 = std::uint64_t;
__extension__ typedef unsigned __int128 u128;

constexpr std::array<int,18> C = {
    8,16,40,42,80,84,85,88,95,120,126,143,145,168,193,240,252,286};
constexpr std::array<int,12> O = {
    10,15,20,30,60,63,132,170,176,190,264,290};
constexpr int CENTER = 50;
constexpr u32 C_FULL = (u32{1} << C.size()) - 1;
constexpr std::uint16_t O_FULL =
    static_cast<std::uint16_t>((std::uint16_t{1} << O.size()) - 1);
constexpr std::size_t O_STATES = std::size_t{1} << O.size();
constexpr u64 EXPECTED_D = 91'205'797'082'400ULL;

struct Event {
    u64 position;
    u32 c_toggle;
    std::uint16_t o_toggle;
    bool center_toggle;
};

struct Cell {
    u64 length;
    u32 c_failed;
    std::uint16_t o_safe;
    bool center_safe;
};

struct Arrangement {
    u64 grid;
    std::size_t raw_events;
    std::size_t grouped_events;
    std::vector<Cell> cells;
};

struct DirectProfile {
    u64 grid;
    u64 mass;
    u64 components;
};

struct LaneStats {
    bool initialized = false;
    u64 profiles = 0;
    u64 nonstrict = 0;
    u64 min_delta = 0;
    u32 min_body = 0;
    u64 min_ties = 0;
    u64 max_cutoff = 0;
    u32 cutoff_body = 0;
    u64 cutoff_ties = 0;
    u64 profile_hash = 1'469'598'103'934'665'603ULL;
};

void require(bool condition,const char* message) {
    if (!condition) throw std::runtime_error(message);
}

u64 checked_lcm(u64 a, u64 b) {
    const u64 g = std::gcd(a,b);
    const u128 value = u128(a/g)*b;
    require(value <= std::numeric_limits<u64>::max(),"lcm overflow");
    return static_cast<u64>(value);
}

std::string c_labels(u32 mask) {
    std::string answer = "{";
    bool first = true;
    for (std::size_t i=0;i<C.size();++i) if (mask&(u32{1}<<i)) {
        if (!first) answer += ',';
        first = false;
        answer += std::to_string(C[i]);
    }
    answer += '}';
    return answer;
}

std::string o_labels(std::uint16_t mask) {
    std::string answer = "{";
    bool first = true;
    for (std::size_t i=0;i<O.size();++i)
        if (mask&(std::uint16_t{1}<<i)) {
            if (!first) answer += ',';
            first = false;
            answer += std::to_string(O[i]);
        }
    answer += '}';
    return answer;
}

void fnv_word(u64& hash, u64 value) {
    constexpr u64 prime = 1'099'511'628'211ULL;
    for (int byte=0;byte<8;++byte) {
        hash ^= (value>>(8*byte))&0xffU;
        hash *= prime;
    }
}

void append_events(std::vector<Event>& events, int speed, u64 grid,
                   int c_index, int o_index, bool center) {
    const u64 denominator = 14ULL*static_cast<u64>(speed);
    require(grid%denominator == 0,"nonintegral event grid");
    const u64 scale = grid/denominator;
    const u32 c_bit = c_index>=0 ? u32{1}<<c_index : 0;
    const std::uint16_t o_bit = o_index>=0
        ? static_cast<std::uint16_t>(std::uint16_t{1}<<o_index) : 0;
    for (int tooth=0;tooth<speed;++tooth) {
        const u64 center_numerator = 14ULL*static_cast<u64>(tooth);
        const u64 left = (center_numerator+denominator-1)%denominator;
        const u64 right = (center_numerator+1)%denominator;
        events.push_back({left*scale,c_bit,o_bit,center});
        events.push_back({right*scale,c_bit,o_bit,center});
    }
}

Arrangement build_arrangement() {
    u64 grid = 1;
    for (int speed:C) grid = checked_lcm(grid,14ULL*speed);
    for (int speed:O) grid = checked_lcm(grid,14ULL*speed);
    grid = checked_lcm(grid,14ULL*CENTER);
    require(grid == EXPECTED_D,"base denominator changed");

    std::vector<Event> raw;
    std::size_t expected_raw = 0;
    for (std::size_t i=0;i<C.size();++i) {
        append_events(raw,C[i],grid,static_cast<int>(i),-1,false);
        expected_raw += 2*C[i];
    }
    for (std::size_t i=0;i<O.size();++i) {
        append_events(raw,O[i],grid,-1,static_cast<int>(i),false);
        expected_raw += 2*O[i];
    }
    append_events(raw,CENTER,grid,-1,-1,true);
    expected_raw += 2*CENTER;
    require(raw.size() == expected_raw,"raw event count changed");
    std::sort(raw.begin(),raw.end(),[](const Event& a,const Event& b) {
        return a.position < b.position;
    });

    std::vector<Cell> cells;
    u32 c_failed = C_FULL;
    std::uint16_t o_safe = 0;
    bool center_safe = false;
    const u32 initial_c = c_failed;
    const std::uint16_t initial_o = o_safe;
    const bool initial_center = center_safe;
    u64 previous = 0;
    std::size_t cursor = 0,groups = 0;
    while (cursor<raw.size()) {
        const u64 position = raw[cursor].position;
        require(position>previous || groups==0,"event order failure");
        if (position>previous)
            cells.push_back({position-previous,c_failed,o_safe,center_safe});
        u32 ct = 0;
        std::uint16_t ot = 0;
        bool qt = false;
        while (cursor<raw.size() && raw[cursor].position==position) {
            ct ^= raw[cursor].c_toggle;
            ot ^= raw[cursor].o_toggle;
            qt ^= raw[cursor].center_toggle;
            ++cursor;
        }
        c_failed ^= ct;
        o_safe ^= ot;
        center_safe ^= qt;
        previous = position;
        ++groups;
    }
    if (previous<grid)
        cells.push_back({grid-previous,c_failed,o_safe,center_safe});
    require(c_failed==initial_c && o_safe==initial_o &&
            center_safe==initial_center,"circular state mismatch");
    require(cells.size()==groups+1,"group/cell count mismatch");
    require(cells.size()==7213,"base cell count changed");
    return {grid,raw.size(),groups,std::move(cells)};
}

void superset_zeta(std::array<u64,O_STATES>& values) {
    for (std::size_t bit=1;bit<O_STATES;bit<<=1)
        for (std::size_t block=0;block<O_STATES;block+=2*bit)
            for (std::size_t offset=0;offset<bit;++offset)
                values[block+offset] += values[block+bit+offset];
}

u64 ceil_ratio(u128 numerator,u128 denominator) {
    require(denominator!=0,"zero ceiling denominator");
    const u128 answer=(numerator+denominator-1)/denominator;
    require(answer<=std::numeric_limits<u64>::max(),"ceiling overflow");
    return static_cast<u64>(answer);
}

DirectProfile direct_profile(const std::vector<int>& speeds,u64 grid) {
    std::vector<std::pair<u64,std::size_t>> events;
    std::size_t expected = 0;
    for (std::size_t index=0;index<speeds.size();++index) {
        const u64 denominator=14ULL*static_cast<u64>(speeds[index]);
        require(grid%denominator==0,"nonintegral direct grid");
        const u64 scale=grid/denominator;
        expected += 2*speeds[index];
        for (int tooth=0;tooth<speeds[index];++tooth) {
            const u64 center_numerator=14ULL*static_cast<u64>(tooth);
            const u64 left=(center_numerator+denominator-1)%denominator;
            const u64 right=(center_numerator+1)%denominator;
            events.emplace_back(left*scale,index);
            events.emplace_back(right*scale,index);
        }
    }
    require(events.size()==expected,"direct event count changed");
    std::sort(events.begin(),events.end());
    std::vector<unsigned char> failing(speeds.size(),1);
    std::size_t active=speeds.size();
    u64 mass=0,components=0,previous=0;
    std::size_t cursor=0;
    while (cursor<events.size()) {
        const u64 position=events[cursor].first;
        if (active==0) mass += position-previous;
        const bool left_safe=active==0;
        while (cursor<events.size() && events[cursor].first==position) {
            const std::size_t index=events[cursor].second;
            if (failing[index]) {failing[index]=0;--active;}
            else {failing[index]=1;++active;}
            ++cursor;
        }
        const bool right_safe=active==0;
        if (right_safe && !left_safe) ++components;
        previous=position;
    }
    if (active==0) mass += grid-previous;
    require(active==speeds.size(),"direct circular state mismatch");
    return {grid,mass,components};
}

std::vector<int> profile_speeds(u32 body,std::uint16_t required) {
    std::vector<int> speeds;
    for (std::size_t i=0;i<C.size();++i)
        if (body&(u32{1}<<i)) speeds.push_back(C[i]);
    for (std::size_t i=0;i<O.size();++i)
        if (required&(std::uint16_t{1}<<i)) speeds.push_back(O[i]);
    speeds.push_back(CENTER);
    return speeds;
}

void run_k(const Arrangement& arrangement,int k) {
    require(k>=5 && k<=9,"petal count outside 5..9");
    const int core_size=9-k;
    std::vector<std::uint16_t> lanes;
    for (std::uint16_t mask=0;mask<=O_FULL;++mask)
        if (std::popcount(mask)==k) lanes.push_back(mask);
    std::array<LaneStats,O_STATES> stats{};

    u64 bodies=0;
    for (u32 body=C_FULL;;body=(body-1)&C_FULL) {
        if (std::popcount(body)==core_size) {
            ++bodies;
            std::array<u64,O_STATES> mass{},starts{},continuations{};
            for (std::size_t i=0;i<arrangement.cells.size();++i) {
                const Cell& current=arrangement.cells[i];
                const Cell& previous=arrangement.cells[
                    (i+arrangement.cells.size()-1)%arrangement.cells.size()];
                const bool current_safe=current.center_safe &&
                    (current.c_failed&body)==0;
                if (!current_safe) continue;
                mass[current.o_safe] += current.length;
                ++starts[current.o_safe];
                const bool previous_safe=previous.center_safe &&
                    (previous.c_failed&body)==0;
                if (previous_safe)
                    ++continuations[current.o_safe&previous.o_safe];
            }
            superset_zeta(mass);
            superset_zeta(starts);
            superset_zeta(continuations);
            for (std::uint16_t required:lanes) {
                LaneStats& lane=stats[required];
                lane.initialized=true;
                ++lane.profiles;
                const u64 components=starts[required]-continuations[required];
                require(components>0,"component count nonpositive");
                const u64 safe_mass=mass[required];
                const u128 lhs=u128(54)*safe_mass;
                const u128 rhs=u128(4)*arrangement.grid;
                const bool positive=lhs>rhs;
                const u64 delta=positive ? static_cast<u64>(lhs-rhs) : 0;
                const u64 cutoff=positive
                    ? ceil_ratio(u128(54)*components*arrangement.grid,
                                 u128(7)*delta) : 0;
                if (!positive) ++lane.nonstrict;
                if (lane.profiles==1 || delta<lane.min_delta) {
                    lane.min_delta=delta;
                    lane.min_body=body;
                    lane.min_ties=1;
                } else if (delta==lane.min_delta) {
                    ++lane.min_ties;
                }
                if (cutoff>lane.max_cutoff) {
                    lane.max_cutoff=cutoff;
                    lane.cutoff_body=body;
                    lane.cutoff_ties=1;
                } else if (cutoff==lane.max_cutoff) {
                    ++lane.cutoff_ties;
                }
                fnv_word(lane.profile_hash,body);
                fnv_word(lane.profile_hash,safe_mass);
                fnv_word(lane.profile_hash,components);
                fnv_word(lane.profile_hash,delta);
                fnv_word(lane.profile_hash,cutoff);
            }
        }
        if (body==0) break;
    }

    u64 total_profiles=0,strict=0,nonstrict=0,global_cutoff=0;
    u64 min_tied_lanes=0,cutoff_tied_lanes=0;
    bool best_set=false,global_min_set=false,global_cutoff_set=false;
    u64 best_min_delta=0,best_cutoff=0,global_min_delta=0;
    std::uint16_t best_required=0,global_min_required=0,global_cutoff_required=0;
    u32 global_min_body=0,global_cutoff_body=0;
    u64 atlas_hash=1'469'598'103'934'665'603ULL;
    for (std::uint16_t required:lanes) {
        const LaneStats& lane=stats[required];
        require(lane.initialized && lane.profiles==bodies,
                "lane profile universe changed");
        total_profiles += lane.profiles;
        if (lane.nonstrict==0) ++strict; else ++nonstrict;
        min_tied_lanes += lane.min_ties>1;
        cutoff_tied_lanes += lane.cutoff_ties>1;
        if (!best_set || lane.min_delta>best_min_delta) {
            best_set=true;
            best_min_delta=lane.min_delta;
            best_required=required;
            best_cutoff=lane.max_cutoff;
        }
        if (!global_min_set || lane.min_delta<global_min_delta) {
            global_min_set=true;
            global_min_delta=lane.min_delta;
            global_min_required=required;
            global_min_body=lane.min_body;
        }
        if (!global_cutoff_set || lane.max_cutoff>global_cutoff) {
            global_cutoff_set=true;
            global_cutoff=lane.max_cutoff;
            global_cutoff_required=required;
            global_cutoff_body=lane.cutoff_body;
        }
        fnv_word(atlas_hash,required);
        fnv_word(atlas_hash,lane.profile_hash);
        fnv_word(atlas_hash,lane.profiles);
        fnv_word(atlas_hash,lane.nonstrict);
        fnv_word(atlas_hash,lane.min_delta);
        fnv_word(atlas_hash,lane.min_body);
        fnv_word(atlas_hash,lane.max_cutoff);
        fnv_word(atlas_hash,lane.cutoff_body);
        std::cout << "LANE k=" << k
                  << " labels=" << o_labels(required)
                  << " profiles=" << lane.profiles
                  << " nonstrict=" << lane.nonstrict
                  << " min_delta=" << lane.min_delta
                  << " min_body=" << c_labels(lane.min_body)
                  << " min_ties=" << lane.min_ties
                  << " max_cutoff=" << lane.max_cutoff
                  << " cutoff_body=" << c_labels(lane.cutoff_body)
                  << " cutoff_ties=" << lane.cutoff_ties
                  << " profile_fnv64=" << lane.profile_hash << '\n';
    }

    const DirectProfile one=direct_profile({1},14);
    require(one.mass==12 && one.components==1,"one-speed control failed");
    const DirectProfile direct_min=direct_profile(
        profile_speeds(global_min_body,global_min_required),arrangement.grid);
    const DirectProfile direct_cutoff=direct_profile(
        profile_speeds(global_cutoff_body,global_cutoff_required),arrangement.grid);
    const u64 replay_min_delta=static_cast<u64>(
        u128(54)*direct_min.mass-u128(4)*arrangement.grid);
    const u64 replay_cutoff=ceil_ratio(
        u128(54)*direct_cutoff.components*arrangement.grid,
        u128(7)*(u128(54)*direct_cutoff.mass-u128(4)*arrangement.grid));
    const u64 direct_cutoff_delta=static_cast<u64>(
        u128(54)*direct_cutoff.mass-u128(4)*arrangement.grid);
    require(replay_min_delta==global_min_delta,"global minimum replay failed");
    require(replay_cutoff==global_cutoff,"global cutoff replay failed");

    std::cout << "SUMMARY k=" << k
              << " core_size=" << core_size
              << " denominator=" << arrangement.grid
              << " raw_events=" << arrangement.raw_events
              << " grouped_events=" << arrangement.grouped_events
              << " cells=" << arrangement.cells.size()
              << " lanes=" << lanes.size()
              << " core_bodies=" << bodies
              << " profiles=" << total_profiles
              << " strict=" << strict
              << " nonstrict=" << nonstrict
              << " global_max_cutoff=" << global_cutoff
              << " global_cutoff_labels=" << o_labels(global_cutoff_required)
              << " global_cutoff_body=" << c_labels(global_cutoff_body)
              << " global_min_delta=" << global_min_delta
              << " global_min_labels=" << o_labels(global_min_required)
              << " global_min_body=" << c_labels(global_min_body)
              << " best=" << o_labels(best_required)
              << " best_min_delta=" << best_min_delta
              << " best_cutoff=" << best_cutoff
              << " atlas_fnv64=" << atlas_hash
              << " min_tied_lanes=" << min_tied_lanes
              << " cutoff_tied_lanes=" << cutoff_tied_lanes
              << " direct_min_mass=" << direct_min.mass
              << " direct_min_components=" << direct_min.components
              << " direct_cutoff_mass=" << direct_cutoff.mass
              << " direct_cutoff_components=" << direct_cutoff.components
              << " direct_cutoff_delta=" << direct_cutoff_delta
              << " direct_replays=3"
              << " status=" << (nonstrict==0 ? "PASS" : "FAIL") << '\n';
}

} // namespace

int main(int argc,char** argv) {
    require(argc==2,"usage: audit <petal-count 5..9>");
    const int k=std::stoi(argv[1]);
    const Arrangement arrangement=build_arrangement();
    run_k(arrangement,k);
}

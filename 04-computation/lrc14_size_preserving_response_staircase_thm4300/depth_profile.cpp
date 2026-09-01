// Scratch-only exact deletion-depth profiler.  For a fixed labelled nine-body
// B it restricts the full literal failure-atom measure to P\B, performs the
// complete 21-bit subset zeta transform, and finds delta_Q(B) over all ranks
// 0..21.  The maintained THM-4296 depth-nine hostile is used as a control.

#define main thm4296_r632_hostile_main
#include "04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/r632_detached_hostile_survivor.cpp"
#undef main

namespace {

struct FullGeometry {
    i64 grid = 0;
    u64 cells = 0;
    i64 pair_ticks = 0;
    std::vector<std::pair<u32, i64>> classes;
};

FullGeometry build_full_geometry(int q, int r) {
    i64 grid = fixed_grid();
    grid = checked_lcm(grid, 14LL * q);
    grid = checked_lcm(grid, 14LL * r);
    std::vector<i64> walls = {0, grid};
    auto add_walls = [&](int speed) {
        const i64 divisor = 14LL * speed;
        require(grid % divisor == 0, "nonintegral wall unit");
        const i64 unit = grid / divisor;
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    for (int speed : kPool) add_walls(speed);
    add_walls(q);
    add_walls(r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    std::map<u32, i64> by_failure;
    FullGeometry geometry;
    geometry.grid = grid;
    geometry.cells = walls.size() - 1;
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        if (!safe_midpoint(q, grid, left, right) ||
            !safe_midpoint(r, grid, left, right))
            continue;
        u32 failure = 0;
        for (unsigned label = 0; label < kPool.size(); ++label)
            if (!safe_midpoint(kPool[label], grid, left, right))
                failure |= u32{1} << label;
        const i64 width = right - left;
        geometry.pair_ticks += width;
        by_failure[failure] += width;
    }
    geometry.classes.assign(by_failure.begin(), by_failure.end());
    return geometry;
}

u64 colex(u32 mask) {
    std::array<std::array<u64, 22>, 31> c{};
    for (unsigned n = 0; n <= 30; ++n) {
        c[n][0] = 1;
        for (unsigned k = 1; k <= 21; ++k)
            c[n][k] = (k > n) ? 0 : ((k == n) ? 1 : c[n - 1][k] + c[n - 1][k - 1]);
    }
    u64 answer = 0;
    unsigned ordinal = 1;
    for (unsigned bit = 0; bit < 30; ++bit)
        if (mask & (u32{1} << bit)) answer += c[bit][ordinal++];
    return answer;
}

struct RankStats {
    u64 candidates = 0;
    u64 active = 0;
    bool have = false;
    i128 maximum_ticks = 0;
    u32 least_maximum = 0;
    bool have_active = false;
    i128 least_active_ticks = 0;
    u32 least_active = 0;
};

void print_fraction(i128 ticks, i128 denominator) {
    const i128 divisor = gcd128(ticks, denominator);
    std::cout << decimal(ticks / divisor) << '/' << decimal(denominator / divisor);
}

void profile_body(const FullGeometry& geometry, u32 body) {
    require(std::popcount(body) == 9 && body < (u32{1} << 30),
            "body must be a labelled nine-set");
    std::array<unsigned char, 21> positions{};
    std::array<signed char, 30> local{};
    local.fill(-1);
    unsigned count = 0;
    for (unsigned bit = 0; bit < 30; ++bit) {
        if (body & (u32{1} << bit)) continue;
        positions[count] = static_cast<unsigned char>(bit);
        local[bit] = static_cast<signed char>(count++);
    }
    require(count == 21, "body complement size changed");

    constexpr u32 states = u32{1} << 21;
    std::vector<i64> mass(states, 0);
    Fnv atom_ledger;
    u64 retained_classes = 0;
    i64 retained_mass = 0;
    for (const auto& [atom, width] : geometry.classes) {
        if (atom & body) continue;
        u32 compact = 0;
        for (unsigned bit = 0; bit < 30; ++bit) {
            if ((atom & (u32{1} << bit)) == 0) continue;
            require(local[bit] >= 0, "atom compression escaped complement");
            compact |= u32{1} << local[bit];
        }
        mass[compact] += width;
        ++retained_classes;
        retained_mass += width;
        atom_ledger.add(atom);
        atom_ledger.add(static_cast<u64>(width));
    }
    for (unsigned bit = 0; bit < 21; ++bit) {
        const u32 flag = u32{1} << bit;
        for (u32 mask = 0; mask < states; ++mask)
            if (mask & flag) mass[mask] += mass[mask ^ flag];
    }

    std::vector<u32> expanded(states, 0);
    std::array<RankStats, 22> rank{};
    const i128 denominator = static_cast<i128>(63) * geometry.grid;
    Fnv profile_ledger;
    for (u32 compact = 0; compact < states; ++compact) {
        if (compact != 0) {
            const u32 low = compact & (u32{0} - compact);
            const unsigned index = std::countr_zero(low);
            expanded[compact] = expanded[compact ^ low] |
                                (u32{1} << positions[index]);
        }
        const u32 repair = expanded[compact];
        const unsigned depth = std::popcount(compact);
        const i128 ticks = static_cast<i128>(63) * mass[compact] -
                           static_cast<i128>(4) * geometry.grid;
        RankStats& stats = rank[depth];
        ++stats.candidates;
        if (!stats.have || ticks > stats.maximum_ticks ||
            (ticks == stats.maximum_ticks && repair < stats.least_maximum)) {
            stats.have = true;
            stats.maximum_ticks = ticks;
            stats.least_maximum = repair;
        }
        if (ticks >= 0) {
            ++stats.active;
            if (!stats.have_active || repair < stats.least_active) {
                stats.have_active = true;
                stats.least_active = repair;
                stats.least_active_ticks = ticks;
            }
        }
        profile_ledger.add(repair);
        profile_ledger.add(static_cast<u64>(mass[compact]));
        add_i128(profile_ledger, ticks);
    }

    unsigned depth = 22;
    for (unsigned k = 0; k <= 21; ++k)
        if (rank[k].active != 0) {
            depth = k;
            break;
        }
    require(depth <= 21, "full deletion complement unexpectedly inactive");
    std::cout << "BODY " << hex8(body) << " LABELS {" << labels(body)
              << "} RETAINED_ATOMS " << retained_classes
              << " RETAINED_ATOM_MASS " << retained_mass
              << " ATOM_FNV " << std::hex << atom_ledger.state << std::dec
              << " DEPTH " << depth
              << " PROFILE_FNV " << std::hex << profile_ledger.state << std::dec
              << '\n';
    if (depth > 0) {
        const RankStats& prior = rank[depth - 1];
        std::cout << "PRIOR RANK " << (depth - 1)
                  << " CANDIDATES " << prior.candidates
                  << " ACTIVE " << prior.active
                  << " MAX_MASK " << hex8(prior.least_maximum)
                  << " MAX_COLEX " << colex(prior.least_maximum)
                  << " MAX_MARGIN_TICKS63 " << decimal(prior.maximum_ticks)
                  << " MAX_SURPLUS ";
        print_fraction(prior.maximum_ticks, denominator);
        std::cout << '\n';
    }
    const RankStats& first = rank[depth];
    std::cout << "FIRST RANK " << depth
              << " CANDIDATES " << first.candidates
              << " ACTIVE " << first.active
              << " LEAST_ACTIVE " << hex8(first.least_active)
              << " LEAST_ACTIVE_COLEX " << colex(first.least_active)
              << " LEAST_ACTIVE_MARGIN_TICKS63 "
              << decimal(first.least_active_ticks) << " LEAST_ACTIVE_SURPLUS ";
    print_fraction(first.least_active_ticks, denominator);
    std::cout << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc >= 4, "usage: depth_profile Q R BODY_HEX...");
        const int q = std::stoi(argv[1]);
        const int r = std::stoi(argv[2]);
        const FullGeometry geometry = build_full_geometry(q, r);
        Fnv geometry_ledger;
        for (const auto& [atom, width] : geometry.classes) {
            geometry_ledger.add(atom);
            geometry_ledger.add(static_cast<u64>(width));
        }
        std::cout << "THM4296_LABELLED_DELETION_DEPTH_PROFILE_V1\n"
                  << "PAIR " << q << ',' << r << " GRID " << geometry.grid
                  << " CELLS " << geometry.cells << " PAIR_TICKS "
                  << geometry.pair_ticks << " FULL_FAILURE_CLASSES "
                  << geometry.classes.size() << " GEOMETRY_FNV " << std::hex
                  << geometry_ledger.state << std::dec << '\n';
        for (int index = 3; index < argc; ++index) {
            const u32 body = static_cast<u32>(std::stoul(argv[index], nullptr, 16));
            profile_body(geometry, body);
        }
        std::cout << "VERDICT PASS FULL_21_BIT_FAILURE_ATOM_ZETA_DEPTHS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "DEPTH_PROFILE_ERROR " << error.what() << '\n';
        return 1;
    }
}

// Detached literal-wall verifier for the post-THM-4277 joint421 scout.
//
// This translation unit imports only the clean-room direct-wall machinery.
// It never imports or calls an endpoint atom, cocycle, primitive prefix, colex
// rank, zeta transform, response atlas, or primary-scout helper.
//
// All claimed-common pairs are checked against all 421 repairs.  For every
// complementary pair, the locked deck is checked in order through the claimed
// first inactive repair, and its literal margin ratio is compared exactly with
// the primary witness ratio.  Fixed-pool walls are built literally once; pair
// walls are rebuilt for every pair.  Contiguous eligible fixed-cell intervals
// are only a lossless batching of direct literal cell-width sums.

#define CARRIER_CLEANROOM_CONTROL_LIBRARY_ONLY
#include "04-computation/lrc14_three_round_learned_carrier_thm4266/cleanroom_resistance_controls.cpp"
#undef CARRIER_CLEANROOM_CONTROL_LIBRARY_ONLY

#include <atomic>
#include <fstream>
#include <iomanip>

namespace {

constexpr std::size_t EXPECTED_COMMON = 148063;
constexpr u64 EXPECTED_COMMON_FNV = UINT64_C(0x465fb756a183167e);
constexpr std::size_t EXPECTED_HOSTILE = 24259;
constexpr u64 EXPECTED_PRIMARY_HOSTILE_FNV = UINT64_C(0xe85cf4266807ae46);
constexpr std::size_t EXPECTED_RESIDUAL = 172322;
constexpr u64 EXPECTED_RESIDUAL_FNV = UINT64_C(0x30b2a7e597ac548c);
constexpr std::size_t EXPECTED_DECK = 421;
constexpr u64 EXPECTED_DECK_FNV = UINT64_C(0x20d63dd42fe8150e);
constexpr u64 EXPECTED_COMMON_TESTS = UINT64_C(62334523);
constexpr u64 EXPECTED_ALL_TESTS = UINT64_C(64116779);
constexpr i64 FIXED_GRID = INT64_C(18241159416480);

struct EdgeLiteral {
    int q = 0;
    int r = 0;
    auto operator<=>(const EdgeLiteral&) const = default;
};

void add_edge_literal(Fnv& ledger, const EdgeLiteral& edge) {
    ledger.add(static_cast<u64>(edge.q));
    ledger.add(static_cast<u64>(edge.r));
}

u64 edge_fnv_literal(const std::vector<EdgeLiteral>& edges) {
    Fnv ledger;
    for (const EdgeLiteral& edge : edges) add_edge_literal(ledger, edge);
    return ledger.state;
}

u64 deck_fnv_literal(const std::vector<u32>& deck) {
    Fnv ledger;
    for (u32 mask : deck) ledger.add(mask);
    return ledger.state;
}

void add_i128_literal(Fnv& ledger, i128 value) {
    const __uint128_t bits = static_cast<__uint128_t>(value);
    ledger.add(static_cast<u64>(bits));
    ledger.add(static_cast<u64>(bits >> 64));
}

i128 parse_signed_i128(const std::string& text) {
    require(!text.empty(), "empty signed integer");
    const bool negative = text.front() == '-';
    const std::size_t begin = negative ? 1 : 0;
    require(begin < text.size(), "sign without digits");
    i128 value = 0;
    for (std::size_t index = begin; index < text.size(); ++index) {
        const char c = text[index];
        require('0' <= c && c <= '9', "bad signed integer digit");
        const int digit = c - '0';
        require(value <= (std::numeric_limits<i128>::max() - digit) / 10,
                "signed integer overflow");
        value = 10 * value + digit;
    }
    return negative ? -value : value;
}

std::vector<u32> read_literal_deck(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open locked joint deck");
    std::vector<u32> deck;
    std::set<u32> distinct;
    std::string word;
    while (input >> word) {
        std::size_t used = 0;
        const u64 wide = std::stoull(word, &used, 16);
        require(used == word.size() && wide < (UINT64_C(1) << 30),
                "bad deck mask token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "deck mask arity/distinctness changed");
        deck.push_back(mask);
    }
    require(deck.size() == EXPECTED_DECK &&
                deck_fnv_literal(deck) == EXPECTED_DECK_FNV,
            "locked joint deck changed");
    return deck;
}

std::vector<EdgeLiteral> read_common_literal(
    const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open common-pair ledger");
    std::vector<EdgeLiteral> edges;
    std::string line;
    while (std::getline(input, line)) {
        require(!line.empty(), "empty common-pair row");
        std::istringstream row(line);
        EdgeLiteral edge;
        char comma = 0;
        char trailing = 0;
        row >> edge.q >> comma >> edge.r;
        require(static_cast<bool>(row) && comma == ',' && edge.q > 0 &&
                    edge.q < edge.r && !(row >> trailing),
                "bad common-pair row");
        edges.push_back(edge);
    }
    require(std::is_sorted(edges.begin(), edges.end()) &&
                std::adjacent_find(edges.begin(), edges.end()) == edges.end() &&
                edges.size() == EXPECTED_COMMON &&
                edge_fnv_literal(edges) == EXPECTED_COMMON_FNV,
            "common-pair ledger changed");
    return edges;
}

struct HostileLiteral {
    EdgeLiteral edge;
    u64 tests = 0;
    u32 mask = 0;
    i128 primary_margin = 0;
    i128 primary_denominator = 0;
};

std::vector<HostileLiteral> read_hostile_literal(
    const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open primary hostile controls");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,tests,first_inactive_mask,margin_num,den",
            "bad primary hostile header");
    std::vector<HostileLiteral> rows;
    Fnv ledger;
    while (std::getline(input, line)) {
        require(!line.empty(), "empty primary hostile row");
        std::array<std::string, 6> field;
        std::size_t begin = 0;
        for (std::size_t index = 0; index < field.size(); ++index) {
            const std::size_t comma = line.find(',', begin);
            require((index + 1 == field.size()) ==
                        (comma == std::string::npos),
                    "bad primary hostile field count");
            field[index] = line.substr(
                begin, comma == std::string::npos ? std::string::npos
                                                  : comma - begin);
            begin = comma == std::string::npos ? line.size() : comma + 1;
        }
        HostileLiteral row;
        row.edge.q = std::stoi(field[0]);
        row.edge.r = std::stoi(field[1]);
        row.tests = std::stoull(field[2]);
        std::size_t used = 0;
        const u64 wide = std::stoull(field[3], &used, 16);
        require(used == field[3].size() && wide < (UINT64_C(1) << 30),
                "bad hostile mask token");
        row.mask = static_cast<u32>(wide);
        row.primary_margin = parse_signed_i128(field[4]);
        row.primary_denominator = parse_signed_i128(field[5]);
        require(row.edge.q > 0 && row.edge.q < row.edge.r &&
                    row.tests >= 1 && row.tests <= EXPECTED_DECK &&
                    std::popcount(row.mask) == 8 &&
                    row.primary_margin < 0 && row.primary_denominator > 0,
                "invalid primary hostile control");
        add_edge_literal(ledger, row.edge);
        ledger.add(row.tests);
        ledger.add(row.mask);
        add_i128_literal(ledger, row.primary_margin);
        add_i128_literal(ledger, row.primary_denominator);
        rows.push_back(row);
    }
    require(std::is_sorted(rows.begin(), rows.end(),
                           [](const HostileLiteral& left,
                              const HostileLiteral& right) {
                               return left.edge < right.edge;
                           }) &&
                std::adjacent_find(
                    rows.begin(), rows.end(),
                    [](const HostileLiteral& left,
                       const HostileLiteral& right) {
                        return left.edge == right.edge;
                    }) == rows.end() &&
                rows.size() == EXPECTED_HOSTILE &&
                ledger.state == EXPECTED_PRIMARY_HOSTILE_FNV,
            "primary hostile-control ledger changed");
    return rows;
}

struct FixedCellLiteral {
    u32 failed = 0;
};

struct IntervalLiteral {
    std::uint16_t left = 0;
    std::uint16_t right = 0;
};

struct RepairLiteral {
    u32 mask = 0;
    std::vector<IntervalLiteral> intervals;
};

struct FixedGeometryLiteral {
    std::vector<i64> walls;
    std::vector<FixedCellLiteral> cells;
    std::vector<RepairLiteral> repairs;
    u64 interval_count = 0;
    u64 ledger = 0;
};

FixedGeometryLiteral build_fixed_literal(const std::vector<u32>& deck) {
    i64 grid = 1;
    for (int speed : POOL) grid = checked_lcm(grid, 14LL * speed);
    require(grid == FIXED_GRID, "literal fixed grid changed");
    FixedGeometryLiteral out;
    out.walls = {0, grid};
    for (int speed : POOL) {
        const i64 unit = grid / (14LL * speed);
        require(unit * (14LL * speed) == grid,
                "nonintegral literal fixed wall unit");
        for (int tooth = 0; tooth < speed; ++tooth) {
            out.walls.push_back((14LL * tooth + 1) * unit);
            out.walls.push_back((14LL * tooth + 13) * unit);
        }
    }
    std::sort(out.walls.begin(), out.walls.end());
    out.walls.erase(std::unique(out.walls.begin(), out.walls.end()),
                    out.walls.end());
    require(out.walls.size() == 7134 && out.walls.front() == 0 &&
                out.walls.back() == grid,
            "literal fixed wall ledger changed");
    out.cells.reserve(out.walls.size() - 1);
    Fnv ledger;
    ledger.add(out.walls.size());
    for (std::size_t index = 1; index < out.walls.size(); ++index) {
        const i64 left = out.walls[index - 1];
        const i64 right = out.walls[index];
        require(left < right, "nonpositive literal fixed cell");
        u32 failed = 0;
        for (unsigned vertex = 0; vertex < POOL.size(); ++vertex) {
            if (!safe_midpoint(POOL[vertex], grid, left, right))
                failed |= u32{1} << vertex;
        }
        out.cells.push_back({failed});
        ledger.add(left);
        ledger.add(right);
        ledger.add(failed);
    }
    out.repairs.reserve(deck.size());
    for (u32 repair : deck) {
        RepairLiteral indexed;
        indexed.mask = repair;
        std::size_t cell = 0;
        while (cell < out.cells.size()) {
            if ((out.cells[cell].failed & ~repair) != 0) {
                ++cell;
                continue;
            }
            const std::size_t left = cell;
            do {
                ++cell;
            } while (cell < out.cells.size() &&
                     (out.cells[cell].failed & ~repair) == 0);
            require(left < UINT16_MAX && cell < UINT16_MAX,
                    "literal interval index exceeds u16");
            indexed.intervals.push_back({
                static_cast<std::uint16_t>(left),
                static_cast<std::uint16_t>(cell)});
        }
        require(!indexed.intervals.empty(), "empty literal repair geometry");
        ledger.add(repair);
        ledger.add(indexed.intervals.size());
        out.interval_count += indexed.intervals.size();
        for (const IntervalLiteral interval : indexed.intervals) {
            ledger.add(interval.left);
            ledger.add(interval.right);
        }
        out.repairs.push_back(std::move(indexed));
    }
    out.ledger = ledger.state;
    require(out.interval_count == 81164,
            "literal repair interval count changed");
    return out;
}

struct LiteralMassGeometry {
    i64 grid = 0;
    i64 pair_ticks = 0;
    u64 joint_cells = 0;
    std::vector<i64> safe_prefix;
};

LiteralMassGeometry build_literal_mass_geometry(
    const EdgeLiteral& edge, const FixedGeometryLiteral& fixed) {
    i64 grid = FIXED_GRID;
    grid = checked_lcm(grid, 14LL * edge.q);
    grid = checked_lcm(grid, 14LL * edge.r);
    require(grid > 0 && grid % FIXED_GRID == 0,
            "literal pair grid does not scale fixed grid");
    const i64 scale = grid / FIXED_GRID;
    std::vector<i64> walls;
    walls.reserve(fixed.walls.size() + 2 * (edge.q + edge.r));
    for (i64 wall : fixed.walls) {
        const i128 scaled = static_cast<i128>(wall) * scale;
        require(scaled <= std::numeric_limits<i64>::max(),
                "scaled fixed wall overflow");
        walls.push_back(static_cast<i64>(scaled));
    }
    auto add_pair_walls = [&](int speed) {
        const i64 unit = grid / (14LL * speed);
        require(unit * (14LL * speed) == grid,
                "nonintegral literal pair wall unit");
        for (int tooth = 0; tooth < speed; ++tooth) {
            walls.push_back((14LL * tooth + 1) * unit);
            walls.push_back((14LL * tooth + 13) * unit);
        }
    };
    add_pair_walls(edge.q);
    add_pair_walls(edge.r);
    std::sort(walls.begin(), walls.end());
    walls.erase(std::unique(walls.begin(), walls.end()), walls.end());
    require(walls.front() == 0 && walls.back() == grid,
            "literal joint boundary lost");

    LiteralMassGeometry out;
    out.grid = grid;
    out.joint_cells = walls.size() - 1;
    out.safe_prefix.assign(fixed.cells.size() + 1, 0);
    std::size_t fixed_cell = 0;
    i64 current_safe = 0;
    for (std::size_t index = 1; index < walls.size(); ++index) {
        const i64 left = walls[index - 1];
        const i64 right = walls[index];
        require(fixed_cell < fixed.cells.size(),
                "literal joint sweep passed fixed cells");
        const i128 scaled_right =
            static_cast<i128>(fixed.walls[fixed_cell + 1]) * scale;
        require(static_cast<i128>(right) <= scaled_right,
                "literal joint cell crosses a fixed wall");
        if (safe_midpoint(edge.q, grid, left, right) &&
            safe_midpoint(edge.r, grid, left, right)) {
            current_safe += right - left;
            out.pair_ticks += right - left;
        }
        if (static_cast<i128>(right) == scaled_right) {
            out.safe_prefix[fixed_cell + 1] =
                out.safe_prefix[fixed_cell] + current_safe;
            current_safe = 0;
            ++fixed_cell;
        }
    }
    require(fixed_cell == fixed.cells.size() && current_safe == 0 &&
                out.safe_prefix.back() == out.pair_ticks,
            "literal joint sweep did not end at fixed boundary");
    return out;
}

i64 repair_mass_literal(const LiteralMassGeometry& geometry,
                        const RepairLiteral& repair) {
    i64 mass = 0;
    for (const IntervalLiteral interval : repair.intervals)
        mass += geometry.safe_prefix[interval.right] -
                geometry.safe_prefix[interval.left];
    return mass;
}

bool rational_less_literal(i128 left_num, i128 left_den,
                           i128 right_num, i128 right_den) {
    require(left_num >= 0 && right_num >= 0 &&
                left_den > 0 && right_den > 0,
            "invalid literal rational comparison");
    bool reversed = false;
    while (true) {
        const i128 left_q = left_num / left_den;
        const i128 right_q = right_num / right_den;
        if (left_q != right_q)
            return reversed ? left_q > right_q : left_q < right_q;
        left_num %= left_den;
        right_num %= right_den;
        if (left_num == 0 || right_num == 0) {
            if (left_num == 0 && right_num == 0) return false;
            const bool raw_less = left_num == 0;
            return reversed ? !raw_less : raw_less;
        }
        std::swap(left_num, left_den);
        std::swap(right_num, right_den);
        reversed = !reversed;
    }
}

bool rational_equal_literal(i128 left_num, i128 left_den,
                            i128 right_num, i128 right_den) {
    return !rational_less_literal(left_num, left_den,
                                  right_num, right_den) &&
           !rational_less_literal(right_num, right_den,
                                  left_num, left_den);
}

bool signed_ratio_equal(i128 left_num, i128 left_den,
                        i128 right_num, i128 right_den) {
    require(left_den > 0 && right_den > 0,
            "nonpositive signed-ratio denominator");
    const i128 left_gcd = gcd128(left_num, left_den);
    const i128 right_gcd = gcd128(right_num, right_den);
    return left_num / left_gcd == right_num / right_gcd &&
           left_den / left_gcd == right_den / right_gcd;
}

struct JobLiteral {
    EdgeLiteral edge;
    bool common = false;
    HostileLiteral hostile;
};

struct PairResultLiteral {
    u64 tests = 0;
    u64 equalities = 0;
    u64 joint_cells = 0;
    i64 grid = 0;
    i64 pair_ticks = 0;
    u32 weakest_mask = 0;
    i128 weakest_margin = 0;
    u32 hostile_mask = 0;
    i64 hostile_mass = 0;
    i128 hostile_margin = 0;
    u64 matrix_ledger = 0;
};

PairResultLiteral verify_literal_job(
    const JobLiteral& job, const FixedGeometryLiteral& fixed) {
    const LiteralMassGeometry geometry =
        build_literal_mass_geometry(job.edge, fixed);
    PairResultLiteral out;
    out.grid = geometry.grid;
    out.pair_ticks = geometry.pair_ticks;
    out.joint_cells = geometry.joint_cells;
    bool weakest_set = false;
    Fnv matrix;
    add_edge_literal(matrix, job.edge);
    matrix.add(job.common);
    matrix.add(geometry.joint_cells);
    matrix.add(geometry.grid);
    matrix.add(geometry.pair_ticks);
    const std::size_t limit =
        job.common ? fixed.repairs.size() : job.hostile.tests;
    for (std::size_t index = 0; index < limit; ++index) {
        const RepairLiteral& repair = fixed.repairs[index];
        const i64 mass = repair_mass_literal(geometry, repair);
        const i128 margin = static_cast<i128>(63) * mass -
                            static_cast<i128>(4) * geometry.grid;
        ++out.tests;
        out.equalities += margin == 0;
        matrix.add(repair.mask);
        matrix.add(mass);
        add_i128_literal(matrix, margin);
        if (job.common) {
            require(margin >= 0,
                    "claimed-common pair has a literal inactive repair");
            if (!weakest_set || margin < out.weakest_margin ||
                (margin == out.weakest_margin &&
                 repair.mask < out.weakest_mask)) {
                weakest_set = true;
                out.weakest_mask = repair.mask;
                out.weakest_margin = margin;
            }
        } else if (index + 1 < limit) {
            require(margin >= 0,
                    "primary hostile witness is not the first inactive repair");
        } else {
            require(margin < 0 && repair.mask == job.hostile.mask,
                    "literal first inactive mask disagrees with primary");
            require(signed_ratio_equal(
                        margin, geometry.grid, job.hostile.primary_margin,
                        job.hostile.primary_denominator),
                    "literal hostile margin ratio disagrees with primary");
            out.hostile_mask = repair.mask;
            out.hostile_mass = mass;
            out.hostile_margin = margin;
        }
    }
    require(out.tests == limit &&
                (job.common ? weakest_set : out.hostile_mask != 0),
            "literal job did not reach its certificate boundary");
    out.matrix_ledger = matrix.state;
    return out;
}

void audit_batch_compression(const EdgeLiteral& edge,
                             const FixedGeometryLiteral& fixed) {
    const LiteralMassGeometry compressed =
        build_literal_mass_geometry(edge, fixed);
    const Geometry direct = build_joint_geometry(edge.q, edge.r);
    require(compressed.grid == direct.grid &&
                compressed.pair_ticks == direct.pair_ticks &&
                compressed.joint_cells == direct.cells.size(),
            "compressed literal geometry disagrees with full joint walls");
    Fnv ledger;
    for (const RepairLiteral& repair : fixed.repairs) {
        i64 direct_mass = 0;
        for (const Cell& cell : direct.cells) {
            if (cell.pair_safe && (cell.failed_pool & ~repair.mask) == 0)
                direct_mass += cell.width;
        }
        const i64 compressed_mass = repair_mass_literal(compressed, repair);
        require(direct_mass == compressed_mass,
                "batched literal mass disagrees with direct cell sum");
        ledger.add(repair.mask);
        ledger.add(direct_mass);
    }
    std::cout << "BATCH_CONTROL PAIR " << edge.q << ',' << edge.r
              << " JOINT_CELLS " << direct.cells.size()
              << " GRID " << direct.grid << " PAIR_TICKS "
              << direct.pair_ticks << " REPAIRS " << fixed.repairs.size()
              << " MASS_FNV " << std::hex << ledger.state << std::dec
              << " PASS\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 5,
                "usage: literal-post4277 JOINT_DECK COMMON_CSV PRIMARY_HOSTILE_CSV LITERAL_HOSTILE_OUT");
        const std::vector<u32> deck = read_literal_deck(argv[1]);
        const std::vector<EdgeLiteral> common = read_common_literal(argv[2]);
        const std::vector<HostileLiteral> hostile =
            read_hostile_literal(argv[3]);
        const FixedGeometryLiteral fixed = build_fixed_literal(deck);

        std::vector<JobLiteral> jobs;
        jobs.reserve(EXPECTED_RESIDUAL);
        std::size_t common_index = 0;
        std::size_t hostile_index = 0;
        while (common_index < common.size() || hostile_index < hostile.size()) {
            if (hostile_index == hostile.size() ||
                (common_index < common.size() &&
                 common[common_index] < hostile[hostile_index].edge)) {
                jobs.push_back({common[common_index], true, {}});
                ++common_index;
            } else {
                require(common_index == common.size() ||
                            hostile[hostile_index].edge < common[common_index],
                        "common and hostile ledgers overlap");
                jobs.push_back({hostile[hostile_index].edge, false,
                                hostile[hostile_index]});
                ++hostile_index;
            }
        }
        require(jobs.size() == EXPECTED_RESIDUAL,
                "literal job universe count changed");
        std::vector<EdgeLiteral> universe;
        universe.reserve(jobs.size());
        for (const JobLiteral& job : jobs) universe.push_back(job.edge);
        require(edge_fnv_literal(universe) == EXPECTED_RESIDUAL_FNV,
                "literal job universe FNV changed");

        std::cout << "JOINT421_LITERAL_POST4277_V1\n";
        std::cout << "DECK " << deck.size() << " FNV " << std::hex
                  << deck_fnv_literal(deck) << std::dec << '\n';
        std::cout << "UNIVERSE " << universe.size() << " FNV " << std::hex
                  << edge_fnv_literal(universe) << " COMMON " << std::dec
                  << common.size() << " COMMON_FNV " << std::hex
                  << edge_fnv_literal(common) << " HOSTILE " << std::dec
                  << hostile.size() << " PRIMARY_HOSTILE_FNV " << std::hex
                  << EXPECTED_PRIMARY_HOSTILE_FNV << std::dec << '\n';
        std::cout << "FIXED_LITERAL_GRID " << FIXED_GRID << " WALLS "
                  << fixed.walls.size() << " CELLS " << fixed.cells.size()
                  << " INTERVALS " << fixed.interval_count << " FNV "
                  << std::hex << fixed.ledger << std::dec << '\n';

        audit_batch_compression({1, 2}, fixed);
        audit_batch_compression({1, 6}, fixed);
        audit_batch_compression({35, 315}, fixed);
        audit_batch_compression({650, 665}, fixed);

        std::vector<PairResultLiteral> results(jobs.size());
        std::atomic<std::size_t> next{0};
        const unsigned hardware = std::thread::hardware_concurrency();
        const unsigned thread_count =
            std::max(1u, std::min(8u, hardware == 0 ? 1u : hardware));
        std::vector<std::thread> workers;
        for (unsigned lane = 0; lane < thread_count; ++lane) {
            workers.emplace_back([&]() {
                while (true) {
                    const std::size_t index = next.fetch_add(1);
                    if (index >= jobs.size()) break;
                    results[index] = verify_literal_job(jobs[index], fixed);
                }
            });
        }
        for (auto& worker : workers) worker.join();

        std::ofstream hostile_out(argv[4]);
        require(static_cast<bool>(hostile_out),
                "cannot create literal hostile output");
        hostile_out << "q,r,tests,first_inactive_mask,literal_mass,literal_margin,literal_grid,primary_margin,primary_den\n";
        u64 common_pairs = 0;
        u64 hostile_pairs = 0;
        u64 common_tests = 0;
        u64 hostile_tests = 0;
        u64 equalities = 0;
        i128 joint_cells = 0;
        i128 grid_sum = 0;
        i128 pair_ticks = 0;
        Fnv row_ledger;
        Fnv hostile_ledger;
        bool global_min_set = false;
        EdgeLiteral global_min_edge;
        u32 global_min_mask = 0;
        i128 global_min_margin = 0;
        i128 global_min_grid = 1;
        for (std::size_t index = 0; index < jobs.size(); ++index) {
            const JobLiteral& job = jobs[index];
            const PairResultLiteral& result = results[index];
            require(result.tests >= 1 && result.tests <= deck.size(),
                    "literal result test count out of range");
            equalities += result.equalities;
            joint_cells += result.joint_cells;
            grid_sum += result.grid;
            pair_ticks += result.pair_ticks;
            add_edge_literal(row_ledger, job.edge);
            row_ledger.add(job.common);
            row_ledger.add(result.tests);
            row_ledger.add(result.equalities);
            row_ledger.add(result.joint_cells);
            row_ledger.add(result.grid);
            row_ledger.add(result.pair_ticks);
            row_ledger.add(result.weakest_mask);
            add_i128_literal(row_ledger, result.weakest_margin);
            row_ledger.add(result.hostile_mask);
            row_ledger.add(result.hostile_mass);
            add_i128_literal(row_ledger, result.hostile_margin);
            row_ledger.add(result.matrix_ledger);
            if (job.common) {
                ++common_pairs;
                common_tests += result.tests;
                const bool less = global_min_set && rational_less_literal(
                    result.weakest_margin, result.grid,
                    global_min_margin, global_min_grid);
                const bool equal = global_min_set && rational_equal_literal(
                    result.weakest_margin, result.grid,
                    global_min_margin, global_min_grid);
                if (!global_min_set || less ||
                    (equal &&
                     std::tie(job.edge.q, job.edge.r, result.weakest_mask) <
                         std::tie(global_min_edge.q, global_min_edge.r,
                                  global_min_mask))) {
                    global_min_set = true;
                    global_min_edge = job.edge;
                    global_min_mask = result.weakest_mask;
                    global_min_margin = result.weakest_margin;
                    global_min_grid = result.grid;
                }
            } else {
                ++hostile_pairs;
                hostile_tests += result.tests;
                add_edge_literal(hostile_ledger, job.edge);
                hostile_ledger.add(result.tests);
                hostile_ledger.add(result.hostile_mask);
                hostile_ledger.add(result.hostile_mass);
                add_i128_literal(hostile_ledger, result.hostile_margin);
                hostile_ledger.add(result.grid);
                add_i128_literal(hostile_ledger,
                                 job.hostile.primary_margin);
                add_i128_literal(hostile_ledger,
                                 job.hostile.primary_denominator);
                hostile_out << job.edge.q << ',' << job.edge.r << ','
                            << result.tests << ',' << std::hex << std::setw(8)
                            << std::setfill('0') << result.hostile_mask
                            << std::dec << std::setfill(' ') << ','
                            << result.hostile_mass << ','
                            << decimal(result.hostile_margin) << ','
                            << result.grid << ','
                            << decimal(job.hostile.primary_margin) << ','
                            << decimal(job.hostile.primary_denominator) << '\n';
            }
        }
        require(hostile_out.good(), "failed writing literal hostile output");
        require(common_pairs == EXPECTED_COMMON &&
                    hostile_pairs == EXPECTED_HOSTILE &&
                    common_tests == EXPECTED_COMMON_TESTS &&
                    common_tests + hostile_tests == EXPECTED_ALL_TESTS &&
                    equalities == 0 && global_min_set &&
                    fixed.ledger == UINT64_C(0x617e42dc22da8fff) &&
                    hostile_ledger.state == UINT64_C(0xddfd5ea485da4c73) &&
                    row_ledger.state == UINT64_C(0x40fcd8d13bbc0851) &&
                    joint_cells == static_cast<i128>(1453859402) &&
                    grid_sum == parse_signed_i128("29393999537086266001920") &&
                    pair_ticks == parse_signed_i128("21595609799092389736288") &&
                    global_min_edge == EdgeLiteral{35, 315} &&
                    global_min_mask == UINT32_C(0x22560801) &&
                    global_min_margin == static_cast<i128>(32377968) &&
                    global_min_grid == FIXED_GRID,
                "literal matrix census changed");

        std::vector<u32> bodies;
        bodies.reserve(EXPECTED_BODIES);
        enumerate_bodies(0, 9, 0, bodies);
        require(bodies.size() == EXPECTED_BODIES,
                "literal body universe changed");
        const BodyAudit body = audit_bodies(bodies, deck);
        require(body.bodies == EXPECTED_BODIES && body.failures == 0 &&
                    body.checks == UINT64_C(405170384) &&
                    body.max_checks == 421 &&
                    body.worst_body == UINT32_C(0x000dad001),
                "locked joint deck misses a nine-body");

        std::cout << "COMMON_MATRIX PAIRS " << common_pairs << " TESTS "
                  << common_tests << " ACTIVE " << common_tests
                  << " INACTIVE 0\n";
        std::cout << "HOSTILE_CONTROLS PAIRS " << hostile_pairs << " TESTS "
                  << hostile_tests << " FIRST_INACTIVE_MATCH "
                  << hostile_pairs << " RATIO_MATCH " << hostile_pairs
                  << " LITERAL_FNV " << std::hex << hostile_ledger.state
                  << std::dec << '\n';
        std::cout << "ALL_TESTS " << common_tests + hostile_tests
                  << " EQUALITIES " << equalities << " JOINT_CELLS "
                  << decimal(joint_cells) << " GRID_SUM "
                  << decimal(grid_sum) << " PAIR_TICKS_SUM "
                  << decimal(pair_ticks) << " ROW_FNV " << std::hex
                  << row_ledger.state << std::dec << '\n';
        std::cout << "MIN_COMMON_GAP PAIR " << global_min_edge.q << ','
                  << global_min_edge.r << " MASK " << std::hex
                  << global_min_mask << std::dec << " LITERAL_MARGIN "
                  << decimal(global_min_margin) << " GRID "
                  << decimal(global_min_grid) << " REDUCED "
                  << fraction(global_min_margin, global_min_grid)
                  << " TARGET_GAP "
                  << fraction(global_min_margin,
                              static_cast<i128>(63) * global_min_grid)
                  << '\n';
        std::cout << "BODY_SCAN BODIES " << body.bodies << " FAILURES "
                  << body.failures << " CHECKS " << body.checks
                  << " MAX_CHECKS " << body.max_checks << " WORST_BODY "
                  << std::hex << body.worst_body << std::dec << '\n';
        std::cout << "THREADS " << thread_count << '\n';
        std::cout << "VERDICT PASS ALL_148063_X_421_DETACHED_LITERAL_WALL_MATRIX_AND_ALL_24259_FIRST_INACTIVE_CONTROLS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "JOINT421_LITERAL_POST4277_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}

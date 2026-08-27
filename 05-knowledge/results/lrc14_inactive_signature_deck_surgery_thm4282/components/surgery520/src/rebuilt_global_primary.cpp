// Scratch generic deterministic endpoint-cocycle scout for a rebuilt deck.
//
// The source universe is reconstructed from the exact post-THM-4271 edge
// ledger by applying (in order) the THM-4276 high-endpoint deletion and the
// THM-4277 rectangle deletion.  At each surviving pair, repairs are tested in
// their locked file order and evaluation stops at the first inactive repair.
//
// For speed, the fixed-pool atom condition F(cell) subset repair is compiled
// once into contiguous cell intervals.  Pair mass is then the sum of exact
// endpoint-cocycle primitive-prefix differences on those intervals.  This is
// identically atom_mass_cegar(build_cocycle_atoms(...), repair): a rank-eight
// repair cannot contain a fixed-pool failure atom of arity greater than eight.

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#endif
#define CARRIER_CEGAR_LIBRARY_ONLY
#include "04-computation/lrc14_three_round_learned_carrier_thm4266/carrier_cegar_descent.cpp"
#undef CARRIER_CEGAR_LIBRARY_ONLY
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include <atomic>
#include <fstream>
#include <map>
#include <mutex>

namespace {

constexpr std::size_t EXPECTED_POST4271 = 174904;
constexpr u64 EXPECTED_POST4271_FNV = UINT64_C(0xb3db855040bcf19e);
constexpr std::size_t EXPECTED_POST4276 = 174741;
constexpr u64 EXPECTED_POST4276_FNV = UINT64_C(0xf13745b05320f83c);
constexpr std::size_t EXPECTED_RECTANGLE = 2419;
constexpr u64 EXPECTED_RECTANGLE_FNV = UINT64_C(0x67b373dc22ac870d);
constexpr std::size_t EXPECTED_POST4277 = 172322;
constexpr u64 EXPECTED_POST4277_FNV = UINT64_C(0x30b2a7e597ac548c);

struct EdgeScout {
    int q = 0;
    int r = 0;
    auto operator<=>(const EdgeScout&) const = default;
};

void add_edge(FnvLocal& ledger, const EdgeScout& edge) {
    ledger.add(static_cast<u64>(edge.q));
    ledger.add(static_cast<u64>(edge.r));
}

u64 edge_fnv(const std::vector<EdgeScout>& edges) {
    FnvLocal ledger;
    for (const EdgeScout& edge : edges) add_edge(ledger, edge);
    return ledger.state;
}

u64 mask_fnv_scout(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

std::vector<EdgeScout> read_post4271(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open post-THM-4271 ledger");
    std::vector<EdgeScout> edges;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::istringstream row(line);
        EdgeScout edge;
        char comma = 0;
        char trailing = 0;
        row >> edge.q >> comma >> edge.r;
        require(static_cast<bool>(row) && comma == ',' && edge.q > 0 &&
                    edge.q < edge.r && !(row >> trailing),
                "bad post-THM-4271 edge row");
        edges.push_back(edge);
    }
    require(std::is_sorted(edges.begin(), edges.end()) &&
                std::adjacent_find(edges.begin(), edges.end()) == edges.end(),
            "post-THM-4271 edge order/distinctness changed");
    require(edges.size() == EXPECTED_POST4271 &&
                edge_fnv(edges) == EXPECTED_POST4271_FNV,
            "post-THM-4271 edge ledger changed");
    return edges;
}

struct Residuals {
    std::vector<EdgeScout> post4276;
    std::vector<EdgeScout> rectangle;
    std::vector<EdgeScout> post4277;
};

Residuals reconstruct_post4277(const std::vector<EdgeScout>& post4271) {
    Residuals out;
    out.post4276.reserve(EXPECTED_POST4276);
    for (const EdgeScout& edge : post4271) {
        const bool retained_boundary =
            (edge.q == 256 || edge.q == 384) && edge.r == 670;
        if (edge.r < 670 || retained_boundary) out.post4276.push_back(edge);
    }
    require(out.post4276.size() == EXPECTED_POST4276 &&
                edge_fnv(out.post4276) == EXPECTED_POST4276_FNV,
            "post-THM-4276 residual changed");
    out.rectangle.reserve(EXPECTED_RECTANGLE);
    out.post4277.reserve(EXPECTED_POST4277);
    for (const EdgeScout& edge : out.post4276) {
        if (450 <= edge.q && edge.q <= 499 &&
            600 <= edge.r && edge.r <= 650) {
            out.rectangle.push_back(edge);
        } else {
            out.post4277.push_back(edge);
        }
    }
    require(out.rectangle.size() == EXPECTED_RECTANGLE &&
                edge_fnv(out.rectangle) == EXPECTED_RECTANGLE_FNV,
            "THM-4277 rectangle ledger changed");
    require(out.post4277.size() == EXPECTED_POST4277 &&
                edge_fnv(out.post4277) == EXPECTED_POST4277_FNV,
            "post-THM-4277 residual changed");
    return out;
}

std::vector<u32> read_deck_scout(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open rebuilt deck");
    std::vector<u32> deck;
    std::set<u32> distinct;
    std::string word;
    while (input >> word) {
        std::size_t used = 0;
        const u64 wide = std::stoull(word, &used, 16);
        require(used == word.size() && wide < (UINT64_C(1) << 30),
                "bad joint repair token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "joint repair arity/distinctness changed");
        deck.push_back(mask);
    }
    require(!deck.empty() && deck.size() <= 1000,
            "rebuilt-deck size outside scratch scout bound");
    return deck;
}

struct IndexedRepair {
    u32 mask = 0;
    std::vector<std::pair<std::uint16_t, std::uint16_t>> intervals;
};

IndexedRepair index_repair(u32 repair, const std::vector<Cell>& cells) {
    IndexedRepair out;
    out.mask = repair;
    std::size_t index = 0;
    while (index < cells.size()) {
        if ((cells[index].failed_pool & ~repair) != 0) {
            ++index;
            continue;
        }
        const std::size_t begin = index;
        do {
            ++index;
        } while (index < cells.size() &&
                 (cells[index].failed_pool & ~repair) == 0);
        require(begin < UINT16_MAX && index < UINT16_MAX,
                "pool wall index exceeds u16");
        out.intervals.push_back({static_cast<std::uint16_t>(begin),
                                 static_cast<std::uint16_t>(index)});
    }
    require(!out.intervals.empty(), "repair has empty fixed-pool geometry");
    return out;
}

// Exact comparison a/b < c/d for nonnegative numerators without forming the
// potentially overflowing cross product.
bool rational_less_nonnegative(i128 a, i128 b, i128 c, i128 d) {
    require(a >= 0 && c >= 0 && b > 0 && d > 0,
            "invalid nonnegative rational comparison");
    bool reversed = false;
    while (true) {
        const i128 aq = a / b;
        const i128 cq = c / d;
        if (aq != cq) return reversed ? aq > cq : aq < cq;
        a %= b;
        c %= d;
        if (a == 0 || c == 0) {
            if (a == 0 && c == 0) return false;
            const bool raw_less = a == 0;
            return reversed ? !raw_less : raw_less;
        }
        std::swap(a, b);
        std::swap(c, d);
        reversed = !reversed;
    }
}

bool rational_equal_nonnegative(i128 a, i128 b, i128 c, i128 d) {
    return !rational_less_nonnegative(a, b, c, d) &&
           !rational_less_nonnegative(c, d, a, b);
}

struct PairScout {
    u64 tests = 0;
    u64 equalities = 0;
    bool common = false;
    u32 first_inactive = 0;
    i128 first_inactive_margin = 0;
    i128 denominator = 0;
    u32 min_common_mask = 0;
    i128 min_common_margin = 0;
    u64 primitive_cells = 0;
    u64 wall_queries = 0;
};

struct WorkerCache {
    std::vector<i128> prefix;
    std::vector<u32> stamp;
    u32 generation = 0;
};

PairScout test_pair(const EdgeScout& edge,
                    const std::vector<IndexedRepair>& repairs,
                    const std::vector<i64>& walls, WorkerCache& cache) {
    const i64 g = std::gcd(edge.q, edge.r);
    const PrimitivePair primitive = build_primitive(edge.q / g, edge.r / g);
    const i128 denominator =
        static_cast<i128>(primitive.grid) * g * COMMON;
    ++cache.generation;
    require(cache.generation != 0, "worker cache generation overflow");
    PairScout out;
    out.denominator = denominator;
    out.primitive_cells = primitive.cells.size();
    bool common_min_set = false;
    auto prefix_at = [&](std::size_t index) -> i128 {
        require(index < walls.size(), "wall lookup out of range");
        if (cache.stamp[index] != cache.generation) {
            cache.prefix[index] =
                actual_safe_prefix_num(primitive, g, walls[index]);
            cache.stamp[index] = cache.generation;
            ++out.wall_queries;
        }
        return cache.prefix[index];
    };
    for (const IndexedRepair& repair : repairs) {
        i128 mass = 0;
        for (const auto [left, right] : repair.intervals)
            mass += prefix_at(right) - prefix_at(left);
        const i128 margin = static_cast<i128>(63) * mass -
                            static_cast<i128>(4) * denominator;
        ++out.tests;
        out.equalities += margin == 0;
        if (margin < 0) {
            out.first_inactive = repair.mask;
            out.first_inactive_margin = margin;
            return out;
        }
        if (!common_min_set || rational_less_nonnegative(
                margin, denominator, out.min_common_margin, denominator) ||
            (margin == out.min_common_margin &&
             repair.mask < out.min_common_mask)) {
            common_min_set = true;
            out.min_common_mask = repair.mask;
            out.min_common_margin = margin;
        }
    }
    require(common_min_set && out.tests == repairs.size(),
            "common pair did not test full deck");
    out.common = true;
    return out;
}

struct LayerScout {
    u64 residual = 0;
    u64 common = 0;
    u64 tests = 0;
    u64 equalities = 0;
    u64 primitive_cells = 0;
    u64 wall_queries = 0;
};

void add_i128_scout(FnvLocal& ledger, i128 value) {
    const __uint128_t bits = static_cast<__uint128_t>(value);
    ledger.add(static_cast<u64>(bits));
    ledger.add(static_cast<u64>(bits >> 64));
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 6 || argc == 8,
                "usage: scout POST4271 JOINT_DECK COMMON_OUT LAYERS_OUT HOSTILE_OUT [BEGIN END]");
        const std::vector<EdgeScout> post4271 = read_post4271(argv[1]);
        const Residuals residuals = reconstruct_post4277(post4271);
        const std::vector<u32> deck = read_deck_scout(argv[2]);
        std::size_t begin = 0;
        std::size_t end = residuals.post4277.size();
        if (argc == 8) {
            begin = std::stoull(argv[6]);
            end = std::stoull(argv[7]);
        }
        require(begin < end && end <= residuals.post4277.size(),
                "invalid residual half-open slice");

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "fixed-pool cell count changed");
        std::vector<i64> walls(cells.size() + 1);
        walls.front() = cells.front().left;
        for (std::size_t index = 0; index < cells.size(); ++index) {
            require(cells[index].left == walls[index],
                    "fixed-pool cells are not contiguous");
            walls[index + 1] = cells[index].right;
        }
        require(walls.front() == 0 && walls.back() == COMMON,
                "fixed-pool wall endpoints changed");
        std::vector<IndexedRepair> repairs;
        repairs.reserve(deck.size());
        u64 intervals = 0;
        FnvLocal geometry_ledger;
        for (u32 mask : deck) {
            repairs.push_back(index_repair(mask, cells));
            geometry_ledger.add(mask);
            geometry_ledger.add(repairs.back().intervals.size());
            intervals += repairs.back().intervals.size();
            for (const auto [left, right] : repairs.back().intervals) {
                geometry_ledger.add(left);
                geometry_ledger.add(right);
            }
        }

        const std::size_t count = end - begin;
        std::vector<PairScout> result(count);
        std::atomic<std::size_t> next{0};
        const unsigned hardware = std::thread::hardware_concurrency();
        const unsigned thread_count =
            std::max(1u, std::min(8u, hardware == 0 ? 1u : hardware));
        std::vector<std::thread> workers;
        std::mutex progress_mutex;
        std::atomic<u64> done{0};
        for (unsigned lane = 0; lane < thread_count; ++lane) {
            workers.emplace_back([&]() {
                WorkerCache cache;
                cache.prefix.resize(walls.size());
                cache.stamp.assign(walls.size(), 0);
                while (true) {
                    const std::size_t local = next.fetch_add(1);
                    if (local >= count) break;
                    result[local] = test_pair(
                        residuals.post4277[begin + local], repairs, walls,
                        cache);
                    const u64 completed = done.fetch_add(1) + 1;
                    if (completed % 8192 == 0 || completed == count) {
                        std::lock_guard<std::mutex> lock(progress_mutex);
                        std::cerr << "PROGRESS COMPLETED " << completed
                                  << " OF " << count << '\n';
                    }
                }
            });
        }
        for (auto& worker : workers) worker.join();

        std::ofstream common_out(argv[3]);
        require(static_cast<bool>(common_out), "cannot create common-pair output");
        std::ofstream hostile_out(argv[5]);
        require(static_cast<bool>(hostile_out),
                "cannot create hostile-control output");
        hostile_out << "q,r,tests,first_inactive_mask,margin_num,den\n";
        std::map<int, LayerScout> layers;
        u64 common = 0;
        u64 tests = 0;
        u64 equalities = 0;
        u64 primitive_cells = 0;
        u64 wall_queries = 0;
        FnvLocal common_ledger;
        FnvLocal hostile_ledger;
        FnvLocal row_ledger;
        bool global_min_set = false;
        EdgeScout global_min_edge;
        u32 global_min_mask = 0;
        i128 global_min_margin = 0;
        i128 global_min_denominator = 1;
        for (std::size_t local = 0; local < count; ++local) {
            const EdgeScout edge = residuals.post4277[begin + local];
            const PairScout& scan = result[local];
            require(scan.tests >= 1 && scan.tests <= deck.size(),
                    "pair test count out of range");
            LayerScout& layer = layers[edge.r];
            ++layer.residual;
            layer.common += scan.common;
            layer.tests += scan.tests;
            layer.equalities += scan.equalities;
            layer.primitive_cells += scan.primitive_cells;
            layer.wall_queries += scan.wall_queries;
            tests += scan.tests;
            equalities += scan.equalities;
            primitive_cells += scan.primitive_cells;
            wall_queries += scan.wall_queries;
            add_edge(row_ledger, edge);
            row_ledger.add(scan.tests);
            row_ledger.add(scan.equalities);
            row_ledger.add(scan.common);
            row_ledger.add(scan.first_inactive);
            add_i128_scout(row_ledger, scan.first_inactive_margin);
            add_i128_scout(row_ledger, scan.denominator);
            row_ledger.add(scan.min_common_mask);
            add_i128_scout(row_ledger, scan.min_common_margin);
            row_ledger.add(scan.primitive_cells);
            row_ledger.add(scan.wall_queries);
            if (!scan.common) {
                add_edge(hostile_ledger, edge);
                hostile_ledger.add(scan.tests);
                hostile_ledger.add(scan.first_inactive);
                add_i128_scout(hostile_ledger, scan.first_inactive_margin);
                add_i128_scout(hostile_ledger, scan.denominator);
                hostile_out << edge.q << ',' << edge.r << ',' << scan.tests
                            << ',' << std::hex << std::setw(8)
                            << std::setfill('0') << scan.first_inactive
                            << std::dec << std::setfill(' ') << ','
                            << decimal(scan.first_inactive_margin) << ','
                            << decimal(scan.denominator) << '\n';
                continue;
            }
            ++common;
            add_edge(common_ledger, edge);
            common_out << edge.q << ',' << edge.r << '\n';
            const bool less = global_min_set && rational_less_nonnegative(
                scan.min_common_margin, scan.denominator,
                global_min_margin, global_min_denominator);
            const bool equal = global_min_set && rational_equal_nonnegative(
                scan.min_common_margin, scan.denominator,
                global_min_margin, global_min_denominator);
            if (!global_min_set || less ||
                (equal && std::tie(edge.q, edge.r, scan.min_common_mask) <
                              std::tie(global_min_edge.q, global_min_edge.r,
                                       global_min_mask))) {
                global_min_set = true;
                global_min_edge = edge;
                global_min_mask = scan.min_common_mask;
                global_min_margin = scan.min_common_margin;
                global_min_denominator = scan.denominator;
            }
        }
        require(common_out.good(), "failed writing common-pair output");
        require(hostile_out.good(), "failed writing hostile-control output");

        std::ofstream layer_out(argv[4]);
        require(static_cast<bool>(layer_out), "cannot create layer output");
        layer_out << "endpoint,residual_pairs,common_active_pairs,repair_tests,equalities,primitive_cells,wall_queries\n";
        for (const auto& [endpoint, layer] : layers) {
            layer_out << endpoint << ',' << layer.residual << ','
                      << layer.common << ',' << layer.tests << ','
                      << layer.equalities << ',' << layer.primitive_cells
                      << ',' << layer.wall_queries << '\n';
        }
        require(layer_out.good(), "failed writing layer output");

        FnvLocal slice_input_ledger;
        for (std::size_t index = begin; index < end; ++index)
            add_edge(slice_input_ledger, residuals.post4277[index]);
        std::cout << "THM4281_REBUILT_POST4277_PRIMARY_ATOM_SCOUT_V1\n";
        std::cout << "POST_THM4271 " << post4271.size() << " FNV "
                  << std::hex << edge_fnv(post4271) << std::dec << '\n';
        std::cout << "POST_THM4276 " << residuals.post4276.size() << " FNV "
                  << std::hex << edge_fnv(residuals.post4276) << std::dec << '\n';
        std::cout << "THM4277_RECTANGLE " << residuals.rectangle.size()
                  << " FNV " << std::hex << edge_fnv(residuals.rectangle)
                  << std::dec << '\n';
        std::cout << "POST_THM4277 " << residuals.post4277.size() << " FNV "
                  << std::hex << edge_fnv(residuals.post4277) << std::dec
                  << '\n';
        std::cout << "JOINT_DECK " << deck.size() << " FNV " << std::hex
                  << mask_fnv_scout(deck) << std::dec << '\n';
        std::cout << "ATOM_SUPPORT FIXED_CELLS " << cells.size()
                  << " WALLS " << walls.size() << " INTERVALS " << intervals
                  << " GEOMETRY_FNV " << std::hex << geometry_ledger.state
                  << std::dec << '\n';
        std::cout << "SLICE BEGIN " << begin << " END " << end << " COUNT "
                  << count << " INPUT_FNV " << std::hex
                  << slice_input_ledger.state << std::dec << '\n';
        std::cout << "SCOUT COMMON_ACTIVE " << common << " COMMON_FNV "
                  << std::hex << common_ledger.state << " ROW_FNV "
                  << row_ledger.state << std::dec << " REPAIR_TESTS " << tests
                  << " EQUALITIES " << equalities << " PRIMITIVE_CELLS "
                  << primitive_cells << " WALL_QUERIES " << wall_queries
                  << '\n';
        std::cout << "HOSTILE_CONTROLS " << count - common << " FNV "
                  << std::hex << hostile_ledger.state << std::dec << '\n';
        if (global_min_set) {
            std::cout << "MIN_COMMON_GAP PAIR " << global_min_edge.q << ','
                      << global_min_edge.r << " MASK " << std::hex
                      << global_min_mask << std::dec << " MARGIN_NUM "
                      << decimal(global_min_margin) << " DEN "
                      << decimal(global_min_denominator) << '\n';
        } else {
            std::cout << "MIN_COMMON_GAP NONE\n";
        }
        std::cout << "LAYERS " << layers.size() << " THREADS " << thread_count
                  << '\n';
        std::cout << "VERDICT PASS FINITE_EXACT_PRIMARY_ATOM_SCOUT_ONLY\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4281_REBUILT_POST4277_PRIMARY_SCOUT_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}

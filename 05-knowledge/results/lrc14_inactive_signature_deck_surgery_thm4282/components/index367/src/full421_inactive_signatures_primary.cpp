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
#include <iomanip>
#include <mutex>

namespace {

constexpr std::size_t EXPECTED_FINAL = 24223;
constexpr u64 EXPECTED_FINAL_FNV = UINT64_C(0x80ec0687d8c7dba7);
constexpr std::size_t EXPECTED_DECK = 421;
constexpr u64 EXPECTED_DECK_FNV = UINT64_C(0x20d63dd42fe8150e);
constexpr std::size_t SIGNATURE_WORDS = (EXPECTED_DECK + 63) / 64;

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

std::vector<u32> read_deck_scout(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open locked joint deck");
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
    require(deck.size() == EXPECTED_DECK &&
                mask_fnv_scout(deck) == EXPECTED_DECK_FNV,
            "locked joint-deck ledger changed");
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

struct WorkerCache {
    std::vector<i128> prefix;
    std::vector<u32> stamp;
    u32 generation = 0;
};

std::vector<EdgeScout> read_final(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open post-THM-4281 residual");
    std::vector<EdgeScout> edges;
    std::string line;
    while (std::getline(input, line)) {
        require(!line.empty(), "empty final residual row");
        std::istringstream row(line);
        EdgeScout edge;
        char comma = 0;
        char trailing = 0;
        row >> edge.q >> comma >> edge.r;
        require(static_cast<bool>(row) && comma == ',' && edge.q > 0 &&
                    edge.q < edge.r && !(row >> trailing),
                "malformed final residual row");
        edges.push_back(edge);
    }
    require(std::is_sorted(edges.begin(), edges.end()) &&
                std::adjacent_find(edges.begin(), edges.end()) == edges.end() &&
                edges.size() == EXPECTED_FINAL &&
                edge_fnv(edges) == EXPECTED_FINAL_FNV,
            "post-THM-4281 residual identity changed");
    return edges;
}

struct SignatureRow {
    std::array<u64, SIGNATURE_WORDS> inactive{};
    u64 inactive_count = 0;
    u64 equalities = 0;
    u64 primitive_cells = 0;
    u64 wall_queries = 0;
};

SignatureRow full_signature(const EdgeScout& edge,
                            const std::vector<IndexedRepair>& repairs,
                            const std::vector<i64>& walls,
                            WorkerCache& cache) {
    const i64 g = std::gcd(edge.q, edge.r);
    const PrimitivePair primitive = build_primitive(edge.q / g, edge.r / g);
    const i128 denominator =
        static_cast<i128>(primitive.grid) * g * COMMON;
    ++cache.generation;
    require(cache.generation != 0, "worker cache generation overflow");
    SignatureRow out;
    out.primitive_cells = primitive.cells.size();
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
    for (std::size_t index = 0; index < repairs.size(); ++index) {
        i128 mass = 0;
        for (const auto [left, right] : repairs[index].intervals)
            mass += prefix_at(right) - prefix_at(left);
        const i128 margin = static_cast<i128>(63) * mass -
                            static_cast<i128>(4) * denominator;
        out.equalities += margin == 0;
        if (margin < 0) {
            out.inactive[index / 64] |= UINT64_C(1) << (index % 64);
            ++out.inactive_count;
        }
    }
    require(out.inactive_count > 0, "final residual unexpectedly common");
    require((out.inactive.back() >> (EXPECTED_DECK % 64)) == 0,
            "signature high padding bits set");
    return out;
}

void add_signature(FnvLocal& ledger, const SignatureRow& row) {
    for (u64 word : row.inactive) ledger.add(word);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 4,
                "usage: full-signatures FINAL_RESIDUAL JOINT_DECK CSV_OUT");
        const std::vector<EdgeScout> edges = read_final(argv[1]);
        const std::vector<u32> deck = read_deck_scout(argv[2]);
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
        for (u32 mask : deck) repairs.push_back(index_repair(mask, cells));

        std::vector<SignatureRow> results(edges.size());
        std::atomic<std::size_t> next{0};
        std::atomic<u64> done{0};
        std::mutex progress_mutex;
        const unsigned hardware = std::thread::hardware_concurrency();
        const unsigned threads =
            std::max(1u, std::min(8u, hardware == 0 ? 1u : hardware));
        std::vector<std::thread> workers;
        for (unsigned lane = 0; lane < threads; ++lane) {
            workers.emplace_back([&]() {
                WorkerCache cache;
                cache.prefix.resize(walls.size());
                cache.stamp.assign(walls.size(), 0);
                while (true) {
                    const std::size_t row = next.fetch_add(1);
                    if (row >= edges.size()) break;
                    results[row] = full_signature(edges[row], repairs, walls, cache);
                    const u64 completed = done.fetch_add(1) + 1;
                    if (completed % 2048 == 0 || completed == edges.size()) {
                        std::lock_guard<std::mutex> lock(progress_mutex);
                        std::cerr << "PROGRESS " << completed << " OF "
                                  << edges.size() << '\n';
                    }
                }
            });
        }
        for (auto& worker : workers) worker.join();

        std::ofstream output(argv[3]);
        require(static_cast<bool>(output), "cannot create signature CSV");
        output << "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6\n";
        FnvLocal signature_ledger;
        FnvLocal row_ledger;
        u64 inactive_total = 0;
        u64 equalities = 0;
        u64 primitive_cells = 0;
        u64 wall_queries = 0;
        for (std::size_t row = 0; row < edges.size(); ++row) {
            const auto& edge = edges[row];
            const auto& result = results[row];
            add_edge(row_ledger, edge);
            row_ledger.add(result.inactive_count);
            add_signature(row_ledger, result);
            add_signature(signature_ledger, result);
            inactive_total += result.inactive_count;
            equalities += result.equalities;
            primitive_cells += result.primitive_cells;
            wall_queries += result.wall_queries;
            output << edge.q << ',' << edge.r << ',' << result.inactive_count;
            for (u64 word : result.inactive) {
                output << ',' << std::hex << std::setw(16)
                       << std::setfill('0') << word << std::dec
                       << std::setfill(' ');
            }
            output << '\n';
        }
        require(output.good(), "failed writing signature CSV");
        std::cout << "THM4281_FULL421_INACTIVE_SIGNATURES_PRIMARY_V1\n"
                  << "FINAL_RESIDUAL " << edges.size() << " FNV "
                  << std::hex << edge_fnv(edges) << std::dec << '\n'
                  << "JOINT_DECK " << deck.size() << " FNV " << std::hex
                  << mask_fnv_scout(deck) << std::dec << '\n'
                  << "SIGNATURE_WORDS " << SIGNATURE_WORDS
                  << " INACTIVE_INCIDENCES " << inactive_total
                  << " EQUALITIES " << equalities
                  << " SIGNATURE_FNV " << std::hex << signature_ledger.state
                  << " ROW_FNV " << row_ledger.state << std::dec << '\n'
                  << "PRIMITIVE_CELLS " << primitive_cells
                  << " WALL_QUERIES " << wall_queries
                  << " THREADS " << threads << '\n'
                  << "VERDICT PASS FINITE_EXACT_FULL_SIGNATURES_PRIMARY\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FULL_SIGNATURE_ERROR " << error.what() << '\n';
        return 1;
    }
}

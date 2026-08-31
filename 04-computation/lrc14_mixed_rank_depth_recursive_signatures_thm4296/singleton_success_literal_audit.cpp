// Detached literal-wall and private-body audit for the exhaustive successful
// singleton-ideal one-replacement census.  This does not import the primary
// endpoint-cocycle candidate engine.

#define SINGLETON_LITERAL_LIBRARY_ONLY
#include "04-computation/lrc14_signature_response_congruence_thm4286/singleton_fibre_literal_verify.cpp"
#undef SINGLETON_LITERAL_LIBRARY_ONLY

#include <atomic>
#include <map>
#include <thread>

namespace {

using PairAudit = std::pair<int, int>;

struct NodeAudit {
    std::size_t index = 0;
    u32 old_mask = 0;
    u32 witness = 0;
    std::size_t expected_rows = 0;
    u64 expected_ideal_fnv = 0;
    std::size_t expected_obligations = 0;
    u64 expected_obligation_fnv = 0;
    std::size_t literal_rows = 0;
    FnvLocal ideal_ledger;
    FnvLocal literal_ledger;
    u64 equalities = 0;
    bool have_minimum = false;
    i128 minimum_margin = 0;
    i128 minimum_grid = 1;
    PairAudit minimum_pair{};
    i128 literal_gap_num = 0;
    i128 literal_gap_den = 1;
    std::vector<std::string> source_fields;
};

std::map<std::size_t, NodeAudit> read_nodes_audit(const std::string& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open witness ledger");
    std::map<std::size_t, NodeAudit> nodes;
    std::string line;
    while (std::getline(in, line)) {
        const auto f = split_csv(line);
        require(f.size() == 17, "bad witness-ledger row");
        NodeAudit node;
        node.index = std::stoul(f[0]);
        node.old_mask = static_cast<u32>(std::stoul(f[1], nullptr, 16));
        node.witness = static_cast<u32>(std::stoul(f[2], nullptr, 16));
        node.expected_rows = std::stoul(f[3]);
        node.expected_ideal_fnv = std::stoull(f[4], nullptr, 16);
        node.expected_obligations = std::stoul(f[5]);
        node.expected_obligation_fnv = std::stoull(f[6], nullptr, 16);
        node.source_fields = f;
        require(node.index < 421 && std::popcount(node.old_mask) == 8 &&
                    std::popcount(node.witness) == 8 &&
                    nodes.emplace(node.index, std::move(node)).second,
                "bad/duplicate node");
    }
    require(nodes.size() == 110, "successful node count changed");
    return nodes;
}

std::set<PairAudit> read_live_audit(const std::string& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open live residual");
    std::set<PairAudit> live;
    std::string line;
    while (std::getline(in, line)) {
        const auto f = split_csv(line);
        require(f.size() == 2 &&
                    live.insert({std::stoi(f[0]), std::stoi(f[1])}).second,
                "bad/duplicate live row");
    }
    require(live.size() == 22647, "live residual count changed");
    return live;
}

std::size_t only_index_audit(const std::array<u64, 7>& words) {
    std::size_t found = 421;
    unsigned count = 0;
    for (std::size_t word = 0; word < words.size(); ++word) {
        count += std::popcount(words[word]);
        if (words[word]) found = 64 * word + std::countr_zero(words[word]);
    }
    return count == 1 ? found : 421;
}

u32 next_audit(u32 value) {
    const u32 low = value & (0u - value);
    const u32 ripple = value + low;
    return ripple | (((value ^ ripple) >> 2) / low);
}

void write_literal_witness_ledger(
    const std::string& path, const std::map<std::size_t, NodeAudit>& nodes) {
    std::ofstream out(
        path, std::ios::out | std::ios::binary | std::ios::trunc);
    require(static_cast<bool>(out), "cannot open literal witness output");
    for (const auto& [index, node] : nodes) {
        const auto& fields = node.source_fields;
        require(fields.size() == 17 && std::stoul(fields[0]) == index,
                "witness source fields changed");
        for (std::size_t field = 0; field < 10; ++field) {
            if (field) out << ',';
            out << fields[field];
        }
        out << ',' << node.minimum_pair.first << ',' << node.minimum_pair.second
            << ',' << decimal(node.literal_gap_num) << '/'
            << decimal(node.literal_gap_den) << ',' << std::hex
            << node.literal_ledger.state << std::dec << ',' << fields[14]
            << ',' << fields[15] << ',' << fields[16] << '\n';
    }
    require(static_cast<bool>(out), "cannot write literal witness output");
}

} // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        if (argc == 2 && std::string(argv[1]) == "--help") {
            std::cout
                << "usage: singleton-success-audit DECK SIGNATURES LIVE "
                   "WITNESSES [--literal-witness-output PATH]\n"
                << "audits the primary witness ledger and optionally writes "
                   "the corrected literal ledger\n";
            return 0;
        }
        require(argc == 5 ||
                    (argc == 7 &&
                     std::string(argv[5]) == "--literal-witness-output"),
                "usage: singleton-success-audit DECK SIGNATURES LIVE "
                "WITNESSES [--literal-witness-output PATH]");
        const std::vector<u32> deck = read_deck(argv[1]);
        const std::vector<SigRow> signatures = read_signatures(argv[2]);
        const std::set<PairAudit> live = read_live_audit(argv[3]);
        std::map<std::size_t, NodeAudit> nodes = read_nodes_audit(argv[4]);
        for (auto& [index, node] : nodes) {
            require(deck[index] == node.old_mask &&
                        std::find(deck.begin(), deck.end(), node.witness) ==
                            deck.end(),
                    "old/witness mask identity changed");
        }

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        std::vector<i64> walls(cells.size() + 1);
        walls.front() = cells.front().left;
        for (std::size_t i = 0; i < cells.size(); ++i) {
            require(cells[i].left == walls[i], "pool cells not contiguous");
            walls[i + 1] = cells[i].right;
        }
        require(walls.front() == 0 && walls.back() == COMMON,
                "pool walls changed");
        std::vector<IndexedRepair> old_geometry;
        old_geometry.reserve(deck.size());
        for (u32 mask : deck) old_geometry.push_back(index_repair(mask, cells));
        std::map<std::size_t, IndexedRepair> witness_geometry;
        for (const auto& [index, node] : nodes)
            witness_geometry.emplace(index, index_repair(node.witness, cells));

        u64 total_literal_rows = 0;
        u64 total_old_tests = 0;
        u64 total_replacement_tests = 0;
        u64 total_equalities = 0;
        FnvLocal global_literal_ledger;
        for (const SigRow& row : signatures) {
            const PairAudit pair_key{row.q, row.r};
            if (!live.contains(pair_key)) continue;
            const std::size_t index = only_index_audit(row.words);
            auto node_it = nodes.find(index);
            if (node_it == nodes.end()) continue;
            NodeAudit& node = node_it->second;
            const LiteralPair pair = build_literal_pair(row.q, row.r);
            std::vector<i128> prefix(walls.size());
            for (std::size_t wall = 0; wall < walls.size(); ++wall)
                prefix[wall] = prefix_at(
                    pair, static_cast<i128>(walls[wall]) * pair.pool_scale);
            std::size_t inactive_count = 0;
            std::size_t inactive_index = deck.size();
            for (std::size_t old = 0; old < old_geometry.size(); ++old) {
                const i128 value =
                    repair_margin(old_geometry[old], prefix, pair.grid);
                ++total_old_tests;
                total_equalities += value == 0;
                node.equalities += value == 0;
                if (value < 0) {
                    ++inactive_count;
                    inactive_index = old;
                }
            }
            const i128 replacement_margin = repair_margin(
                witness_geometry.at(index), prefix, pair.grid);
            ++total_replacement_tests;
            total_equalities += replacement_margin == 0;
            node.equalities += replacement_margin == 0;
            require(inactive_count == 1 && inactive_index == index &&
                        replacement_margin > 0,
                    "literal singleton/replacement sign failed");
            if (!node.have_minimum || fraction_less_exact(
                    replacement_margin, pair.grid,
                    node.minimum_margin, node.minimum_grid)) {
                node.have_minimum = true;
                node.minimum_margin = replacement_margin;
                node.minimum_grid = pair.grid;
                node.minimum_pair = pair_key;
            }
            ++node.literal_rows;
            node.ideal_ledger.add(row.q);
            node.ideal_ledger.add(row.r);
            node.literal_ledger.add(index);
            node.literal_ledger.add(node.old_mask);
            node.literal_ledger.add(node.witness);
            node.literal_ledger.add(row.q);
            node.literal_ledger.add(row.r);
            add_i128(node.literal_ledger, replacement_margin);
            add_i128(node.literal_ledger, pair.grid);
            global_literal_ledger.add(index);
            global_literal_ledger.add(row.q);
            global_literal_ledger.add(row.r);
            add_i128(global_literal_ledger, replacement_margin);
            add_i128(global_literal_ledger, pair.grid);
            ++total_literal_rows;
        }

        for (auto& [index, node] : nodes) {
            require(node.literal_rows == node.expected_rows &&
                        node.ideal_ledger.state == node.expected_ideal_fnv &&
                        node.equalities == 0 && node.have_minimum,
                    "node literal/ideal summary changed");
            const i128 divisor = gcd_i128(
                node.minimum_margin,
                static_cast<i128>(63) * node.minimum_grid);
            node.literal_gap_num = node.minimum_margin / divisor;
            node.literal_gap_den =
                (static_cast<i128>(63) * node.minimum_grid) / divisor;
        }
        require(total_literal_rows == 1219 && total_equalities == 0,
                "global literal census changed");

        constexpr std::size_t body_count = 14307150;
        std::vector<u32> bodies(body_count);
        u32 body = (u32{1} << 9) - 1;
        for (std::size_t ordinal = 0; ordinal < bodies.size(); ++ordinal) {
            bodies[ordinal] = body;
            if (ordinal + 1 < bodies.size()) body = next_audit(body);
        }
        require(bodies.back() == UINT32_C(0x3fe00000),
                "body enumeration endpoint changed");
        std::vector<std::int16_t> singleton_owner(body_count, -1);
        std::atomic<u64> uncovered{0};
        constexpr unsigned thread_count = 8;
        std::vector<std::thread> workers;
        for (unsigned lane = 0; lane < thread_count; ++lane) {
            workers.emplace_back([&, lane]() {
                const std::size_t begin = body_count * lane / thread_count;
                const std::size_t end = body_count * (lane + 1) / thread_count;
                for (std::size_t ordinal = begin; ordinal < end; ++ordinal) {
                    unsigned hits = 0;
                    std::size_t owner = deck.size();
                    for (std::size_t index = 0; index < deck.size(); ++index) {
                        if ((bodies[ordinal] & deck[index]) != 0) continue;
                        owner = index;
                        if (++hits == 2) break;
                    }
                    if (hits == 0) ++uncovered;
                    else if (hits == 1)
                        singleton_owner[ordinal] = static_cast<std::int16_t>(owner);
                }
            });
        }
        for (auto& worker : workers) worker.join();
        require(uncovered.load() == 0, "original deck fails body cover");

        std::array<u64, 421> private_counts{};
        std::array<FnvLocal, 421> private_ledgers{};
        u64 singleton_bodies = 0;
        for (std::size_t ordinal = 0; ordinal < body_count; ++ordinal) {
            const int owner = singleton_owner[ordinal];
            if (owner < 0) continue;
            ++singleton_bodies;
            ++private_counts[owner];
            private_ledgers[owner].add(bodies[ordinal]);
            auto node_it = nodes.find(owner);
            if (node_it != nodes.end())
                require((bodies[ordinal] & node_it->second.witness) == 0,
                        "witness misses private body");
        }
        require(singleton_bodies == 3512 && private_counts[318] == 0,
                "private-body census changed");
        for (const auto& [index, node] : nodes)
            require(private_counts[index] == node.expected_obligations &&
                        private_ledgers[index].state ==
                            node.expected_obligation_fnv,
                    "node private-body ledger changed");

        if (argc == 7) write_literal_witness_ledger(argv[6], nodes);

        FnvLocal node_summary;
        std::cout << "SINGLETON_SUCCESS_LITERAL_AUDIT_V1\n"
                  << "NODES " << nodes.size() << " ROWS "
                  << total_literal_rows << " OLD_TESTS " << total_old_tests
                  << " REPLACEMENT_TESTS " << total_replacement_tests
                  << " EQUALITIES " << total_equalities
                  << " LITERAL_FNV " << std::hex
                  << global_literal_ledger.state << std::dec
                  << '\n'
                  << "BODY_SCAN " << body_count << " THREADS "
                  << thread_count << " UNCOVERED " << uncovered.load()
                  << " SINGLETON_BODIES " << singleton_bodies << '\n';
        for (const auto& [index, node] : nodes) {
            node_summary.add(index);
            node_summary.add(node.expected_rows);
            node_summary.add(node.old_mask);
            node_summary.add(node.witness);
            node_summary.add(node.ideal_ledger.state);
            node_summary.add(private_counts[index]);
            node_summary.add(private_ledgers[index].state);
            node_summary.add(node.literal_ledger.state);
            std::cout << "NODE J " << index << " ROWS "
                      << node.literal_rows << " IDEAL_FNV " << std::hex
                      << node.ideal_ledger.state << " OLD " << std::setw(8)
                      << std::setfill('0') << node.old_mask << " WITNESS "
                      << std::setw(8) << node.witness << std::dec
                      << std::setfill(' ') << " OBLIGATIONS "
                      << private_counts[index] << " OBLIGATION_FNV "
                      << std::hex << private_ledgers[index].state
                      << " LITERAL_FNV " << node.literal_ledger.state
                      << std::dec << " WEAKEST " << node.minimum_pair.first
                      << ',' << node.minimum_pair.second << " GAP "
                      << decimal(node.literal_gap_num) << '/'
                      << decimal(node.literal_gap_den) << '\n';
        }
        std::cout << "NODE_SUMMARY_FNV " << std::hex << node_summary.state
                  << std::dec << '\n'
                  << "BODY_COVER_LOGIC ORIGINAL_COVER_PLUS_ALL_PRIVATE_BODIES_HIT\n"
                  << "VERDICT PASS DETACHED_LITERAL_AND_BODY_RESPONSE_AUDIT\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SINGLETON_SUCCESS_LITERAL_AUDIT_ERROR "
                  << error.what() << '\n';
        return 1;
    }
}

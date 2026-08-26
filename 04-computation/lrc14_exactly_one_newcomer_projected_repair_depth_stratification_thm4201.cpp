#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#endif
#define main thm4188_primary_main
#include "lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.cpp"
#undef main
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#include <unordered_set>

namespace {

u32 probe_mask(std::initializer_list<int> values) {
    u32 result = 0;
    for (int value : values) {
        const auto found = std::find(POOL.begin(), POOL.end(), value);
        require(found != POOL.end(), "probe label absent from pool");
        result |= u32{1} << static_cast<int>(found - POOL.begin());
    }
    return result;
}

std::vector<u32> probe_minimal_projected_edges(u32 anchor,
                                                u32 allowed,
                                                const EdgeLayer& layer) {
    std::vector<u32> projected;
    projected.reserve(layer.edges.size());
    for (u32 edge : layer.edges) {
        if ((edge & anchor) != 0) continue;
        const u32 image = edge & allowed;
        if (image == 0) return {0};
        projected.push_back(image);
    }
    std::sort(projected.begin(), projected.end(), [](u32 left, u32 right) {
        const int left_size = std::popcount(left);
        const int right_size = std::popcount(right);
        if (left_size != right_size) return left_size < right_size;
        return left < right;
    });
    projected.erase(std::unique(projected.begin(), projected.end()), projected.end());
    std::vector<u32> minimal;
    for (u32 edge : projected) {
        bool redundant = false;
        for (u32 old : minimal) {
            if ((old & edge) == old) {
                redundant = true;
                break;
            }
        }
        if (!redundant) minimal.push_back(edge);
    }
    return minimal;
}

bool probe_cover_search(const std::vector<u32>& edges,
                        u32 chosen,
                        int remaining,
                        std::unordered_set<u32>& dead,
                        u64& nodes,
                        u32& witness) {
    ++nodes;
    u32 uncovered = 0;
    bool found_uncovered = false;
    for (u32 edge : edges) {
        if ((edge & chosen) == 0) {
            if (edge == 0) return false;
            uncovered = edge;
            found_uncovered = true;
            break;
        }
    }
    if (!found_uncovered) {
        witness = chosen;
        return true;
    }
    if (remaining == 0) return false;
    if (dead.find(chosen) != dead.end()) return false;
    u32 options = uncovered;
    while (options != 0) {
        const u32 bit = options & (~options + 1u);
        if (probe_cover_search(edges, chosen | bit, remaining - 1,
                               dead, nodes, witness)) {
            return true;
        }
        options ^= bit;
    }
    dead.insert(chosen);
    return false;
}

i128 probe_body_delta(const std::vector<Cell>& cells, int q, u32 body) {
    i64 previous = safe_prefix(0, q);
    i128 mass = 0;
    for (const Cell& cell : cells) {
        const i64 current = safe_prefix(cell.right, q);
        const i64 contribution = current - previous;
        if ((cell.failed_pool & body) == 0) mass += contribution;
        previous = current;
    }
    return static_cast<i128>(9) * mass - static_cast<i128>(8) * q * COMMON;
}

}  // namespace

int main() {
    try {
        const std::vector<Cell> cells = build_pool_cells();
        const std::array<u32, 3> anchors = {
            probe_mask({40, 143, 252}),
            probe_mask({80, 143, 252}),
            probe_mask({143, 240, 252})};
        const u32 forbidden_originals = probe_mask({120, 126});
        const u32 universe = (u32{1} << 30) - 1;
        Fnv1a64 ledger;
        u64 d5_rows = 0;
        u64 d5_cover_rows = 0;
        u64 d6_rows = 0;
        std::cout << "PROBE exactly_one_original_all_q_native_repairs\n";
        for (int q : EXPECTED_Q7_FAILURES) {
            const AtomMass mass = build_atom_mass(cells, q, 7);
            const EdgeLayer layer5 = build_layer(mass, q, 5);
            const EdgeLayer layer6 = build_layer(mass, q, 6);
            for (u32 anchor : anchors) {
                const u32 allowed = universe & ~anchor & ~forbidden_originals;
                bool closed = false;
                for (const EdgeLayer* layer : {&layer5, &layer6}) {
                    const std::vector<u32> edges =
                        probe_minimal_projected_edges(anchor, allowed, *layer);
                    std::unordered_set<u32> dead;
                    u64 nodes = 0;
                    u32 witness = 0;
                    const bool has_cover = probe_cover_search(
                        edges, 0, 7, dead, nodes, witness);
                    if (layer->arity == 5) {
                        ++d5_rows;
                        if (has_cover) ++d5_cover_rows;
                    } else {
                        ++d6_rows;
                    }
                    ledger.add_u64_le(static_cast<u64>(q));
                    ledger.add_u64_le(anchor);
                    ledger.add_u64_le(static_cast<u64>(layer->arity));
                    ledger.add_u64_le(layer->edges.size());
                    ledger.add_u64_le(edges.size());
                    ledger.add_u64_le(static_cast<u64>(has_cover));
                    ledger.add_u64_le(nodes);
                    ledger.add_u64_le(witness);
                    std::cout << "ROW q " << q
                              << " anchor " << labels(anchor)
                              << " D" << layer->arity
                              << " raw " << layer->edges.size()
                              << " projected_minimal " << edges.size()
                              << " has_cover_le7 " << has_cover
                              << " nodes " << nodes;
                    if (has_cover) {
                        std::unordered_set<u32> dead6;
                        u64 nodes6 = 0;
                        u32 witness6 = 0;
                        const bool has_cover6 = probe_cover_search(
                            edges, 0, 6, dead6, nodes6, witness6);
                        std::cout << " witness " << labels(witness)
                                  << " has_cover_le6 " << has_cover6
                                  << " nodes_le6 " << nodes6
                                  << " body_delta "
                                  << decimal(probe_body_delta(cells, q, anchor | witness));
                    } else if (layer->arity == 6) {
                        bool found_deeper = false;
                        for (int deeper_budget = 8;
                             deeper_budget <= 12 && !found_deeper;
                             ++deeper_budget) {
                            std::unordered_set<u32> deeper_dead;
                            u64 deeper_nodes = 0;
                            u32 deeper_witness = 0;
                            found_deeper = probe_cover_search(
                                edges, 0, deeper_budget, deeper_dead,
                                deeper_nodes, deeper_witness);
                            std::cout << " has_cover_le" << deeper_budget << ' '
                                      << found_deeper
                                      << " nodes_le" << deeper_budget << ' '
                                      << deeper_nodes;
                            if (found_deeper) {
                                std::cout << " witness" << deeper_budget << ' '
                                          << labels(deeper_witness);
                            }
                        }
                    }
                    std::cout << '\n';
                    if (!has_cover) {
                        closed = true;
                        break;
                    }
                }
                require(closed, "one-original anchor not closed by native E5/E6");
            }
        }
        require(d5_rows == 69 && d5_cover_rows == 4 && d6_rows == 4,
                "primary native staircase split changed");
        std::cout << "SUMMARY D5_ROWS " << d5_rows
                  << " D5_ZERO_COVER_ROWS " << (d5_rows - d5_cover_rows)
                  << " D5_COVER_ROWS " << d5_cover_rows
                  << " D6_ROWS " << d6_rows
                  << " D6_ZERO_COVER_ROWS " << d6_rows << '\n';
        std::cout << "PRIMARY_SEMANTIC_LEDGER_FNV1A64_LE "
                  << hex64(ledger.value()) << '\n';
        std::cout << "VERDICT PASS\n";
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}

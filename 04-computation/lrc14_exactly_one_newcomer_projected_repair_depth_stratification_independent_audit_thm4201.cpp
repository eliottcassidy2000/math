#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#endif
#define main thm4188_joint_main
#include "lrc14_all_newcomer_zero_original_joint_wall_audit_thm4188.cpp"
#undef main
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace {

struct ExhaustiveCoverCount {
    u64 candidates = 0;
    u64 covers = 0;
    u32 first = 0;
    std::vector<u32> witnesses;
    std::size_t projected_unique_edges = 0;
};

std::vector<u32> independent_projected_unique_edges(u32 anchor,
                                                     u32 allowed,
                                                     const Layer& layer) {
    std::vector<u32> edges;
    edges.reserve(layer.edges.size());
    for (u32 edge : layer.edges) {
        if ((edge & anchor) == 0) edges.push_back(edge & allowed);
    }
    std::sort(edges.begin(), edges.end(), [](u32 left, u32 right) {
        const int left_size = std::popcount(left);
        const int right_size = std::popcount(right);
        if (left_size != right_size) return left_size < right_size;
        return left < right;
    });
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

ExhaustiveCoverCount exhaustive_projected_seven_covers(u32 anchor,
                                                        u32 allowed,
                                                        const Layer& layer) {
    const std::vector<u32> edges =
        independent_projected_unique_edges(anchor, allowed, layer);

    std::vector<int> vertices;
    for (int vertex = 0; vertex < 30; ++vertex) {
        if ((allowed & (u32{1} << vertex)) != 0) vertices.push_back(vertex);
    }
    require(vertices.size() == 25, "exactly-one optional ground changed");

    ExhaustiveCoverCount result;
    result.projected_unique_edges = edges.size();
    for_each_k_subset(25, 7, [&](u32 local) {
        ++result.candidates;
        u32 candidate = 0;
        u32 bits = local;
        while (bits != 0) {
            const u32 bit = bits & (~bits + 1u);
            candidate |= u32{1} << vertices[std::countr_zero(bit)];
            bits ^= bit;
        }
        bool cover = true;
        for (u32 edge : edges) {
            if ((edge & candidate) == 0) {
                cover = false;
                break;
            }
        }
        if (cover) {
            ++result.covers;
            if (result.first == 0) result.first = candidate;
            result.witnesses.push_back(candidate);
        }
    });
    require(result.candidates == 480700,
            "exactly-one seven-body candidate count changed");
    return result;
}

struct IndependentTau {
    int tau = 0;
    u32 witness = 0;
    u64 nodes = 0;
};

IndependentTau independent_projected_tau(u32 anchor,
                                         u32 allowed,
                                         const Layer& layer,
                                         int maximum_budget) {
    const std::vector<u32> edges =
        independent_projected_unique_edges(anchor, allowed, layer);
    require(std::find(edges.begin(), edges.end(), u32{0}) == edges.end(),
            "empty projected edge has infinite transversal number");
    IndependentTau result;
    for (int budget = 0; budget <= maximum_budget; ++budget) {
        std::unordered_set<u32> dead;
        u64 nodes = 0;
        const bool covered = cover_search(edges, 0, budget, dead, nodes);
        result.nodes += nodes;
        if (covered) {
            result.tau = budget;
            // The inherited solver is decision-only.  The exact value, not a
            // canonical witness, is the independent sharpening sidecar.
            return result;
        }
    }
    return result;
}

}  // namespace

int main() {
    try {
        const FixedGeometry fixed = build_fixed_geometry();
        const Support support = build_support(fixed);
        const std::array<u32, 3> anchors = {
            label_mask({40, 143, 252}),
            label_mask({80, 143, 252}),
            label_mask({143, 240, 252})};
        const u32 forbidden_originals = label_mask({120, 126});
        const u32 universe = (u32{1} << 30) - 1;
        Ledger ledger;
        u64 rows_d5 = 0;
        u64 rows_d6 = 0;
        u64 d5_cover_rows = 0;
        std::cout << "AUDIT exactly_one_original_all_q_joint_wall_exhaustive\n";
        for (int q : EXPECTED_Q7) {
            const JointMass mass = build_joint_mass(fixed, support, q);
            const Layer layer5 = build_layer(mass, support, 5);
            const Layer layer6 = build_layer(mass, support, 6);
            for (u32 anchor : anchors) {
                const u32 allowed = universe & ~anchor & ~forbidden_originals;
                const ExhaustiveCoverCount d5 =
                    exhaustive_projected_seven_covers(anchor, allowed, layer5);
                ++rows_d5;
                ledger.add(static_cast<u64>(q));
                ledger.add(anchor);
                ledger.add(5);
                ledger.add(layer5.edges.size());
                ledger.add(d5.projected_unique_edges);
                ledger.add(d5.covers);
                ledger.add(d5.first);
                for (u32 witness : d5.witnesses) ledger.add(witness);
                std::cout << "ROW q " << q
                          << " anchor " << labels(anchor)
                          << " D5 raw " << layer5.edges.size()
                          << " projected_unique " << d5.projected_unique_edges
                          << " exact_seven_covers " << d5.covers;
                if (d5.covers != 0) {
                    ++d5_cover_rows;
                    std::cout << " witnesses";
                    for (u32 witness : d5.witnesses) {
                        const i64 direct_mass = direct_joint_mass(
                            fixed, q, universe & ~(anchor | witness));
                        const i128 direct_delta = static_cast<i128>(63) * direct_mass -
                                                  static_cast<i128>(4) * mass.denominator;
                        std::cout << " [" << labels(witness) << ';'
                                  << decimal(direct_delta) << ']';
                    }
                    const ExhaustiveCoverCount d6 =
                        exhaustive_projected_seven_covers(anchor, allowed, layer6);
                    const IndependentTau tau6 =
                        independent_projected_tau(anchor, allowed, layer6, 10);
                    ++rows_d6;
                    ledger.add(static_cast<u64>(q));
                    ledger.add(anchor);
                    ledger.add(6);
                    ledger.add(layer6.edges.size());
                    ledger.add(d6.projected_unique_edges);
                    ledger.add(d6.covers);
                    ledger.add(d6.first);
                    ledger.add(static_cast<u64>(tau6.tau));
                    ledger.add(tau6.nodes);
                    std::cout << " D6_raw " << layer6.edges.size()
                              << " D6_projected_unique " << d6.projected_unique_edges
                              << " D6_exact_seven_covers " << d6.covers
                              << " D6_tau " << tau6.tau
                              << " D6_tau_search_nodes " << tau6.nodes;
                    require(d6.covers == 0,
                            "native D6 retained an exactly-one seven-cover");
                    const int expected_tau = (q == 105) ? 9 : 10;
                    require(tau6.tau == expected_tau,
                            "native D6 projected transversal number changed");
                }
                std::cout << '\n';
            }
        }
        require(rows_d5 == 69, "native D5 row count changed");
        require(d5_cover_rows == 4 && rows_d6 == 4,
                "native D5/D6 staircase split changed");
        std::cout << "SUMMARY D5_ROWS " << rows_d5
                  << " D5_ZERO_COVER_ROWS " << (rows_d5 - d5_cover_rows)
                  << " D5_COVER_ROWS " << d5_cover_rows
                  << " D6_ROWS " << rows_d6
                  << " D6_ZERO_COVER_ROWS " << rows_d6 << '\n';
        std::cout << "SEMANTIC_LEDGER_FNV1A64_LE " << hex64(ledger.value()) << '\n';
        std::cout << "VERDICT PASS\n";
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}

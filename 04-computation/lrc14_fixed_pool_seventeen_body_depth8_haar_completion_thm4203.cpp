// Primary exact pool-wall audit for THM-4203.
// This file is maintained-path ready: after moving it into 04-computation,
// the first include branch is selected without changing the source bytes.

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#endif
#define main thm4188_primary_main
#if __has_include("lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.cpp")
#include "lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.cpp"
#else
#include "../04-computation/lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.cpp"
#endif
#undef main
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

#include <unordered_set>

namespace {

constexpr u64 EXPECTED_CANDIDATE_LEDGER = UINT64_C(0x68ff828dae998c28);

u32 candidate_mask(std::initializer_list<int> values) {
    u32 result = 0;
    for (int value : values) {
        const auto found = std::find(POOL.begin(), POOL.end(), value);
        require(found != POOL.end(), "candidate label absent from pool");
        result |= u32{1} << static_cast<int>(found - POOL.begin());
    }
    return result;
}

bool candidate_is_cover(u32 body, const std::vector<u32>& edges) {
    return std::all_of(edges.begin(), edges.end(), [&](u32 edge) {
        return (body & edge) != 0;
    });
}

i64 candidate_fixed_mass(const std::vector<Cell>& cells, u32 body) {
    i64 result = 0;
    for (const Cell& cell : cells) {
        if ((cell.failed_pool & body) == 0) result += cell.right - cell.left;
    }
    return result;
}

struct CandidateCoverResult {
    bool found = false;
    u32 witness = 0;
    u64 nodes = 0;
    u64 dead_states = 0;
};

bool candidate_cover_search(const std::vector<u32>& edges,
                            u32 chosen,
                            int remaining,
                            std::unordered_set<u32>& dead,
                            CandidateCoverResult& result) {
    ++result.nodes;
    u32 uncovered = 0;
    for (u32 edge : edges) {
        if ((edge & chosen) == 0) {
            uncovered = edge;
            break;
        }
    }
    if (uncovered == 0) {
        result.found = true;
        result.witness = chosen;
        return true;
    }
    if (remaining == 0 || dead.contains(chosen)) return false;
    while (uncovered != 0) {
        const u32 bit = uncovered & (~uncovered + 1u);
        if (candidate_cover_search(edges, chosen | bit, remaining - 1,
                                   dead, result)) return true;
        uncovered ^= bit;
    }
    dead.insert(chosen);
    return false;
}

CandidateCoverResult candidate_search(const std::vector<u32>& edges,
                                      int budget) {
    std::unordered_set<u32> dead;
    CandidateCoverResult result;
    candidate_cover_search(edges, 0, budget, dead, result);
    result.dead_states = dead.size();
    return result;
}

void candidate_solver_controls() {
    const std::vector<u32> positive = {
        (u32{1} << 0) | (u32{1} << 1),
        (u32{1} << 1) | (u32{1} << 2)};
    const std::vector<u32> negative = {
        (u32{1} << 0) | (u32{1} << 1),
        (u32{1} << 2) | (u32{1} << 3)};
    const CandidateCoverResult yes = candidate_search(positive, 1);
    const CandidateCoverResult no = candidate_search(negative, 1);
    require(yes.found && std::popcount(yes.witness) == 1,
            "positive cover control failed");
    require(!no.found, "negative cover control failed");
}

struct CandidateExhaustiveEleven {
    u64 bodies = 0;
    u64 edge_checks = 0;
    u64 maximum_checks = 0;
    u32 closest_body = 0;
    u32 missed_edge = 0;
    u32 cover = 0;
};

struct CandidateCoverCensus {
    const std::vector<u32>& edges;
    const std::vector<Cell>& cells;
    std::unordered_set<u32> visited;
    std::vector<u32> covers;
    u64 nodes = 0;
    u64 dead_leaves = 0;
    int minimum_size = 31;
    i128 minimum_delta = 0;
    u32 minimum_body = 0;
    bool have_minimum = false;

    CandidateCoverCensus(const std::vector<u32>& edge_bank,
                         const std::vector<Cell>& cell_bank)
        : edges(edge_bank), cells(cell_bank) {}

    i128 body_delta(u32 body) const {
        return static_cast<i128>(63) * candidate_fixed_mass(cells, body) -
               static_cast<i128>(4) * COMMON;
    }

    void search(u32 chosen, int remaining) {
        if (!visited.insert(chosen).second) return;
        ++nodes;
        u32 uncovered = 0;
        for (u32 edge : edges) {
            if ((edge & chosen) == 0) {
                uncovered = edge;
                break;
            }
        }
        if (uncovered == 0) {
            covers.push_back(chosen);
            const int size = std::popcount(chosen);
            minimum_size = std::min(minimum_size, size);
            const i128 value = body_delta(chosen);
            if (!have_minimum || value < minimum_delta) {
                have_minimum = true;
                minimum_delta = value;
                minimum_body = chosen;
            }
            return;
        }
        if (remaining == 0) {
            ++dead_leaves;
            return;
        }
        while (uncovered != 0) {
            const u32 bit = uncovered & (~uncovered + 1u);
            search(chosen | bit, remaining - 1);
            uncovered ^= bit;
        }
    }
};

void candidate_extend_by(u32 body, int next_vertex, int needed,
                         std::vector<u32>& output) {
    if (needed == 0) {
        output.push_back(body);
        return;
    }
    for (int vertex = next_vertex; vertex <= 30 - needed; ++vertex) {
        const u32 bit = u32{1} << vertex;
        if ((body & bit) == 0) {
            candidate_extend_by(body | bit, vertex + 1, needed - 1, output);
        }
    }
}

CandidateExhaustiveEleven candidate_exhaust_all_eleven(
        std::vector<u32> edges) {
    auto mix = [](u64 value) {
        value += UINT64_C(0x9e3779b97f4a7c15);
        value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
        return value ^ (value >> 31);
    };
    std::sort(edges.begin(), edges.end(), [&](u32 left, u32 right) {
        const u64 left_hash = mix(left);
        const u64 right_hash = mix(right);
        return left_hash != right_hash ? left_hash < right_hash : left < right;
    });

    CandidateExhaustiveEleven result;
    u32 body = (u32{1} << 11) - 1;
    const u32 limit = u32{1} << 30;
    while (body < limit) {
        u64 checked = 0;
        u32 missed = 0;
        for (u32 edge : edges) {
            ++checked;
            if ((edge & body) == 0) {
                missed = edge;
                break;
            }
        }
        ++result.bodies;
        result.edge_checks += checked;
        if (checked > result.maximum_checks) {
            result.maximum_checks = checked;
            result.closest_body = body;
            result.missed_edge = missed;
        }
        if (missed == 0) {
            result.cover = body;
            break;
        }
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    return result;
}

}  // namespace

int main() {
    try {
        candidate_solver_controls();
        const std::vector<Cell> cells = build_pool_cells();
        const AtomMass mass50 = build_atom_mass(cells, 50, 8);
        EdgeLayer e7 = build_layer(mass50, 50, 7);
        EdgeLayer e8 = build_layer(mass50, 50, 8);
        std::sort(e7.edges.begin(), e7.edges.end());
        std::sort(e8.edges.begin(), e8.edges.end());
        require(e7.edges.size() == UINT64_C(821737) && e7.equalities == 0,
                "q50 E7 layer changed");
        require(e8.edges.size() == UINT64_C(4088855) && e8.equalities == 0,
                "q50 E8 layer changed");

        const u32 e7_cover = candidate_mask(
            {8,80,85,88,95,143,145,168,193,240,252});
        const u32 e8_miss = candidate_mask(
            {16,84,120,176,190,264,286,290});
        const u32 e8_cover14 = candidate_mask(
            {8,10,16,80,85,88,95,143,145,168,176,193,240,252});
        require(std::popcount(e7_cover) == 11,
                "E7 hostile has wrong cardinality");
        require(candidate_is_cover(e7_cover, e7.edges),
                "explicit E7 eleven-cover changed");
        require((e7_cover & e8_miss) == 0,
                "E8 hostile miss ceased to be disjoint");
        require(std::binary_search(e8.edges.begin(), e8.edges.end(), e8_miss),
                "explicit E8 missed repair changed");
        require(std::popcount(e8_cover14) == 14 &&
                    candidate_is_cover(e8_cover14, e8.edges),
                "explicit E8 fourteen-cover changed");
        const i128 body_delta =
            static_cast<i128>(63) * candidate_fixed_mass(cells, e7_cover) -
            static_cast<i128>(4) * COMMON;
        require(body_delta == static_cast<i128>(139755964331592LL),
                "direct E7-hostile body margin changed");
        const i128 cover14_body_delta =
            static_cast<i128>(63) * candidate_fixed_mass(cells, e8_cover14) -
            static_cast<i128>(4) * COMMON;
        require(cover14_body_delta == static_cast<i128>(77740560077796LL),
                "direct E8-cover body margin changed");

        CandidateCoverCensus cover_census(e8.edges, cells);
        cover_census.visited.reserve(UINT64_C(3000000));
        cover_census.search(0, 17);
        std::sort(cover_census.covers.begin(), cover_census.covers.end());
        const u64 terminal_size14 = static_cast<u64>(std::count_if(
            cover_census.covers.begin(), cover_census.covers.end(),
            [](u32 body) { return std::popcount(body) == 14; }));
        const u64 terminal_size15 = static_cast<u64>(std::count_if(
            cover_census.covers.begin(), cover_census.covers.end(),
            [](u32 body) { return std::popcount(body) == 15; }));
        const u64 terminal_size16 = static_cast<u64>(std::count_if(
            cover_census.covers.begin(), cover_census.covers.end(),
            [](u32 body) { return std::popcount(body) == 16; }));
        const u64 terminal_size17 = static_cast<u64>(std::count_if(
            cover_census.covers.begin(), cover_census.covers.end(),
            [](u32 body) { return std::popcount(body) == 17; }));
        require(cover_census.nodes == UINT64_C(4991151) &&
                    cover_census.dead_leaves == UINT64_C(1428692) &&
                    cover_census.covers.size() == 57410 &&
                    terminal_size14 == 35 && terminal_size15 == 817 &&
                    terminal_size16 == 8296 && terminal_size17 == 48262 &&
                    cover_census.minimum_size == 14,
                "q50 E8 complete cover census changed");
        require(cover_census.minimum_body == candidate_mask(
                    {8,10,16,42,60,80,85,88,95,143,145,168,176,193,240,252,264}) &&
                    cover_census.minimum_delta ==
                        static_cast<i128>(34853530580568LL),
                "q50 E8 direct cover minimum changed");
        require(cover_census.minimum_delta > 0,
                "an E8 cover lost its positive direct Haar margin");

        // The recursive search stops as soon as a branch first becomes a
        // cover.  A terminal cover need not be inclusion-minimal because the
        // branch may already contain a redundant vertex.  Nevertheless every
        // bounded cover contains a terminal search cover: always branch on a
        // vertex it contributes to the first uncovered edge.  Close the full
        // terminal list upward through size seventeen and deduplicate.
        std::vector<u32> all_cover_bodies = cover_census.covers;
        for (u32 body : cover_census.covers) {
            const int size = std::popcount(body);
            for (int target = size + 1; target <= 17; ++target) {
                candidate_extend_by(body, 0, target - size, all_cover_bodies);
            }
        }
        std::sort(all_cover_bodies.begin(), all_cover_bodies.end());
        all_cover_bodies.erase(
            std::unique(all_cover_bodies.begin(), all_cover_bodies.end()),
            all_cover_bodies.end());
        const u64 all_size14 = static_cast<u64>(std::count_if(
            all_cover_bodies.begin(), all_cover_bodies.end(),
            [](u32 body) { return std::popcount(body) == 14; }));
        const u64 all_size15 = static_cast<u64>(std::count_if(
            all_cover_bodies.begin(), all_cover_bodies.end(),
            [](u32 body) { return std::popcount(body) == 15; }));
        const u64 all_size16 = static_cast<u64>(std::count_if(
            all_cover_bodies.begin(), all_cover_bodies.end(),
            [](u32 body) { return std::popcount(body) == 16; }));
        const u64 all_size17 = static_cast<u64>(std::count_if(
            all_cover_bodies.begin(), all_cover_bodies.end(),
            [](u32 body) { return std::popcount(body) == 17; }));
        i128 all_minimum_delta = 0;
        u32 all_minimum_body = 0;
        bool have_all_minimum = false;
        for (u32 body : all_cover_bodies) {
            const i128 value = cover_census.body_delta(body);
            if (!have_all_minimum || value < all_minimum_delta) {
                have_all_minimum = true;
                all_minimum_delta = value;
                all_minimum_body = body;
            }
        }
        require(all_cover_bodies.size() == 66468 && all_size14 == 35 &&
                    all_size15 == 876 && all_size16 == 9320 &&
                    all_size17 == 56237,
                "expanded E8 transversal census changed");
        require(all_minimum_body == candidate_mask(
                    {8,10,16,42,60,80,85,88,95,143,145,168,176,193,240,252,264}) &&
                    all_minimum_delta ==
                        static_cast<i128>(34853530580568LL),
                "expanded E8 direct minimum changed");
        require(all_minimum_delta > 0,
                "an E8 transversal through size seventeen failed directly");
        const CandidateExhaustiveEleven brute =
            candidate_exhaust_all_eleven(e8.edges);
        require(brute.bodies == UINT64_C(54627300) && brute.cover == 0,
                "literal C(30,11) audit changed");

        Fnv1a64 terminal_cover_fingerprint;
        for (u32 body : cover_census.covers) {
            terminal_cover_fingerprint.add_u64_le(body);
            terminal_cover_fingerprint.add_u64_le(
                static_cast<u64>(cover_census.body_delta(body)));
        }
        Fnv1a64 all_cover_fingerprint;
        Fnv1a64 exact17_cover_fingerprint;
        for (u32 body : all_cover_bodies) {
            all_cover_fingerprint.add_u64_le(body);
            all_cover_fingerprint.add_u64_le(
                static_cast<u64>(cover_census.body_delta(body)));
            if (std::popcount(body) == 17) {
                exact17_cover_fingerprint.add_u64_le(body);
                exact17_cover_fingerprint.add_u64_le(
                    static_cast<u64>(cover_census.body_delta(body)));
            }
        }

        Fnv1a64 ledger;
        for (u64 value : {
                static_cast<u64>(e7.edges.size()), e7.equalities,
                static_cast<u64>(e8.edges.size()), e8.equalities,
                static_cast<u64>(e7_cover), static_cast<u64>(e8_miss),
                static_cast<u64>(body_delta), cover_census.nodes,
                cover_census.dead_leaves,
                static_cast<u64>(cover_census.covers.size()),
                static_cast<u64>(cover_census.minimum_size),
                static_cast<u64>(cover_census.minimum_body),
                static_cast<u64>(cover_census.minimum_delta),
                terminal_size14, terminal_size15, terminal_size16,
                terminal_size17, terminal_cover_fingerprint.value(),
                static_cast<u64>(all_cover_bodies.size()), all_size14,
                all_size15, all_size16, all_size17,
                static_cast<u64>(all_minimum_body),
                static_cast<u64>(all_minimum_delta),
                exact17_cover_fingerprint.value(), all_cover_fingerprint.value(),
                brute.bodies, brute.edge_checks,
                brute.maximum_checks, static_cast<u64>(brute.closest_body),
                static_cast<u64>(brute.missed_edge),
                static_cast<u64>(e8_cover14),
                static_cast<u64>(cover14_body_delta)}) {
            ledger.add_u64_le(value);
        }

        std::cout << "AUDIT complete_full_pool_depth8_primary\n";
        std::cout << "Q50 E7 " << e7.edges.size()
                  << " EQUALITIES " << e7.equalities
                  << " COVER11 " << labels(e7_cover)
                  << " DIRECT_63M_MINUS_4D " << decimal(body_delta) << '\n';
        std::cout << "Q50 E8 " << e8.edges.size()
                  << " EQUALITIES " << e8.equalities
                  << " HOSTILE_MISSED_EDGE " << labels(e8_miss)
                  << " EDGE_DELTA " << decimal(delta(e8_miss, mass50, 50))
                  << '\n';
        std::cout << "E8 TERMINAL_COVER_CENSUS_BUDGET17 NODES "
                  << cover_census.nodes
                  << " DEAD_LEAVES " << cover_census.dead_leaves
                  << " COVERS " << cover_census.covers.size()
                  << " SIZE14 " << terminal_size14
                  << " SIZE15 " << terminal_size15
                  << " SIZE16 " << terminal_size16
                  << " SIZE17 " << terminal_size17
                  << " MINIMUM_SIZE " << cover_census.minimum_size
                  << " MIN_DIRECT_63M_MINUS_4D "
                  << decimal(cover_census.minimum_delta)
                  << " MIN_BODY " << labels(cover_census.minimum_body)
                  << " TERMINAL_COVER_MARGIN_FNV1A64_LE "
                  << hex64(terminal_cover_fingerprint.value()) << '\n';
        std::cout << "E8 ALL_COVERS_THROUGH17 "
                  << all_cover_bodies.size()
                  << " SIZE14 " << all_size14
                  << " SIZE15 " << all_size15
                  << " SIZE16 " << all_size16
                  << " SIZE17 " << all_size17
                  << " MIN_DIRECT_63M_MINUS_4D "
                  << decimal(all_minimum_delta)
                  << " MIN_BODY " << labels(all_minimum_body)
                  << " ALL_COVER_MARGIN_FNV1A64_LE "
                  << hex64(all_cover_fingerprint.value())
                  << " EXACT17_COVER_MARGIN_FNV1A64_LE "
                  << hex64(exact17_cover_fingerprint.value()) << '\n';
        std::cout << "E8 COVER14 " << labels(e8_cover14)
                  << " DIRECT_63M_MINUS_4D "
                  << decimal(cover14_body_delta) << '\n';
        std::cout << "BRUTE_C30_11 BODIES " << brute.bodies
                  << " COVER NONE EDGE_CHECKS " << brute.edge_checks
                  << " MAX_CHECKS " << brute.maximum_checks
                  << " CLOSEST_BODY " << labels(brute.closest_body)
                  << " MISSED_EDGE " << labels(brute.missed_edge) << '\n';
        std::cout << "BODY_COUNTS C30_11 54627300 C30_12 86493225"
                     " C30_13 119759850 C30_14 145422675"
                     " C30_15 155117520 C30_16 145422675"
                     " C30_17 119759850\n";
        std::cout << "JOINT_11_BODY_COUNT_WITH_THM4191"
                     " 54627300+30045015=84672315\n";
        std::cout << "SEMANTIC_LEDGER_FNV1A64_LE "
                  << hex64(ledger.value()) << '\n';
        if (EXPECTED_CANDIDATE_LEDGER != 0) {
            require(ledger.value() == EXPECTED_CANDIDATE_LEDGER,
                    "primary candidate semantic ledger changed");
        }
        std::cout << "VERDICT PASS ALL_POOL_BODIES_LE17_SAFE\n";
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}

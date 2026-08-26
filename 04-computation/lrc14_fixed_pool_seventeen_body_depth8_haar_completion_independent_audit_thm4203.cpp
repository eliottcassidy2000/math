// Independent explicit joint-wall audit for THM-4203.
// This file is maintained-path ready: after moving it into 04-computation,
// the first include branch is selected without changing the source bytes.

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#endif
#define main thm4188_independent_main
#if __has_include("lrc14_all_newcomer_zero_original_joint_wall_audit_thm4188.cpp")
#include "lrc14_all_newcomer_zero_original_joint_wall_audit_thm4188.cpp"
#else
#include "../04-computation/lrc14_all_newcomer_zero_original_joint_wall_audit_thm4188.cpp"
#endif
#undef main
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace {

constexpr u64 EXPECTED_CANDIDATE_LEDGER = UINT64_C(0xd9e1bd77881a64a7);

Support candidate_support_nine(const FixedGeometry& fixed) {
    Support support;
    for (const FixedCell& cell : fixed.cells) {
        if (std::popcount(cell.failed) <= 9) support.masks.push_back(cell.failed);
    }
    std::sort(support.masks.begin(), support.masks.end());
    support.masks.erase(std::unique(support.masks.begin(), support.masks.end()),
                        support.masks.end());
    support.position.reserve(2 * support.masks.size());
    for (u32 index = 0; index < support.masks.size(); ++index) {
        support.position.emplace(support.masks[index], index);
    }
    return support;
}

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

i64 candidate_fixed_mass(const FixedGeometry& fixed, u32 body) {
    i64 result = 0;
    for (const FixedCell& cell : fixed.cells) {
        if ((cell.failed & body) == 0) result += cell.right - cell.left;
    }
    return result;
}

struct CandidateCoverCensus {
    const std::vector<u32>& edges;
    const FixedGeometry& fixed;
    std::unordered_set<u32> visited;
    std::vector<u32> covers;
    u64 nodes = 0;
    u64 dead_leaves = 0;
    int minimum_size = 31;
    i128 minimum_delta = 0;
    u32 minimum_body = 0;
    bool have_minimum = false;

    CandidateCoverCensus(const std::vector<u32>& edge_bank,
                         const FixedGeometry& fixed_geometry)
        : edges(edge_bank), fixed(fixed_geometry) {}

    i128 body_delta(u32 body) const {
        return static_cast<i128>(63) * candidate_fixed_mass(fixed, body) -
               static_cast<i128>(4) * POOL_DENOMINATOR;
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

}  // namespace

int main() {
    try {
        cover_solver_controls();
        const FixedGeometry fixed = build_fixed_geometry();
        const Support support = candidate_support_nine(fixed);
        const JointMass mass50 = build_joint_mass(fixed, support, 50);
        Layer e7 = build_layer(mass50, support, 7);
        Layer e8 = build_layer(mass50, support, 8);
        std::sort(e7.edges.begin(), e7.edges.end());
        std::sort(e8.edges.begin(), e8.edges.end());
        require(e7.edges.size() == UINT64_C(821737) && e7.equalities == 0,
                "joint q50 E7 layer changed");
        require(e8.edges.size() == UINT64_C(4088855) && e8.equalities == 0,
                "joint q50 E8 layer changed");

        const u32 e7_cover = candidate_mask(
            {8,80,85,88,95,143,145,168,193,240,252});
        const u32 e8_miss = candidate_mask(
            {16,84,120,176,190,264,286,290});
        const u32 e8_cover14 = candidate_mask(
            {8,10,16,80,85,88,95,143,145,168,176,193,240,252});
        require(candidate_is_cover(e7_cover, e7.edges),
                "joint E7 hostile cover changed");
        require((e7_cover & e8_miss) == 0 &&
                    std::binary_search(e8.edges.begin(), e8.edges.end(), e8_miss),
                "joint E8 hostile repair changed");
        require(std::popcount(e8_cover14) == 14 &&
                    candidate_is_cover(e8_cover14, e8.edges),
                "joint E8 fourteen-cover changed");
        const i128 body_delta =
            static_cast<i128>(63) * candidate_fixed_mass(fixed, e7_cover) -
            static_cast<i128>(4) * POOL_DENOMINATOR;
        require(body_delta == static_cast<i128>(139755964331592LL),
                "joint direct body margin changed");
        const i128 cover14_body_delta =
            static_cast<i128>(63) * candidate_fixed_mass(fixed, e8_cover14) -
            static_cast<i128>(4) * POOL_DENOMINATOR;
        require(cover14_body_delta == static_cast<i128>(77740560077796LL),
                "joint direct E8-cover body margin changed");
        require(direct_joint_mass(fixed, 50, e8_miss) ==
                    subset_mass(e8_miss, mass50, support),
                "literal joint-wall repair mass disagreed");

        CandidateCoverCensus cover_census(e8.edges, fixed);
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
                "joint E8 complete cover census changed");
        require(cover_census.minimum_body == candidate_mask(
                    {8,10,16,42,60,80,85,88,95,143,145,168,176,193,240,252,264}) &&
                    cover_census.minimum_delta ==
                        static_cast<i128>(34853530580568LL),
                "joint E8 cover minimum changed");
        require(cover_census.minimum_delta > 0,
                "an E8 cover lost its direct positive Haar margin");

        // A terminal search cover need not be inclusion-minimal: a branch may
        // choose a redundant vertex before it first covers every edge.  Every
        // bounded cover still contains some terminal search cover by following
        // only its vertices through the first-uncovered-edge recursion.  Close
        // the terminal list upward through size seventeen and deduplicate.
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
                "expanded joint E8 transversal census changed");
        require(all_minimum_body == candidate_mask(
                    {8,10,16,42,60,80,85,88,95,143,145,168,176,193,240,252,264}) &&
                    all_minimum_delta ==
                        static_cast<i128>(34853530580568LL),
                "expanded joint E8 direct minimum changed");
        require(all_minimum_delta > 0,
                "a joint E8 transversal through size seventeen failed directly");

        std::vector<int> outside_cover14;
        for (int vertex = 0; vertex < 30; ++vertex) {
            if ((e8_cover14 & (u32{1} << vertex)) == 0) {
                outside_cover14.push_back(vertex);
            }
        }
        u64 e9_candidates = 0;
        u64 e9_lawful = 0;
        u64 e9_equalities = 0;
        u32 e9_first = 0;
        u32 e9_best = 0;
        i128 e9_first_delta = 0;
        i128 e9_best_delta = std::numeric_limits<i64>::min();
        auto e9_callback = [&](u32 deletion) {
            ++e9_candidates;
            const i128 value = threshold_delta(deletion, mass50, support);
            if (value == 0) ++e9_equalities;
            if (value < 0) return;
            ++e9_lawful;
            if (e9_first == 0) {
                e9_first = deletion;
                e9_first_delta = value;
            }
            if (e9_best == 0 || value > e9_best_delta) {
                e9_best = deletion;
                e9_best_delta = value;
            }
        };
        choose_from_vertices(outside_cover14, 9, 0, 0, e9_callback);
        require(e9_candidates == UINT64_C(11440) &&
                    e9_lawful == UINT64_C(108) && e9_equalities == 0,
                "joint E9 repair census changed");
        const u32 e9_root_control =
            candidate_mask({15,30,42,60,63,84,120,132,264});
        require(threshold_delta(e9_root_control, mass50, support) ==
                    static_cast<i128>(4603247032482LL),
                "joint E9 root control changed");
        require(e9_first == candidate_mask({15,30,40,60,63,120,170,190,290}) &&
                    e9_first_delta == static_cast<i128>(285809592222LL),
                "joint E9 lexicographic first repair changed");
        require(e9_best == candidate_mask({40,42,63,120,132,170,190,264,290}) &&
                    e9_best_delta == static_cast<i128>(15738299433612LL),
                "joint E9 best repair changed");
        require(direct_joint_mass(fixed, 50, e9_first) ==
                    subset_mass(e9_first, mass50, support) &&
                    direct_joint_mass(fixed, 50, e9_best) ==
                    subset_mass(e9_best, mass50, support),
                "literal joint-wall E9 repair mass disagreed");

        Ledger terminal_cover_fingerprint;
        for (u32 body : cover_census.covers) {
            terminal_cover_fingerprint.add(body);
            terminal_cover_fingerprint.add(
                static_cast<u64>(cover_census.body_delta(body)));
        }
        Ledger all_cover_fingerprint;
        Ledger exact17_cover_fingerprint;
        for (u32 body : all_cover_bodies) {
            all_cover_fingerprint.add(body);
            all_cover_fingerprint.add(
                static_cast<u64>(cover_census.body_delta(body)));
            if (std::popcount(body) == 17) {
                exact17_cover_fingerprint.add(body);
                exact17_cover_fingerprint.add(
                    static_cast<u64>(cover_census.body_delta(body)));
            }
        }

        Ledger ledger;
        for (u64 value : {
                mass50.joint_cells, static_cast<u64>(support.masks.size()),
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
                static_cast<u64>(e8_cover14),
                static_cast<u64>(cover14_body_delta), e9_candidates, e9_lawful,
                e9_equalities, static_cast<u64>(e9_first),
                static_cast<u64>(e9_first_delta), static_cast<u64>(e9_best),
                static_cast<u64>(e9_best_delta)}) {
            ledger.add(value);
        }

        std::cout << "AUDIT fixed_pool_depth8_independent_joint_wall\n";
        std::cout << "JOINT_CELLS " << mass50.joint_cells
                  << " SUPPORT " << support.masks.size() << '\n';
        std::cout << "Q50 E7 " << e7.edges.size()
                  << " EQUALITIES " << e7.equalities
                  << " COVER11 " << labels(e7_cover)
                  << " DIRECT_63M_MINUS_4D " << decimal(body_delta) << '\n';
        std::cout << "Q50 E8 " << e8.edges.size()
                  << " EQUALITIES " << e8.equalities
                  << " HOSTILE_MISSED_EDGE " << labels(e8_miss)
                  << " EDGE_DELTA "
                  << decimal(threshold_delta(e8_miss, mass50, support)) << '\n';
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
        std::cout << "E9_DISJOINT_COVER14 CANDIDATES " << e9_candidates
                  << " LAWFUL " << e9_lawful
                  << " EQUALITIES " << e9_equalities
                  << " FIRST " << labels(e9_first)
                  << " JOINT_DELTA " << decimal(e9_first_delta)
                  << " BEST " << labels(e9_best)
                  << " JOINT_DELTA " << decimal(e9_best_delta) << '\n';
        std::cout << "BODY_COUNTS C30_11 54627300 C30_12 86493225"
                     " C30_13 119759850 C30_14 145422675"
                     " C30_15 155117520 C30_16 145422675"
                     " C30_17 119759850\n";
        std::cout << "SEMANTIC_LEDGER_FNV1A64_LE "
                  << hex64(ledger.value()) << '\n';
        if (EXPECTED_CANDIDATE_LEDGER != 0) {
            require(ledger.value() == EXPECTED_CANDIDATE_LEDGER,
                    "independent candidate semantic ledger changed");
        }
        std::cout << "VERDICT PASS ALL_POOL_BODIES_LE17_SAFE"
                     " TAU_E8_Q50_EQ_14 E9_REPAIRS_108\n";
    } catch (const std::exception& error) {
        std::cerr << "ERROR " << error.what() << '\n';
        return 1;
    }
}

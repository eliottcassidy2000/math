// Complete labelled response quotient for the eighteen failed bodies at
// (256,667), with an exact quotient cover and obligation-packing lower bound.

#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

namespace {

constexpr std::size_t OBLIGATIONS_667 = 18;
constexpr u32 FULL_RESPONSE_667 = (u32{1} << OBLIGATIONS_667) - 1;
constexpr std::array<u32, OBLIGATIONS_667> FAILURES_667 = {
    UINT32_C(0x0514a488), UINT32_C(0x0516a401),
    UINT32_C(0x0516a408), UINT32_C(0x07142409),
    UINT32_C(0x07142488), UINT32_C(0x07162401),
    UINT32_C(0x07162408), UINT32_C(0x07162480),
    UINT32_C(0x071e2008), UINT32_C(0x071e2400),
    UINT32_C(0x0d14a401), UINT32_C(0x0d14a408),
    UINT32_C(0x0d542408), UINT32_C(0x0f142401),
    UINT32_C(0x0f142408), UINT32_C(0x0f1c2008),
    UINT32_C(0x151aa400), UINT32_C(0x27162400)};

}  // namespace

int main() {
    try {
        std::cout.setf(std::ios::unitbuf);
        init_choose8_local();
        FnvLocal failure_ledger;
        u32 union_mask = 0;
        for (u32 body : FAILURES_667) {
            require(std::popcount(body) == 9, "failed body arity changed");
            failure_ledger.add(body);
            union_mask |= body;
        }
        require(failure_ledger.state == UINT64_C(0xce3d9a9b4b429af8) &&
                    union_mask == UINT32_C(0x3f5ea489) &&
                    std::popcount(union_mask) == 17,
                "endpoint-667 failure universe changed");

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        const ActiveUniverse active = build_active_universe(cells, 256, 667);
        std::vector<u32> response(EXPECTED_REPAIRS, 0);
        u64 incidences = 0;
        for (std::size_t obligation = 0; obligation < FAILURES_667.size();
             ++obligation) {
            enumerate_disjoint_repairs(
                FAILURES_667[obligation], [&](u32, u64 rank) {
                    if (!active.active[rank]) return;
                    response[rank] |= u32{1} << obligation;
                    ++incidences;
                });
        }

        constexpr std::size_t PATTERNS = std::size_t{1} << OBLIGATIONS_667;
        std::vector<u64> multiplicity(PATTERNS, 0);
        std::vector<u32> least(PATTERNS, 0);
        FnvLocal nonempty_ledger;
        u64 nonempty = 0;
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            const u32 pattern = response[rank];
            ++multiplicity[pattern];
            const u32 repair = unrank_colex8(rank);
            if (least[pattern] == 0 || repair < least[pattern])
                least[pattern] = repair;
            if (pattern != 0) {
                ++nonempty;
                nonempty_ledger.add(repair);
                nonempty_ledger.add(pattern);
            }
        }
        require(std::accumulate(multiplicity.begin(), multiplicity.end(),
                                u64{0}) == EXPECTED_REPAIRS,
                "response multiplicity census changed");

        std::vector<u32> present;
        FnvLocal class_ledger;
        for (u32 pattern = 1; pattern <= FULL_RESPONSE_667; ++pattern) {
            if (multiplicity[pattern] == 0) continue;
            present.push_back(pattern);
            class_ledger.add(pattern);
            class_ledger.add(multiplicity[pattern]);
            class_ledger.add(least[pattern]);
        }

        std::vector<unsigned char> has_superset(PATTERNS, 0);
        for (u32 pattern : present) has_superset[pattern] = 1;
        for (std::size_t bit = 0; bit < OBLIGATIONS_667; ++bit) {
            for (u32 pattern = 0; pattern <= FULL_RESPONSE_667; ++pattern) {
                if ((pattern & (u32{1} << bit)) == 0)
                    has_superset[pattern] |=
                        has_superset[pattern | (u32{1} << bit)];
            }
        }
        std::vector<u32> maximal;
        FnvLocal maximal_ledger;
        for (u32 pattern : present) {
            bool strict = false;
            for (std::size_t bit = 0; bit < OBLIGATIONS_667; ++bit) {
                if ((pattern & (u32{1} << bit)) == 0 &&
                    has_superset[pattern | (u32{1} << bit)]) {
                    strict = true;
                    break;
                }
            }
            if (!strict) {
                maximal.push_back(pattern);
                maximal_ledger.add(pattern);
                maximal_ledger.add(multiplicity[pattern]);
                maximal_ledger.add(least[pattern]);
            }
        }
        require(!maximal.empty(), "no maximal response classes");

        constexpr int INF = 1000;
        std::vector<int> distance(PATTERNS, INF);
        std::vector<u32> previous(PATTERNS, 0);
        std::vector<u32> previous_pattern(PATTERNS, 0);
        distance[0] = 0;
        for (u32 covered = 0; covered <= FULL_RESPONSE_667; ++covered) {
            if (distance[covered] == INF) continue;
            for (u32 pattern : maximal) {
                const u32 next = covered | pattern;
                if (distance[next] > distance[covered] + 1) {
                    distance[next] = distance[covered] + 1;
                    previous[next] = covered;
                    previous_pattern[next] = pattern;
                }
            }
        }
        require(distance[FULL_RESPONSE_667] < INF,
                "response quotient has no cover");
        std::vector<u32> cover_patterns;
        for (u32 state = FULL_RESPONSE_667; state != 0;
             state = previous[state]) {
            require(previous_pattern[state] != 0,
                    "cover predecessor missing");
            cover_patterns.push_back(previous_pattern[state]);
        }
        std::vector<u32> cover_masks;
        for (u32 pattern : cover_patterns) cover_masks.push_back(least[pattern]);
        std::sort(cover_masks.begin(), cover_masks.end());
        require(cover_masks.size() ==
                    static_cast<std::size_t>(distance[FULL_RESPONSE_667]),
                "cover recovery changed");

        std::array<u32, OBLIGATIONS_667> conflicts{};
        for (u32 pattern : present) {
            for (std::size_t left = 0; left < OBLIGATIONS_667; ++left) {
                if ((pattern & (u32{1} << left)) == 0) continue;
                conflicts[left] |= pattern & ~(u32{1} << left);
            }
        }
        u32 packing = 0;
        for (u32 subset = 1; subset <= FULL_RESPONSE_667; ++subset) {
            bool valid = true;
            for (std::size_t bit = 0; bit < OBLIGATIONS_667; ++bit) {
                if ((subset & (u32{1} << bit)) != 0 &&
                    (subset & conflicts[bit]) != 0) {
                    valid = false;
                    break;
                }
            }
            if (valid && (std::popcount(subset) > std::popcount(packing) ||
                          (std::popcount(subset) == std::popcount(packing) &&
                           subset < packing)))
                packing = subset;
        }
        require(std::popcount(packing) <= distance[FULL_RESPONSE_667],
                "packing exceeds quotient cover");

        FnvLocal cover_ledger;
        for (u32 mask : cover_masks) cover_ledger.add(mask);
        FnvLocal packing_ledger;
        for (std::size_t bit = 0; bit < OBLIGATIONS_667; ++bit) {
            if ((packing & (u32{1} << bit)) != 0) {
                packing_ledger.add(bit);
                packing_ledger.add(FAILURES_667[bit]);
            }
        }
        require(active.count == UINT64_C(1431416) &&
                    active.fnv == UINT64_C(0xa244bee8f3ec231c) &&
                    incidences == UINT64_C(7093) && nonempty == 3696 &&
                    multiplicity[0] == UINT64_C(5849229) &&
                    nonempty_ledger.state == UINT64_C(0x04759f17896f1ef7) &&
                    present.size() == 118 &&
                    class_ledger.state == UINT64_C(0xd60a7fa42a3422b5) &&
                    maximal.size() == 5 &&
                    maximal_ledger.state == UINT64_C(0x9ca2982b7a0df199) &&
                    distance[FULL_RESPONSE_667] == 2 &&
                    cover_masks == std::vector<u32>{UINT32_C(0x02a05206),
                                                    UINT32_C(0x10e05240)} &&
                    cover_ledger.state == UINT64_C(0x297bb05f4bfdfa99) &&
                    packing == UINT32_C(0x00012000) &&
                    packing_ledger.state == UINT64_C(0x7efa8b70a985ed41),
                "frozen endpoint667 response quotient changed");

        std::cout << "ENDPOINT667_FULL_RESPONSE_ATLAS_V1\n"
                  << "PAIR 256,667 FAILURES " << FAILURES_667.size()
                  << " FAILURE_FNV " << std::hex << failure_ledger.state
                  << " UNION " << union_mask << std::dec
                  << " UNION_SIZE 17\n"
                  << "REPAIRS " << EXPECTED_REPAIRS << " ACTIVE "
                  << active.count << " ACTIVE_FNV " << std::hex << active.fnv
                  << std::dec << " ZETA_OPERATIONS "
                  << active.zeta_operations << '\n'
                  << "COMPLEMENT_CHECKS "
                  << FAILURES_667.size() * DISJOINT_REPAIRS_PER_BODY
                  << " ACTIVE_INCIDENCES " << incidences
                  << " NONEMPTY_REPAIRS " << nonempty
                  << " EMPTY_REPAIRS " << multiplicity[0]
                  << " RESPONSE_FNV " << std::hex
                  << nonempty_ledger.state << std::dec << '\n'
                  << "NONEMPTY_CLASSES " << present.size()
                  << " CLASS_FNV " << std::hex << class_ledger.state
                  << std::dec << " MAXIMAL_CLASSES " << maximal.size()
                  << " MAXIMAL_FNV " << std::hex << maximal_ledger.state
                  << std::dec << '\n';
        for (u32 pattern : maximal) {
            std::cout << "MAXIMAL_PATTERN " << std::hex << pattern
                      << " LEAST " << least[pattern] << std::dec
                      << " MULTIPLICITY " << multiplicity[pattern]
                      << " COVER " << std::popcount(pattern) << '\n';
        }
        std::cout << "MINIMUM_AUGMENTATION "
                  << distance[FULL_RESPONSE_667] << " COVER_MASKS ";
        for (std::size_t index = 0; index < cover_masks.size(); ++index) {
            if (index != 0) std::cout << ',';
            std::cout << std::hex << cover_masks[index] << std::dec;
        }
        std::cout << " COVER_FNV " << std::hex << cover_ledger.state
                  << std::dec << '\n'
                  << "PACKING_SIZE " << std::popcount(packing)
                  << " PACKING_BITS " << std::hex << packing
                  << " PACKING_FNV " << packing_ledger.state << std::dec
                  << " OBLIGATION_INDICES ";
        bool first = true;
        for (std::size_t bit = 0; bit < OBLIGATIONS_667; ++bit) {
            if ((packing & (u32{1} << bit)) == 0) continue;
            if (!first) std::cout << ',';
            first = false;
            std::cout << bit;
        }
        std::cout << '\n'
                  << "LOWER_BOUND " << std::popcount(packing)
                  << " BY_RESPONSE_PACKING MATCHES_COVER "
                  << (std::popcount(packing) == distance[FULL_RESPONSE_667]
                          ? "YES" : "NO")
                  << '\n'
                  << "VERDICT PASS EXACT_COMPLETE_RESPONSE_QUOTIENT_AND_COVER\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT667_ATLAS_ERROR " << error.what() << '\n';
        return 1;
    }
}

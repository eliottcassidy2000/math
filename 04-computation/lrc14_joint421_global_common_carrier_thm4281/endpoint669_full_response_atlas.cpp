// Complete labelled response quotient for the four failed bodies at
// (256,669), plus a detached union-complement one-atom minimum certificate.

#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

namespace {

constexpr std::array<u32, 4> FAILURES_669 = {
    UINT32_C(0x0f047001), UINT32_C(0x1f302001),
    UINT32_C(0x1f803001), UINT32_C(0x37003088)};

i128 response_atom_mass(const AtomData& atoms, u32 repair) {
    i128 mass = 0;
    for (const auto& [failure, value] : atoms.mass) {
        if ((failure & ~repair) == 0) mass += value;
    }
    return mass;
}

template <class Callback>
void choose_allowed_rec(const std::vector<int>& positions, int start,
                        int chosen, u32 mask, Callback& callback) {
    if (chosen == 8) {
        callback(mask);
        return;
    }
    const int needed = 8 - chosen;
    for (int index = start;
         index <= static_cast<int>(positions.size()) - needed; ++index) {
        choose_allowed_rec(positions, index + 1, chosen + 1,
                           mask | (u32{1} << positions[index]), callback);
    }
}

}  // namespace

int main() {
    try {
        std::cout.setf(std::ios::unitbuf);
        init_choose8_local();
        FnvLocal failure_ledger;
        u32 union_mask = 0;
        for (u32 body : FAILURES_669) {
            require(std::popcount(body) == 9, "failed body arity changed");
            failure_ledger.add(body);
            union_mask |= body;
        }
        require(failure_ledger.state == UINT64_C(0x6a4e5c7582d44240) &&
                    union_mask == UINT32_C(0x3fb47089) &&
                    std::popcount(union_mask) == 16,
                "endpoint-669 failure universe changed");

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        const ActiveUniverse active = build_active_universe(cells, 256, 669);
        std::vector<unsigned char> response(EXPECTED_REPAIRS, 0);
        u64 incidences = 0;
        for (std::size_t obligation = 0; obligation < FAILURES_669.size();
             ++obligation) {
            enumerate_disjoint_repairs(
                FAILURES_669[obligation], [&](u32, u64 rank) {
                    if (!active.active[rank]) return;
                    response[rank] |= static_cast<unsigned char>(1u << obligation);
                    ++incidences;
                });
        }

        std::array<u64, 16> multiplicity{};
        std::array<u32, 16> least{};
        FnvLocal nonempty_ledger;
        u64 nonempty = 0;
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            const unsigned pattern = response[rank];
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

        const i64 g = std::gcd(256, 669);
        const PrimitivePair primitive = build_primitive(256, 669);
        const AtomData atoms = build_cocycle_atoms(cells, primitive, g);
        const i128 denominator =
            static_cast<i128>(primitive.grid) * g * COMMON;
        std::vector<int> complement;
        for (int bit = 0; bit < 30; ++bit) {
            if ((union_mask & (u32{1} << bit)) == 0) complement.push_back(bit);
        }
        require(complement.size() == 14,
                "union-complement size changed");
        u64 complement_repairs = 0;
        u64 complement_active = 0;
        u32 least_active = 0;
        u32 closest_repair = 0;
        i128 closest_margin = std::numeric_limits<i128>::min();
        FnvLocal complement_active_ledger;
        auto test_union_complement = [&](u32 repair) {
            ++complement_repairs;
            require(std::popcount(repair) == 8 &&
                        (repair & union_mask) == 0,
                    "bad union-complement repair");
            const i128 margin = static_cast<i128>(63) *
                                    response_atom_mass(atoms, repair) -
                                static_cast<i128>(4) * denominator;
            if (closest_repair == 0 || margin > closest_margin ||
                (margin == closest_margin && repair < closest_repair)) {
                closest_repair = repair;
                closest_margin = margin;
            }
            if (margin < 0) return;
            ++complement_active;
            complement_active_ledger.add(repair);
            if (least_active == 0 || repair < least_active)
                least_active = repair;
        };
        choose_allowed_rec(complement, 0, 0, 0, test_union_complement);
        require(complement_repairs == UINT64_C(3003),
                "union-complement repair census changed");
        require(complement_active == multiplicity[15] &&
                    least_active == least[15],
                "one-atom/full-quotient cross-check failed");

        std::vector<unsigned> present_patterns;
        for (unsigned pattern = 1; pattern < 16; ++pattern)
            if (multiplicity[pattern] != 0) present_patterns.push_back(pattern);
        std::vector<u32> minimum_witness;
        for (std::size_t target = 1;
             target <= present_patterns.size() && minimum_witness.empty();
             ++target) {
            std::vector<unsigned> chosen;
            std::vector<u32> best_masks;
            std::function<void(std::size_t, unsigned)> visit =
                [&](std::size_t start, unsigned covered) {
                    if (chosen.size() == target) {
                        if (covered != 15) return;
                        std::vector<u32> masks;
                        for (unsigned pattern : chosen)
                            masks.push_back(least[pattern]);
                        std::sort(masks.begin(), masks.end());
                        if (best_masks.empty() || masks < best_masks)
                            best_masks = std::move(masks);
                        return;
                    }
                    const std::size_t needed = target - chosen.size();
                    for (std::size_t index = start;
                         index + needed <= present_patterns.size(); ++index) {
                        chosen.push_back(present_patterns[index]);
                        visit(index + 1, covered | present_patterns[index]);
                        chosen.pop_back();
                    }
                };
            visit(0, 0);
            minimum_witness = std::move(best_masks);
        }
        require(!minimum_witness.empty(), "response quotient has no cover");

        FnvLocal class_ledger;
        u64 nonempty_classes = 0;
        std::cout << "ENDPOINT669_FULL_RESPONSE_ATLAS_V1\n"
                  << "PAIR 256,669 FAILURES 4 FAILURE_FNV " << std::hex
                  << failure_ledger.state << " UNION " << union_mask
                  << std::dec << " UNION_SIZE 16\n"
                  << "REPAIRS " << EXPECTED_REPAIRS << " ACTIVE "
                  << active.count << " ACTIVE_FNV " << std::hex << active.fnv
                  << std::dec << " ZETA_OPERATIONS "
                  << active.zeta_operations << '\n'
                  << "COMPLEMENT_CHECKS "
                  << FAILURES_669.size() * DISJOINT_REPAIRS_PER_BODY
                  << " ACTIVE_INCIDENCES " << incidences
                  << " NONEMPTY_REPAIRS " << nonempty
                  << " EMPTY_REPAIRS " << multiplicity[0]
                  << " RESPONSE_FNV " << std::hex
                  << nonempty_ledger.state << std::dec << '\n';
        for (unsigned pattern = 1; pattern < 16; ++pattern) {
            if (multiplicity[pattern] == 0) continue;
            ++nonempty_classes;
            class_ledger.add(pattern);
            class_ledger.add(multiplicity[pattern]);
            class_ledger.add(least[pattern]);
            std::cout << "PATTERN " << std::hex << pattern << " LEAST "
                      << least[pattern] << std::dec << " MULTIPLICITY "
                      << multiplicity[pattern] << " COVER "
                      << std::popcount(pattern) << '\n';
        }
        std::cout << "NONEMPTY_CLASSES " << nonempty_classes
                  << " CLASS_FNV " << std::hex << class_ledger.state
                  << std::dec << '\n'
                  << "UNION_COMPLEMENT CANDIDATES " << complement_repairs
                  << " ACTIVE " << complement_active << " ACTIVE_FNV "
                  << std::hex << complement_active_ledger.state
                  << " LEAST " << least_active << " CLOSEST "
                  << closest_repair << std::dec << " CLOSEST_MARGIN_NUM "
                  << decimal(closest_margin)
                  << " DEN " << decimal(denominator) << '\n'
                  << "ONE_ATOM "
                  << (complement_active == 0 ? "NO_GO" : "EXISTS")
                  << "\nMINIMUM_AUGMENTATION " << minimum_witness.size()
                  << " WITNESSES ";
        for (std::size_t index = 0; index < minimum_witness.size(); ++index) {
            if (index != 0) std::cout << ',';
            std::cout << std::hex << minimum_witness[index] << std::dec;
        }
        std::cout << " ZERO_IMPOSSIBLE INHERITED_FAILURES_4"
                  << " LOWER_BY_COMPLETE_RESPONSE_QUOTIENT\n"
                  << "VERDICT PASS EXACT_COMPLETE_RESPONSE_QUOTIENT_AND_MINIMUM\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT669_ATLAS_ERROR " << error.what() << '\n';
        return 1;
    }
}

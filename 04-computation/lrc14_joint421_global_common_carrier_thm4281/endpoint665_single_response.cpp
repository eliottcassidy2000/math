// Complete one-obligation response scan for the sole endpoint-665 hostile row.

#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

namespace {

constexpr u32 FAILURE_665 = UINT32_C(0x151a3400);

i128 mass_665(const AtomData& atoms, u32 repair) {
    i128 mass = 0;
    for (const auto& [failure, value] : atoms.mass) {
        if ((failure & ~repair) == 0) mass += value;
    }
    return mass;
}

}  // namespace

int main() {
    try {
        std::cout.setf(std::ios::unitbuf);
        init_choose8_local();
        require(std::popcount(FAILURE_665) == 9,
                "endpoint665 failure arity changed");
        FnvLocal failure_ledger;
        failure_ledger.add(FAILURE_665);
        require(failure_ledger.state == UINT64_C(0x183de0540d550fac),
                "endpoint665 failure ledger changed");

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        const ActiveUniverse active = build_active_universe(cells, 520, 665);
        const i64 g = std::gcd(520, 665);
        const PrimitivePair primitive = build_primitive(520 / g, 665 / g);
        const AtomData atoms = build_cocycle_atoms(cells, primitive, g);
        const i128 denominator =
            static_cast<i128>(primitive.grid) * g * COMMON;

        std::vector<u32> responders;
        u32 best = 0;
        i128 best_margin = std::numeric_limits<i128>::min();
        enumerate_disjoint_repairs(FAILURE_665, [&](u32 repair, u64 rank) {
            const i128 margin = static_cast<i128>(63) * mass_665(atoms, repair) -
                                static_cast<i128>(4) * denominator;
            require((margin >= 0) == static_cast<bool>(active.active[rank]),
                    "direct/active-universe disagreement");
            if (best == 0 || margin > best_margin ||
                (margin == best_margin && repair < best)) {
                best = repair;
                best_margin = margin;
            }
            if (margin >= 0) responders.push_back(repair);
        });
        require(!responders.empty(), "one-obligation response is empty");
        std::sort(responders.begin(), responders.end());
        require(std::adjacent_find(responders.begin(), responders.end()) ==
                    responders.end(),
                "duplicate response repair");
        FnvLocal responder_ledger;
        for (u32 repair : responders) responder_ledger.add(repair);
        const u32 least = responders.front();
        const i128 least_margin =
            static_cast<i128>(63) * mass_665(atoms, least) -
            static_cast<i128>(4) * denominator;
        require(least_margin >= 0 && (least & FAILURE_665) == 0,
                "least response witness changed");
        require(active.count == UINT64_C(1924688) &&
                    active.fnv == UINT64_C(0x07ecf882b7e3a319) &&
                    responders.size() == 1281 &&
                    responder_ledger.state == UINT64_C(0x98d0e1a04abf4e2f) &&
                    least == UINT32_C(0x0004409f) &&
                    least_margin == static_cast<i128>(1490243161299984000LL) &&
                    best == UINT32_C(0x2284c003) &&
                    best_margin == static_cast<i128>(5088108568163100480LL) &&
                    denominator == static_cast<i128>(17661820193412595200ULL),
                "frozen endpoint665 response ledger changed");

        std::cout << "ENDPOINT665_SINGLE_RESPONSE_V1\n"
                  << "PAIR 520,665 FAILURE " << std::hex << FAILURE_665
                  << " FAILURE_FNV " << failure_ledger.state << std::dec
                  << "\nREPAIRS " << EXPECTED_REPAIRS << " ACTIVE "
                  << active.count << " ACTIVE_FNV " << std::hex << active.fnv
                  << std::dec << " ZETA_OPERATIONS "
                  << active.zeta_operations << '\n'
                  << "DISJOINT_CANDIDATES " << DISJOINT_REPAIRS_PER_BODY
                  << " ACTIVE_RESPONDERS " << responders.size()
                  << " RESPONDER_FNV " << std::hex << responder_ledger.state
                  << " LEAST " << least << std::dec
                  << " LEAST_MARGIN_NUM " << decimal(least_margin)
                  << " DEN " << decimal(denominator) << '\n'
                  << "BEST " << std::hex << best << std::dec
                  << " BEST_MARGIN_NUM " << decimal(best_margin)
                  << " DEN " << decimal(denominator) << '\n'
                  << "MINIMUM_AUGMENTATION 1 WITNESS " << std::hex << least
                  << std::dec
                  << " ZERO_IMPOSSIBLE INHERITED_FAILURES_1\n"
                  << "VERDICT PASS EXACT_COMPLETE_SINGLE_RESPONSE_AND_MINIMUM\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT665_RESPONSE_ERROR " << error.what() << '\n';
        return 1;
    }
}

#include "04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/carrier_scan_support.cpp"

#include <iomanip>

namespace {

i128 direct_mass632(const AtomData& atoms, u32 mask) {
    i128 mass = 0;
    for (const auto& [atom, value] : atoms.mass)
        if ((atom & ~mask) == 0) mass += value;
    return mass;
}

void add_signed_i128(FnvLocal& f, i128 value) {
    const __int128_t signed_value = value;
    const __uint128_t bits = static_cast<__uint128_t>(signed_value);
    f.add(static_cast<u64>(bits));
    f.add(static_cast<u64>(bits >> 64));
}

}  // namespace

int main() {
    try {
        init_choose8_local();
        constexpr u32 body = UINT32_C(0x1d106401);
        require(std::popcount(body) == 9, "obstruction body rank changed");
        const std::vector<Cell> cells = build_pool_cells();
        const ActiveUniverse active = build_active_universe(cells, 256, 632);
        const i64 g = std::gcd(256, 632);
        const PrimitivePair primitive = build_primitive(256 / g, 632 / g);
        const AtomData atoms = build_cocycle_atoms(cells, primitive, g);
        const i128 denominator = static_cast<i128>(primitive.grid) * g * COMMON;

        u64 candidates = 0;
        u64 active_candidates = 0;
        u64 equalities = 0;
        FnvLocal mask_ledger;
        FnvLocal margin_ledger;
        i128 maximum_margin = std::numeric_limits<i128>::min();
        i128 minimum_margin = std::numeric_limits<i128>::max();
        u32 maximum_mask = 0;
        u32 minimum_mask = 0;
        enumerate_disjoint_repairs(body, [&](u32 mask, u64 rank) {
            ++candidates;
            mask_ledger.add(mask);
            const i128 margin = static_cast<i128>(63) * direct_mass632(atoms, mask) -
                                static_cast<i128>(4) * denominator;
            const bool direct_active = margin >= 0;
            require(direct_active == static_cast<bool>(active.active[rank]),
                    "direct/full activity mismatch");
            active_candidates += direct_active;
            equalities += margin == 0;
            if (margin > maximum_margin ||
                (margin == maximum_margin && (maximum_mask == 0 || mask < maximum_mask))) {
                maximum_margin = margin;
                maximum_mask = mask;
            }
            if (margin < minimum_margin ||
                (margin == minimum_margin && (minimum_mask == 0 || mask < minimum_mask))) {
                minimum_margin = margin;
                minimum_mask = mask;
            }
            margin_ledger.add(mask);
            add_signed_i128(margin_ledger, margin);
        });
        require(candidates == DISJOINT_REPAIRS_PER_BODY && active_candidates == 0 &&
                    equalities == 0 && maximum_margin < 0,
                "zero-response obstruction failed");
        std::cout << "LRC14_632_ZERO_RESPONSE_V1\n"
                  << "PAIR 256,632 BODY " << std::hex << std::setw(8)
                  << std::setfill('0') << body << std::dec << std::setfill(' ')
                  << '\n'
                  << "ACTIVE_UNIVERSE " << active.count << " FNV " << std::hex
                  << active.fnv << std::dec << '\n'
                  << "DISJOINT_RANK8 " << candidates << " MASK_FNV " << std::hex
                  << mask_ledger.state << std::dec << " ACTIVE "
                  << active_candidates << " EQUALITIES " << equalities << '\n'
                  << "MARGIN_FNV " << std::hex << margin_ledger.state << std::dec
                  << " MAXIMUM " << decimal(maximum_margin) << " MAX_MASK "
                  << std::hex << std::setw(8) << std::setfill('0')
                  << maximum_mask << std::dec << std::setfill(' ')
                  << " MINIMUM " << decimal(minimum_margin) << " MIN_MASK "
                  << std::hex << std::setw(8) << std::setfill('0')
                  << minimum_mask << std::dec << std::setfill(' ') << '\n'
                  << "VERDICT PASS FINITE_EXACT_NO_RANK8_RESPONSE\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "LRC14_632_ZERO_RESPONSE_ERROR " << e.what() << '\n';
        return 1;
    }
}

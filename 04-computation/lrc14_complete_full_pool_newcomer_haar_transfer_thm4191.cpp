// Primary exact full-pool transversal audit for THM-4191.

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#endif
#define main thm4188_canonical_main
#include "lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.cpp"
#undef main
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace {

constexpr u64 EXPECTED_GLOBAL_LEDGER = UINT64_C(0x53664fb38b90d1c0);

struct ExhaustiveResult {
    u64 candidates = 0;
    u64 edge_checks = 0;
    u64 maximum = 0;
    u32 closest = 0;
    u32 missed = 0;
    u32 cover = 0;
    u64 fingerprint = FNV1A64_OFFSET;
};

ExhaustiveResult exhaust_all_ten(std::vector<u32> edges) {
    auto mix = [](u64 value) {
        value += UINT64_C(0x9e3779b97f4a7c15);
        value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
        return value ^ (value >> 31);
    };
    std::sort(edges.begin(), edges.end(), [&](u32 left, u32 right) {
        const u64 left_hash = mix(left);
        const u64 right_hash = mix(right);
        if (left_hash != right_hash) return left_hash < right_hash;
        return left < right;
    });
    ExhaustiveResult result;
    Fnv1a64 digest;
    u32 body = (u32{1} << 10) - 1;
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
        ++result.candidates;
        result.edge_checks += checked;
        if (checked > result.maximum) {
            result.maximum = checked;
            result.closest = body;
            result.missed = missed;
        }
        digest.add_u64_le(body);
        digest.add_u64_le(missed);
        digest.add_u64_le(checked);
        if (missed == 0) {
            result.cover = body;
            break;
        }
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    result.fingerprint = digest.value();
    return result;
}

bool is_cover(u32 candidate, const std::vector<u32>& edges) {
    return std::all_of(edges.begin(), edges.end(), [&](u32 edge) {
        return (candidate & edge) != 0;
    });
}

u32 mask_of_values(std::initializer_list<int> values) {
    u32 mask = 0;
    for (int value : values) {
        const auto found = std::find(POOL.begin(), POOL.end(), value);
        require(found != POOL.end(), "control label missing");
        mask |= u32{1} << std::distance(POOL.begin(), found);
    }
    return mask;
}

}  // namespace

int main() {
    const auto cells = build_pool_cells();
    Fnv1a64 ledger;
    std::cout << "AUDIT full_pool_resonance_transversal_primary\n";

    const AtomMass mass50 = build_atom_mass(cells, 50, 7);
    const EdgeLayer q50e6 = build_layer(mass50, 50, 6);
    const EdgeLayer q50e7 = build_layer(mass50, 50, 7);
    require(std::find(POOL.begin(), POOL.end(), 50) == POOL.end(),
            "q50 ceased to be an outside-pool newcomer");
    require(q50e6.edges.size() == 85324 && q50e6.equalities == 0,
            "q50 E6 control layer changed");
    require(q50e7.edges.size() == 821737 && q50e7.equalities == 0,
            "q50 E7 base layer changed");
    const u32 control6 = mask_of_values({80,85,88,95,145,168,193,240,252,286});
    require(is_cover(control6, q50e6.edges), "q50 E6 explicit cover changed");
    const ExhaustiveResult base = exhaust_all_ten(q50e7.edges);
    require(base.candidates == UINT64_C(30045015),
            "q50 full body universe changed");
    require(base.cover == 0, "q50 E7 has a ten-cover");
    require(UINT64_C(8436285) + 3 * UINT64_C(4686825) +
                3 * UINT64_C(2220075) + UINT64_C(888030) ==
            UINT64_C(30045015), "original-count decomposition changed");
    ledger.add_u64_le(50);
    ledger.add_u64_le(7);
    ledger.add_u64_le(q50e7.edges.size());
    ledger.add_u64_le(q50e7.equalities);
    ledger.add_u64_le(base.candidates);
    ledger.add_u64_le(base.edge_checks);
    ledger.add_u64_le(base.maximum);
    ledger.add_u64_le(base.closest);
    ledger.add_u64_le(base.missed);
    ledger.add_u64_le(base.fingerprint);
    std::cout << "CONTROL_Q50_E6 COVER " << labels(control6)
              << " SIZE " << std::popcount(control6) << '\n';
    std::cout << "BODY_DECOMPOSITION 8436285+3*4686825+3*2220075+888030=30045015\n";
    std::cout << "BASE_Q 50 E7 " << q50e7.edges.size()
              << " EQUALITIES " << q50e7.equalities
              << " BODIES " << base.candidates
              << " COVER NONE EDGE_CHECKS " << base.edge_checks
              << " MAX_CHECKS " << base.maximum
              << " CLOSEST_BODY " << labels(base.closest)
              << " MISSED_EDGE " << labels(base.missed)
              << " ROW_FNV1A64_LE " << hex64(base.fingerprint) << '\n';

    for (int q : EXPECTED_Q7_FAILURES) {
        require(std::find(POOL.begin(), POOL.end(), q) == POOL.end(),
                "resonance is not an outside-pool newcomer");
        const AtomMass mass = build_atom_mass(cells, q, 6);
        const EdgeLayer e6 = build_layer(mass, q, 6);
        const ExhaustiveResult result = exhaust_all_ten(e6.edges);
        require(result.cover == 0, "resonance E6 has ten-cover on full pool");
        require(result.candidates == UINT64_C(30045015),
                "full C(30,10) universe not exhausted");
        ledger.add_u64_le(q);
        ledger.add_u64_le(e6.edges.size());
        ledger.add_u64_le(e6.equalities);
        ledger.add_u64_le(result.candidates);
        ledger.add_u64_le(result.edge_checks);
        ledger.add_u64_le(result.maximum);
        ledger.add_u64_le(result.closest);
        ledger.add_u64_le(result.missed);
        ledger.add_u64_le(result.fingerprint);
        std::cout << "Q " << q << " E6 " << e6.edges.size()
                  << " EQUALITIES " << e6.equalities
                  << " BODIES " << result.candidates
                  << " COVER NONE EDGE_CHECKS " << result.edge_checks
                  << " MAX_CHECKS " << result.maximum
                  << " CLOSEST_BODY " << labels(result.closest)
                  << " MISSED_EDGE " << labels(result.missed)
                  << " ROW_FNV1A64_LE " << hex64(result.fingerprint) << '\n';
    }
    require(ledger.value() == EXPECTED_GLOBAL_LEDGER,
            "primary semantic ledger changed");
    std::cout << "GLOBAL_FNV1A64_LE " << hex64(ledger.value()) << '\n';
    std::cout << "VERDICT FULL_POOL_TAU_GT_10_BASE_AND_RESONANCES\n";
}

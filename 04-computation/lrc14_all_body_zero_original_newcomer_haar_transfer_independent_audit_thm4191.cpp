// Independent exact audit for THM-4191.
//
// Geometry comes from the THM-4188 independent joint-wall implementation:
// every newcomer wall is explicitly refined with the fixed-pool walls.  The
// new cover path uses recursive subset generation rather than Gosper masks.

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#endif
#define main thm4188_independent_canonical_main
#include "lrc14_all_newcomer_zero_original_joint_wall_audit_thm4188.cpp"
#undef main
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace {

constexpr u32 ORIGINAL_MASK_PROBE =
    (u32{1} << 15) | (u32{1} << 16) | (u32{1} << 18);
constexpr u32 U_MASK_PROBE = ((u32{1} << 30) - 1) & ~ORIGINAL_MASK_PROBE;
constexpr u32 HOSTILE_Q50_E6_COVER_PROBE = UINT32_C(0x17187400);

bool hits_every_edge_probe(u32 body, const std::vector<u32>& edges) {
    for (u32 edge : edges) {
        if ((body & edge) == 0) return false;
    }
    return true;
}

u32 lift_local_probe(u32 local) {
    u32 global = 0;
    int local_index = 0;
    for (int global_index = 0; global_index < 30; ++global_index) {
        if ((ORIGINAL_MASK_PROBE & (u32{1} << global_index)) != 0) continue;
        if ((local & (u32{1} << local_index)) != 0) {
            global |= u32{1} << global_index;
        }
        ++local_index;
    }
    require(local_index == 27, "independent local/global lift changed");
    return global;
}

struct IndependentExhaustion {
    u64 candidates = 0;
    u64 edge_checks = 0;
    u64 max_edge_checks = 0;
    u32 closest_body = 0;
    u32 closest_miss = 0;
    u32 cover = 0;
    u64 fingerprint = FNV_OFFSET;
};

IndependentExhaustion independent_exhaust(std::vector<u32> edges) {
    std::stable_sort(edges.begin(), edges.end(), [](u32 left, u32 right) {
        const int ls = std::popcount(left & U_MASK_PROBE);
        const int rs = std::popcount(right & U_MASK_PROBE);
        if (ls != rs) return ls < rs;
        return left < right;
    });
    IndependentExhaustion result;
    Ledger ledger;
    auto inspect = [&](u32 local) {
        const u32 body = lift_local_probe(local);
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
        if (checked > result.max_edge_checks) {
            result.max_edge_checks = checked;
            result.closest_body = body;
            result.closest_miss = missed;
        }
        ledger.add(body);
        ledger.add(missed);
        ledger.add(checked);
        if (missed == 0) result.cover = body;
    };
    for_each_k_subset(27, 10, inspect);
    result.fingerprint = ledger.value();
    return result;
}

}  // namespace

int main() {
    const FixedGeometry fixed = build_fixed_geometry();
    const Support support = build_support(fixed);
    Ledger global;
    std::cout << "AUDIT independent_joint_wall_all_23_resonance_E6_transversal\n";

    const JointMass mass50 = build_joint_mass(fixed, support, 50);
    const Layer hostile6 = build_layer(mass50, support, 6);
    require(std::popcount(HOSTILE_Q50_E6_COVER_PROBE) == 10,
            "independent q50 E6 hostile cardinality changed");
    require((HOSTILE_Q50_E6_COVER_PROBE & ORIGINAL_MASK_PROBE) == 0,
            "independent q50 E6 hostile left zero-original ground");
    require(hits_every_edge_probe(HOSTILE_Q50_E6_COVER_PROBE,
                                  hostile6.edges),
            "independent q50 E6 hostile ceased to cover");
    std::cout << "HOSTILE_Q50 E6 " << hostile6.edges.size()
              << " TEN_COVER " << labels(HOSTILE_Q50_E6_COVER_PROBE) << '\n';
    const Layer base7 = build_layer(mass50, support, 7);
    const IndependentExhaustion base_result = independent_exhaust(base7.edges);
    require(base_result.candidates == UINT64_C(8436285),
            "independent q50 body universe changed");
    require(base_result.cover == 0, "independent q50 E7 has ten-cover");
    global.add(50);
    global.add(7);
    global.add(base7.edges.size());
    global.add(base_result.fingerprint);
    std::cout << "BASE_Q 50 JOINT_DENOMINATOR " << mass50.denominator
              << " JOINT_CELLS " << mass50.joint_cells
              << " E7 " << base7.edges.size()
              << " EQUALITIES " << base7.equalities
              << " TEN_BODIES " << base_result.candidates
              << " COVER NONE"
              << " EDGE_CHECKS " << base_result.edge_checks
              << " MAX_CHECKS " << base_result.max_edge_checks
              << " CLOSEST_BODY " << labels(base_result.closest_body)
              << " MISSED_EDGE " << labels(base_result.closest_miss)
              << " ROW_FNV1A64_LE " << hex64(base_result.fingerprint) << '\n';

    for (int q : EXPECTED_Q7) {
        const JointMass mass = build_joint_mass(fixed, support, q);
        const Layer layer = build_layer(mass, support, 6);
        const IndependentExhaustion result = independent_exhaust(layer.edges);
        require(result.candidates == UINT64_C(8436285),
                "independent ten-body universe changed");
        require(result.cover == 0, "independent native E6 has ten-cover");
        require(result.closest_miss != 0, "independent closest miss absent");

        global.add(static_cast<u64>(q));
        global.add(layer.edges.size());
        global.add(layer.equalities);
        global.add(result.candidates);
        global.add(result.edge_checks);
        global.add(result.max_edge_checks);
        global.add(result.closest_body);
        global.add(result.closest_miss);
        global.add(result.fingerprint);

        std::cout << "Q " << q
                  << " JOINT_DENOMINATOR " << mass.denominator
                  << " JOINT_CELLS " << mass.joint_cells
                  << " E6 " << layer.edges.size()
                  << " EQUALITIES " << layer.equalities
                  << " TEN_BODIES " << result.candidates
                  << " COVER NONE"
                  << " EDGE_CHECKS " << result.edge_checks
                  << " MAX_CHECKS " << result.max_edge_checks
                  << " CLOSEST_BODY " << labels(result.closest_body)
                  << " MISSED_EDGE " << labels(result.closest_miss)
                  << " ROW_FNV1A64_LE " << hex64(result.fingerprint) << '\n';
    }
    require(global.value() == UINT64_C(0x18bb68e547da2a9e),
            "THM-4191 independent semantic ledger changed");
    std::cout << "GLOBAL_FNV1A64_LE " << hex64(global.value()) << '\n';
    std::cout << "VERDICT INDEPENDENT_ALL_RESONANCES_NATIVE_E6_TAU_GT_10\n";
}

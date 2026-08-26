// Independent exact full-pool transversal audit for THM-4191.

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

constexpr u64 EXPECTED_GLOBAL_LEDGER_FULL = UINT64_C(0xe9ae70b15ef96dfb);

struct ExhaustiveFull {
    u64 candidates = 0;
    u64 checks = 0;
    u64 maximum = 0;
    u32 closest = 0;
    u32 missed = 0;
    u32 cover = 0;
    u64 fingerprint = FNV_OFFSET;
};

u64 mix_edge(u64 value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

ExhaustiveFull exhaust_full(std::vector<u32> edges) {
    std::stable_sort(edges.begin(), edges.end(), [](u32 left, u32 right) {
        const u64 lh = mix_edge(left);
        const u64 rh = mix_edge(right);
        if (lh != rh) return lh < rh;
        return left < right;
    });
    ExhaustiveFull result;
    Ledger digest;
    auto inspect = [&](u32 body) {
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
        result.checks += checked;
        if (checked > result.maximum) {
            result.maximum = checked;
            result.closest = body;
            result.missed = missed;
        }
        digest.add(body);
        digest.add(missed);
        digest.add(checked);
        if (missed == 0) result.cover = body;
    };
    for_each_k_subset(30, 10, inspect);
    result.fingerprint = digest.value();
    return result;
}

bool covers(u32 body, const std::vector<u32>& edges) {
    return std::all_of(edges.begin(), edges.end(), [&](u32 edge) {
        return (body & edge) != 0;
    });
}

}  // namespace

int main() {
    const FixedGeometry fixed = build_fixed_geometry();
    const Support support = build_support(fixed);
    Ledger global;
    std::cout << "AUDIT full_pool_resonance_transversal_independent_joint_wall\n";

    const JointMass mass50 = build_joint_mass(fixed, support, 50);
    const Layer q50e6 = build_layer(mass50, support, 6);
    const Layer q50e7 = build_layer(mass50, support, 7);
    require(!is_pool_label(50), "independent q50 is a pool label");
    require(q50e6.edges.size() == 85324 && q50e6.equalities == 0,
            "independent q50 E6 changed");
    require(q50e7.edges.size() == 821737 && q50e7.equalities == 0,
            "independent q50 E7 changed");
    const u32 control = label_mask({80,85,88,95,145,168,193,240,252,286});
    require(covers(control, q50e6.edges), "independent q50 E6 control changed");
    const ExhaustiveFull base = exhaust_full(q50e7.edges);
    require(base.candidates == UINT64_C(30045015),
            "independent q50 body universe changed");
    require(base.cover == 0, "independent q50 E7 has a ten-cover");
    require(UINT64_C(8436285) + 3 * UINT64_C(4686825) +
                3 * UINT64_C(2220075) + UINT64_C(888030) ==
            UINT64_C(30045015), "independent decomposition changed");
    global.add(50);
    global.add(7);
    global.add(q50e7.edges.size());
    global.add(base.candidates);
    global.add(base.checks);
    global.add(base.maximum);
    global.add(base.closest);
    global.add(base.missed);
    global.add(base.fingerprint);
    std::cout << "CONTROL_Q50_E6 COVER " << labels(control) << '\n';
    std::cout << "BODY_DECOMPOSITION 8436285+3*4686825+3*2220075+888030=30045015\n";
    std::cout << "BASE_Q 50 JOINT_DENOMINATOR " << mass50.denominator
              << " JOINT_CELLS " << mass50.joint_cells
              << " E7 " << q50e7.edges.size()
              << " EQUALITIES " << q50e7.equalities
              << " BODIES " << base.candidates
              << " COVER NONE EDGE_CHECKS " << base.checks
              << " MAX_CHECKS " << base.maximum
              << " CLOSEST_BODY " << labels(base.closest)
              << " MISSED_EDGE " << labels(base.missed)
              << " ROW_FNV1A64_LE " << hex64(base.fingerprint) << '\n';

    for (int q : EXPECTED_Q7) {
        require(!is_pool_label(q), "independent resonance is a pool label");
        const JointMass mass = build_joint_mass(fixed, support, q);
        const Layer e6 = build_layer(mass, support, 6);
        const ExhaustiveFull result = exhaust_full(e6.edges);
        require(result.candidates == UINT64_C(30045015),
                "independent resonance body universe changed");
        require(result.cover == 0, "independent resonance E6 has a ten-cover");
        global.add(q);
        global.add(e6.edges.size());
        global.add(e6.equalities);
        global.add(result.candidates);
        global.add(result.checks);
        global.add(result.maximum);
        global.add(result.closest);
        global.add(result.missed);
        global.add(result.fingerprint);
        std::cout << "Q " << q
                  << " JOINT_DENOMINATOR " << mass.denominator
                  << " JOINT_CELLS " << mass.joint_cells
                  << " E6 " << e6.edges.size()
                  << " EQUALITIES " << e6.equalities
                  << " BODIES " << result.candidates
                  << " COVER NONE EDGE_CHECKS " << result.checks
                  << " MAX_CHECKS " << result.maximum
                  << " CLOSEST_BODY " << labels(result.closest)
                  << " MISSED_EDGE " << labels(result.missed)
                  << " ROW_FNV1A64_LE " << hex64(result.fingerprint) << '\n';
    }
    require(global.value() == EXPECTED_GLOBAL_LEDGER_FULL,
            "independent full-pool semantic ledger changed");
    std::cout << "GLOBAL_FNV1A64_LE " << hex64(global.value()) << '\n';
    std::cout << "VERDICT INDEPENDENT_FULL_POOL_TAU_GT_10_ALL_Q_CLASSES\n";
}

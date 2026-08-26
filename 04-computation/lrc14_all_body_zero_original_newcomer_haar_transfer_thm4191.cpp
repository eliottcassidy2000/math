// Primary exact audit for THM-4191.
//
// Include the canonical THM-4188 primary implementation into this translation
// unit so its exact pool-cell integration and edge construction are reused
// byte-for-byte.  The new part is independent of its anchor hierarchy: it
// exhausts all C(27,10) zero-original bodies against each native E6(q).

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

constexpr u32 ORIGINAL_MASK =
    (u32{1} << 15) | (u32{1} << 16) | (u32{1} << 18);
constexpr u32 U_MASK = ((u32{1} << 30) - 1) & ~ORIGINAL_MASK;
constexpr u32 HOSTILE_Q50_E6_COVER = UINT32_C(0x17187400);

bool hits_every_edge(u32 body, const std::vector<u32>& edges) {
    return std::all_of(edges.begin(), edges.end(), [body](u32 edge) {
        return (body & edge) != 0;
    });
}

u32 lift_local_mask(u32 local) {
    u32 global = 0;
    int local_index = 0;
    for (int global_index = 0; global_index < 30; ++global_index) {
        if ((ORIGINAL_MASK & (u32{1} << global_index)) != 0) continue;
        if ((local & (u32{1} << local_index)) != 0) {
            global |= u32{1} << global_index;
        }
        ++local_index;
    }
    require(local_index == 27, "local/global lift changed");
    return global;
}

struct Exhaustion {
    u64 candidates = 0;
    u64 edge_checks = 0;
    u64 max_edge_checks = 0;
    u32 closest_body = 0;
    u32 closest_missed_edge = 0;
    u32 cover = 0;
    u64 fingerprint = FNV1A64_OFFSET;
};

Exhaustion exhaust_ten_covers(std::vector<u32> edges) {
    std::sort(edges.begin(), edges.end(), [](u32 left, u32 right) {
        const int left_size = std::popcount(left & U_MASK);
        const int right_size = std::popcount(right & U_MASK);
        if (left_size != right_size) return left_size < right_size;
        return left < right;
    });

    Exhaustion result;
    Fnv1a64 fingerprint;
    u32 local = (u32{1} << 10) - 1;
    const u32 local_limit = u32{1} << 27;
    while (local < local_limit) {
        const u32 body = lift_local_mask(local);
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
            result.closest_missed_edge = missed;
        }
        fingerprint.add_u64_le(body);
        fingerprint.add_u64_le(missed);
        fingerprint.add_u64_le(checked);
        if (missed == 0) {
            result.cover = body;
            break;
        }
        const u32 next = next_combination(local);
        if (next <= local) break;
        local = next;
    }
    result.fingerprint = fingerprint.value();
    return result;
}

}  // namespace

int main() {
    const auto cells = build_pool_cells();
    Fnv1a64 global;
    std::cout << "AUDIT all_23_resonance_native_E6_projected_transversal\n";
    std::cout << "BODY_UNIVERSE C(27,10)=8436285 ORIGINAL_MASK "
              << labels(ORIGINAL_MASK) << '\n';

    const AtomMass mass50 = build_atom_mass(cells, 50, 7);
    const EdgeLayer hostile6 = build_layer(mass50, 50, 6);
    require(std::popcount(HOSTILE_Q50_E6_COVER) == 10,
            "q50 E6 hostile cardinality changed");
    require((HOSTILE_Q50_E6_COVER & ORIGINAL_MASK) == 0,
            "q50 E6 hostile left the zero-original ground set");
    require(hits_every_edge(HOSTILE_Q50_E6_COVER, hostile6.edges),
            "q50 E6 hostile ceased to be a cover");
    std::cout << "HOSTILE_Q50 E6 " << hostile6.edges.size()
              << " TEN_COVER " << labels(HOSTILE_Q50_E6_COVER) << '\n';
    const EdgeLayer base7 = build_layer(mass50, 50, 7);
    const Exhaustion base_result = exhaust_ten_covers(base7.edges);
    require(base_result.candidates == UINT64_C(8436285),
            "q50 ten-body universe not exhausted");
    require(base_result.cover == 0, "q50 E7 has a ten-cover");
    global.add_u64_le(50);
    global.add_u64_le(7);
    global.add_u64_le(base7.edges.size());
    global.add_u64_le(base_result.fingerprint);
    std::cout << "BASE_Q 50 E7 " << base7.edges.size()
              << " EQUALITIES " << base7.equalities
              << " TEN_BODIES " << base_result.candidates
              << " COVER NONE"
              << " EDGE_CHECKS " << base_result.edge_checks
              << " MAX_CHECKS " << base_result.max_edge_checks
              << " CLOSEST_BODY " << labels(base_result.closest_body)
              << " MISSED_EDGE " << labels(base_result.closest_missed_edge)
              << " ROW_FNV1A64_LE " << hex64(base_result.fingerprint) << '\n';

    for (int q : EXPECTED_Q7_FAILURES) {
        const AtomMass mass = build_atom_mass(cells, q, 6);
        EdgeLayer layer = build_layer(mass, q, 6);
        const Exhaustion result = exhaust_ten_covers(layer.edges);
        require(result.candidates == UINT64_C(8436285),
                "ten-body universe not exhausted");
        require(result.cover == 0, "native E6 has a ten-cover");
        require(result.closest_missed_edge != 0, "closest miss missing");

        global.add_u64_le(static_cast<u64>(q));
        global.add_u64_le(layer.edges.size());
        global.add_u64_le(layer.equalities);
        global.add_u64_le(result.candidates);
        global.add_u64_le(result.edge_checks);
        global.add_u64_le(result.max_edge_checks);
        global.add_u64_le(result.closest_body);
        global.add_u64_le(result.closest_missed_edge);
        global.add_u64_le(result.fingerprint);

        std::cout << "Q " << q
                  << " E6 " << layer.edges.size()
                  << " EQUALITIES " << layer.equalities
                  << " TEN_BODIES " << result.candidates
                  << " COVER NONE"
                  << " EDGE_CHECKS " << result.edge_checks
                  << " MAX_CHECKS " << result.max_edge_checks
                  << " CLOSEST_BODY " << labels(result.closest_body)
                  << " MISSED_EDGE " << labels(result.closest_missed_edge)
                  << " ROW_FNV1A64_LE " << hex64(result.fingerprint) << '\n';
    }
    require(global.value() == UINT64_C(0x14e15a5a00ad0764),
            "THM-4191 primary semantic ledger changed");
    std::cout << "GLOBAL_FNV1A64_LE " << hex64(global.value()) << '\n';
    std::cout << "VERDICT ALL_RESONANCES_NATIVE_E6_TAU_GT_10\n";
}

// FINITE-EXACT scout: sharpen the q=50 global E7 transversal boundary after
// THM-4191. This imports the canonical exact pool integrator as its geometry
// engine and checks explicit upper witnesses against the complete E7 deck.

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

u32 mask(std::initializer_list<int> values) {
    u32 result = 0;
    for (int value : values) {
        const auto it = std::find(POOL.begin(), POOL.end(), value);
        require(it != POOL.end(), "probe label absent from pool");
        result |= u32{1} << std::distance(POOL.begin(), it);
    }
    return result;
}

bool covers(u32 candidate, const std::vector<u32>& edges) {
    return std::all_of(edges.begin(), edges.end(), [&](u32 edge) {
        return (candidate & edge) != 0;
    });
}

void inspect(const char* name, u32 candidate, const std::vector<u32>& edges,
             u64 expected_misses, u32 expected_common) {
    u64 misses = 0;
    u32 intersection = (u32{1} << 30) - 1;
    u32 union_mask = 0;
    u32 first = 0;
    for (u32 edge : edges) {
        if ((candidate & edge) != 0) continue;
        if (misses == 0) first = edge;
        ++misses;
        intersection &= edge;
        union_mask |= edge;
    }
    require(misses == expected_misses, "unexpected missed-edge count");
    require(intersection == expected_common, "unexpected common vertex set");
    std::cout << name << " SIZE " << std::popcount(candidate)
              << " MISSES " << misses
              << " COMMON " << labels(intersection)
              << " UNION " << labels(union_mask)
              << " FIRST " << labels(first) << '\n';
    if (intersection != 0) {
        const u32 vertex = intersection & (~intersection + 1);
        const u32 extended = candidate | vertex;
        require(std::popcount(extended) == std::popcount(candidate) + 1,
                "extension did not add one vertex");
        require(covers(extended, edges), "common-vertex extension is not cover");
        std::cout << name << " EXTENDED_COVER " << labels(extended)
                  << " SIZE " << std::popcount(extended) << '\n';
    }
}

}  // namespace

int main() {
    const auto cells = build_pool_cells();
    const AtomMass mass50 = build_atom_mass(cells, 50, 7);
    const EdgeLayer e7 = build_layer(mass50, 50, 7);
    require(e7.edges.size() == 821737 && e7.equalities == 0,
            "q50 E7 control changed");

    inspect("DEPTH6_HOSTILE_TEN_BODY",
            mask({80,85,88,95,145,168,193,240,252,286}), e7.edges,
            278, mask({8}));
    inspect("FULL_POOL_CLOSEST_TEN_BODY",
            mask({80,85,88,95,143,145,168,193,240,252}), e7.edges,
            275, mask({8}));
}

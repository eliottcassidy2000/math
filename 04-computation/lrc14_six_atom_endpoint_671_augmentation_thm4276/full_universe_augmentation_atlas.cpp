// Exact full-universe failed-body augmentation probe for the two THM-4271
// boundary rows. This is a temporary research executable.

#define CASCADE_LIBRARY_ONLY
#include "04-computation/lrc14_three_round_learned_carrier_thm4266/cascade_pair_exhaustive_primary.cpp"
#undef CASCADE_LIBRARY_ONLY

#include <array>
#include <map>

namespace {

constexpr std::array<u32, 26> FAIL_256 = {{
    0x06166401,0x07067001,0x07106409,0x07126088,0x07126401,0x07162401,
    0x07163400,0x07166008,0x0d107401,0x0d10e401,0x0d146401,0x0d186401,
    0x0d246401,0x0d506401,0x0d906401,0x0f106401,0x0f142401,0x0f142408,
    0x0f143400,0x0f146008,0x15923400,0x17162400,0x17922008,0x1d106401,
    0x1d902401,0x1f142400
}};
constexpr u32 FAIL_384 = 0x0d186401;

struct PairMasses {
    std::vector<i128> mass;
    i128 denominator = 0;
};

PairMasses build_pair_masses(const std::vector<Cell>& cells, i64 q, i64 r) {
    const i64 g = std::gcd(q, r);
    const PrimitivePair primitive = build_primitive(q / g, r / g);
    const AtomData atoms = build_cocycle_atoms(cells, primitive, g);
    PairMasses out;
    out.mass.assign(EXPECTED_REPAIRS, 0);
    u64 operations = 0;
    for (const auto& [mask, value] : atoms.mass) {
        add_supersets_pair(mask, 8 - std::popcount(mask), 0, 0,
                           value, out.mass, operations);
    }
    require(operations == UINT64_C(152170690), "zeta count changed");
    out.denominator = static_cast<i128>(primitive.grid) * g * COMMON;
    return out;
}

}  // namespace

int main() {
    try {
        init_choose8_local();
        const std::vector<Cell> cells = build_pool_cells();
        const PairMasses a = build_pair_masses(cells, 256, 671);
        const PairMasses b = build_pair_masses(cells, 384, 671);
        std::map<u32, u32> least_by_pattern;
        u64 active_a = 0;
        u64 active_b = 0;
        u64 joint_candidates = 0;
        u32 repair = (u32{1} << 8) - 1;
        const u32 limit = u32{1} << 30;
        u64 rank = 0;
        while (repair != 0 && repair < limit) {
            const bool aa = static_cast<i128>(63) * a.mass[rank] >=
                            static_cast<i128>(4) * a.denominator;
            const bool ab = static_cast<i128>(63) * b.mass[rank] >=
                            static_cast<i128>(4) * b.denominator;
            active_a += aa;
            active_b += ab;
            u32 pattern = 0;
            if (aa) {
                for (unsigned i = 0; i < FAIL_256.size(); ++i) {
                    if ((repair & FAIL_256[i]) == 0) pattern |= u32{1} << i;
                }
            }
            if (ab && (repair & FAIL_384) == 0) pattern |= u32{1} << 26;
            if (pattern != 0) {
                ++joint_candidates;
                auto [it, inserted] = least_by_pattern.emplace(pattern, repair);
                if (!inserted && repair < it->second) it->second = repair;
            }
            ++rank;
            const u32 next = next_combination(repair);
            if (next <= repair) break;
            repair = next;
        }
        require(rank == EXPECTED_REPAIRS, "repair enumeration changed");
        require(active_a == 1721339 && active_b == 2444056,
                "active census changed");
        std::cout << "ROUND5_FULL_UNIVERSE_AUGMENTATION_V1\n";
        std::cout << "ACTIVE_256 " << active_a << " ACTIVE_384 " << active_b
                  << " COVERING_CANDIDATES " << joint_candidates
                  << " DISTINCT_PATTERNS " << least_by_pattern.size() << "\n";
        for (const auto& [pattern, least] : least_by_pattern) {
            std::cout << "PATTERN " << std::hex << pattern << " LEAST "
                      << least << std::dec << " COVER "
                      << std::popcount(pattern) << "\n";
        }
        std::cout << "VERDICT PASS FULL_REPAIR_PATTERN_ATLAS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ROUND5_AUGMENTATION_ERROR " << error.what() << '\n';
        return 1;
    }
}

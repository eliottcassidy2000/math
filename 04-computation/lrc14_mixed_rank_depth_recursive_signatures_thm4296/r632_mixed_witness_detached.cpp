// Detached literal-wall certificate for the exact mixed-rank repair number.

#define main r632_detached_hostile_main
#include "r632_detached_hostile_survivor.cpp"
#undef main

namespace {

constexpr std::array<std::size_t, 3> kMixedPacking = {2, 22, 63};
constexpr std::array<u32, 3> kMixedWitness = {
    0x00619324, 0x201813a4, 0x21888126};

template <class Callback>
void mixed_choose9_rec(const std::vector<unsigned char>& positions,
                       std::size_t start, unsigned chosen, u32 mask,
                       Callback& callback) {
    if (chosen == 9) {
        callback(mask);
        return;
    }
    const std::size_t needed = 9 - chosen;
    for (std::size_t index = start;
         index + needed <= positions.size(); ++index)
        mixed_choose9_rec(positions, index + 1, chosen + 1,
                          mask | (u32{1} << positions[index]), callback);
}

template <class Callback>
u64 mixed_enumerate9(u32 forbidden, Callback callback) {
    std::vector<unsigned char> positions;
    for (unsigned bit = 0; bit < 30; ++bit)
        if ((forbidden & (u32{1} << bit)) == 0)
            positions.push_back(static_cast<unsigned char>(bit));
    if (positions.size() < 9) return 0;
    u64 count = 0;
    auto counted = [&](u32 mask) {
        require(std::popcount(mask) == 9 && (mask & forbidden) == 0,
                "rank-nine mixed enumeration escaped");
        callback(mask);
        ++count;
    };
    mixed_choose9_rec(positions, 0, 0, 0, counted);
    return count;
}

}  // namespace


int main(int argc, char** argv) {
    try {
        require(argc == 3,
                "usage: detached-mixed FAILURES RESPONSE_CSV");
        const std::vector<u32> failures = read_failures(argv[1]);
        const Geometry geometry = build_geometry();
        Fnv packing_ledger;
        u64 rank8_checks = 0;
        u64 rank9_checks = 0;
        std::cout << "R632_MIXED89_WITNESS_DETACHED_V1\n"
                  << "PAIR 256,632 FAILURES 66 FNV " << std::hex
                  << kFailuresFnv << std::dec << '\n';
        for (std::size_t left = 0; left < kMixedPacking.size(); ++left) {
            for (std::size_t right = left + 1;
                 right < kMixedPacking.size(); ++right) {
                const std::size_t left_index = kMixedPacking[left];
                const std::size_t right_index = kMixedPacking[right];
                const u32 united = failures[left_index] | failures[right_index];
                Audit rank8;
                const u64 count8 = enumerate_complement(
                    united, [&](u32 mask) { consume(rank8, geometry, mask); });
                Audit rank9;
                const u64 count9 = mixed_enumerate9(
                    united, [&](u32 mask) { consume(rank9, geometry, mask); });
                require(count8 == choose(30 - std::popcount(united), 8) &&
                            count9 == choose(30 - std::popcount(united), 9) &&
                            rank8.active == 0 && rank9.active == 0 &&
                            rank8.maximum.ticks < 0 &&
                            rank9.maximum.ticks < 0,
                        "mixed packing pair has a common active responder");
                rank8_checks += count8;
                rank9_checks += count9;
                packing_ledger.add(left_index);
                packing_ledger.add(right_index);
                packing_ledger.add(failures[left_index]);
                packing_ledger.add(failures[right_index]);
                packing_ledger.add(count8);
                add_i128(packing_ledger, rank8.maximum.ticks);
                packing_ledger.add(rank8.least_maximum);
                packing_ledger.add(count9);
                add_i128(packing_ledger, rank9.maximum.ticks);
                packing_ledger.add(rank9.least_maximum);
                std::cout << "PACKING_PAIR " << left_index << ','
                          << right_index << " BODIES "
                          << hex8(failures[left_index]) << ','
                          << hex8(failures[right_index])
                          << " RANK8_CANDIDATES " << count8
                          << " RANK8_MAX_MARGIN "
                          << decimal(rank8.maximum.ticks) << " RANK8_MASK "
                          << hex8(rank8.least_maximum)
                          << " RANK9_CANDIDATES " << count9
                          << " RANK9_MAX_MARGIN "
                          << decimal(rank9.maximum.ticks) << " RANK9_MASK "
                          << hex8(rank9.least_maximum) << '\n';
            }
        }
        std::ofstream responses(argv[2]);
        require(static_cast<bool>(responses), "cannot create mixed responses");
        responses << "mask_hex,rank,body_ordinal,body_hex,margin_ticks63,grid\n";
        std::vector<unsigned char> covered(failures.size(), 0);
        Fnv witness_ledger;
        Fnv response_ledger;
        u64 incidences = 0;
        for (u32 mask : kMixedWitness) {
            require(std::popcount(mask) == 9, "mixed witness rank changed");
            const Margin value = margin(geometry, mask);
            require(value.ticks > 0, "mixed witness is not strictly active");
            witness_ledger.add(mask);
            witness_ledger.add(static_cast<u64>(value.mass));
            add_i128(witness_ledger, value.ticks);
            u64 hits = 0;
            u64 packing_hits = 0;
            for (std::size_t index = 0; index < failures.size(); ++index) {
                if ((mask & failures[index]) != 0) continue;
                covered[index] = 1;
                ++hits;
                ++incidences;
                packing_hits +=
                    std::find(kMixedPacking.begin(), kMixedPacking.end(), index)
                    != kMixedPacking.end();
                response_ledger.add(mask);
                response_ledger.add(index);
                response_ledger.add(failures[index]);
                responses << hex8(mask) << ",9," << index << ','
                          << hex8(failures[index]) << ','
                          << decimal(value.ticks) << ',' << geometry.grid
                          << '\n';
            }
            require(packing_hits == 1,
                    "mixed witness pays wrong packing multiplicity");
            std::cout << "WITNESS MASK " << hex8(mask) << " RANK 9 LABELS {"
                      << labels(mask) << "} MASS " << value.mass
                      << " MARGIN_TICKS63 " << decimal(value.ticks)
                      << " HITS " << hits << " PACKING_HITS "
                      << packing_hits << '\n';
        }
        require(responses.good(), "mixed response write failed");
        require(std::accumulate(covered.begin(), covered.end(), u64{0}) ==
                    failures.size(),
                "mixed witness does not cover all failures");
        std::cout << "PACKING SIZE 3 RANK8_CHECKS " << rank8_checks
                  << " RANK9_CHECKS " << rank9_checks << " LEDGER_FNV "
                  << std::hex << packing_ledger.state << std::dec << '\n'
                  << "MIXED_MINIMUM 3 LOWER three-obligation-mixed-activity-"
                     "packing UPPER explicit-three-rank9-cover\n"
                  << "COVERED 66 INCIDENCES " << incidences
                  << " WITNESS_FNV " << std::hex << witness_ledger.state
                  << " RESPONSE_FNV " << response_ledger.state << std::dec
                  << '\n'
                  << "SCOPE IMPORT_FREE_LITERAL_WALL_FIXED_PAIR_LABELLED_"
                     "RANK8_OR_RANK9_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_MIXED89_MINIMUM_THREE\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "R632_MIXED_DETACHED_ERROR " << error.what() << '\n';
        return 1;
    }
}

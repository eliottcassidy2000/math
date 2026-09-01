// Independent complement-subset audit of the unique endpoint-593 failure's
// complete rank-eight/rank-nine response family.

#define ENDPOINT595_REPAIRED_RAW_MAIN endpoint593_complement_hidden_main
#include "04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4309/strict_repair_raw_audit.cpp"
#undef ENDPOINT595_REPAIRED_RAW_MAIN

int main() {
    try {
        constexpr u32 body = UINT32_C(0x34087401);
        constexpr u32 universe = (u32{1} << 30) - 1;
        const u32 complement = universe ^ body;
        require(std::popcount(body) == 9 && std::popcount(complement) == 21,
                "failure/complement rank changed");
        const Geometry geometry = build_geometry(96, 593);
        std::array<u64, 2> events{}, responders{};
        std::array<std::vector<u32>, 2> response_masks;
        std::array<Fnv, 2> ledgers;
        std::array<u32, 2> least{UINT32_MAX, UINT32_MAX};
        for (u32 mask = complement;; mask = (mask - 1) & complement) {
            const int rank = std::popcount(mask);
            if (rank == 8 || rank == 9) {
                const std::size_t slot = rank - 8;
                ++events[slot];
                if (margin(geometry, mask).ticks >= 0) {
                    ++responders[slot];
                    response_masks[slot].push_back(mask);
                    least[slot] = std::min(least[slot], mask);
                }
            }
            if (mask == 0) break;
        }
        // The complement-subset walk is descending.  Canonical response
        // identities are hashed in increasing raw-mask order, independently
        // of the enumeration used to discover them.
        for (auto& masks : response_masks) {
            std::sort(masks.begin(), masks.end());
        }
        for (std::size_t slot = 0; slot < response_masks.size(); ++slot) {
            for (u32 mask : response_masks[slot]) ledgers[slot].add(mask);
        }
        require(events[0] == 203490 && events[1] == 293930,
                "complement subset counts changed");
        require(responders[0] == 1636 &&
                    ledgers[0].state == UINT64_C(0x56f82f5dc11db83b) &&
                    least[0] == UINT32_C(0x0134012c),
                "rank-eight response identity changed");
        require(responders[1] == 16209 &&
                    ledgers[1].state == UINT64_C(0x0f615d806860553f) &&
                    least[1] == UINT32_C(0x0036092c),
                "rank-nine response identity changed");
        std::cout << "LRC14_ENDPOINT593_COMPLEMENT_RESPONSE_AUDIT_V1\n"
                  << "FAILURE 96,593,34087401 COMPLEMENT_RANK 21\n";
        for (int rank = 8; rank <= 9; ++rank) {
            const std::size_t slot = rank - 8;
            std::cout << "RANK " << rank << " COMPLEMENT_SUBSETS "
                      << events[slot] << " ACTIVE_RESPONDERS "
                      << responders[slot] << " FNV " << std::hex
                      << ledgers[slot].state << " LEAST " << hex8(least[slot])
                      << std::dec << '\n';
        }
        std::cout << "SCOPE INDEPENDENT_COMPLEMENT_GENERATION_FIXED_SINGLE_"
                     "FAILURE_NO_PHYSICAL_ENTRY_NO_LRC14\nVERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT593_COMPLEMENT_ERROR " << error.what() << '\n';
        return 1;
    }
}

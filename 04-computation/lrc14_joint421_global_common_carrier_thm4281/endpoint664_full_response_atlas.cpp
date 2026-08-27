// Complete labelled response quotient for the two failed bodies at (256,664).

#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

namespace {

constexpr std::array<u32, 2> FAILURES_664 = {
    UINT32_C(0x0d143408), UINT32_C(0x2d143400)};

}  // namespace

int main() {
    try {
        std::cout.setf(std::ios::unitbuf);
        init_choose8_local();
        FnvLocal failure_ledger;
        u32 union_mask = 0;
        for (u32 body : FAILURES_664) {
            require(std::popcount(body) == 9, "failed body arity changed");
            failure_ledger.add(body);
            union_mask |= body;
        }
        require(failure_ledger.state == UINT64_C(0xe9872402af5e2795) &&
                    union_mask == UINT32_C(0x2d143408) &&
                    std::popcount(union_mask) == 10,
                "endpoint664 failure universe changed");

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        const ActiveUniverse active = build_active_universe(cells, 256, 664);
        std::vector<unsigned char> response(EXPECTED_REPAIRS, 0);
        u64 incidences = 0;
        for (std::size_t obligation = 0; obligation < FAILURES_664.size();
             ++obligation) {
            enumerate_disjoint_repairs(
                FAILURES_664[obligation], [&](u32, u64 rank) {
                    if (!active.active[rank]) return;
                    response[rank] |=
                        static_cast<unsigned char>(1u << obligation);
                    ++incidences;
                });
        }

        std::array<u64, 4> multiplicity{};
        std::array<u32, 4> least{};
        FnvLocal response_ledger;
        u64 nonempty = 0;
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            const unsigned pattern = response[rank];
            ++multiplicity[pattern];
            const u32 repair = unrank_colex8(rank);
            if (least[pattern] == 0 || repair < least[pattern])
                least[pattern] = repair;
            if (pattern != 0) {
                ++nonempty;
                response_ledger.add(repair);
                response_ledger.add(pattern);
            }
        }
        require(std::accumulate(multiplicity.begin(), multiplicity.end(),
                                u64{0}) == EXPECTED_REPAIRS,
                "response census changed");
        std::vector<u32> witnesses;
        if (multiplicity[3] != 0) {
            witnesses.push_back(least[3]);
        } else {
            require(multiplicity[1] != 0 && multiplicity[2] != 0,
                    "two-obligation quotient is not coverable");
            witnesses = {least[1], least[2]};
            std::sort(witnesses.begin(), witnesses.end());
        }
        FnvLocal witness_ledger;
        unsigned union_response = 0;
        for (u32 witness : witnesses) {
            witness_ledger.add(witness);
            for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
                if (unrank_colex8(rank) == witness) {
                    union_response |= response[rank];
                    break;
                }
            }
        }
        require(union_response == 3, "selected witnesses miss a failure");

        FnvLocal class_ledger;
        u64 classes = 0;
        for (unsigned pattern = 1; pattern < 4; ++pattern) {
            if (multiplicity[pattern] == 0) continue;
            ++classes;
            class_ledger.add(pattern);
            class_ledger.add(multiplicity[pattern]);
            class_ledger.add(least[pattern]);
        }
        require(active.count == UINT64_C(1637771) &&
                    active.fnv == UINT64_C(0x1f612f6504d891bd) &&
                    incidences == UINT64_C(1461) && nonempty == 1354 &&
                    multiplicity == std::array<u64, 4>{
                        UINT64_C(5851571), UINT64_C(254),
                        UINT64_C(993), UINT64_C(107)} &&
                    response_ledger.state == UINT64_C(0x74100225fdfc1b26) &&
                    least[1] == UINT32_C(0x20884283) &&
                    least[2] == UINT32_C(0x0002ca89) &&
                    least[3] == UINT32_C(0x00c0c125) && classes == 3 &&
                    class_ledger.state == UINT64_C(0xbf65c791751a6556) &&
                    witnesses == std::vector<u32>{UINT32_C(0x00c0c125)} &&
                    witness_ledger.state == UINT64_C(0x80ed30f8cee1cf5b),
                "frozen endpoint664 response quotient changed");
        std::cout << "ENDPOINT664_FULL_RESPONSE_ATLAS_V1\n"
                  << "PAIR 256,664 FAILURES 2 FAILURE_FNV " << std::hex
                  << failure_ledger.state << " UNION " << union_mask
                  << std::dec << " UNION_SIZE 10\n"
                  << "REPAIRS " << EXPECTED_REPAIRS << " ACTIVE "
                  << active.count << " ACTIVE_FNV " << std::hex << active.fnv
                  << std::dec << " ZETA_OPERATIONS "
                  << active.zeta_operations << '\n'
                  << "COMPLEMENT_CHECKS "
                  << FAILURES_664.size() * DISJOINT_REPAIRS_PER_BODY
                  << " ACTIVE_INCIDENCES " << incidences
                  << " NONEMPTY_REPAIRS " << nonempty
                  << " EMPTY_REPAIRS " << multiplicity[0]
                  << " RESPONSE_FNV " << std::hex << response_ledger.state
                  << std::dec << '\n';
        for (unsigned pattern = 1; pattern < 4; ++pattern) {
            if (multiplicity[pattern] == 0) continue;
            std::cout << "PATTERN " << std::hex << pattern << " LEAST "
                      << least[pattern] << std::dec << " MULTIPLICITY "
                      << multiplicity[pattern] << " COVER "
                      << std::popcount(pattern) << '\n';
        }
        std::cout << "NONEMPTY_CLASSES " << classes << " CLASS_FNV "
                  << std::hex << class_ledger.state << std::dec << '\n'
                  << "MINIMUM_AUGMENTATION " << witnesses.size()
                  << " WITNESSES ";
        for (std::size_t index = 0; index < witnesses.size(); ++index) {
            if (index != 0) std::cout << ',';
            std::cout << std::hex << witnesses[index] << std::dec;
        }
        std::cout << " WITNESS_FNV " << std::hex << witness_ledger.state
                  << std::dec << " ZERO_IMPOSSIBLE INHERITED_FAILURES_2"
                  << (witnesses.size() == 1
                          ? " ONE_WITNESS_FULL_RESPONSE\n"
                          : " ONE_IMPOSSIBLE_NO_FULL_RESPONSE\n")
                  << "VERDICT PASS EXACT_COMPLETE_RESPONSE_QUOTIENT_AND_MINIMUM\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT664_ATLAS_ERROR " << error.what() << '\n';
        return 1;
    }
}

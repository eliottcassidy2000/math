// Complete single-obligation response census and strict-inactive capacity for
// the endpoint-593 failure of the THM-4309 3,925-mask carrier.

#define ENDPOINT595_REPAIRED_RAW_MAIN endpoint593_response_hidden_main
#include "04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4309/strict_repair_raw_audit.cpp"
#undef ENDPOINT595_REPAIRED_RAW_MAIN

namespace {
constexpr u64 kFinalDeleteFnv593 = UINT64_C(0xff4c932f9a7adac8);
constexpr u64 kFinalCarrierFnv593 = UINT64_C(0x6fbd0bffcf0ed78b);
constexpr u64 kTop593Fnv593 = UINT64_C(0x5424c07fa724011f);
constexpr u64 kFull432Fnv593 = UINT64_C(0xa7ed492c64d1c0d8);
constexpr u32 kFailureBody593 = UINT32_C(0x34087401);

std::vector<u32> read_flexible_masks593(const std::filesystem::path& path,
                                        std::size_t expected,
                                        u64 expected_fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open mask ledger");
    std::vector<u32> masks;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask_agent(token);
        require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                    distinct.insert(mask).second,
                "mask ledger rank/distinctness changed");
        masks.push_back(mask);
    }
    require(input.eof() && masks.size() == expected &&
                masks_fnv_agent(masks) == expected_fnv,
            "mask ledger identity changed");
    return masks;
}

std::vector<PairAgent> read_pair_audit_rows593(
    const std::filesystem::path& path, std::size_t expected) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open pair-audit ledger");
    std::string line;
    require(std::getline(input, line) &&
                line.starts_with("q,r,active,active_fnv,"),
            "pair-audit header changed");
    std::vector<PairAgent> rows;
    std::set<std::pair<int, int>> distinct;
    while (std::getline(input, line)) {
        require(!line.empty(), "blank pair-audit row");
        const std::size_t first = line.find(',');
        const std::size_t second = line.find(',', first + 1);
        require(first != std::string::npos && second != std::string::npos,
                "malformed pair-audit row");
        const int q = std::stoi(line.substr(0, first));
        const int r = std::stoi(line.substr(first + 1, second - first - 1));
        require(q > 0 && q < r && distinct.insert({q, r}).second,
                "pair-audit row invalid/duplicate");
        rows.push_back(PairAgent{q, r});
    }
    require(input.eof() && rows.size() == expected,
            "pair-audit row count changed");
    return rows;
}

std::vector<u32> reconstruct_c3925_593(
    const std::filesystem::path& base,
    const std::filesystem::path& add45,
    const std::filesystem::path& suffix9,
    const std::filesystem::path& repairs76,
    const std::filesystem::path& additions4,
    const std::filesystem::path& deletions73,
    const std::filesystem::path& additions10,
    const std::filesystem::path& final_delete) {
    const auto repairs = read_mixed617(repairs76, 76, kRepairFnv);
    const auto old_additions = read_mixed617(additions4, 4, kOldAdditionFnv);
    const auto old_deletions = read_mixed617(deletions73, 73, kOldDeleteFnv);
    const auto new_additions = read_mixed617(additions10, 10, kNewAdditionFnv);
    const auto final_deletions = read_flexible_masks593(
        final_delete, 5104, kFinalDeleteFnv593);
    std::vector<u32> augmented = build_mixed_carrier(base, add45, suffix9);
    std::set<u32> distinct(augmented.begin(), augmented.end());
    for (u32 mask : repairs) {
        require(distinct.insert(mask).second, "repair overlap");
        augmented.push_back(mask);
    }
    require(augmented.size() == 9088 &&
                masks_fnv_agent(augmented) == kAugmentedFnv,
            "augmented carrier changed");
    const std::set<u32> old_deleted(old_deletions.begin(), old_deletions.end());
    std::vector<u32> c596;
    for (u32 mask : augmented)
        if (!old_deleted.contains(mask)) c596.push_back(mask);
    for (u32 mask : old_additions) c596.push_back(mask);
    require(c596.size() == 9019 && masks_fnv_agent(c596) == kC596Fnv,
            "C596 identity changed");
    const std::set<u32> deleted(final_deletions.begin(), final_deletions.end());
    std::vector<u32> carrier;
    for (u32 mask : c596)
        if (!deleted.contains(mask)) carrier.push_back(mask);
    distinct = std::set<u32>(carrier.begin(), carrier.end());
    for (u32 mask : new_additions) {
        require(distinct.insert(mask).second, "addition10 overlap");
        carrier.push_back(mask);
    }
    require(carrier.size() == 3925 &&
                masks_fnv_agent(carrier) == kFinalCarrierFnv593,
            "C3925 identity changed");
    return carrier;
}
}  // namespace

#ifndef ENDPOINT593_RESPONSE_CAPACITY_MAIN
#define ENDPOINT593_RESPONSE_CAPACITY_MAIN main
#endif

int ENDPOINT593_RESPONSE_CAPACITY_MAIN(int argc, char** argv) {
    try {
        require(argc == 17,
                "usage: response JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE "
                "PREFIX391_CSV TOP594_CSV TOP593 ADDITION_OUT INACTIVE_OUT");
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> carrier = reconstruct_c3925_593(
            argv[2], argv[3], argv[4], argv[7], argv[8], argv[9], argv[10],
            argv[11]);
        (void)read_pairs_repaired(argv[5], 22647, kUniverseFnv);
        (void)read_pairs_repaired(argv[6], 1624, kOldUnionFnv);
        const auto prefix391 = read_pair_audit_rows593(argv[12], 391);
        const auto top594 = read_pair_audit_rows593(argv[13], 25);
        const auto top593_set = read_pairs_repaired(argv[14], 16, kTop593Fnv593);
        std::set<std::pair<int, int>> target_set;
        std::set<std::pair<int, int>> prior_set;
        for (PairAgent pair : prefix391)
            require(target_set.insert({pair.q, pair.r}).second &&
                        prior_set.insert({pair.q, pair.r}).second,
                    "prefix row overlap");
        for (PairAgent pair : top594)
            require(target_set.insert({pair.q, pair.r}).second &&
                        prior_set.insert({pair.q, pair.r}).second,
                    "top594 row overlap");
        for (const auto& pair : top593_set)
            require(target_set.insert(pair).second, "top593 row overlap");
        require(target_set.size() == 432, "full target row count changed");
        require(prior_set.size() == 416, "prior target row count changed");
        std::vector<Geometry> geometry;
        std::vector<bool> is_prior;
        for (const auto& [q, r] : target_set) {
            geometry.push_back(build_geometry(q, r));
            is_prior.push_back(prior_set.contains({q, r}));
        }
        Fnv target_ledger;
        for (const auto& [q, r] : target_set) {
            target_ledger.add(q);
            target_ledger.add(r);
        }
        require(target_ledger.state == kFull432Fnv593,
                "full target identity changed");

        const Geometry failure_geometry = build_geometry(96, 593);
        std::array<u64, 2> all_count{}, active_count{}, response_count{};
        std::array<Fnv, 2> response_ledger;
        std::array<u32, 2> least_response{UINT32_MAX, UINT32_MAX};
        const u32 limit = u32{1} << 30;
        for (int rank = 8; rank <= 9; ++rank) {
            const std::size_t slot = rank - 8;
            for (u32 mask = (u32{1} << rank) - 1; mask < limit;
                 mask = next_combination(mask)) {
                ++all_count[slot];
                if (margin(failure_geometry, mask).ticks < 0) continue;
                ++active_count[slot];
                if ((mask & kFailureBody593) != 0) continue;
                ++response_count[slot];
                response_ledger[slot].add(mask);
                least_response[slot] = std::min(least_response[slot], mask);
                require(std::find(carrier.begin(), carrier.end(), mask) ==
                            carrier.end(),
                        "active response already in failing carrier");
            }
        }
        require(response_count[0] + response_count[1] > 0,
                "single failure has no rank8/9 response");
        const u32 addition = std::min(least_response[0], least_response[1]);
        require((addition & kFailureBody593) == 0 &&
                    margin(failure_geometry, addition).ticks >= 0,
                "selected response invalid");

        std::vector<u32> inactive_prefix416;
        std::vector<u32> inactive_full432;
        Fnv sign_ledger;
        u64 equality_cells = 0;
        for (u32 mask : carrier) {
            bool active416 = false;
            bool active432 = false;
            for (std::size_t index = 0; index < geometry.size(); ++index) {
                const Margin exact = margin(geometry[index], mask);
                equality_cells += exact.ticks == 0;
                const bool active = exact.ticks >= 0;
                sign_ledger.add(mask);
                sign_ledger.add(active);
                if (is_prior[index]) active416 |= active;
                active432 |= active;
            }
            if (!active416) inactive_prefix416.push_back(mask);
            if (!active432) inactive_full432.push_back(mask);
        }

        std::ofstream addition_out(argv[15]);
        std::ofstream inactive_out(argv[16]);
        require(static_cast<bool>(addition_out) &&
                    static_cast<bool>(inactive_out),
                "cannot create response/capacity ledgers");
        addition_out << hex8(addition) << '\n';
        for (u32 mask : inactive_full432) inactive_out << hex8(mask) << '\n';
        require(addition_out.good() && inactive_out.good(),
                "response/capacity ledger write failed");

        std::cout << "LRC14_ENDPOINT593_SINGLE_RESPONSE_CAPACITY_V1\n"
                  << "C3925 " << carrier.size() << " FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec << " TARGET_ROWS "
                  << target_set.size() << " TARGET_FNV " << std::hex
                  << target_ledger.state << std::dec << '\n'
                  << "FAILURE 96,593," << hex8(kFailureBody593)
                  << " OBLIGATIONS 1 EXACT_RESPONSE_MIN 1\n";
        for (int rank = 8; rank <= 9; ++rank) {
            const std::size_t slot = rank - 8;
            std::cout << "RANK " << rank << " ALL " << all_count[slot]
                      << " ACTIVE " << active_count[slot] << " RESPONDERS "
                      << response_count[slot] << " RESPONDER_FNV " << std::hex
                      << response_ledger[slot].state << " LEAST "
                      << hex8(least_response[slot]) << std::dec << '\n';
        }
        std::cout << "SELECTED_ADDITION " << hex8(addition) << " RANK "
                  << std::popcount(addition) << " ADDITION_FNV " << std::hex
                  << masks_fnv_agent(std::vector<u32>{addition}) << std::dec
                  << '\n'
                  << "STRICT_INACTIVE_PREFIX416 " << inactive_prefix416.size()
                  << " FNV " << std::hex
                  << masks_fnv_agent(inactive_prefix416) << std::dec << '\n'
                  << "STRICT_INACTIVE_FULL432 " << inactive_full432.size()
                  << " FNV " << std::hex
                  << masks_fnv_agent(inactive_full432) << std::dec << '\n'
                  << "SIGN_CELLS " << carrier.size() * target_set.size()
                  << " SIGN_FNV " << std::hex << sign_ledger.state << std::dec
                  << " EQUALITIES " << equality_cells << '\n'
                  << "SCOPE COMPLETE_RANK8_RANK9_SINGLE_FAILURE_RESPONSE_"
                     "AND_STRICT_INACTIVITY_FIXED_432_ROWS_NO_PHYSICAL_ENTRY_"
                     "NO_LRC14\nVERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT593_RESPONSE_CAPACITY_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}

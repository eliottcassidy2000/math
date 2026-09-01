// Generic full-391-row raw replay for a final active-deletion ledger.  This
// deliberately proves simultaneous safety by replay rather than inferring it
// from individual redundancy.

#define ENDPOINT595_REPAIRED_RAW_MAIN endpoint595_baseline_hidden_main
#include "strict_repair_raw_audit.cpp"
#undef ENDPOINT595_REPAIRED_RAW_MAIN

namespace {
std::vector<u32> read_flexible_deletions(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open final deletion ledger");
    std::vector<u32> masks;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask_agent(token);
        require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                    distinct.insert(mask).second,
                "final deletion invalid/duplicate");
        masks.push_back(mask);
    }
    require(input.eof() && !masks.empty(), "empty final deletion ledger");
    return masks;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 15,
                "usage: full_raw JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION REPAIRS76 OLD_ADD4 OLD_DELETE73 NEW_ADD10 "
                "FINAL_DELETE PAIR_CSV FAILURE_CSV WORKERS");
        const unsigned worker_count = std::stoul(argv[14]);
        require(worker_count >= 1 && worker_count <= 64, "workers escaped 1..64");
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const auto repairs = read_mixed617(argv[7], 76, kRepairFnv);
        const auto old_additions = read_mixed617(argv[8], 4, kOldAdditionFnv);
        const auto old_deletions = read_mixed617(argv[9], 73, kOldDeleteFnv);
        const auto new_additions = read_mixed617(argv[10], 10, kNewAdditionFnv);
        const auto final_deletions = read_flexible_deletions(argv[11]);
        Fnv deletion_ledger;
        for (u32 mask : final_deletions) deletion_ledger.add(mask);

        std::vector<u32> augmented = build_mixed_carrier(argv[2], argv[3], argv[4]);
        std::set<u32> distinct(augmented.begin(), augmented.end());
        for (u32 mask : repairs) {
            require(distinct.insert(mask).second, "repair overlap");
            augmented.push_back(mask);
        }
        require(augmented.size() == 9088 && masks_fnv_agent(augmented) == kAugmentedFnv,
                "augmented carrier changed");
        const std::set<u32> old_deleted(old_deletions.begin(), old_deletions.end());
        std::vector<u32> c596;
        for (u32 mask : augmented)
            if (!old_deleted.contains(mask)) c596.push_back(mask);
        for (u32 mask : old_additions) c596.push_back(mask);
        require(c596.size() == 9019 && masks_fnv_agent(c596) == kC596Fnv,
                "C596 carrier changed");
        const std::set<u32> final_deleted(final_deletions.begin(),
                                          final_deletions.end());
        for (u32 mask : final_deletions)
            require(std::find(c596.begin(), c596.end(), mask) != c596.end() &&
                        !joint_set.contains(mask),
                    "final deletion absent/joint");
        std::vector<u32> carrier;
        for (u32 mask : c596)
            if (!final_deleted.contains(mask)) carrier.push_back(mask);
        distinct = std::set<u32>(carrier.begin(), carrier.end());
        for (u32 mask : new_additions) {
            require(distinct.insert(mask).second, "new addition overlap");
            carrier.push_back(mask);
        }
        const u64 rank8 = std::count_if(carrier.begin(), carrier.end(),
                                       [](u32 mask) { return std::popcount(mask) == 8; });

        const auto universe = read_pairs_repaired(argv[5], 22647, kUniverseFnv);
        const auto old_union = read_pairs_repaired(argv[6], 1624, kOldUnionFnv);
        std::vector<PairAgent> band;
        std::set<PairKey> band_set;
        for (const auto& [q, r] : universe)
            if (r >= 596 || (r == 595 && !old_union.contains({q, r}))) {
                band.push_back(PairAgent{q, r});
                band_set.insert({q, r});
            }
        require(band.size() == 391, "full preservation band changed");

        std::vector<PairAudit617> audits(band.size());
        std::atomic<std::size_t> cursor{0};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        auto worker = [&]() {
            try {
                while (true) {
                    const std::size_t index = cursor.fetch_add(1);
                    if (index >= band.size()) return;
                    audits[index] = audit_pair617(band[index], joint, carrier, joint_set);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
        };
        std::vector<std::thread> workers;
        for (unsigned index = 0; index < worker_count; ++index)
            workers.emplace_back(worker);
        for (auto& thread : workers) thread.join();
        if (worker_error) std::rethrow_exception(worker_error);

        std::ofstream pair_out(argv[12]);
        std::ofstream failure_out(argv[13]);
        require(static_cast<bool>(pair_out) && static_cast<bool>(failure_out),
                "cannot create full raw ledgers");
        pair_out << "q,r,active,active_fnv,active_joint,active_nonjoint,exposed,"
                    "exposed_fnv,minimum_hits,maximum_hits,failures,failure_fnv\n";
        failure_out << "q,r,body_hex\n";
        Fnv pair_ledger;
        u64 total_exposed = 0, total_hits = 0, total_failures = 0;
        for (std::size_t index = 0; index < band.size(); ++index) {
            const auto pair = band[index];
            const auto& audit = audits[index];
            total_exposed += audit.exposed;
            total_hits += audit.hit_incidences;
            total_failures += audit.failures;
            pair_ledger.add(pair.q); pair_ledger.add(pair.r);
            pair_ledger.add(audit.active); pair_ledger.add(audit.active_fnv);
            pair_ledger.add(audit.active_joint); pair_ledger.add(audit.active_nonjoint);
            pair_ledger.add(audit.exposed); pair_ledger.add(audit.exposed_fnv);
            pair_ledger.add(audit.minimum_hits); pair_ledger.add(audit.maximum_hits);
            pair_ledger.add(audit.failures); pair_ledger.add(audit.failure_fnv);
            pair_out << pair.q << ',' << pair.r << ',' << audit.active << ','
                     << std::hex << audit.active_fnv << std::dec << ','
                     << audit.active_joint << ',' << audit.active_nonjoint << ','
                     << audit.exposed << ',' << std::hex << audit.exposed_fnv
                     << std::dec << ',' << audit.minimum_hits << ','
                     << audit.maximum_hits << ',' << audit.failures << ','
                     << std::hex << audit.failure_fnv << std::dec << '\n';
            for (u32 body : audit.failure_bodies)
                failure_out << pair.q << ',' << pair.r << ',' << hex8(body) << '\n';
        }
        require(pair_out.good() && failure_out.good(),
                "full-prefix final ledgers failed to write");
        std::cout << "LRC14_ENDPOINT595_FINAL_FULL_PREFIX_RAW_AUDIT_V1\n"
                  << "C596 9019 FNV " << std::hex << kC596Fnv << std::dec
                  << " DELETE " << final_deletions.size() << " DELETE_FNV "
                  << std::hex << deletion_ledger.state << " ADD10_FNV "
                  << kNewAdditionFnv << std::dec << '\n'
                  << "FINAL_CARRIER " << carrier.size() << " FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec << " RANK8 " << rank8
                  << " RANK9 " << carrier.size() - rank8 << '\n'
                  << "ROWS " << band.size() << " ROW_FNV " << std::hex
                  << pair_fnv_repaired(band_set) << std::dec << " WORKERS "
                  << worker_count << " BODY_UNIVERSE_PER_ROW " << kBodyCount617
                  << " TOTAL_BODY_TESTS " << kBodyCount617 * band.size() << '\n'
                  << "SUMMARY TOTAL_EXPOSED " << total_exposed
                  << " TOTAL_HIT_INCIDENCES " << total_hits << " FAILURES "
                  << total_failures << " PAIR_LEDGER_FNV " << std::hex
                  << pair_ledger.state << std::dec << '\n'
                  << "SCOPE DIRECT_SIMULTANEOUS_FULL_PREFIX_RAW_REPLAY_"
                     "NO_INDIVIDUAL_REDUNDANCY_INFERENCE_NO_PHYSICAL_ENTRY_"
                     "NO_LRC14\nVERDICT "
                  << (total_failures == 0 ? "PASS" : "FAIL") << '\n';
        return total_failures == 0 ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "FINAL_FULL_PREFIX_RAW_ERROR " << error.what() << '\n';
        return 1;
    }
}

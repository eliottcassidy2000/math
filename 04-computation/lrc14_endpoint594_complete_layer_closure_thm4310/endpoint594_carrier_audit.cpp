// Direct raw replay of the THM-4309 3,925-mask carrier on the exact
// endpoint-594 layer, promoted for THM-4310.

#define ENDPOINT595_REPAIRED_RAW_MAIN endpoint594_strict_hidden_main
#include "04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4309/strict_repair_raw_audit.cpp"
#undef ENDPOINT595_REPAIRED_RAW_MAIN

namespace {
constexpr u64 kTop594Fnv = UINT64_C(0xcce015c81f7121d9);
constexpr u64 kFinalDeleteFnvScout = UINT64_C(0xff4c932f9a7adac8);
constexpr u64 kFinalCarrierFnvScout = UINT64_C(0x6fbd0bffcf0ed78b);

std::vector<u32> read_final_deletions(const std::filesystem::path& path) {
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

std::vector<PairAgent> read_target_rows(const std::filesystem::path& path) {
    const auto row_set = read_pairs_repaired(path, 25, kTop594Fnv);
    std::vector<PairAgent> rows;
    for (const auto& [q, r] : row_set) {
        require(r == 594, "target escaped endpoint 594");
        rows.push_back(PairAgent{q, r});
    }
    return rows;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 16 || argc == 17,
                "usage: audit JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION REPAIRS76 OLD_ADD4 OLD_DELETE73 ADD10 "
                "FINAL_DELETE TARGET25 PAIR_CSV FAILURE_CSV WORKERS "
                "[HOSTILE_DELETE_ONE]");
        const unsigned worker_count = std::stoul(argv[15]);
        require(worker_count >= 1 && worker_count <= 64,
                "workers escaped 1..64");

        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const auto repairs = read_mixed617(argv[7], 76, kRepairFnv);
        const auto old_additions = read_mixed617(argv[8], 4, kOldAdditionFnv);
        const auto old_deletions = read_mixed617(argv[9], 73, kOldDeleteFnv);
        const auto new_additions = read_mixed617(argv[10], 10, kNewAdditionFnv);
        const auto final_deletions = read_final_deletions(argv[11]);

        Fnv deletion_ledger;
        for (u32 mask : final_deletions) deletion_ledger.add(mask);
        require(final_deletions.size() == 5104 &&
                    deletion_ledger.state == kFinalDeleteFnvScout,
                "final deletion identity changed");

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
        const u64 base_rank8 = std::count_if(
            carrier.begin(), carrier.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        require(carrier.size() == 3925 && base_rank8 == 3858 &&
                    masks_fnv_agent(carrier) == kFinalCarrierFnvScout,
                "final carrier identity changed");
        u32 hostile_delete = 0;
        if (argc == 17) {
            hostile_delete = parse_mask_agent(argv[16]);
            require(!joint_set.contains(hostile_delete),
                    "hostile deletion entered joint deck");
            const auto old_size = carrier.size();
            std::erase(carrier, hostile_delete);
            require(carrier.size() + 1 == old_size,
                    "hostile deletion absent/nonunique");
        }
        const u64 rank8 = std::count_if(carrier.begin(), carrier.end(),
                                       [](u32 mask) { return std::popcount(mask) == 8; });

        // Read the complete inherited universe as an identity/control even
        // while the exact 25-row carrier target is supplied separately.
        (void)read_pairs_repaired(argv[5], 22647, kUniverseFnv);
        (void)read_pairs_repaired(argv[6], 1624, kOldUnionFnv);
        const auto rows = read_target_rows(argv[12]);

        std::vector<PairAudit617> audits(rows.size());
        std::atomic<std::size_t> cursor{0};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        auto worker = [&]() {
            try {
                while (true) {
                    const std::size_t index = cursor.fetch_add(1);
                    if (index >= rows.size()) return;
                    audits[index] = audit_pair617(rows[index], joint, carrier, joint_set);
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

        std::ofstream pair_out(argv[13]);
        std::ofstream failure_out(argv[14]);
        require(static_cast<bool>(pair_out) && static_cast<bool>(failure_out),
                "cannot create audit ledgers");
        pair_out << "q,r,active,active_fnv,active_joint,active_nonjoint,exposed,"
                    "exposed_fnv,minimum_hits,maximum_hits,failures,failure_fnv\n";
        failure_out << "q,r,body_hex\n";
        Fnv pair_ledger, obligation_ledger;
        u64 total_exposed = 0, total_hits = 0, total_failures = 0;
        u64 failed_rows = 0;
        for (std::size_t index = 0; index < rows.size(); ++index) {
            const auto pair = rows[index];
            const auto& audit = audits[index];
            total_exposed += audit.exposed;
            total_hits += audit.hit_incidences;
            total_failures += audit.failures;
            failed_rows += audit.failures != 0;
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
            for (u32 body : audit.failure_bodies) {
                failure_out << pair.q << ',' << pair.r << ',' << hex8(body) << '\n';
                obligation_ledger.add(pair.q);
                obligation_ledger.add(pair.r);
                obligation_ledger.add(body);
            }
        }
        require(pair_out.good() && failure_out.good(),
                "audit ledgers failed to write");

        std::cout << "LRC14_ENDPOINT594_C3925_RAW_AUDIT_V1\n"
                  << "C596 9019 FNV " << std::hex << kC596Fnv << std::dec
                  << " DELETE " << final_deletions.size() << " DELETE_FNV "
                  << std::hex << deletion_ledger.state << " ADD10_FNV "
                  << kNewAdditionFnv << std::dec << '\n'
                  << "BASE_C595 3925 FNV " << std::hex << kFinalCarrierFnvScout
                  << std::dec << " RANK8 " << base_rank8 << " RANK9 67 "
                     "JOINT_RETAINED 421\n"
                  << "AUDIT_CARRIER " << carrier.size() << " FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec << " RANK8 " << rank8
                  << " RANK9 " << carrier.size() - rank8 << " HOSTILE_DELETE "
                  << (hostile_delete == 0 ? std::string("NONE")
                                          : hex8(hostile_delete)) << '\n'
                  << "ROWS " << rows.size() << " ENDPOINT 594 ROW_FNV "
                  << std::hex << kTop594Fnv << std::dec << " WORKERS "
                  << worker_count << " BODY_UNIVERSE_PER_ROW " << kBodyCount617
                  << " TOTAL_BODY_TESTS " << kBodyCount617 * rows.size() << '\n';
        for (std::size_t index = 0; index < rows.size(); ++index) {
            const auto pair = rows[index];
            const auto& audit = audits[index];
            std::cout << "PAIR " << pair.q << ',' << pair.r << " ACTIVE "
                      << audit.active << " ACTIVE_FNV " << std::hex
                      << audit.active_fnv << std::dec << " ACTIVE_JOINT "
                      << audit.active_joint << " ACTIVE_NONJOINT "
                      << audit.active_nonjoint << " EXPOSED " << audit.exposed
                      << " HIT_RANGE " << audit.minimum_hits << ".."
                      << audit.maximum_hits << " HIT_INCIDENCES "
                      << audit.hit_incidences << " FAILURES " << audit.failures
                      << " FAILURE_FNV " << std::hex << audit.failure_fnv
                      << std::dec << '\n';
        }
        std::cout << "SUMMARY TOTAL_EXPOSED " << total_exposed
                  << " TOTAL_HIT_INCIDENCES " << total_hits << " FAILURES "
                  << total_failures << " FAILED_ROWS " << failed_rows
                  << " OBLIGATION_FNV " << std::hex << obligation_ledger.state
                  << " PAIR_LEDGER_FNV " << pair_ledger.state << std::dec << '\n'
                  << "SCOPE FINITE_EXACT_FIXED_POOL_ENDPOINT594_LAYER_"
                     "C3925_CARRIER_ONLY_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT " << (total_failures == 0 ? "PASS" : "HOSTILE_FAIL")
                  << '\n';
        return total_failures == 0 ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT594_C3925_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}

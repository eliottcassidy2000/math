// Independent raw-body replay of the 9,027-mask endpoint-595 repair carrier.
// This path consumes the frozen additions/deletions but not the response atlas.

#define ENDPOINT617_RAW_VERIFY_MAIN endpoint595_repaired_hidden_main
#include "04-computation/lrc14_size_preserving_response_staircase_thm4300/endpoint617_exchanged_carrier_raw_verify.cpp"
#undef ENDPOINT617_RAW_VERIFY_MAIN

#include <atomic>
#include <fstream>
#include <mutex>
#include <thread>

namespace {
constexpr u64 kRepairFnv = UINT64_C(0x64ce5f9d1ec8c4c2);
constexpr u64 kOldAdditionFnv = UINT64_C(0xdc0eebaebf688c65);
constexpr u64 kOldDeleteFnv = UINT64_C(0x9240b264ab65aa62);
constexpr u64 kNewAdditionFnv = UINT64_C(0x6740cc137170afc5);
constexpr u64 kNewDeleteFnv = UINT64_C(0x3b5ca775eedae38b);
constexpr u64 kAugmentedFnv = UINT64_C(0x55e8588798885ae5);
constexpr u64 kC596Fnv = UINT64_C(0x892fef44a9e6b37e);
constexpr u64 kUniverseFnv = UINT64_C(0xdf5374d4aca67677);
constexpr u64 kOldUnionFnv = UINT64_C(0x11414a33ab91fef6);
constexpr u64 kTop28Fnv = UINT64_C(0x47981ce64825ef2a);
constexpr u64 kC595Fnv = UINT64_C(0x194ce45fffbe78b0);
constexpr u64 kRepairedPairLedgerFnv = UINT64_C(0xfbaa57899a6a4d29);

using PairKey = std::pair<int, int>;

std::set<PairKey> read_pairs_repaired(const std::filesystem::path& path,
                                      std::size_t expected, u64 expected_fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open pair set");
    std::set<PairKey> rows;
    Fnv ledger;
    std::string line;
    while (std::getline(input, line)) {
        const auto comma = line.find(',');
        require(comma != std::string::npos, "malformed pair row");
        const PairKey pair{std::stoi(line.substr(0, comma)),
                           std::stoi(line.substr(comma + 1))};
        require(pair.first > 0 && pair.first < pair.second &&
                    rows.insert(pair).second,
                "pair invalid/duplicate");
        ledger.add(pair.first);
        ledger.add(pair.second);
    }
    require(input.eof() && rows.size() == expected && ledger.state == expected_fnv,
            "pair set identity changed");
    return rows;
}

u64 pair_fnv_repaired(const std::set<PairKey>& rows) {
    Fnv ledger;
    for (const auto& [q, r] : rows) {
        ledger.add(q);
        ledger.add(r);
    }
    return ledger.state;
}
}  // namespace

#ifndef ENDPOINT595_REPAIRED_RAW_MAIN
#define ENDPOINT595_REPAIRED_RAW_MAIN main
#endif

int ENDPOINT595_REPAIRED_RAW_MAIN(int argc, char** argv) {
    try {
        require(argc == 15,
                "usage: raw JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE OLD_UNION "
                "REPAIRS76 OLD_ADD4 OLD_DELETE73 NEW_ADD10 NEW_DELETE2 "
                "PAIR_CSV FAILURE_CSV WORKERS");
        const unsigned worker_count = std::stoul(argv[14]);
        require(worker_count >= 1 && worker_count <= 28, "workers escaped 1..28");
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const auto repairs = read_mixed617(argv[7], 76, kRepairFnv);
        const auto old_additions = read_mixed617(argv[8], 4, kOldAdditionFnv);
        const auto old_deletions = read_mixed617(argv[9], 73, kOldDeleteFnv);
        const auto new_additions = read_mixed617(argv[10], 10, kNewAdditionFnv);
        const auto new_deletions = read_mixed617(argv[11], 2, kNewDeleteFnv);

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
        const std::set<u32> new_deleted(new_deletions.begin(), new_deletions.end());
        for (u32 mask : new_deletions)
            require(std::find(c596.begin(), c596.end(), mask) != c596.end() &&
                        !joint_set.contains(mask),
                    "new deletion absent/joint");
        std::vector<u32> carrier;
        for (u32 mask : c596)
            if (!new_deleted.contains(mask)) carrier.push_back(mask);
        distinct = std::set<u32>(carrier.begin(), carrier.end());
        for (u32 mask : new_additions) {
            require(distinct.insert(mask).second, "new addition overlap");
            carrier.push_back(mask);
        }
        require(carrier.size() == 9027 && masks_fnv_agent(carrier) == kC595Fnv,
                "repaired carrier identity changed");
        const u64 rank8 = std::count_if(carrier.begin(), carrier.end(),
                                       [](u32 m) { return std::popcount(m) == 8; });
        require(rank8 == 8959, "repaired carrier rank census changed");

        const auto universe = read_pairs_repaired(argv[5], 22647, kUniverseFnv);
        const auto old_union = read_pairs_repaired(argv[6], 1624, kOldUnionFnv);
        std::set<PairKey> residual;
        std::set_difference(universe.begin(), universe.end(), old_union.begin(),
                            old_union.end(), std::inserter(residual, residual.end()));
        std::vector<PairAgent> top;
        std::set<PairKey> top_set;
        for (const auto& [q, r] : residual)
            if (r == 595) {
                top.push_back(PairAgent{q, r});
                top_set.insert({q, r});
            }
        require(top.size() == 28 && pair_fnv_repaired(top_set) == kTop28Fnv,
                "endpoint-595 target layer changed");

        std::vector<PairAudit617> audits(top.size());
        std::atomic<std::size_t> cursor{0};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        auto worker = [&]() {
            try {
                while (true) {
                    const std::size_t index = cursor.fetch_add(1);
                    if (index >= top.size()) return;
                    audits[index] = audit_pair617(top[index], joint, carrier, joint_set);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
        };
        std::vector<std::thread> workers;
        for (unsigned i = 0; i < worker_count; ++i) workers.emplace_back(worker);
        for (auto& thread : workers) thread.join();
        if (worker_error) std::rethrow_exception(worker_error);

        std::ofstream pair_out(argv[12]);
        std::ofstream failure_out(argv[13]);
        require(static_cast<bool>(pair_out) && static_cast<bool>(failure_out),
                "cannot create raw ledgers");
        pair_out << "q,r,active,active_fnv,active_joint,active_nonjoint,exposed,"
                    "exposed_fnv,minimum_hits,maximum_hits,failures,failure_fnv\n";
        failure_out << "q,r,body_hex\n";
        Fnv pair_ledger;
        u64 total_exposed = 0;
        u64 total_hits = 0;
        u64 total_failures = 0;
        for (std::size_t index = 0; index < top.size(); ++index) {
            const PairAgent pair = top[index];
            const PairAudit617& audit = audits[index];
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
        require(pair_out.good() && failure_out.good() && total_failures == 0 &&
                    total_exposed == 303583 &&
                    total_hits == UINT64_C(10749547) &&
                    pair_ledger.state == kRepairedPairLedgerFnv,
                "repaired endpoint-595 raw replay failed");

        std::cout << "LRC14_ENDPOINT595_REPAIRED_RAW_AUDIT_V1\n"
                  << "C596 9019 FNV " << std::hex << kC596Fnv << std::dec
                  << " DELETE2_FNV " << std::hex << kNewDeleteFnv
                  << " ADD10_FNV " << kNewAdditionFnv << std::dec << '\n'
                  << "C595 9027 FNV " << std::hex << masks_fnv_agent(carrier)
                  << std::dec << " RANK8 " << rank8 << " RANK9 "
                  << carrier.size() - rank8 << " JOINT_RETAINED 421\n"
                  << "ROWS 28 ENDPOINT 595 FNV " << std::hex
                  << pair_fnv_repaired(top_set) << std::dec << " WORKERS "
                  << worker_count << " BODY_UNIVERSE_PER_ROW " << kBodyCount617
                  << " TOTAL_BODY_TESTS " << kBodyCount617 * top.size() << '\n';
        for (std::size_t index = 0; index < top.size(); ++index) {
            const auto pair = top[index];
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
                  << total_failures << " PAIR_LEDGER_FNV " << std::hex
                  << pair_ledger.state << std::dec << '\n'
                  << "SCOPE INDEPENDENT_RAW_BODY_REPLAY_FIXED_POOL_"
                     "ENDPOINT595_LAYER_ONLY_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS ALL_28_ENDPOINT595_ROWS_CLOSE\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT595_REPAIRED_RAW_ERROR " << error.what() << '\n';
        return 1;
    }
}

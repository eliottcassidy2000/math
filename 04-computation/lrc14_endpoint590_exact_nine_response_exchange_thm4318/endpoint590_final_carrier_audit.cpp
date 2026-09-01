// Frontier-only exact replay of THM-4313/THM-4314's unchanged final
// size-3925 carrier on the canonical 13-row endpoint-590 layer. This is a
// finite fixed-pool scout;
// it asserts neither a physical entry nor LRC(14).

#define ENDPOINT593_RESPONSE_CAPACITY_MAIN endpoint593_capacity_hidden_for_endpoint590
#include "04-computation/lrc14_endpoint593_response_exchange_thm4311/endpoint593_response_capacity.cpp"
#undef ENDPOINT593_RESPONSE_CAPACITY_MAIN

#include <atomic>
#include <mutex>
#include <thread>

namespace {
constexpr u64 kAddition593Fnv = UINT64_C(0x60873ef7a2b4ab90);
constexpr u64 kDeletion593Fnv = UINT64_C(0x4c14214a64ec202c);
constexpr u64 kExchangedCarrierFnv = UINT64_C(0xc9e5faef52ca5707);
constexpr u64 kCover43Fnv = UINT64_C(0xca3cb80f471f2e7e);
constexpr u64 kDeletion43Fnv = UINT64_C(0x0dd14ef0fe3eec62);
constexpr u64 kTop590Fnv = UINT64_C(0x44aa8a793d162cf9);
constexpr u64 kFinalCarrier590Fnv = UINT64_C(0xa0d08a38c10bdab7);

std::vector<u32> read_cover590(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open THM4313 cover");
    std::string line;
    require(std::getline(input, line) && line.starts_with("mask_hex,"),
            "THM4313 cover header changed");
    std::vector<u32> result;
    std::set<u32> distinct;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto comma = line.find(',');
        require(comma != std::string::npos, "malformed THM4313 cover row");
        const u32 mask = parse_mask_agent(line.substr(0, comma));
        require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                    distinct.insert(mask).second,
                "THM4313 cover rank/distinctness changed");
        result.push_back(mask);
    }
    require(input.eof() && result.size() == 43,
            "THM4313 cover count changed");
    return result;
}

std::vector<u32> read_deletions590(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open THM4313 deletion ledger");
    std::vector<u32> result;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask_agent(token);
        require(distinct.insert(mask).second, "duplicate THM4313 deletion");
        result.push_back(mask);
    }
    require(input.eof() && result.size() == 43,
            "THM4313 deletion count changed");
    return result;
}
}

#ifndef ENDPOINT590_FINAL_CARRIER_AUDIT_MAIN
#define ENDPOINT590_FINAL_CARRIER_AUDIT_MAIN main
#endif

int ENDPOINT590_FINAL_CARRIER_AUDIT_MAIN(int argc, char** argv) {
    try {
        require(argc == 20,
                "usage: endpoint590 JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE "
                "ADDITION593 DELETE593 TOP590 COVER43 DELETE43 PAIR_CSV "
                "FAILURE_CSV WORKERS");
        const unsigned worker_count = std::stoul(argv[19]);
        require(worker_count >= 1 && worker_count <= 64,
                "workers escaped 1..64");

        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> old_carrier = reconstruct_c3925_593(
            argv[2], argv[3], argv[4], argv[7], argv[8], argv[9], argv[10],
            argv[11]);
        (void)read_pairs_repaired(argv[5], 22647, kUniverseFnv);
        (void)read_pairs_repaired(argv[6], 1624, kOldUnionFnv);

        const auto add593 = read_flexible_masks593(argv[12], 1, kAddition593Fnv);
        const auto del593 = read_flexible_masks593(argv[13], 1, kDeletion593Fnv);
        std::vector<u32> carrier593;
        for (u32 mask : old_carrier)
            if (mask != del593.front()) carrier593.push_back(mask);
        require(carrier593.size() + 1 == old_carrier.size(),
                "THM4311 deletion absent/nonunique");
        carrier593.push_back(add593.front());
        require(carrier593.size() == 3925 &&
                    masks_fnv_agent(carrier593) == kExchangedCarrierFnv,
                "THM4311 carrier identity changed");
        const std::unordered_set<u32> carrier593_set(carrier593.begin(),
                                                     carrier593.end());

        const auto additions = read_cover590(argv[15]);
        const auto deletions = read_deletions590(argv[16]);
        require(additions.size() == 43 &&
                    masks_fnv_agent(additions) == kCover43Fnv,
                "THM4313 cover identity changed");
        require(deletions.size() == 43 &&
                    masks_fnv_agent(deletions) == kDeletion43Fnv,
                "THM4313 deletion identity changed");
        const std::unordered_set<u32> deletion_set(deletions.begin(),
                                                   deletions.end());
        std::vector<u32> carrier;
        for (u32 mask : carrier593) {
            if (deletion_set.contains(mask)) {
                require(!joint_set.contains(mask),
                        "THM4313 deletion entered joint deck");
            } else {
                carrier.push_back(mask);
            }
        }
        require(carrier.size() + deletions.size() == carrier593.size(),
                "THM4313 deletion absent from carrier");
        for (u32 mask : additions) {
            require(!joint_set.contains(mask) && !carrier593_set.contains(mask),
                    "THM4313 addition overlaps inherited carrier/joint deck");
            carrier.push_back(mask);
        }
        const u64 rank8 = std::count_if(
            carrier.begin(), carrier.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        require(carrier.size() == 3925 && rank8 == 3818 &&
                    masks_fnv_agent(carrier) == kFinalCarrier590Fnv &&
                    std::set<u32>(carrier.begin(), carrier.end()).size() ==
                        carrier.size(),
                "THM4313 final carrier identity changed");
        for (u32 mask : joint)
            require(std::find(carrier.begin(), carrier.end(), mask) !=
                        carrier.end(),
                    "THM4313 carrier lost protected joint mask");

        const auto row_set = read_pairs_repaired(argv[14], 13, kTop590Fnv);
        std::vector<PairAgent> rows;
        Fnv row_ledger;
        for (const auto& [q, r] : row_set) {
            require(r == 590, "frontier row escaped endpoint590");
            rows.push_back(PairAgent{q, r});
            row_ledger.add(q); row_ledger.add(r);
        }
        require(rows.size() == 13 && row_ledger.state == kTop590Fnv,
                "endpoint590 row identity changed");

        std::vector<PairAudit617> audits(rows.size());
        std::atomic<std::size_t> cursor{0};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        auto worker = [&]() {
            try {
                while (true) {
                    const std::size_t index = cursor.fetch_add(1);
                    if (index >= rows.size()) return;
                    audits[index] =
                        audit_pair617(rows[index], joint, carrier, joint_set);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
        };
        std::vector<std::thread> workers;
        for (unsigned index = 0; index < worker_count; ++index)
            workers.emplace_back(worker);
        for (std::thread& thread : workers) thread.join();
        if (worker_error) std::rethrow_exception(worker_error);

        std::ofstream pair_out(argv[17]);
        std::ofstream failure_out(argv[18]);
        require(pair_out && failure_out, "cannot create endpoint590 ledgers");
        pair_out << "q,r,active,active_fnv,active_joint,active_nonjoint,exposed,"
                    "exposed_fnv,minimum_hits,maximum_hits,failures,failure_fnv\n";
        failure_out << "q,r,body_hex\n";
        Fnv pair_ledger, failure_ledger;
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
                failure_out << pair.q << ',' << pair.r << ',' << hex8(body)
                            << '\n';
                failure_ledger.add(pair.q); failure_ledger.add(pair.r);
                failure_ledger.add(body);
            }
        }
        require(pair_out.good() && failure_out.good(),
                "endpoint590 ledger write failed");
        require(total_exposed == 154023 && total_hits == 4256322 &&
                    total_failures == 100 && failed_rows == 1,
                "endpoint590 aggregate census changed");
        require(failure_ledger.state == UINT64_C(0x8d19cba1e86e53b5) &&
                    pair_ledger.state == UINT64_C(0x0c0eb3343c5a35bf),
                "endpoint590 failure/pair ledger identity changed");

        std::cout << "LRC14_ENDPOINT590_UNCHANGED_CARRIER_AUDIT_V1\n"
                  << "CARRIER " << carrier.size() << " FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec << " RANK8 " << rank8
                  << " RANK9 " << carrier.size() - rank8 << " JOINT_RETAINED "
                  << joint.size() << '\n'
                  << "ROWS " << rows.size() << " ENDPOINT 590 ROW_FNV "
                  << std::hex << row_ledger.state << std::dec << " WORKERS "
                  << worker_count << " BODY_TESTS "
                  << kBodyCount617 * rows.size() << '\n';
        for (std::size_t index = 0; index < rows.size(); ++index) {
            const auto pair = rows[index];
            const auto& audit = audits[index];
            std::cout << "PAIR " << pair.q << ',' << pair.r << " ACTIVE "
                      << audit.active << " EXPOSED " << audit.exposed
                      << " HIT_RANGE " << audit.minimum_hits << ".."
                      << audit.maximum_hits << " FAILURES " << audit.failures
                      << " FAILURE_FNV " << std::hex << audit.failure_fnv
                      << std::dec << '\n';
        }
        std::cout << "SUMMARY EXPOSED " << total_exposed
                  << " HIT_INCIDENCES " << total_hits << " FAILURES "
                  << total_failures << " FAILED_ROWS " << failed_rows
                  << " FAILURE_FNV " << std::hex << failure_ledger.state
                  << " PAIR_LEDGER_FNV " << pair_ledger.state << std::dec << '\n'
                  << "SCOPE FINITE_EXACT_FIXED_POOL_ENDPOINT590_LAYER_"
                     "THM4313_THM4314_UNCHANGED_CARRIER_ONLY_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT "
                  << (total_failures == 0 ? "PASS" : "HOSTILE_FAIL") << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT590_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}

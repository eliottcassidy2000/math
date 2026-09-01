// Generic exact scout of the post-endpoint-590 nine-for-nine exchanged
// carrier on a frozen single-endpoint typed frontier. Scratch only.

#define ENDPOINT590_QUOTIENT_MAIN endpoint590_quotient_hidden_for_generic_endpoint
#include "04-computation/lrc14_endpoint590_exact_nine_response_exchange_thm4318/endpoint590_protected_deletion_quotient.cpp"
#undef ENDPOINT590_QUOTIENT_MAIN

#include <atomic>
#include <mutex>
#include <thread>

namespace {
constexpr u64 kAddition9Fnv = UINT64_C(0xd1cf49e4b811b958);
constexpr u64 kDeletion9Fnv = UINT64_C(0x3546eb56552b4cde);
constexpr u64 kGenericEndpointCarrierFnv = UINT64_C(0xeeae5518d84ccac5);

std::vector<u32> read_delete9_generic(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open endpoint590 delete9");
    std::vector<u32> result;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask_agent(token);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "endpoint590 delete9 rank/distinctness changed");
        result.push_back(mask);
    }
    require(input.eof() && result.size() == 9,
            "endpoint590 delete9 count changed");
    return result;
}
}

#ifndef GENERIC_ENDPOINT_EXCHANGED_CARRIER_AUDIT_MAIN
#define GENERIC_ENDPOINT_EXCHANGED_CARRIER_AUDIT_MAIN main
#endif

int GENERIC_ENDPOINT_EXCHANGED_CARRIER_AUDIT_MAIN(int argc, char** argv) {
    try {
        require(argc == 25,
                "usage: generic_endpoint JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE "
                "ADDITION593 DELETE593 COVER43 DELETE43 TOP COVER9 "
                "DELETE9 ENDPOINT ROWS ROW_FNV_HEX PAIR_OUT FAILURE_OUT WORKERS");
        const int expected_endpoint = std::stoi(argv[19]);
        const std::size_t expected_rows = std::stoull(argv[20]);
        const u64 expected_row_fnv = std::stoull(argv[21], nullptr, 16);
        const unsigned worker_count = std::stoul(argv[24]);
        require(expected_endpoint > 0 && expected_rows > 0,
                "invalid endpoint/row count");
        require(worker_count >= 1 && worker_count <= 64,
                "workers escaped 1..64");
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> base =
            reconstruct_final_carrier_quotient(argv, joint);
        const std::unordered_set<u32> base_set(base.begin(), base.end());
        const std::vector<u32> additions = read_cover_quotient(argv[17]);
        const std::vector<u32> deletions = read_delete9_generic(argv[18]);
        require(masks_fnv_agent(additions) == kAddition9Fnv &&
                    masks_fnv_agent(deletions) == kDeletion9Fnv,
                "endpoint590 exchange ledger identity changed");
        const std::unordered_set<u32> deletion_set(deletions.begin(),
                                                   deletions.end());
        std::vector<u32> carrier;
        for (u32 mask : base) {
            if (deletion_set.contains(mask)) {
                require(!joint_set.contains(mask), "exchange deleted joint mask");
            } else {
                carrier.push_back(mask);
            }
        }
        require(carrier.size() + deletions.size() == base.size(),
                "exchange deletion absent from base carrier");
        for (u32 mask : additions) {
            require(!base_set.contains(mask) && !joint_set.contains(mask),
                    "exchange addition overlaps base/joint deck");
            carrier.push_back(mask);
        }
        const u64 rank8 = std::count_if(carrier.begin(), carrier.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        require(carrier.size() == 3925 && rank8 == 3809 &&
                    masks_fnv_agent(carrier) == kGenericEndpointCarrierFnv &&
                    std::set<u32>(carrier.begin(), carrier.end()).size() ==
                        carrier.size(),
                "endpoint590 exchanged carrier identity changed");
        for (u32 mask : joint)
            require(std::find(carrier.begin(), carrier.end(), mask) !=
                        carrier.end(),
                    "endpoint590 exchange lost protected joint mask");

        const auto row_set = read_pairs_repaired(
            argv[16], expected_rows, expected_row_fnv);
        std::vector<PairAgent> rows;
        Fnv row_ledger;
        for (const auto& [q, r] : row_set) {
            require(r == expected_endpoint, "frontier row escaped endpoint");
            rows.push_back({q, r});
            row_ledger.add(q); row_ledger.add(r);
        }

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
        for (auto& thread : workers) thread.join();
        if (worker_error) std::rethrow_exception(worker_error);

        std::ofstream pair_out(argv[22]);
        std::ofstream failure_out(argv[23]);
        require(pair_out && failure_out, "cannot create endpoint ledgers");
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
                "endpoint ledger write failed");

        std::cout << "LRC14_GENERIC_ENDPOINT_EXCHANGED_CARRIER_AUDIT_V1\n"
                  << "CARRIER " << carrier.size() << " FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec << " RANK8 " << rank8
                  << " RANK9 " << carrier.size() - rank8
                  << " JOINT_RETAINED " << joint.size() << '\n'
                  << "ROWS " << rows.size() << " ENDPOINT "
                  << expected_endpoint << " ROW_FNV "
                  << std::hex << row_ledger.state << std::dec << " WORKERS "
                  << worker_count << " BODY_TESTS "
                  << rows.size() * kBodyCount617 << '\n';
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
                  << "SCOPE FINITE_EXACT_FIXED_POOL_SINGLE_ENDPOINT_LAYER_"
                     "POST590_EXCHANGED_CARRIER_ONLY_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT "
                  << (total_failures == 0 ? "PASS" : "HOSTILE_FAIL") << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "GENERIC_ENDPOINT_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}

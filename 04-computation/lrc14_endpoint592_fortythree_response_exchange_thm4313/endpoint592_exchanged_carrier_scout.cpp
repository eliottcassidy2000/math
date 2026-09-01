// Frontier-only hostile replay of the endpoint-593 exchanged carrier on the
// exact 35-row endpoint-592 typed frontier.

#define ENDPOINT593_RESPONSE_CAPACITY_MAIN endpoint593_response_hidden_for_endpoint592
#include "04-computation/lrc14_endpoint593_response_exchange_thm4311/endpoint593_response_capacity.cpp"
#undef ENDPOINT593_RESPONSE_CAPACITY_MAIN

#include <atomic>
#include <mutex>
#include <thread>

namespace {
constexpr u64 kTop592FnvScout = UINT64_C(0x3eb23833c35b9266);
constexpr u64 kExchangedCarrierFnvScout = UINT64_C(0xc9e5faef52ca5707);
constexpr u64 kAddition1FnvScout = UINT64_C(0x60873ef7a2b4ab90);
constexpr u64 kDelete1FnvScout = UINT64_C(0x4c14214a64ec202c);
}

int main(int argc, char** argv) {
    try {
        require(argc == 18,
                "usage: scout JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE "
                "ADDITION1 DELETE1 TOP592 PAIR_CSV FAILURE_CSV WORKERS");
        const unsigned worker_count = std::stoul(argv[17]);
        require(worker_count >= 1 && worker_count <= 64,
                "workers escaped 1..64");

        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> carrier = reconstruct_c3925_593(
            argv[2], argv[3], argv[4], argv[7], argv[8], argv[9], argv[10],
            argv[11]);
        const auto addition = read_flexible_masks593(
            argv[12], 1, kAddition1FnvScout);
        const auto deletion = read_flexible_masks593(
            argv[13], 1, kDelete1FnvScout);
        require(!joint_set.contains(addition.front()) &&
                    !joint_set.contains(deletion.front()),
                "exchange entered protected joint deck");

        std::vector<u32> exchanged;
        for (u32 mask : carrier)
            if (mask != deletion.front()) exchanged.push_back(mask);
        require(exchanged.size() + 1 == carrier.size(),
                "deletion absent/nonunique");
        require(std::find(exchanged.begin(), exchanged.end(), addition.front()) ==
                    exchanged.end(),
                "addition overlaps exchanged carrier");
        exchanged.push_back(addition.front());
        const u64 rank8 = std::count_if(
            exchanged.begin(), exchanged.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        require(exchanged.size() == 3925 && rank8 == 3857 &&
                    masks_fnv_agent(exchanged) == kExchangedCarrierFnvScout,
                "exchanged carrier identity changed");
        for (u32 mask : joint)
            require(std::find(exchanged.begin(), exchanged.end(), mask) !=
                        exchanged.end(),
                    "exchange lost protected joint mask");

        (void)read_pairs_repaired(argv[5], 22647, kUniverseFnv);
        (void)read_pairs_repaired(argv[6], 1624, kOldUnionFnv);
        const auto row_set = read_pairs_repaired(argv[14], 35, kTop592FnvScout);
        std::vector<PairAgent> rows;
        for (const auto& [q, r] : row_set) {
            require(r == 592, "target escaped endpoint 592");
            rows.push_back(PairAgent{q, r});
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
                        audit_pair617(rows[index], joint, exchanged, joint_set);
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

        std::ofstream pair_out(argv[15]);
        std::ofstream failure_out(argv[16]);
        require(pair_out && failure_out, "cannot create frontier ledgers");
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
                failure_out << pair.q << ',' << pair.r << ',' << hex8(body)
                            << '\n';
                obligation_ledger.add(pair.q);
                obligation_ledger.add(pair.r);
                obligation_ledger.add(body);
            }
        }
        require(pair_out.good() && failure_out.good(),
                "frontier ledger write failed");

        std::cout << "LRC14_ENDPOINT592_EXCHANGED_CARRIER_SCOUT_V1\n"
                  << "CARRIER " << exchanged.size() << " FNV " << std::hex
                  << masks_fnv_agent(exchanged) << std::dec << " RANK8 "
                  << rank8 << " RANK9 " << exchanged.size() - rank8
                  << " JOINT_RETAINED " << joint.size() << '\n'
                  << "ROWS " << rows.size() << " ENDPOINT 592 ROW_FNV "
                  << std::hex << kTop592FnvScout << std::dec << " WORKERS "
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
                  << " OBLIGATION_FNV " << std::hex << obligation_ledger.state
                  << " PAIR_LEDGER_FNV " << pair_ledger.state << std::dec << '\n'
                  << "SCOPE FINITE_EXACT_FIXED_POOL_ENDPOINT592_LAYER_"
                     "EXCHANGED_CARRIER_ONLY_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT " << (total_failures == 0 ? "PASS" : "HOSTILE_FAIL")
                  << '\n';
        return total_failures == 0 ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT592_SCOUT_ERROR " << error.what() << '\n';
        return 1;
    }
}

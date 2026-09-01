// Complete direct raw replay of a proposed nine-addition/nine-deletion
// endpoint-590 exchange on all 480 inherited rows plus 13 endpoint-590 rows.
// Scratch only; no physical-entry claim.

#define ENDPOINT590_QUOTIENT_MAIN endpoint590_quotient_hidden_for_exchange493
#include "endpoint590_protected_deletion_quotient.cpp"
#undef ENDPOINT590_QUOTIENT_MAIN

namespace {

constexpr u64 kTarget493FnvExchange = UINT64_C(0x1fef91ec25d074e5);
constexpr u64 kSelectedCover9FnvExchange = UINT64_C(0xd1cf49e4b811b958);
constexpr u64 kSelectedDelete9FnvExchange = UINT64_C(0x3546eb56552b4cde);
constexpr u64 kSelectedCarrierFnvExchange = UINT64_C(0xeeae5518d84ccac5);

std::vector<u32> read_delete9_exchange(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open endpoint590 deletions");
    std::vector<u32> result;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask_agent(token);
        require(distinct.insert(mask).second, "duplicate endpoint590 deletion");
        result.push_back(mask);
    }
    require(input.eof() && result.size() == 9,
            "endpoint590 deletion size changed");
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 24,
                "usage: exchange493 JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE "
                "ADDITION593 DELETE593 COVER43 DELETE43 FULL467_PAIR TOP591 "
                "TOP590 COVER9 DELETE9 PAIR_OUT FAILURE_OUT WORKERS");
        const unsigned worker_count = std::stoul(argv[23]);
        require(worker_count >= 1 && worker_count <= 64,
                "workers escaped 1..64");
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> carrier =
            reconstruct_final_carrier_quotient(argv, joint);
        const std::unordered_set<u32> carrier_set(carrier.begin(), carrier.end());
        const std::vector<u32> additions = read_cover_quotient(argv[19]);
        const std::vector<u32> deletions = read_delete9_exchange(argv[20]);
        require(masks_fnv_agent(additions) == kSelectedCover9FnvExchange &&
                    masks_fnv_agent(deletions) == kSelectedDelete9FnvExchange,
                "selected endpoint590 exchange ledgers changed");
        const std::set<u32> deletion_set(deletions.begin(), deletions.end());
        std::vector<u32> exchanged;
        for (u32 mask : carrier) {
            if (deletion_set.contains(mask)) {
                require(!joint_set.contains(mask), "exchange deleted joint mask");
            } else {
                exchanged.push_back(mask);
            }
        }
        require(exchanged.size() + deletions.size() == carrier.size(),
                "exchange deletion absent from carrier");
        for (u32 mask : additions) {
            require(!carrier_set.contains(mask) && !joint_set.contains(mask),
                    "exchange addition overlaps carrier/joint deck");
            exchanged.push_back(mask);
        }
        require(exchanged.size() == 3925 &&
                    std::set<u32>(exchanged.begin(), exchanged.end()).size() ==
                        exchanged.size(),
                "exchanged carrier size/distinctness changed");
        require(masks_fnv_agent(exchanged) == kSelectedCarrierFnvExchange,
                "selected endpoint590 carrier identity changed");
        for (u32 mask : joint)
            require(std::find(exchanged.begin(), exchanged.end(), mask) !=
                        exchanged.end(),
                    "exchange lost protected joint mask");

        std::vector<PairAgent> full467;
        (void)read_full467_minimum_one(argv[16], full467);
        const auto top591 = read_pairs_repaired(argv[17], 13, kTop591Fnv590);
        const auto top590 =
            read_pairs_repaired(argv[18], 13, kTop590FnvQuotient);
        std::set<std::pair<int, int>> target_set;
        for (PairAgent pair : full467)
            require(target_set.insert({pair.q, pair.r}).second,
                    "full467 target overlap");
        for (const auto& pair : top591)
            require(target_set.insert(pair).second, "endpoint591 target overlap");
        for (const auto& pair : top590)
            require(target_set.insert(pair).second, "endpoint590 target overlap");
        require(target_set.size() == 493, "target493 count changed");
        std::vector<PairAgent> rows;
        Fnv target_ledger;
        for (const auto& [q, r] : target_set) {
            rows.push_back({q, r});
            target_ledger.add(q); target_ledger.add(r);
        }
        require(target_ledger.state == kTarget493FnvExchange,
                "target493 identity changed");

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

        std::ofstream pair_out(argv[21]);
        std::ofstream failure_out(argv[22]);
        require(pair_out && failure_out, "cannot create exchange audit outputs");
        pair_out << "q,r,active,active_fnv,active_joint,active_nonjoint,exposed,"
                    "exposed_fnv,minimum_hits,maximum_hits,failures,failure_fnv\n";
        failure_out << "q,r,body_hex\n";
        Fnv pair_ledger, failure_ledger;
        u64 total_exposed = 0, total_hits = 0, total_failures = 0;
        std::size_t failed_rows = 0;
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
                "exchange audit output write failed");

        const u64 rank8 = std::count_if(exchanged.begin(), exchanged.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        std::cout << "LRC14_ENDPOINT590_EXCHANGE493_AUDIT_V1\n"
                  << "BASE_CARRIER 3925 FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec << " ADDITIONS 9 FNV "
                  << std::hex << masks_fnv_agent(additions) << std::dec
                  << " DELETIONS 9 FNV " << std::hex
                  << masks_fnv_agent(deletions) << std::dec << '\n'
                  << "EXCHANGED_CARRIER 3925 FNV " << std::hex
                  << masks_fnv_agent(exchanged) << std::dec << " RANK8 "
                  << rank8 << " RANK9 " << exchanged.size() - rank8
                  << " JOINT_RETAINED 421\n"
                  << "ROWS 493 ROW_FNV " << std::hex << target_ledger.state
                  << std::dec << " BODY_TESTS " << rows.size() * kBodyCount617
                  << '\n'
                  << "SUMMARY EXPOSED " << total_exposed << " HIT_INCIDENCES "
                  << total_hits << " FAILURES " << total_failures
                  << " FAILED_ROWS " << failed_rows << " FAILURE_FNV "
                  << std::hex << failure_ledger.state << " PAIR_LEDGER_FNV "
                  << pair_ledger.state << std::dec << '\n'
                  << "SCOPE DIRECT_SIMULTANEOUS_FIXED_NINE_ADD_NINE_DELETE_"
                     "TARGET493_RAW_REPLAY_PROTECTED_JOINT_NO_PHYSICAL_ENTRY_"
                     "NO_LRC14\nVERDICT "
                  << (total_failures == 0 ? "PASS" : "HOSTILE_FAIL") << '\n';
        return total_failures == 0 ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT590_EXCHANGE493_ERROR " << error.what() << '\n';
        return 1;
    }
}

// Complete raw audit for a proposed same-cardinality endpoint-592 exchange.
// The additions and deletions are external ledgers of equal cardinality.

#define ENDPOINT593_RESPONSE_CAPACITY_MAIN endpoint593_capacity_hidden_for_592_exchange
#include "04-computation/lrc14_endpoint593_response_exchange_thm4311/endpoint593_response_capacity.cpp"
#undef ENDPOINT593_RESPONSE_CAPACITY_MAIN

#include <atomic>
#include <mutex>
#include <thread>

namespace {
constexpr u64 kAddition593Fnv = UINT64_C(0x60873ef7a2b4ab90);
constexpr u64 kDeletion593Fnv = UINT64_C(0x4c14214a64ec202c);
constexpr u64 kExchangedCarrierFnv = UINT64_C(0xc9e5faef52ca5707);
constexpr u64 kTop592Fnv = UINT64_C(0x3eb23833c35b9266);
constexpr u64 kFull467Fnv = UINT64_C(0x2d6aa988098aa5eb);
constexpr u64 kCover43Fnv = UINT64_C(0xca3cb80f471f2e7e);
constexpr u64 kDeletion43Fnv = UINT64_C(0x0dd14ef0fe3eec62);
constexpr u64 kFinalCarrierFnv = UINT64_C(0xa0d08a38c10bdab7);

std::vector<u32> read_cover_masks_exchange(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open endpoint592 cover");
    std::string line;
    require(std::getline(input, line) && line.starts_with("mask_hex,"),
            "endpoint592 cover header changed");
    std::vector<u32> result;
    std::set<u32> distinct;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto comma = line.find(',');
        require(comma != std::string::npos, "malformed endpoint592 cover row");
        const u32 mask = parse_mask_agent(line.substr(0, comma));
        require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                    distinct.insert(mask).second,
                "cover rank/distinctness changed");
        result.push_back(mask);
    }
    require(input.eof() && !result.empty(), "empty/malformed endpoint592 cover");
    return result;
}

std::vector<u32> read_deletions_exchange(const std::filesystem::path& path,
                                         std::size_t expected) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open endpoint592 deletions");
    std::vector<u32> result;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask_agent(token);
        require(distinct.insert(mask).second, "duplicate endpoint592 deletion");
        result.push_back(mask);
    }
    require(input.eof() && result.size() == expected,
            "endpoint592 deletion count changed");
    return result;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 23,
                "usage: audit JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE OLD_UNION "
                "REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE PREFIX391 TOP594 "
                "TOP593 ADDITION593 DELETE593 TOP592 COVER_CSV DELETE_TXT "
                "PAIR_CSV FAILURE_CSV WORKERS");
        const unsigned worker_count = std::stoul(argv[22]);
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
        const auto add593 = read_flexible_masks593(argv[15], 1, kAddition593Fnv);
        const auto del593 = read_flexible_masks593(argv[16], 1, kDeletion593Fnv);
        std::vector<u32> carrier;
        for (u32 mask : old_carrier)
            if (mask != del593.front()) carrier.push_back(mask);
        require(carrier.size() + 1 == old_carrier.size(),
                "THM4311 deletion absent/nonunique");
        carrier.push_back(add593.front());
        require(carrier.size() == 3925 &&
                    masks_fnv_agent(carrier) == kExchangedCarrierFnv,
                "THM4311 exchanged carrier changed");
        const std::unordered_set<u32> carrier_set(carrier.begin(), carrier.end());

        const auto additions = read_cover_masks_exchange(argv[18]);
        const auto deletions = read_deletions_exchange(argv[19], additions.size());
        require(additions.size() == 43 &&
                    masks_fnv_agent(additions) == kCover43Fnv,
                "frozen endpoint592 cover identity changed");
        require(deletions.size() == 43 &&
                    masks_fnv_agent(deletions) == kDeletion43Fnv,
                "frozen endpoint592 deletion identity changed");
        const std::unordered_set<u32> deletion_set(deletions.begin(), deletions.end());
        std::vector<u32> exchanged;
        for (u32 mask : carrier) {
            if (deletion_set.contains(mask)) {
                require(!joint_set.contains(mask), "deletion entered joint deck");
            } else {
                exchanged.push_back(mask);
            }
        }
        require(exchanged.size() + deletions.size() == carrier.size(),
                "deletion absent from carrier");
        for (u32 mask : additions) {
            require(!joint_set.contains(mask) && !carrier_set.contains(mask),
                    "addition overlaps carrier/joint deck");
            exchanged.push_back(mask);
        }
        require(exchanged.size() == carrier.size() &&
                    std::set<u32>(exchanged.begin(), exchanged.end()).size() ==
                        exchanged.size(),
                "exchange size/distinctness changed");
        const u64 exchanged_rank8 = std::count_if(
            exchanged.begin(), exchanged.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        require(exchanged_rank8 == 3818 &&
                    masks_fnv_agent(exchanged) == kFinalCarrierFnv,
                "frozen final carrier rank/FNV changed");
        for (u32 mask : joint)
            require(std::find(exchanged.begin(), exchanged.end(), mask) !=
                        exchanged.end(),
                    "exchange lost protected joint mask");

        const auto prefix391 = read_pair_audit_rows593(argv[12], 391);
        const auto top594 = read_pair_audit_rows593(argv[13], 25);
        const auto top593 = read_pairs_repaired(argv[14], 16, kTop593Fnv593);
        const auto top592 = read_pairs_repaired(argv[17], 35, kTop592Fnv);
        std::set<std::pair<int, int>> target_set;
        for (PairAgent pair : prefix391)
            require(target_set.insert({pair.q, pair.r}).second,
                    "prefix391 overlap");
        for (PairAgent pair : top594)
            require(target_set.insert({pair.q, pair.r}).second, "top594 overlap");
        for (const auto& pair : top593)
            require(target_set.insert(pair).second, "top593 overlap");
        for (const auto& pair : top592)
            require(target_set.insert(pair).second, "top592 overlap");
        require(target_set.size() == 467, "full target count changed");
        std::vector<PairAgent> rows;
        Fnv target_ledger;
        for (const auto& [q, r] : target_set) {
            rows.push_back(PairAgent{q, r});
            target_ledger.add(q); target_ledger.add(r);
        }
        require(target_ledger.state == kFull467Fnv,
                "full467 target identity changed");

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
        for (std::thread& thread : workers) thread.join();
        if (worker_error) std::rethrow_exception(worker_error);

        std::ofstream pair_out(argv[20]);
        std::ofstream failure_out(argv[21]);
        require(pair_out && failure_out, "cannot create raw audit ledgers");
        pair_out << "q,r,active,active_fnv,active_joint,active_nonjoint,exposed,"
                    "exposed_fnv,minimum_hits,maximum_hits,failures,failure_fnv\n";
        failure_out << "q,r,body_hex\n";
        Fnv pair_ledger, failure_ledger;
        u64 total_exposed = 0, total_hits = 0, total_failures = 0;
        for (std::size_t index = 0; index < rows.size(); ++index) {
            const auto pair = rows[index];
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
            for (u32 body : audit.failure_bodies) {
                failure_out << pair.q << ',' << pair.r << ',' << hex8(body)
                            << '\n';
                failure_ledger.add(pair.q); failure_ledger.add(pair.r);
                failure_ledger.add(body);
            }
        }
        require(pair_out.good() && failure_out.good(),
                "raw audit ledger write failed");

        std::cout << "LRC14_ENDPOINT592_FULL467_EXCHANGE_RAW_AUDIT_V1\n"
                  << "BASE_CARRIER " << carrier.size() << " FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec << " ADDITIONS "
                  << additions.size() << " ADDITION_FNV " << std::hex
                  << masks_fnv_agent(additions) << std::dec << " DELETIONS "
                  << deletions.size() << " DELETION_FNV " << std::hex
                  << masks_fnv_agent(deletions) << std::dec << '\n'
                  << "EXCHANGED_CARRIER " << exchanged.size() << " FNV "
                  << std::hex << masks_fnv_agent(exchanged) << std::dec
                  << " RANK8 " << exchanged_rank8 << " RANK9 "
                  << exchanged.size() - exchanged_rank8
                  << " JOINT_RETAINED " << joint.size() << '\n'
                  << "ROWS " << rows.size() << " ROW_FNV " << std::hex
                  << target_ledger.state << std::dec << " BODY_TESTS "
                  << kBodyCount617 * rows.size() << '\n'
                  << "SUMMARY EXPOSED " << total_exposed << " HIT_INCIDENCES "
                  << total_hits << " FAILURES " << total_failures
                  << " FAILURE_FNV " << std::hex << failure_ledger.state
                  << " PAIR_LEDGER_FNV " << pair_ledger.state << std::dec << '\n'
                  << "SCOPE DIRECT_SIMULTANEOUS_FULL467_RAW_REPLAY_FIXED_POOL_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\nVERDICT "
                  << (total_failures == 0 ? "PASS" : "HOSTILE_FAIL") << '\n';
        return total_failures == 0 ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT592_EXCHANGE_ERROR " << error.what() << '\n';
        return 1;
    }
}

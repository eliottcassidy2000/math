// Primary finite-exact audit for THM-4303.  Reconstruct THM-4302's C_596,
// derive the current residual top layer from the frozen universe and typed
// union, and replay every labelled nine-body on every one of those rows.

#define ENDPOINT617_RAW_VERIFY_MAIN thm4303_raw_hidden_main
#include "../lrc14_size_preserving_response_staircase_thm4300/endpoint617_exchanged_carrier_raw_verify.cpp"
#undef ENDPOINT617_RAW_VERIFY_MAIN

#include <atomic>
#include <fstream>
#include <map>
#include <mutex>
#include <thread>

namespace {

constexpr u64 kRepairFnv = UINT64_C(0x64ce5f9d1ec8c4c2);
constexpr u64 kAdditionFnv = UINT64_C(0xdc0eebaebf688c65);
constexpr u64 kDeleteFnv = UINT64_C(0x9240b264ab65aa62);
constexpr u64 kAugmentedFnv = UINT64_C(0x55e8588798885ae5);
constexpr u64 kCarrierFnv = UINT64_C(0x892fef44a9e6b37e);
constexpr u64 kOldUnionFnv = UINT64_C(0x11414a33ab91fef6);
constexpr u64 kThm4302UnionFnv = UINT64_C(0xb1c8ecf1dd4a71c5);
constexpr u64 kThm4302ResidualFnv = UINT64_C(0x7da11cd038486887);
constexpr u64 kTop28Fnv = UINT64_C(0x47981ce64825ef2a);

using PairKey = std::pair<int, int>;

std::set<PairKey> read_pairs4303(const std::filesystem::path& path,
                                 std::size_t expected_count,
                                 u64 expected_fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open pair ledger");
    std::set<PairKey> rows;
    Fnv ledger;
    std::string line;
    PairKey previous{0, 0};
    bool have_previous = false;
    while (std::getline(input, line)) {
        require(!line.empty(), "blank pair row");
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos &&
                    line.find(',', comma + 1) == std::string::npos,
                "malformed pair row");
        const PairKey pair{std::stoi(line.substr(0, comma)),
                           std::stoi(line.substr(comma + 1))};
        require(pair.first > 0 && pair.first < pair.second,
                "invalid pair row");
        require(!have_previous || previous < pair,
                "pair ledger order changed");
        require(rows.insert(pair).second, "duplicate pair row");
        previous = pair;
        have_previous = true;
        ledger.add(pair.first);
        ledger.add(pair.second);
    }
    require(input.eof() && rows.size() == expected_count &&
                ledger.state == expected_fnv,
            "pair ledger identity changed");
    return rows;
}

u64 pair_fnv4303(const auto& rows) {
    Fnv ledger;
    for (const auto& pair : rows) {
        ledger.add(pair.first);
        ledger.add(pair.second);
    }
    return ledger.state;
}

std::set<PairKey> set_union4303(const std::set<PairKey>& left,
                                const std::set<PairKey>& right) {
    std::set<PairKey> answer = left;
    answer.insert(right.begin(), right.end());
    return answer;
}

std::set<PairKey> set_difference4303(const std::set<PairKey>& left,
                                     const std::set<PairKey>& right) {
    std::set<PairKey> answer;
    std::set_difference(left.begin(), left.end(), right.begin(), right.end(),
                        std::inserter(answer, answer.end()));
    return answer;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 13,
                "usage: raw_audit JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION1624 REPAIRS76 ADDITIONS4 DELETE73 PAIR_CSV_OUT "
                "FAILURE_CSV_OUT WORKERS");
        const unsigned worker_count = std::stoul(argv[12]);
        require(worker_count >= 1 && worker_count <= 28,
                "worker count escaped 1..28");

        const std::vector<u32> joint = read_masks_agent(
            argv[1], kJointCountAgent, kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> repairs =
            read_mixed617(argv[7], 76, kRepairFnv);
        const std::vector<u32> additions =
            read_mixed617(argv[8], 4, kAdditionFnv);
        const std::vector<u32> deletions =
            read_mixed617(argv[9], 73, kDeleteFnv);

        std::vector<u32> augmented =
            build_mixed_carrier(argv[2], argv[3], argv[4]);
        std::set<u32> distinct(augmented.begin(), augmented.end());
        for (u32 repair : repairs) {
            require(distinct.insert(repair).second, "repair overlap");
            augmented.push_back(repair);
        }
        require(augmented.size() == 9088 &&
                    masks_fnv_agent(augmented) == kAugmentedFnv,
                "THM-4302 augmented carrier changed");
        const std::set<u32> deletion_set(deletions.begin(), deletions.end());
        std::vector<u32> carrier;
        for (u32 mask : augmented)
            if (!deletion_set.contains(mask)) carrier.push_back(mask);
        for (u32 addition : additions) {
            require(distinct.insert(addition).second, "addition overlap");
            carrier.push_back(addition);
        }
        require(carrier.size() == 9019 &&
                    masks_fnv_agent(carrier) == kCarrierFnv,
                "THM-4302 C596 identity changed");
        const u64 rank8 = std::count_if(
            carrier.begin(), carrier.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        require(rank8 == 8961, "THM-4302 carrier ranks changed");
        for (u32 mask : joint)
            require(std::find(carrier.begin(), carrier.end(), mask) !=
                        carrier.end(),
                    "joint coordinate absent from carrier");

        const std::set<PairKey> universe =
            read_pairs4303(argv[5], 22647, kResidualFnvAgent);
        const std::set<PairKey> old_union =
            read_pairs4303(argv[6], 1624, kOldUnionFnv);
        require(std::includes(universe.begin(), universe.end(),
                              old_union.begin(), old_union.end()),
                "old typed union escaped universe");
        std::set<PairKey> k596;
        for (PairKey pair : universe)
            if (pair.second >= 596) k596.insert(pair);
        const std::set<PairKey> thm4302_union =
            set_union4303(old_union, k596);
        const std::set<PairKey> thm4302_residual =
            set_difference4303(universe, thm4302_union);
        require(thm4302_union.size() == 1633 &&
                    pair_fnv4303(thm4302_union) == kThm4302UnionFnv &&
                    thm4302_residual.size() == 21014 &&
                    pair_fnv4303(thm4302_residual) == kThm4302ResidualFnv,
                "THM-4302 typed partition changed");
        int maximum_endpoint = 0;
        for (PairKey pair : thm4302_residual)
            maximum_endpoint = std::max(maximum_endpoint, pair.second);
        std::vector<PairAgent> top;
        std::set<PairKey> top_keys;
        for (PairKey pair : thm4302_residual)
            if (pair.second == maximum_endpoint) {
                top.push_back(PairAgent{pair.first, pair.second});
                top_keys.insert(pair);
            }
        require(maximum_endpoint == 595 && top.size() == 28 &&
                    pair_fnv4303(top_keys) == kTop28Fnv,
                "derived top layer changed");

        std::vector<PairAudit617> audits(top.size());
        std::atomic<std::size_t> cursor{0};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        auto worker = [&]() {
            try {
                while (true) {
                    const std::size_t index = cursor.fetch_add(1);
                    if (index >= top.size()) return;
                    audits[index] =
                        audit_pair617(top[index], joint, carrier, joint_set);
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

        const std::map<PairKey, std::pair<u64, u64>> expected_failures = {
            {{96, 595}, {116, UINT64_C(0xfedacdbff3f31981)}},
            {{100, 595}, {13, UINT64_C(0x3ac9ac8b4b9ad93f)}},
            {{210, 595}, {16, UINT64_C(0xa6a226f12c168d3a)}}};
        std::ofstream pair_out(argv[10]);
        std::ofstream failure_out(argv[11]);
        require(static_cast<bool>(pair_out) && static_cast<bool>(failure_out),
                "cannot create raw ledgers");
        pair_out << "q,r,active,active_fnv,active_joint,active_nonjoint,"
                    "exposed,exposed_fnv,minimum_hits,maximum_hits,failures,"
                    "failure_fnv\n";
        failure_out << "q,r,body_hex\n";

        u64 closed_rows = 0;
        u64 failed_rows = 0;
        u64 total_exposed = 0;
        u64 total_hit_incidences = 0;
        u64 total_failures = 0;
        Fnv pair_ledger;
        Fnv global_failure_ledger;
        for (std::size_t index = 0; index < top.size(); ++index) {
            const PairAgent pair = top[index];
            const PairAudit617& audit = audits[index];
            const auto expected = expected_failures.find({pair.q, pair.r});
            if (expected == expected_failures.end()) {
                require(audit.failures == 0,
                        "unexpected endpoint-595 failing row");
                ++closed_rows;
            } else {
                require(audit.failures == expected->second.first &&
                            audit.failure_fnv == expected->second.second,
                        "endpoint-595 failure identity changed");
                ++failed_rows;
            }
            total_exposed += audit.exposed;
            total_hit_incidences += audit.hit_incidences;
            total_failures += audit.failures;
            pair_ledger.add(pair.q);
            pair_ledger.add(pair.r);
            pair_ledger.add(audit.active);
            pair_ledger.add(audit.active_fnv);
            pair_ledger.add(audit.active_joint);
            pair_ledger.add(audit.active_nonjoint);
            pair_ledger.add(audit.exposed);
            pair_ledger.add(audit.exposed_fnv);
            pair_ledger.add(audit.minimum_hits);
            pair_ledger.add(audit.maximum_hits);
            pair_ledger.add(audit.failures);
            pair_ledger.add(audit.failure_fnv);
            pair_out << pair.q << ',' << pair.r << ',' << audit.active << ','
                     << std::hex << audit.active_fnv << std::dec << ','
                     << audit.active_joint << ',' << audit.active_nonjoint
                     << ',' << audit.exposed << ',' << std::hex
                     << audit.exposed_fnv << std::dec << ','
                     << audit.minimum_hits << ',' << audit.maximum_hits << ','
                     << audit.failures << ',' << std::hex << audit.failure_fnv
                     << std::dec << '\n';
            for (u32 body : audit.failure_bodies) {
                failure_out << pair.q << ',' << pair.r << ',' << hex8(body)
                            << '\n';
                global_failure_ledger.add(pair.q);
                global_failure_ledger.add(pair.r);
                global_failure_ledger.add(body);
            }
        }
        require(pair_out.good() && failure_out.good(),
                "raw ledger write failed");
        require(closed_rows == 25 && failed_rows == 3 &&
                    total_failures == 145,
                "endpoint-595 closure partition changed");

        std::cout << "THM4303_ENDPOINT595_PRIMARY_RAW_AUDIT_V1\n"
                  << "CARRIER 9019 FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec << " RANK8 "
                  << rank8 << " RANK9 " << carrier.size() - rank8
                  << " JOINT_RETAINED 421\n"
                  << "ROW_SOURCE DERIVED_TOP_LAYER_OF_THM4302_RESIDUAL "
                     "ROWS "
                  << top.size() << " ENDPOINT " << maximum_endpoint
                  << " FNV " << std::hex << pair_fnv4303(top_keys)
                  << std::dec << '\n'
                  << "WORKERS " << worker_count
                  << " BODY_UNIVERSE_PER_ROW " << kBodyCount617
                  << " TOTAL_BODY_TESTS " << kBodyCount617 * top.size()
                  << '\n';
        for (std::size_t index = 0; index < top.size(); ++index) {
            const PairAgent pair = top[index];
            const PairAudit617& audit = audits[index];
            std::cout << "PAIR " << pair.q << ',' << pair.r << " ACTIVE "
                      << audit.active << " ACTIVE_FNV " << std::hex
                      << audit.active_fnv << std::dec << " ACTIVE_JOINT "
                      << audit.active_joint << " ACTIVE_NONJOINT "
                      << audit.active_nonjoint << " EXPOSED " << audit.exposed
                      << " HIT_RANGE " << audit.minimum_hits << ".."
                      << audit.maximum_hits << " HIT_INCIDENCES "
                      << audit.hit_incidences << " FAILURES "
                      << audit.failures << " FAILURE_FNV " << std::hex
                      << audit.failure_fnv << std::dec << '\n';
        }
        std::cout << "SUMMARY CLOSED_ROWS " << closed_rows
                  << " FAILED_ROWS " << failed_rows << " TOTAL_EXPOSED "
                  << total_exposed << " TOTAL_HIT_INCIDENCES "
                  << total_hit_incidences << " TOTAL_FAILURES "
                  << total_failures << " PAIR_LEDGER_FNV " << std::hex
                  << pair_ledger.state << " GLOBAL_FAILURE_FNV "
                  << global_failure_ledger.state << std::dec << '\n'
                  << "FAILURE_PAIRS 96,595:116:fedacdbff3f31981 "
                     "100,595:13:3ac9ac8b4b9ad93f "
                     "210,595:16:a6a226f12c168d3a\n"
                  << "SCOPE FINITE_EXACT_FIXED_C596_CARRIER_TOP595_ONLY_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACTLY_25_OF_28_TOP_ROWS_CLOSE\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4303_PRIMARY_RAW_AUDIT_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}

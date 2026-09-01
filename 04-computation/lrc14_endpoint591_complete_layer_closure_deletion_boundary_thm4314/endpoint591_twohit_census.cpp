// Exact singleton-deletion robustness and two-hit boundary for the thirteen
// endpoint-591 rows closed by the THM-4313 final carrier.

#define ENDPOINT593_RESPONSE_CAPACITY_MAIN endpoint593_hidden_for_twohit
#include "04-computation/lrc14_endpoint593_response_exchange_thm4311/endpoint593_response_capacity.cpp"
#undef ENDPOINT593_RESPONSE_CAPACITY_MAIN

namespace {

constexpr u64 kAddition593Fnv = UINT64_C(0x60873ef7a2b4ab90);
constexpr u64 kDeletion593Fnv = UINT64_C(0x4c14214a64ec202c);
constexpr u64 kExchangedCarrierFnv = UINT64_C(0xc9e5faef52ca5707);
constexpr u64 kCover43Fnv = UINT64_C(0xca3cb80f471f2e7e);
constexpr u64 kDeletion43Fnv = UINT64_C(0x0dd14ef0fe3eec62);
constexpr u64 kTop591Fnv = UINT64_C(0xfc332c0697c671c7);
constexpr u64 kFinalCarrier591Fnv = UINT64_C(0xa0d08a38c10bdab7);

std::vector<u32> read_cover591(const std::filesystem::path& path) {
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

std::vector<u32> read_deletions591(const std::filesystem::path& path) {
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

struct TwoHitObligation {
    u32 body = 0;
    u32 first = 0;
    u32 second = 0;
};

struct TwoHitRow {
    u64 exposed = 0;
    u64 below_two = 0;
    std::vector<TwoHitObligation> obligations;
};

TwoHitRow census_twohit(PairAgent pair, const std::vector<u32>& joint,
                        const std::vector<u32>& carrier,
                        const std::unordered_set<u32>& joint_set) {
    const Geometry geometry = build_geometry(pair.q, pair.r);
    std::vector<u32> active_joint, active_nonjoint;
    for (u32 mask : joint)
        if (margin(geometry, mask).ticks >= 0) active_joint.push_back(mask);
    for (u32 mask : carrier)
        if (!joint_set.contains(mask) && margin(geometry, mask).ticks >= 0)
            active_nonjoint.push_back(mask);
    TwoHitRow result;
    u64 bodies = 0;
    const u32 limit = u32{1} << 30;
    for (u32 body = (u32{1} << 9) - 1; body < limit;
         body = next_combination(body)) {
        ++bodies;
        bool joint_hit = false;
        for (u32 mask : active_joint)
            if ((mask & body) == 0) {
                joint_hit = true;
                break;
            }
        if (joint_hit) continue;
        ++result.exposed;
        std::array<u32, 3> witnesses{};
        unsigned hits = 0;
        for (u32 mask : active_nonjoint) {
            if ((mask & body) != 0) continue;
            witnesses[hits++] = mask;
            if (hits == 3) break;
        }
        if (hits < 2) {
            ++result.below_two;
        } else if (hits == 2) {
            if (witnesses[1] < witnesses[0])
                std::swap(witnesses[0], witnesses[1]);
            result.obligations.push_back(
                {body, witnesses[0], witnesses[1]});
        }
    }
    require(bodies == kBodyCount617, "body universe changed");
    return result;
}

}  // namespace

#ifndef ENDPOINT591_TWOHIT_MAIN
#define ENDPOINT591_TWOHIT_MAIN main
#endif

int ENDPOINT591_TWOHIT_MAIN(int argc, char** argv) {
    try {
        require(argc == 20,
                "usage: twohit JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE "
                "ADDITION593 DELETE593 TOP591 COVER43 DELETE43 CENSUS_CSV "
                "PAIR_CSV WORKERS");
        const unsigned worker_count = std::stoul(argv[19]);
        require(worker_count >= 1 && worker_count <= 64,
                "workers escaped range");
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> old_carrier = reconstruct_c3925_593(
            argv[2], argv[3], argv[4], argv[7], argv[8], argv[9], argv[10],
            argv[11]);
        (void)read_pairs_repaired(argv[5], 22647, kUniverseFnv);
        (void)read_pairs_repaired(argv[6], 1624, kOldUnionFnv);
        const auto add593 = read_flexible_masks593(
            argv[12], 1, kAddition593Fnv);
        const auto del593 = read_flexible_masks593(
            argv[13], 1, kDeletion593Fnv);
        std::vector<u32> carrier593;
        for (u32 mask : old_carrier)
            if (mask != del593.front()) carrier593.push_back(mask);
        carrier593.push_back(add593.front());
        require(carrier593.size() == 3925 &&
                    masks_fnv_agent(carrier593) == kExchangedCarrierFnv,
                "THM4311 carrier changed");
        const auto additions = read_cover591(argv[15]);
        const auto deletions = read_deletions591(argv[16]);
        require(masks_fnv_agent(additions) == kCover43Fnv &&
                    masks_fnv_agent(deletions) == kDeletion43Fnv,
                "THM4313 exchange changed");
        const std::set<u32> deletion_set(deletions.begin(), deletions.end());
        std::vector<u32> carrier;
        for (u32 mask : carrier593)
            if (!deletion_set.contains(mask)) carrier.push_back(mask);
        carrier.insert(carrier.end(), additions.begin(), additions.end());
        require(carrier.size() == 3925 &&
                    masks_fnv_agent(carrier) == kFinalCarrier591Fnv,
                "final carrier changed");

        const auto row_set = read_pairs_repaired(argv[14], 13, kTop591Fnv);
        std::vector<PairAgent> rows;
        for (const auto& [q, r] : row_set) rows.push_back({q, r});
        std::vector<TwoHitRow> audits(rows.size());
        std::atomic<std::size_t> cursor{0};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        auto worker = [&]() {
            try {
                while (true) {
                    const std::size_t index = cursor.fetch_add(1);
                    if (index >= rows.size()) return;
                    audits[index] =
                        census_twohit(rows[index], joint, carrier, joint_set);
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

        std::ofstream census_out(argv[17]);
        std::ofstream pair_out(argv[18]);
        require(census_out && pair_out, "cannot create census outputs");
        census_out << "q,r,body_hex,first_mask_hex,second_mask_hex\n";
        pair_out << "first_mask_hex,second_mask_hex,obligations,rows,body_fnv\n";
        std::map<std::pair<u32, u32>,
                 std::vector<std::tuple<int, int, u32>>> by_pair;
        Fnv census_ledger;
        u64 total_exposed = 0, total_twohit = 0, total_below_two = 0;
        for (std::size_t index = 0; index < rows.size(); ++index) {
            total_exposed += audits[index].exposed;
            total_below_two += audits[index].below_two;
            total_twohit += audits[index].obligations.size();
            for (const auto& obligation : audits[index].obligations) {
                census_out << rows[index].q << ',' << rows[index].r << ','
                           << hex8(obligation.body) << ','
                           << hex8(obligation.first) << ','
                           << hex8(obligation.second) << '\n';
                census_ledger.add(rows[index].q);
                census_ledger.add(rows[index].r);
                census_ledger.add(obligation.body);
                census_ledger.add(obligation.first);
                census_ledger.add(obligation.second);
                by_pair[{obligation.first, obligation.second}].push_back(
                    {rows[index].q, rows[index].r, obligation.body});
            }
        }
        std::size_t maximum_pair_load = 0;
        std::pair<u32, u32> least_maximum{};
        for (const auto& [pair, obligations] : by_pair) {
            Fnv body_ledger;
            std::set<std::pair<int, int>> pair_rows;
            for (const auto& [q, r, body] : obligations) {
                body_ledger.add(q); body_ledger.add(r); body_ledger.add(body);
                pair_rows.insert({q, r});
            }
            pair_out << hex8(pair.first) << ',' << hex8(pair.second) << ','
                     << obligations.size() << ',' << pair_rows.size() << ','
                     << std::hex << body_ledger.state << std::dec << '\n';
            if (obligations.size() > maximum_pair_load) {
                maximum_pair_load = obligations.size();
                least_maximum = pair;
            }
        }
        require(census_out.good() && pair_out.good(), "census write failed");
        require(total_exposed == 28791 && total_below_two == 0,
                "singleton robustness changed");
        std::cout << "LRC14_ENDPOINT591_TWOHIT_CENSUS_V1\n"
                  << "CARRIER 3925 FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec << " ROWS 13 EXPOSED "
                  << total_exposed << " BELOW_TWO " << total_below_two << '\n';
        for (std::size_t index = 0; index < rows.size(); ++index)
            std::cout << "ROW " << rows[index].q << ",591 EXPOSED "
                      << audits[index].exposed << " TWOHIT "
                      << audits[index].obligations.size() << '\n';
        std::cout << "TOTAL_TWOHIT " << total_twohit << " CENSUS_FNV "
                  << std::hex << census_ledger.state << std::dec
                  << " DISTINCT_WITNESS_PAIRS " << by_pair.size()
                  << " MAX_PAIR_LOAD " << maximum_pair_load
                  << " LEAST_MAX_PAIR " << hex8(least_maximum.first) << ','
                  << hex8(least_maximum.second) << '\n'
                  << "CONSEQUENCE ALL_NONJOINT_SINGLE_DELETIONS_SAFE_WITH_"
                     "JOINT_DECK_RETAINED_EACH_LISTED_NONJOINT_WITNESS_PAIR_"
                     "IS_A_HOSTILE_DOUBLE_DELETION\n"
                  << "SCOPE FIXED_CARRIER_ENDPOINT591_LAYER_PROTECTED_JOINT_"
                     "DELETION_MECHANISM_ONLY_NO_LRC14\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT591_TWOHIT_ERROR " << error.what() << '\n';
        return 1;
    }
}

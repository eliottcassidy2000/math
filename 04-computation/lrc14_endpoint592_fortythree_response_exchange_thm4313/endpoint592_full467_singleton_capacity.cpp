// Exact one-deletion capacity scout for an endpoint-592 response cover.
//
// The fixed carrier is THM-4311's size-3925 exchange.  We adjoin every mask
// in a supplied endpoint-592 cover, retain the 421 joint masks, and inspect
// the complete 432 inherited rows plus the 35 endpoint-592 rows.  The output
// is an exact singleton/private-obligation quotient only; it deliberately
// makes no simultaneous-deletion inference.

#define ENDPOINT593_RESPONSE_CAPACITY_MAIN endpoint593_capacity_hidden_for_592_capacity
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
constexpr u64 kAugmented43Fnv = UINT64_C(0xb8221545b2d668d0);

std::vector<u32> read_cover_masks(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open endpoint592 cover");
    std::string line;
    require(std::getline(input, line) && line.starts_with("mask_hex,"),
            "endpoint592 cover header changed");
    std::vector<u32> masks;
    std::set<u32> distinct;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto comma = line.find(',');
        require(comma != std::string::npos, "malformed endpoint592 cover row");
        const u32 mask = parse_mask_agent(line.substr(0, comma));
        require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                    distinct.insert(mask).second,
                "cover masks escaped ranks 8/9 or repeated");
        masks.push_back(mask);
    }
    require(input.eof() && !masks.empty(), "empty/malformed endpoint592 cover");
    return masks;
}

struct SingletonRow592 {
    PairAgent pair;
    u64 body_tests = 0;
    u64 exposed = 0;
    std::vector<std::pair<u32, u32>> singleton;  // body, sole witness
};

struct InheritedSingleton592 {
    PairAgent pair;
    u32 body = 0;
    u32 witness = 0;
};

std::vector<InheritedSingleton592> read_inherited_singletons592(
    const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open inherited singleton ledger");
    std::string line;
    require(std::getline(input, line) &&
                line == "q,r,body_hex,witness_hex,witness_rank",
            "inherited singleton header changed");
    std::vector<InheritedSingleton592> result;
    Fnv ledger;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        std::string body_token, witness_token;
        unsigned rank = 0;
        InheritedSingleton592 item;
        fields >> item.pair.q >> item.pair.r >> body_token >> witness_token >> rank;
        item.body = parse_mask_agent(body_token);
        item.witness = parse_mask_agent(witness_token);
        require(fields && std::popcount(item.body) == 9 &&
                    std::popcount(item.witness) == rank,
                "malformed inherited singleton row");
        result.push_back(item);
        ledger.add(item.pair.q); ledger.add(item.pair.r);
        ledger.add(item.body); ledger.add(item.witness);
    }
    require(input.eof() && result.size() == 1520 &&
                ledger.state == UINT64_C(0x39fedd8ceb347304),
            "inherited singleton identity changed");
    return result;
}

SingletonRow592 audit_singletons592(
    PairAgent pair, const std::vector<u32>& joint,
    const std::unordered_set<u32>& joint_set,
    const std::vector<u32>& augmented) {
    const Geometry geometry = build_geometry(pair.q, pair.r);
    std::vector<u32> active_joint;
    std::vector<u32> active_nonjoint;
    for (u32 mask : joint)
        if (margin(geometry, mask).ticks >= 0) active_joint.push_back(mask);
    for (u32 mask : augmented)
        if (!joint_set.contains(mask) && margin(geometry, mask).ticks >= 0)
            active_nonjoint.push_back(mask);

    SingletonRow592 result;
    result.pair = pair;
    const u32 limit = u32{1} << 30;
    for (u32 body = (u32{1} << 9) - 1; body < limit;
         body = next_combination(body)) {
        ++result.body_tests;
        bool joint_hit = false;
        for (u32 mask : active_joint)
            if ((mask & body) == 0) {
                joint_hit = true;
                break;
            }
        if (joint_hit) continue;
        ++result.exposed;
        unsigned hits = 0;
        u32 sole = 0;
        for (u32 mask : active_nonjoint)
            if ((mask & body) == 0) {
                ++hits;
                sole = mask;
                if (hits == 2) break;
            }
        require(hits != 0, "augmented carrier failed full467 target");
        if (hits == 1) result.singleton.push_back({body, sole});
    }
    require(result.body_tests == kBodyCount617, "body universe changed");
    return result;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 24,
                "usage: quotient JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE "
                "PREFIX391 TOP594 TOP593 ADDITION593 DELETE593 TOP592 "
                "COVER_CSV INHERITED_SINGLETON432 SINGLETON_CSV PROTECTED_TXT "
                "SAFE_STATS_CSV WORKERS");
        const unsigned worker_count = std::stoul(argv[23]);
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
        require(!joint_set.contains(add593.front()) &&
                    !joint_set.contains(del593.front()),
                "THM4311 exchange entered joint deck");

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

        const std::vector<u32> additions = read_cover_masks(argv[18]);
        require(additions.size() == 43 &&
                    masks_fnv_agent(additions) == kCover43Fnv,
                "frozen endpoint592 cover identity changed");
        std::vector<u32> augmented = carrier;
        for (u32 mask : additions) {
            require(!joint_set.contains(mask) && !carrier_set.contains(mask),
                    "endpoint592 addition overlaps carrier/joint deck");
            augmented.push_back(mask);
        }
        require(std::set<u32>(augmented.begin(), augmented.end()).size() ==
                    augmented.size(),
                "augmented carrier repeated a mask");
        require(augmented.size() == 3968 &&
                    masks_fnv_agent(augmented) == kAugmented43Fnv,
                "frozen augmented carrier identity changed");

        const auto prefix391 = read_pair_audit_rows593(argv[12], 391);
        const auto top594 = read_pair_audit_rows593(argv[13], 25);
        const auto top593 = read_pairs_repaired(argv[14], 16, kTop593Fnv593);
        const auto top592 = read_pairs_repaired(argv[17], 35, kTop592Fnv);
        std::set<std::pair<int, int>> target_set;
        for (PairAgent pair : prefix391)
            require(target_set.insert({pair.q, pair.r}).second,
                    "prefix391 overlap");
        for (PairAgent pair : top594)
            require(target_set.insert({pair.q, pair.r}).second,
                    "top594 overlap");
        for (const auto& pair : top593)
            require(target_set.insert(pair).second, "top593 overlap");
        for (const auto& pair : top592)
            require(target_set.insert(pair).second, "top592 overlap");
        require(target_set.size() == 467, "full target count changed");
        std::vector<PairAgent> rows;
        std::vector<Geometry> geometries;
        std::map<std::pair<int, int>, std::size_t> target_row_index;
        Fnv target_ledger;
        for (const auto& [q, r] : target_set) {
            target_row_index.emplace(std::pair{q, r}, rows.size());
            rows.push_back(PairAgent{q, r});
            geometries.push_back(build_geometry(q, r));
            target_ledger.add(q);
            target_ledger.add(r);
        }
        require(target_ledger.state == kFull467Fnv,
                "full467 target identity changed");

        // Adding masks cannot create a singleton on an inherited row.  Thus
        // the exact inherited singleton set is obtained by filtering the
        // THM-4311 ledger; only the 35 new frontier rows require raw scans.
        std::vector<PairAgent> frontier_rows;
        for (const auto& [q, r] : top592)
            frontier_rows.push_back(PairAgent{q, r});
        std::sort(frontier_rows.begin(), frontier_rows.end(),
                  [](PairAgent a, PairAgent b) {
                      return std::tie(a.q, a.r) < std::tie(b.q, b.r);
                  });
        std::vector<SingletonRow592> audits(frontier_rows.size());
        std::atomic<std::size_t> cursor{0};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        auto worker = [&]() {
            try {
                while (true) {
                    const std::size_t index = cursor.fetch_add(1);
                    if (index >= frontier_rows.size()) return;
                    audits[index] = audit_singletons592(
                        frontier_rows[index], joint, joint_set, augmented);
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

        std::ofstream singleton_out(argv[20]);
        std::ofstream protected_out(argv[21]);
        std::ofstream safe_stats_out(argv[22]);
        require(singleton_out && protected_out && safe_stats_out,
                "cannot create capacity ledgers");
        singleton_out << "q,r,body_hex,witness_hex,witness_rank,witness_origin\n";
        protected_out << "mask_hex,private_obligations,private_rows,origin\n";
        safe_stats_out << "mask_hex,rank,active_rows\n";

        const std::unordered_set<u32> additions_set(additions.begin(), additions.end());
        std::map<u32, std::pair<u64, std::set<std::size_t>>> private_by_mask;
        Fnv singleton_ledger;
        u64 body_tests = 0, exposed = 0, private_obligations = 0;
        u64 inherited_survivors = 0, inherited_resolved = 0;
        const auto inherited = read_inherited_singletons592(argv[19]);
        for (const auto& item : inherited) {
            const auto found =
                target_row_index.find({item.pair.q, item.pair.r});
            require(found != target_row_index.end() && item.pair.r != 592,
                    "inherited singleton escaped prior432 target");
            require(carrier_set.contains(item.witness) &&
                        !joint_set.contains(item.witness),
                    "inherited singleton witness escaped carrier nonjoint deck");
            bool resolved = false;
            for (u32 addition : additions) {
                if (margin(geometries[found->second], addition).ticks >= 0 &&
                    (addition & item.body) == 0) {
                    resolved = true;
                    break;
                }
            }
            if (resolved) {
                ++inherited_resolved;
                continue;
            }
            ++inherited_survivors;
            ++private_obligations;
            auto& entry = private_by_mask[item.witness];
            ++entry.first;
            entry.second.insert(found->second);
            singleton_ledger.add(item.pair.q); singleton_ledger.add(item.pair.r);
            singleton_ledger.add(item.body); singleton_ledger.add(item.witness);
            singleton_out << item.pair.q << ',' << item.pair.r << ','
                          << hex8(item.body) << ',' << hex8(item.witness) << ','
                          << std::popcount(item.witness) << ",carrier\n";
        }
        for (std::size_t audit_index = 0; audit_index < audits.size();
             ++audit_index) {
            const auto& audit = audits[audit_index];
            body_tests += audit.body_tests;
            exposed += audit.exposed;
            const auto target_index =
                target_row_index.find({audit.pair.q, audit.pair.r});
            require(target_index != target_row_index.end(),
                    "frontier singleton row escaped target");
            for (const auto& [body, witness] : audit.singleton) {
                ++private_obligations;
                auto& item = private_by_mask[witness];
                ++item.first;
                item.second.insert(target_index->second);
                singleton_ledger.add(audit.pair.q);
                singleton_ledger.add(audit.pair.r);
                singleton_ledger.add(body);
                singleton_ledger.add(witness);
                singleton_out << audit.pair.q << ',' << audit.pair.r << ','
                              << hex8(body) << ',' << hex8(witness) << ','
                              << std::popcount(witness) << ','
                              << (additions_set.contains(witness) ? "addition" :
                                  (carrier_set.contains(witness) ? "carrier" :
                                   "unexpected"))
                              << '\n';
            }
        }
        Fnv protected_ledger;
        std::set<u32> protected_set;
        u64 protected_carrier = 0, protected_addition = 0;
        for (const auto& [mask, item] : private_by_mask) {
            protected_set.insert(mask);
            protected_ledger.add(mask);
            const bool is_addition = additions_set.contains(mask);
            protected_addition += is_addition;
            protected_carrier += !is_addition && carrier_set.contains(mask);
            protected_out << hex8(mask) << ',' << item.first << ','
                          << item.second.size() << ','
                          << (is_addition ? "addition" : "carrier") << '\n';
        }

        struct SafeStat { u32 mask; unsigned rank; unsigned active_rows; };
        std::vector<SafeStat> safe;
        for (u32 mask : carrier) {
            if (joint_set.contains(mask) || protected_set.contains(mask)) continue;
            unsigned active_rows = 0;
            for (const Geometry& geometry : geometries)
                active_rows += margin(geometry, mask).ticks >= 0;
            safe.push_back({mask, static_cast<unsigned>(std::popcount(mask)),
                            active_rows});
        }
        std::sort(safe.begin(), safe.end(), [](const SafeStat& left,
                                               const SafeStat& right) {
            return std::tie(left.active_rows, left.rank, left.mask) <
                   std::tie(right.active_rows, right.rank, right.mask);
        });
        Fnv safe_ledger;
        for (const auto& item : safe) {
            safe_ledger.add(item.mask);
            safe_stats_out << hex8(item.mask) << ',' << item.rank << ','
                           << item.active_rows << '\n';
        }
        require(singleton_out.good() && protected_out.good() &&
                    safe_stats_out.good(),
                "capacity ledger write failed");

        require(body_tests == UINT64_C(500750250) && exposed == 709299,
                "endpoint592 raw quotient census changed");
        require(inherited_survivors == 668 && inherited_resolved == 852,
                "inherited singleton filter changed");
        require(private_obligations == 1510 &&
                    singleton_ledger.state == UINT64_C(0x4dfcadbd1d5c0c31),
                "full467 singleton ledger changed");
        require(private_by_mask.size() == 412 &&
                    protected_ledger.state == UINT64_C(0x6b6559e5bf1a529d) &&
                    protected_carrier == 369 && protected_addition == 43,
                "full467 protected-mask quotient changed");
        require(safe.size() == 3135 &&
                    safe_ledger.state == UINT64_C(0xbd2d6c6eabf2203b),
                "full467 singleton-safe carrier ledger changed");

        const u64 nonjoint_total = augmented.size() - joint.size();
        std::cout << "LRC14_ENDPOINT592_FULL467_SINGLETON_CAPACITY_V1\n"
                  << "BASE_EXCHANGED_CARRIER " << carrier.size() << " FNV "
                  << std::hex << masks_fnv_agent(carrier) << std::dec
                  << " COVER_ADDITIONS " << additions.size() << " COVER_FNV "
                  << std::hex << masks_fnv_agent(additions) << std::dec
                  << " AUGMENTED " << augmented.size() << " AUGMENTED_FNV "
                  << std::hex << masks_fnv_agent(augmented) << std::dec << '\n'
                  << "TARGET_ROWS " << rows.size() << " TARGET_FNV "
                  << std::hex << target_ledger.state << std::dec
                  << " RAW_FRONTIER_ROWS " << frontier_rows.size()
                  << " RAW_BODY_TESTS " << body_tests << " RAW_EXPOSED "
                  << exposed
                  << '\n'
                  << "INHERITED_SINGLETONS " << inherited.size()
                  << " INHERITED_SURVIVORS " << inherited_survivors
                  << " INHERITED_RESOLVED_BY_COVER " << inherited_resolved
                  << '\n'
                  << "PRIVATE_OBLIGATIONS " << private_obligations
                  << " SINGLETON_FNV " << std::hex << singleton_ledger.state
                  << std::dec << " PROTECTED_TOTAL " << private_by_mask.size()
                  << " PROTECTED_FNV " << std::hex << protected_ledger.state
                  << std::dec << " PROTECTED_CARRIER " << protected_carrier
                  << " PROTECTED_ADDITIONS " << protected_addition << '\n'
                  << "NONJOINT_TOTAL " << nonjoint_total
                  << " OLD_NONJOINT " << carrier.size() - joint.size()
                  << " OLD_SINGLETON_SAFE " << safe.size() << " SAFE_FNV "
                  << std::hex << safe_ledger.state << std::dec << '\n'
                  << "SCOPE EXACT_SINGLE_DELETION_CAPACITY_ONLY_FIXED_FULL467_"
                     "TARGET_JOINT_DECK_RETAINED_NO_SIMULTANEOUS_INFERENCE_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\nVERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT592_CAPACITY_ERROR " << error.what() << '\n';
        return 1;
    }
}

// Exact protected-joint singleton-deletion quotient after adjoining a
// nine-mask endpoint-590 response cover.  For the inherited full467 layer,
// only rows whose frozen final-carrier nonjoint minimum is one need raw scans;
// adding masks cannot create a new singleton.  Endpoint591/590 use their
// complete frozen low-multiplicity ledgers. Scratch only.

#define ENDPOINT591_TWOHIT_MAIN endpoint591_hidden_for_endpoint590_quotient
#include "04-computation/lrc14_endpoint591_complete_layer_closure_deletion_boundary_thm4314/endpoint591_twohit_census.cpp"
#undef ENDPOINT591_TWOHIT_MAIN

namespace {

constexpr u64 kFull467Fnv590 = UINT64_C(0x2d6aa988098aa5eb);
constexpr u64 kTop591Fnv590 = UINT64_C(0xfc332c0697c671c7);
constexpr u64 kTop590FnvQuotient = UINT64_C(0x44aa8a793d162cf9);
constexpr u64 kLow591Fnv590 = UINT64_C(0x8dad3ce6f3e587c2);
constexpr u64 kLow590FnvQuotient = UINT64_C(0xa5fc0bfbe9f4cc13);

struct LowItemQuotient {
    PairAgent pair{};
    u32 body = 0;
    unsigned hits = 0;
    u32 first = 0;
    u32 second = 0;
};

struct PrivateItemQuotient {
    PairAgent pair{};
    u32 body = 0;
    u32 witness = 0;
    std::string source;
};

struct ScannedRowQuotient {
    PairAgent pair{};
    u64 bodies = 0;
    u64 exposed = 0;
    std::vector<PrivateItemQuotient> singleton;
};

std::vector<u32> read_cover_quotient(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open endpoint590 cover");
    std::string line;
    require(std::getline(input, line) && line.starts_with("mask_hex,"),
            "endpoint590 cover header changed");
    std::vector<u32> result;
    std::set<u32> distinct;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto comma = line.find(',');
        require(comma != std::string::npos, "malformed endpoint590 cover row");
        const u32 mask = parse_mask_agent(line.substr(0, comma));
        require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                    distinct.insert(mask).second,
                "cover mask rank/distinctness changed");
        result.push_back(mask);
    }
    require(input.eof() && result.size() == 9,
            "endpoint590 cover size changed");
    return result;
}

std::vector<u32> reconstruct_final_carrier_quotient(char** argv,
                                                     const std::vector<u32>& joint) {
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
    carrier593.push_back(add593.front());
    require(carrier593.size() == 3925 &&
                masks_fnv_agent(carrier593) == kExchangedCarrierFnv,
            "THM4311 carrier changed");
    const auto additions = read_cover591(argv[14]);
    const auto deletions = read_deletions591(argv[15]);
    require(masks_fnv_agent(additions) == kCover43Fnv &&
                masks_fnv_agent(deletions) == kDeletion43Fnv,
            "THM4313 ledgers changed");
    const std::set<u32> deletion_set(deletions.begin(), deletions.end());
    std::vector<u32> carrier;
    for (u32 mask : carrier593)
        if (!deletion_set.contains(mask)) carrier.push_back(mask);
    carrier.insert(carrier.end(), additions.begin(), additions.end());
    require(carrier.size() == 3925 &&
                masks_fnv_agent(carrier) == kFinalCarrier591Fnv,
            "THM4313 final carrier changed");
    for (u32 mask : joint)
        require(std::find(carrier.begin(), carrier.end(), mask) != carrier.end(),
                "joint deck changed");
    return carrier;
}

std::vector<PairAgent> read_full467_minimum_one(
    const std::filesystem::path& path, std::vector<PairAgent>& all_rows) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open full467 pair audit");
    std::string line;
    require(std::getline(input, line) && line.starts_with("q,r,active,"),
            "full467 pair header changed");
    std::vector<PairAgent> scan_rows;
    std::set<std::pair<int, int>> distinct;
    Fnv row_ledger;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        PairAgent pair;
        u64 active = 0, active_joint = 0, active_nonjoint = 0, exposed = 0;
        u64 minimum = 0, maximum = 0, failures = 0;
        std::string active_fnv, exposed_fnv, failure_fnv;
        fields >> pair.q >> pair.r >> active >> active_fnv >> active_joint
               >> active_nonjoint >> exposed >> exposed_fnv >> minimum
               >> maximum >> failures >> failure_fnv;
        require(fields && failures == 0 && minimum >= 1 &&
                    distinct.insert({pair.q, pair.r}).second,
                "malformed/hostile full467 pair row");
        all_rows.push_back(pair);
        row_ledger.add(pair.q); row_ledger.add(pair.r);
        if (minimum == 1) scan_rows.push_back(pair);
    }
    require(input.eof() && all_rows.size() == 467 && scan_rows.size() == 42 &&
                row_ledger.state == kFull467Fnv590,
            "full467 row/minimum identity changed");
    return scan_rows;
}

std::vector<LowItemQuotient> read_low_quotient(
    const std::filesystem::path& path, u64 expected_fnv,
    std::size_t expected_rows, int endpoint) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open low-witness ledger");
    std::string line;
    require(std::getline(input, line) &&
                line == "q,r,body_hex,hits,first_mask_hex,second_mask_hex",
            "low-witness header changed");
    std::vector<LowItemQuotient> result;
    Fnv ledger;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::vector<std::string> fields;
        std::size_t start = 0;
        while (true) {
            const auto comma = line.find(',', start);
            fields.push_back(line.substr(start, comma - start));
            if (comma == std::string::npos) break;
            start = comma + 1;
        }
        require(fields.size() == 6, "malformed low-witness row");
        LowItemQuotient item;
        item.pair.q = std::stoi(fields[0]);
        item.pair.r = std::stoi(fields[1]);
        item.body = parse_mask_agent(fields[2]);
        item.hits = std::stoul(fields[3]);
        item.first = fields[4].empty() ? 0 : parse_mask_agent(fields[4]);
        item.second = fields[5].empty() ? 0 : parse_mask_agent(fields[5]);
        require(item.pair.r == endpoint && std::popcount(item.body) == 9 &&
                    item.hits <= 2,
                "low-witness item escaped endpoint/range");
        result.push_back(item);
        ledger.add(item.pair.q); ledger.add(item.pair.r);
        ledger.add(item.body); ledger.add(item.hits);
        ledger.add(item.first); ledger.add(item.second);
    }
    require(input.eof() && result.size() == expected_rows &&
                ledger.state == expected_fnv,
            "low-witness ledger identity changed");
    return result;
}

bool addition_resolves_quotient(const LowItemQuotient& item,
                                const std::vector<u32>& additions) {
    const Geometry geometry = build_geometry(item.pair.q, item.pair.r);
    for (u32 mask : additions)
        if (margin(geometry, mask).ticks >= 0 && (mask & item.body) == 0)
            return true;
    return false;
}

ScannedRowQuotient scan_private_quotient(
    PairAgent pair, const std::vector<u32>& joint,
    const std::unordered_set<u32>& joint_set,
    const std::vector<u32>& augmented) {
    const Geometry geometry = build_geometry(pair.q, pair.r);
    std::vector<u32> active_joint, active_nonjoint;
    for (u32 mask : joint)
        if (margin(geometry, mask).ticks >= 0) active_joint.push_back(mask);
    for (u32 mask : augmented)
        if (!joint_set.contains(mask) && margin(geometry, mask).ticks >= 0)
            active_nonjoint.push_back(mask);
    ScannedRowQuotient result;
    result.pair = pair;
    const u32 limit = u32{1} << 30;
    for (u32 body = (u32{1} << 9) - 1; body < limit;
         body = next_combination(body)) {
        ++result.bodies;
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
                sole = mask;
                if (++hits == 2) break;
            }
        require(hits != 0, "augmented carrier failed inherited row");
        if (hits == 1)
            result.singleton.push_back({pair, body, sole, "full467_raw"});
    }
    require(result.bodies == kBodyCount617, "body universe changed");
    return result;
}

}  // namespace

#ifndef ENDPOINT590_QUOTIENT_MAIN
#define ENDPOINT590_QUOTIENT_MAIN main
#endif

int ENDPOINT590_QUOTIENT_MAIN(int argc, char** argv) {
    try {
        require(argc == 27,
                "usage: quotient590 JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE "
                "ADDITION593 DELETE593 COVER43 DELETE43 FULL467_PAIR TOP591 "
                "LOW591 TOP590 LOW590 COVER9 SINGLETON_OUT PROTECTED_OUT "
                "SAFE_OUT DELETE9_OUT WORKERS");
        const unsigned worker_count = std::stoul(argv[26]);
        require(worker_count >= 1 && worker_count <= 64,
                "workers escaped 1..64");
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> carrier =
            reconstruct_final_carrier_quotient(argv, joint);
        const std::unordered_set<u32> carrier_set(carrier.begin(), carrier.end());
        const std::vector<u32> additions = read_cover_quotient(argv[21]);
        std::vector<u32> augmented = carrier;
        for (u32 mask : additions) {
            require(!carrier_set.contains(mask) && !joint_set.contains(mask),
                    "endpoint590 addition overlaps carrier/joint deck");
            augmented.push_back(mask);
        }
        require(augmented.size() == 3934 &&
                    std::set<u32>(augmented.begin(), augmented.end()).size() ==
                        augmented.size(),
                "augmented carrier size/distinctness changed");

        std::vector<PairAgent> full467;
        const auto scan_rows = read_full467_minimum_one(argv[16], full467);
        const auto top591_pairs =
            read_pairs_repaired(argv[17], 13, kTop591Fnv590);
        const auto top590_pairs =
            read_pairs_repaired(argv[19], 13, kTop590FnvQuotient);
        std::vector<PairAgent> rows480 = full467;
        for (const auto& [q, r] : top591_pairs) rows480.push_back({q, r});
        std::sort(rows480.begin(), rows480.end(), [](PairAgent left,
                                                     PairAgent right) {
            return std::tie(left.q, left.r) < std::tie(right.q, right.r);
        });
        require(std::adjacent_find(rows480.begin(), rows480.end(),
                    [](PairAgent left, PairAgent right) {
                        return std::tie(left.q, left.r) ==
                               std::tie(right.q, right.r);
                    }) == rows480.end() && rows480.size() == 480,
                "inherited480 overlap/count changed");
        Fnv rows480_ledger;
        for (PairAgent pair : rows480) {
            rows480_ledger.add(pair.q); rows480_ledger.add(pair.r);
        }
        std::vector<PairAgent> rows493 = rows480;
        for (const auto& [q, r] : top590_pairs) rows493.push_back({q, r});
        std::sort(rows493.begin(), rows493.end(), [](PairAgent left,
                                                     PairAgent right) {
            return std::tie(left.q, left.r) < std::tie(right.q, right.r);
        });
        require(std::adjacent_find(rows493.begin(), rows493.end(),
                    [](PairAgent left, PairAgent right) {
                        return std::tie(left.q, left.r) ==
                               std::tie(right.q, right.r);
                    }) == rows493.end() && rows493.size() == 493,
                "target493 overlap/count changed");
        Fnv rows493_ledger;
        for (PairAgent pair : rows493) {
            rows493_ledger.add(pair.q); rows493_ledger.add(pair.r);
        }

        std::vector<ScannedRowQuotient> audits(scan_rows.size());
        std::atomic<std::size_t> cursor{0};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        auto worker = [&]() {
            try {
                while (true) {
                    const std::size_t index = cursor.fetch_add(1);
                    if (index >= scan_rows.size()) return;
                    audits[index] = scan_private_quotient(
                        scan_rows[index], joint, joint_set, augmented);
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

        std::vector<PrivateItemQuotient> private_items;
        u64 raw_bodies = 0, raw_exposed = 0;
        for (const auto& audit : audits) {
            raw_bodies += audit.bodies;
            raw_exposed += audit.exposed;
            private_items.insert(private_items.end(), audit.singleton.begin(),
                                 audit.singleton.end());
        }
        require(raw_bodies == scan_rows.size() * kBodyCount617,
                "raw quotient body total changed");

        const auto low591 = read_low_quotient(argv[18], kLow591Fnv590, 28, 591);
        const auto low590 = read_low_quotient(argv[20], kLow590FnvQuotient,
                                              2730, 590);
        u64 resolved591 = 0, retained591 = 0;
        for (const auto& item : low591) {
            if (item.hits != 1) continue;
            require(carrier_set.contains(item.first),
                    "endpoint591 sole witness escaped carrier");
            if (addition_resolves_quotient(item, additions)) ++resolved591;
            else {
                ++retained591;
                if (!joint_set.contains(item.first))
                    private_items.push_back(
                        {item.pair, item.body, item.first, "endpoint591_low"});
            }
        }
        u64 covered_failures590 = 0, resolved_one590 = 0, retained_one590 = 0;
        u64 retained_joint590 = 0, retained_nonjoint590 = 0;
        for (const auto& item : low590) {
            const bool resolves = addition_resolves_quotient(item, additions);
            if (item.hits == 0) {
                require(resolves, "cover failed frozen endpoint590 body");
                ++covered_failures590;
            } else if (item.hits == 1) {
                require(carrier_set.contains(item.first),
                        "endpoint590 sole witness escaped carrier");
                if (resolves) ++resolved_one590;
                else {
                    ++retained_one590;
                    if (joint_set.contains(item.first)) {
                        ++retained_joint590;
                    } else {
                        ++retained_nonjoint590;
                        private_items.push_back(
                            {item.pair, item.body, item.first,
                             "endpoint590_low"});
                    }
                }
            }
        }
        require(covered_failures590 == 100,
                "cover did not close exact endpoint590 failure set");

        std::sort(private_items.begin(), private_items.end(),
                  [](const auto& left, const auto& right) {
            return std::tie(left.pair.q, left.pair.r, left.body, left.witness) <
                   std::tie(right.pair.q, right.pair.r, right.body, right.witness);
        });
        std::ofstream singleton_out(argv[22]);
        std::ofstream protected_out(argv[23]);
        std::ofstream safe_out(argv[24]);
        std::ofstream deletion_out(argv[25]);
        require(singleton_out && protected_out && safe_out && deletion_out,
                "cannot create quotient outputs");
        singleton_out << "q,r,body_hex,witness_hex,witness_rank,source\n";
        Fnv singleton_ledger;
        std::map<u32, std::pair<u64, std::set<std::pair<int, int>>>> private_by_mask;
        for (const auto& item : private_items) {
            require(carrier_set.contains(item.witness) &&
                        !joint_set.contains(item.witness),
                    "protected nonjoint private witness escaped old carrier");
            singleton_out << item.pair.q << ',' << item.pair.r << ','
                          << hex8(item.body) << ',' << hex8(item.witness) << ','
                          << std::popcount(item.witness) << ',' << item.source
                          << '\n';
            singleton_ledger.add(item.pair.q); singleton_ledger.add(item.pair.r);
            singleton_ledger.add(item.body); singleton_ledger.add(item.witness);
            auto& entry = private_by_mask[item.witness];
            ++entry.first;
            entry.second.insert({item.pair.q, item.pair.r});
        }

        protected_out << "mask_hex,rank,private_obligations,private_rows\n";
        Fnv protected_ledger;
        for (const auto& [mask, data] : private_by_mask) {
            protected_ledger.add(mask);
            protected_out << hex8(mask) << ',' << std::popcount(mask) << ','
                          << data.first << ',' << data.second.size() << '\n';
        }

        std::vector<Geometry> geometry480, geometry493;
        for (PairAgent pair : rows480)
            geometry480.push_back(build_geometry(pair.q, pair.r));
        for (PairAgent pair : rows493)
            geometry493.push_back(build_geometry(pair.q, pair.r));
        struct SafeStatQuotient {
            u32 mask = 0;
            unsigned rank = 0;
            unsigned active480 = 0;
            unsigned active493 = 0;
        };
        std::vector<SafeStatQuotient> safe;
        for (u32 mask : carrier) {
            if (joint_set.contains(mask) || private_by_mask.contains(mask))
                continue;
            SafeStatQuotient item;
            item.mask = mask;
            item.rank = std::popcount(mask);
            for (const Geometry& geometry : geometry480)
                item.active480 += margin(geometry, mask).ticks >= 0;
            for (const Geometry& geometry : geometry493)
                item.active493 += margin(geometry, mask).ticks >= 0;
            safe.push_back(item);
        }
        std::sort(safe.begin(), safe.end(), [](const auto& left,
                                               const auto& right) {
            return std::tie(left.active493, left.active480, left.rank, left.mask) <
                   std::tie(right.active493, right.active480, right.rank, right.mask);
        });
        require(safe.size() >= 9, "fewer than nine singleton-safe old masks");
        safe_out << "mask_hex,rank,active_rows480,active_rows493\n";
        Fnv safe_ledger;
        for (const auto& item : safe) {
            safe_ledger.add(item.mask);
            safe_out << hex8(item.mask) << ',' << item.rank << ','
                     << item.active480 << ',' << item.active493 << '\n';
        }
        std::vector<u32> deletion;
        for (std::size_t index = 0; index < 9; ++index) {
            deletion.push_back(safe[index].mask);
            deletion_out << hex8(safe[index].mask) << '\n';
        }
        const u64 deletion_fnv = masks_fnv_agent(deletion);
        require(singleton_out.good() && protected_out.good() && safe_out.good()
                    && deletion_out.good(),
                "quotient output write failed");

        std::cout << "LRC14_ENDPOINT590_PROTECTED_DELETION_QUOTIENT_V1\n"
                  << "BASE_CARRIER 3925 FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec
                  << " COVER9_FNV " << std::hex
                  << masks_fnv_agent(additions) << std::dec
                  << " AUGMENTED " << augmented.size() << " AUGMENTED_FNV "
                  << std::hex << masks_fnv_agent(augmented) << std::dec << '\n'
                  << "INHERITED_ROWS 480 ROW_FNV " << std::hex
                  << rows480_ledger.state << std::dec << " TARGET_ROWS 493 ROW_FNV "
                  << std::hex << rows493_ledger.state << std::dec << '\n'
                  << "FULL467_RAW_ROWS " << scan_rows.size()
                  << " RAW_BODY_TESTS " << raw_bodies << " RAW_EXPOSED "
                  << raw_exposed << '\n'
                  << "ENDPOINT591_SOLE 2 RESOLVED " << resolved591
                  << " RETAINED " << retained591 << '\n'
                  << "ENDPOINT590_FAILURES_COVERED " << covered_failures590
                  << " ONE_WITNESS 612 RESOLVED " << resolved_one590
                  << " RETAINED " << retained_one590 << " RETAINED_JOINT "
                  << retained_joint590 << " RETAINED_NONJOINT "
                  << retained_nonjoint590 << '\n'
                  << "PRIVATE_NONJOINT_OBLIGATIONS " << private_items.size()
                  << " SINGLETON_FNV " << std::hex << singleton_ledger.state
                  << std::dec << " PROTECTED_OLD_MASKS "
                  << private_by_mask.size() << " PROTECTED_FNV " << std::hex
                  << protected_ledger.state << std::dec << '\n'
                  << "OLD_NONJOINT " << carrier.size() - joint.size()
                  << " SINGLETON_SAFE_OLD " << safe.size() << " SAFE_FNV "
                  << std::hex << safe_ledger.state << std::dec << '\n'
                  << "LOW_ACTIVITY_DELETE9_FNV " << std::hex << deletion_fnv
                  << std::dec << " MASKS";
        for (u32 mask : deletion) std::cout << ' ' << hex8(mask);
        std::cout << '\n'
                  << "SCOPE EXACT_SINGLE_NONJOINT_DELETION_QUOTIENT_WITH_"
                     "ALL_421_JOINT_MASKS_RETAINED_FIXED_AUGMENTED_TARGET493_"
                     "NO_SIMULTANEOUS_DELETION_INFERENCE_NO_PHYSICAL_ENTRY_"
                     "NO_LRC14\nVERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT590_QUOTIENT_ERROR " << error.what() << '\n';
        return 1;
    }
}

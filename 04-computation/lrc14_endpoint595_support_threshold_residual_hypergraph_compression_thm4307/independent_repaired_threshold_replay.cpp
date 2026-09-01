// Independent verifier for the repaired endpoint-595 support-threshold
// deletion.  It reconstructs D350, audits the hostile residual, verifies the
// exact H37 cover/packing, checks Dfinal=D350\H37, and replays every body whose
// witness set is changed by the final simultaneous deletion.

#define ENDPOINT617_RAW_VERIFY_MAIN threshold_replay_hidden_main
#include "04-computation/lrc14_size_preserving_response_staircase_thm4300/endpoint617_exchanged_carrier_raw_verify.cpp"
#undef ENDPOINT617_RAW_VERIFY_MAIN

#include <atomic>
#include <fstream>
#include <mutex>
#include <thread>

namespace {

constexpr u64 kRepairFnvThr = UINT64_C(0x64ce5f9d1ec8c4c2);
constexpr u64 kAddition4FnvThr = UINT64_C(0xdc0eebaebf688c65);
constexpr u64 kDelete73FnvThr = UINT64_C(0x9240b264ab65aa62);
constexpr u64 kAddition10FnvThr = UINT64_C(0x6740cc137170afc5);
constexpr u64 kAugmentedFnvThr = UINT64_C(0x55e8588798885ae5);
constexpr u64 kC596FnvThr = UINT64_C(0x892fef44a9e6b37e);
constexpr u64 kOldUnionFnvThr = UINT64_C(0x11414a33ab91fef6);

using PairKeyThr = std::pair<int, int>;

std::set<PairKeyThr> read_pair_set_thr(const std::filesystem::path& path,
                                       std::size_t expected,
                                       u64 expected_fnv) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open pair set");
    std::set<PairKeyThr> rows;
    Fnv ledger;
    std::string line;
    PairKeyThr previous{0, 0};
    bool have_previous = false;
    while (std::getline(input, line)) {
        require(!line.empty(), "blank pair row");
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos &&
                    line.find(',', comma + 1) == std::string::npos,
                "malformed pair row");
        const PairKeyThr pair{std::stoi(line.substr(0, comma)),
                              std::stoi(line.substr(comma + 1))};
        require(pair.first > 0 && pair.first < pair.second,
                "invalid pair row");
        require(!have_previous || previous < pair, "pair order changed");
        require(rows.insert(pair).second, "duplicate pair row");
        previous = pair;
        have_previous = true;
        ledger.add(pair.first);
        ledger.add(pair.second);
    }
    require(input.eof() && rows.size() == expected &&
                ledger.state == expected_fnv,
            "pair-set identity changed");
    return rows;
}

std::vector<u32> reconstruct_c596_thr(
    const std::filesystem::path& base,
    const std::filesystem::path& add45,
    const std::filesystem::path& suffix9,
    const std::filesystem::path& repairs76,
    const std::filesystem::path& additions4,
    const std::filesystem::path& deletions73) {
    std::vector<u32> augmented = build_mixed_carrier(base, add45, suffix9);
    const std::vector<u32> repairs =
        read_mixed617(repairs76, 76, kRepairFnvThr);
    const std::vector<u32> additions =
        read_mixed617(additions4, 4, kAddition4FnvThr);
    const std::vector<u32> deletions =
        read_mixed617(deletions73, 73, kDelete73FnvThr);
    std::set<u32> distinct(augmented.begin(), augmented.end());
    for (u32 mask : repairs) {
        require(distinct.insert(mask).second, "repair overlap");
        augmented.push_back(mask);
    }
    require(augmented.size() == 9088 &&
                masks_fnv_agent(augmented) == kAugmentedFnvThr,
            "augmented identity changed");
    const std::set<u32> deleted(deletions.begin(), deletions.end());
    std::vector<u32> carrier;
    for (u32 mask : augmented)
        if (!deleted.contains(mask)) carrier.push_back(mask);
    for (u32 mask : additions) {
        require(distinct.insert(mask).second, "addition4 overlap");
        carrier.push_back(mask);
    }
    require(carrier.size() == 9019 &&
                masks_fnv_agent(carrier) == kC596FnvThr,
            "C596 identity changed");
    return carrier;
}

std::vector<u32> read_threshold_ledger(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open threshold ledger");
    std::vector<u32> masks;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask_agent(token);
        require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                    distinct.insert(mask).second,
                "threshold ledger rank/distinctness changed");
        masks.push_back(mask);
    }
    require(input.eof(), "threshold ledger parse failed");
    return masks;
}

struct LocalAuditThr {
    PairAgent pair;
    std::vector<u32> removed_active;
    u64 tested = 0;
    u64 failures = 0;
    u64 failure_fnv = kEmptyFnv;
    std::vector<u32> failure_bodies;
    u64 active_survivors = 0;
    u64 active_survivor_fnv = kEmptyFnv;
};

LocalAuditThr audit_local_threshold(
    PairAgent pair, std::vector<u32> removed_active,
    const std::vector<u32>& joint, const std::unordered_set<u32>& joint_set,
    const std::vector<u32>& final_carrier) {
    const Geometry geometry = build_geometry(pair.q, pair.r);
    std::vector<u32> active_joint;
    std::vector<u32> active_nonjoint;
    Fnv active_ledger;
    for (u32 mask : joint)
        if (margin(geometry, mask).ticks >= 0) active_joint.push_back(mask);
    for (u32 mask : final_carrier) {
        if (margin(geometry, mask).ticks < 0) continue;
        active_ledger.add(mask);
        if (!joint_set.contains(mask)) active_nonjoint.push_back(mask);
    }
    for (u32 mask : removed_active)
        require(margin(geometry, mask).ticks >= 0,
                "row group contains inactive mask");

    LocalAuditThr audit;
    audit.pair = pair;
    audit.removed_active = std::move(removed_active);
    audit.active_survivors = active_joint.size() + active_nonjoint.size();
    audit.active_survivor_fnv = active_ledger.state;
    Fnv failure_ledger;
    const u32 limit = u32{1} << 30;
    for (u32 body = (u32{1} << 9) - 1; body < limit;
         body = next_combination(body)) {
        bool removed_hit = false;
        for (u32 mask : audit.removed_active)
            if ((mask & body) == 0) {
                removed_hit = true;
                break;
            }
        if (!removed_hit) continue;
        ++audit.tested;
        bool survivor_hit = false;
        for (u32 mask : active_joint)
            if ((mask & body) == 0) {
                survivor_hit = true;
                break;
            }
        if (!survivor_hit)
            for (u32 mask : active_nonjoint)
                if ((mask & body) == 0) {
                    survivor_hit = true;
                    break;
                }
        if (!survivor_hit) {
            ++audit.failures;
            failure_ledger.add(body);
            audit.failure_bodies.push_back(body);
        }
    }
    audit.failure_fnv = failure_ledger.state;
    return audit;
}

struct ObligationThr {
    PairAgent pair;
    u32 body = 0;
    bool operator==(const ObligationThr& other) const {
        return pair.q == other.pair.q && pair.r == other.pair.r &&
               body == other.body;
    }
};

std::vector<ObligationThr> read_obligations_thr(
    const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open hostile failure ledger");
    std::string line;
    require(std::getline(input, line) && line == "q,r,body_hex",
            "hostile failure header changed");
    std::vector<ObligationThr> obligations;
    std::set<std::tuple<int, int, u32>> distinct;
    std::tuple<int, int, u32> previous{};
    bool have_previous = false;
    while (std::getline(input, line)) {
        require(!line.empty(), "blank hostile failure row");
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        int q = 0, r = 0;
        std::string token;
        fields >> q >> r >> token;
        const u32 body = parse_mask_agent(token);
        const auto key = std::tuple{q, r, body};
        require(fields && q > 0 && q < r && std::popcount(body) == 9 &&
                    distinct.insert(key).second &&
                    (!have_previous || previous < key),
                "hostile failure row invalid/out of order");
        obligations.push_back({PairAgent{q, r}, body});
        previous = key;
        have_previous = true;
    }
    require(input.eof() && obligations.size() == 84,
            "hostile failure count changed");
    return obligations;
}

std::vector<std::size_t> read_packing_thr(
    const std::filesystem::path& path, std::size_t obligation_count) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open packing ledger");
    std::vector<std::size_t> support;
    std::set<std::size_t> distinct;
    std::size_t index = 0;
    while (input >> index) {
        require(index < obligation_count && distinct.insert(index).second,
                "packing index invalid/duplicate");
        support.push_back(index);
    }
    require(input.eof() && support.size() == 37,
            "packing support count changed");
    return support;
}

struct ResponseThr {
    std::array<u64, 2> word{};
    bool operator<(const ResponseThr& other) const { return word < other.word; }
};

struct TypeThr {
    u64 count = 0;
    u32 least = UINT32_MAX;
};

std::vector<u32> build_repaired_carrier_thr(
    const std::vector<u32>& c596, const std::vector<u32>& additions10,
    const std::vector<u32>& deletions) {
    const std::set<u32> deletion_set(deletions.begin(), deletions.end());
    std::vector<u32> carrier;
    for (u32 mask : c596)
        if (!deletion_set.contains(mask)) carrier.push_back(mask);
    std::set<u32> distinct(carrier.begin(), carrier.end());
    for (u32 mask : additions10) {
        require(distinct.insert(mask).second, "addition10 overlap");
        carrier.push_back(mask);
    }
    return carrier;
}

struct DiffAuditThr {
    std::size_t affected_rows = 0;
    u64 support_incidences = 0;
    u64 tested = 0;
    u64 failures = 0;
    u64 row_fnv = kEmptyFnv;
    std::vector<ObligationThr> failure_bodies;
};

DiffAuditThr run_diff_audit_thr(
    const std::vector<u32>& deletions,
    const std::vector<std::vector<std::size_t>>& support,
    const std::map<u32, std::size_t>& c596_index,
    const std::vector<PairAgent>& target, const std::vector<u32>& joint,
    const std::unordered_set<u32>& joint_set,
    const std::vector<u32>& carrier, unsigned worker_count) {
    std::map<std::size_t, std::vector<u32>> row_groups;
    u64 support_incidences = 0;
    for (u32 mask : deletions) {
        const std::size_t mask_index = c596_index.at(mask);
        for (std::size_t row_index : support[mask_index]) {
            row_groups[row_index].push_back(mask);
            ++support_incidences;
        }
    }
    std::vector<std::pair<std::size_t, std::vector<u32>>> jobs;
    for (const auto& [row_index, group] : row_groups)
        jobs.push_back({row_index, group});
    std::vector<LocalAuditThr> audits(jobs.size());
    std::atomic<std::size_t> cursor{0};
    std::exception_ptr worker_error;
    std::mutex error_mutex;
    auto worker = [&]() {
        try {
            while (true) {
                const std::size_t index = cursor.fetch_add(1);
                if (index >= jobs.size()) return;
                audits[index] = audit_local_threshold(
                    target[jobs[index].first], jobs[index].second, joint,
                    joint_set, carrier);
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

    DiffAuditThr summary;
    summary.affected_rows = audits.size();
    summary.support_incidences = support_incidences;
    Fnv row_ledger;
    for (const LocalAuditThr& audit : audits) {
        Fnv group_ledger;
        for (u32 mask : audit.removed_active) group_ledger.add(mask);
        summary.tested += audit.tested;
        summary.failures += audit.failures;
        row_ledger.add(audit.pair.q);
        row_ledger.add(audit.pair.r);
        row_ledger.add(audit.removed_active.size());
        row_ledger.add(group_ledger.state);
        row_ledger.add(audit.tested);
        row_ledger.add(audit.failures);
        row_ledger.add(audit.failure_fnv);
        row_ledger.add(audit.active_survivors);
        row_ledger.add(audit.active_survivor_fnv);
        for (u32 body : audit.failure_bodies)
            summary.failure_bodies.push_back({audit.pair, body});
    }
    summary.row_fnv = row_ledger.state;
    std::sort(summary.failure_bodies.begin(), summary.failure_bodies.end(),
              [](const ObligationThr& left, const ObligationThr& right) {
                  return std::tuple{left.pair.q, left.pair.r, left.body} <
                         std::tuple{right.pair.q, right.pair.r, right.body};
              });
    return summary;
}

}  // namespace

#ifndef THRESHOLD_DELETION_REPLAY_MAIN
#define THRESHOLD_DELETION_REPLAY_MAIN obsolete_threshold_deletion_replay_main
#endif

int THRESHOLD_DELETION_REPLAY_MAIN(int argc, char** argv) {
    try {
        require(argc == 14,
                "usage: replay JOINT BASE8951 ADD45 SUFFIX9 RESIDUAL OLD_UNION "
                "REPAIRS76 ADDITIONS4 DELETE73 ADDITIONS10 THRESHOLD "
                "DELETE_LEDGER WORKERS");
        const unsigned threshold = std::stoul(argv[11]);
        const unsigned worker_count = std::stoul(argv[13]);
        require(threshold <= 391 && worker_count >= 1 && worker_count <= 64,
                "threshold/workers escaped bounds");
        const std::vector<u32> joint = read_masks_agent(
            argv[1], kJointCountAgent, kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> c596 = reconstruct_c596_thr(
            argv[2], argv[3], argv[4], argv[7], argv[8], argv[9]);
        const std::vector<u32> additions10 =
            read_mixed617(argv[10], 10, kAddition10FnvThr);
        const std::vector<u32> frozen_deletions = read_threshold_ledger(argv[12]);

        const std::set<PairKeyThr> old_union =
            read_pair_set_thr(argv[6], 1624, kOldUnionFnvThr);
        const std::vector<PairAgent> raw_band = read_band_agent(argv[5], 595);
        require(raw_band.size() == 394, "raw r>=595 band changed");
        std::vector<PairAgent> target;
        for (PairAgent pair : raw_band)
            if (pair.r >= 596 || !old_union.contains({pair.q, pair.r}))
                target.push_back(pair);
        require(target.size() == 391, "target changed");
        std::vector<Geometry> geometries;
        geometries.reserve(target.size());
        for (PairAgent pair : target)
            geometries.push_back(build_geometry(pair.q, pair.r));

        std::vector<std::vector<std::size_t>> support(c596.size());
        Fnv sign_ledger;
        Fnv support_ledger;
        u64 equality_cells = 0;
        for (std::size_t mask_index = 0; mask_index < c596.size(); ++mask_index) {
            Fnv one_support;
            for (std::size_t row_index = 0; row_index < target.size(); ++row_index) {
                const Margin exact = margin(geometries[row_index],
                                            c596[mask_index]);
                const bool active = exact.ticks >= 0;
                equality_cells += exact.ticks == 0;
                sign_ledger.add(c596[mask_index]);
                sign_ledger.add(active);
                if (!active) continue;
                support[mask_index].push_back(row_index);
                one_support.add(target[row_index].q);
                one_support.add(target[row_index].r);
            }
            support_ledger.add(c596[mask_index]);
            support_ledger.add(support[mask_index].size());
            support_ledger.add(one_support.state);
        }
        require(equality_cells == 0, "activity equality appeared");

        std::vector<u32> expected_deletions;
        std::map<u32, std::size_t> c596_index;
        for (std::size_t index = 0; index < c596.size(); ++index) {
            c596_index.emplace(c596[index], index);
            // The joint deck is a protected structural subset.  The frozen
            // threshold family ranges only over nonjoint C596 masks.
            if (!joint_set.contains(c596[index]) &&
                support[index].size() <= threshold)
                expected_deletions.push_back(c596[index]);
        }
        require(frozen_deletions == expected_deletions,
                "frozen ledger is not the complete carrier-order threshold set");
        const std::unordered_set<u32> deletion_set(frozen_deletions.begin(),
                                                    frozen_deletions.end());
        for (u32 mask : frozen_deletions)
            require(!joint_set.contains(mask), "threshold deleted joint mask");

        std::vector<u32> final_carrier;
        for (u32 mask : c596)
            if (!deletion_set.contains(mask)) final_carrier.push_back(mask);
        std::set<u32> final_distinct(final_carrier.begin(), final_carrier.end());
        for (u32 mask : additions10) {
            require(final_distinct.insert(mask).second,
                    "addition10 overlap after threshold deletion");
            final_carrier.push_back(mask);
        }
        require(final_carrier.size() == 9029 - frozen_deletions.size(),
                "final size changed");
        const u64 rank8 = std::count_if(
            final_carrier.begin(), final_carrier.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        for (u32 mask : joint)
            require(std::find(final_carrier.begin(), final_carrier.end(), mask) !=
                        final_carrier.end(),
                    "final carrier lost joint mask");

        std::map<std::size_t, std::vector<u32>> row_groups;
        u64 support_incidences = 0;
        for (u32 mask : frozen_deletions) {
            const std::size_t mask_index = c596_index.at(mask);
            for (std::size_t row_index : support[mask_index]) {
                row_groups[row_index].push_back(mask);
                ++support_incidences;
            }
        }
        std::vector<LocalAuditThr> audits(row_groups.size());
        std::vector<std::pair<std::size_t, std::vector<u32>>> jobs;
        for (const auto& [row_index, group] : row_groups)
            jobs.push_back({row_index, group});
        std::atomic<std::size_t> cursor{0};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        auto worker = [&]() {
            try {
                while (true) {
                    const std::size_t index = cursor.fetch_add(1);
                    if (index >= jobs.size()) return;
                    audits[index] = audit_local_threshold(
                        target[jobs[index].first], jobs[index].second, joint,
                        joint_set, final_carrier);
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

        Fnv deletion_ledger;
        for (u32 mask : frozen_deletions) deletion_ledger.add(mask);
        std::cout << "LRC14_THRESHOLD_DELETION_REPLAY_V1\n"
                  << "THRESHOLD " << threshold << " TARGET_ROWS "
                  << target.size() << " C596 " << c596.size() << " C596_FNV "
                  << std::hex << masks_fnv_agent(c596) << std::dec << '\n'
                  << "SUPPORT_SIGN_CELLS " << c596.size() * target.size()
                  << " SIGN_FNV " << std::hex << sign_ledger.state << std::dec
                  << " EQUALITIES " << equality_cells << " SUPPORT_LEDGER_FNV "
                  << std::hex << support_ledger.state << std::dec << '\n'
                  << "DELETE_COUNT " << frozen_deletions.size() << " DELETE_FNV "
                  << std::hex << deletion_ledger.state << std::dec
                  << " SUPPORT_INCIDENCES " << support_incidences
                  << " AFFECTED_ROWS " << row_groups.size() << '\n'
                  << "FINAL_CARRIER " << final_carrier.size() << " FNV "
                  << std::hex << masks_fnv_agent(final_carrier) << std::dec
                  << " RANK8 " << rank8 << " RANK9 "
                  << final_carrier.size() - rank8 << " JOINT_RETAINED "
                  << joint.size() << '\n';

        u64 total_tested = 0;
        u64 total_failures = 0;
        Fnv row_ledger;
        for (const LocalAuditThr& audit : audits) {
            Fnv group_ledger;
            for (u32 mask : audit.removed_active) group_ledger.add(mask);
            total_tested += audit.tested;
            total_failures += audit.failures;
            row_ledger.add(audit.pair.q);
            row_ledger.add(audit.pair.r);
            row_ledger.add(audit.removed_active.size());
            row_ledger.add(group_ledger.state);
            row_ledger.add(audit.tested);
            row_ledger.add(audit.failures);
            row_ledger.add(audit.failure_fnv);
            row_ledger.add(audit.active_survivors);
            row_ledger.add(audit.active_survivor_fnv);
            std::cout << "ROW " << audit.pair.q << ',' << audit.pair.r
                      << " REMOVED_ACTIVE " << audit.removed_active.size()
                      << " GROUP_FNV " << std::hex << group_ledger.state
                      << std::dec << " TESTED " << audit.tested << " FAILURES "
                      << audit.failures << " FAILURE_FNV " << std::hex
                      << audit.failure_fnv << std::dec << " ACTIVE_SURVIVORS "
                      << audit.active_survivors << " ACTIVE_FNV " << std::hex
                      << audit.active_survivor_fnv << std::dec << '\n';
        }
        require(total_failures == 0, "threshold deletion is unsafe");
        std::cout << "SUMMARY AFFECTED_ROWS " << audits.size()
                  << " UNAFFECTED_ROWS " << target.size() - audits.size()
                  << " TESTED_RELEVANT_BODIES " << total_tested
                  << " FAILURES " << total_failures << " ROW_LEDGER_FNV "
                  << std::hex << row_ledger.state << std::dec << '\n'
                  << "SCOPE COMPLETE_CARRIER_ORDER_NONJOINT_SUPPORT_THRESHOLD_"
                     "ROW_LOCAL_RELEVANT_BODY_AUDIT_FIXED_391_ROWS_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THRESHOLD_REPLAY_ERROR " << error.what() << '\n';
        return 1;
    }
}

int main(int argc, char** argv) {
    try {
        require(argc == 18,
                "usage: repaired_replay JOINT BASE8951 ADD45 SUFFIX9 "
                "RESIDUAL OLD_UNION REPAIRS76 ADDITIONS4 DELETE73 "
                "ADDITIONS10 THRESHOLD D350 RETAIN37 DFINAL "
                "HOSTILE_FAILURE84 PACKING37 WORKERS");
        const unsigned threshold = std::stoul(argv[11]);
        const unsigned worker_count = std::stoul(argv[17]);
        require(threshold == 350 && worker_count >= 1 && worker_count <= 64,
                "threshold/workers escaped frozen bounds");

        const std::vector<u32> joint = read_masks_agent(
            argv[1], kJointCountAgent, kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> c596 = reconstruct_c596_thr(
            argv[2], argv[3], argv[4], argv[7], argv[8], argv[9]);
        const std::vector<u32> additions10 =
            read_mixed617(argv[10], 10, kAddition10FnvThr);
        const std::vector<u32> d350 = read_threshold_ledger(argv[12]);
        const std::vector<u32> retain37 = read_threshold_ledger(argv[13]);
        const std::vector<u32> dfinal = read_threshold_ledger(argv[14]);
        const std::vector<ObligationThr> hostile =
            read_obligations_thr(argv[15]);
        const std::vector<std::size_t> packing =
            read_packing_thr(argv[16], hostile.size());

        const std::set<PairKeyThr> old_union =
            read_pair_set_thr(argv[6], 1624, kOldUnionFnvThr);
        const std::vector<PairAgent> raw_band = read_band_agent(argv[5], 595);
        require(raw_band.size() == 394, "raw r>=595 band changed");
        std::vector<PairAgent> target;
        std::set<PairKeyThr> target_set;
        for (PairAgent pair : raw_band)
            if (pair.r >= 596 || !old_union.contains({pair.q, pair.r})) {
                target.push_back(pair);
                target_set.insert({pair.q, pair.r});
            }
        require(target.size() == 391, "target changed");
        Fnv target_ledger;
        for (const auto& [q, r] : target_set) {
            target_ledger.add(q);
            target_ledger.add(r);
        }
        require(target_ledger.state == UINT64_C(0xc732a1532c12b9f6),
                "target ledger changed");

        std::vector<Geometry> target_geometry;
        for (PairAgent pair : target)
            target_geometry.push_back(build_geometry(pair.q, pair.r));
        std::vector<std::vector<std::size_t>> support(c596.size());
        std::map<u32, std::size_t> c596_index;
        Fnv sign_ledger;
        Fnv support_ledger;
        u64 equality_cells = 0;
        for (std::size_t mask_index = 0; mask_index < c596.size(); ++mask_index) {
            c596_index.emplace(c596[mask_index], mask_index);
            Fnv one_support;
            for (std::size_t row_index = 0; row_index < target.size(); ++row_index) {
                const Margin exact = margin(target_geometry[row_index],
                                            c596[mask_index]);
                const bool active = exact.ticks >= 0;
                equality_cells += exact.ticks == 0;
                sign_ledger.add(c596[mask_index]);
                sign_ledger.add(active);
                if (!active) continue;
                support[mask_index].push_back(row_index);
                one_support.add(target[row_index].q);
                one_support.add(target[row_index].r);
            }
            support_ledger.add(c596[mask_index]);
            support_ledger.add(support[mask_index].size());
            support_ledger.add(one_support.state);
        }
        require(equality_cells == 0 &&
                    sign_ledger.state == UINT64_C(0x0853edf378527886) &&
                    support_ledger.state == UINT64_C(0xc93aae41a6892f37),
                "support matrix identity changed");

        std::vector<u32> expected_d350;
        for (std::size_t index = 0; index < c596.size(); ++index)
            if (!joint_set.contains(c596[index]) &&
                support[index].size() <= threshold)
                expected_d350.push_back(c596[index]);
        require(d350 == expected_d350 && d350.size() == 5141 &&
                    masks_fnv_agent(d350) == UINT64_C(0x03921cf597ee9863),
                "D350 is not the exact carrier-order threshold set");
        const std::set<u32> d350_set(d350.begin(), d350.end());
        const std::set<u32> retain_set(retain37.begin(), retain37.end());
        require(retain37.size() == 37 && retain_set.size() == 37 &&
                    std::includes(d350_set.begin(), d350_set.end(),
                                  retain_set.begin(), retain_set.end()) &&
                    masks_fnv_agent(retain37) ==
                        UINT64_C(0x531da886ae2e455a),
                "H37 identity/subset changed");
        std::vector<u32> expected_dfinal;
        for (u32 mask : d350)
            if (!retain_set.contains(mask)) expected_dfinal.push_back(mask);
        require(dfinal == expected_dfinal && dfinal.size() == 5104 &&
                    masks_fnv_agent(dfinal) ==
                        UINT64_C(0xff4c932f9a7adac8),
                "Dfinal != D350\\H37");

        const std::vector<u32> hostile_carrier =
            build_repaired_carrier_thr(c596, additions10, d350);
        const std::vector<u32> final_carrier =
            build_repaired_carrier_thr(c596, additions10, dfinal);
        const u64 hostile_rank8 = std::count_if(
            hostile_carrier.begin(), hostile_carrier.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        const u64 final_rank8 = std::count_if(
            final_carrier.begin(), final_carrier.end(),
            [](u32 mask) { return std::popcount(mask) == 8; });
        require(hostile_carrier.size() == 3888 && hostile_rank8 == 3826 &&
                    masks_fnv_agent(hostile_carrier) ==
                        UINT64_C(0x4d1082144fd457d8),
                "hostile carrier identity changed");
        require(final_carrier.size() == 3925 && final_rank8 == 3858 &&
                    masks_fnv_agent(final_carrier) ==
                        UINT64_C(0x6fbd0bffcf0ed78b),
                "final carrier identity changed");
        for (u32 mask : joint)
            require(std::find(final_carrier.begin(), final_carrier.end(), mask) !=
                        final_carrier.end(),
                    "final carrier lost protected joint mask");

        const DiffAuditThr hostile_audit = run_diff_audit_thr(
            d350, support, c596_index, target, joint, joint_set,
            hostile_carrier, worker_count);
        require(hostile_audit.failure_bodies == hostile &&
                    hostile_audit.failures == 84,
                "independent hostile residual differs from frozen 84");
        Fnv hostile_failure_ledger;
        for (const ObligationThr& obligation : hostile) {
            hostile_failure_ledger.add(obligation.pair.q);
            hostile_failure_ledger.add(obligation.pair.r);
            hostile_failure_ledger.add(obligation.body);
        }
        require(hostile_failure_ledger.state ==
                    UINT64_C(0x3a4207b0cb910c10),
                "hostile failure ledger identity changed");

        std::vector<Geometry> hostile_geometry;
        for (const ObligationThr& obligation : hostile)
            hostile_geometry.push_back(
                build_geometry(obligation.pair.q, obligation.pair.r));
        std::map<ResponseThr, TypeThr> types;
        std::map<u32, ResponseThr> response_by_mask;
        Fnv responder_ledger;
        u64 responders = 0;
        for (u32 mask : d350) {
            ResponseThr response;
            for (std::size_t index = 0; index < hostile.size(); ++index)
                if ((mask & hostile[index].body) == 0 &&
                    margin(hostile_geometry[index], mask).ticks >= 0)
                    response.word[index / 64] |=
                        UINT64_C(1) << (index % 64);
            response_by_mask.emplace(mask, response);
            if (response.word[0] == 0 && response.word[1] == 0) continue;
            ++responders;
            responder_ledger.add(mask);
            TypeThr& type = types[response];
            ++type.count;
            type.least = std::min(type.least, mask);
        }
        Fnv type_ledger;
        for (const auto& [response, type] : types) {
            type_ledger.add(response.word[0]);
            type_ledger.add(response.word[1]);
            type_ledger.add(type.count);
            type_ledger.add(type.least);
        }
        require(responders == 70 &&
                    responder_ledger.state ==
                        UINT64_C(0xb18cf54b71b4a093) &&
                    types.size() == 53 &&
                    type_ledger.state == UINT64_C(0x578ebd478ecee5a),
                "residual response quotient identity changed");

        ResponseThr cover_union;
        Fnv retain_response_ledger;
        for (u32 mask : retain37) {
            const ResponseThr response = response_by_mask.at(mask);
            cover_union.word[0] |= response.word[0];
            cover_union.word[1] |= response.word[1];
            retain_response_ledger.add(mask);
            retain_response_ledger.add(response.word[0]);
            retain_response_ledger.add(response.word[1]);
        }
        require(cover_union.word[0] == UINT64_MAX &&
                    cover_union.word[1] == (UINT64_C(1) << 20) - 1,
                "H37 does not cover all 84 residual obligations");
        for (const auto& [response, type] : types) {
            unsigned packed_hits = 0;
            for (std::size_t index : packing)
                packed_hits +=
                    (response.word[index / 64] >> (index % 64)) & 1U;
            require(packed_hits <= 1,
                    "packing lower certificate violated by response type");
        }
        Fnv packing_ledger;
        for (std::size_t index : packing) packing_ledger.add(index);
        require(retain_response_ledger.state ==
                    UINT64_C(0x061eaf6e2f6fbc8f),
                "H37 mask/response identity changed");

        const DiffAuditThr final_audit = run_diff_audit_thr(
            dfinal, support, c596_index, target, joint, joint_set,
            final_carrier, worker_count);
        require(final_audit.failures == 0 &&
                    final_audit.failure_bodies.empty(),
                "final relevant-body replay found a failure");

        std::cout
            << "LRC14_INDEPENDENT_REPAIRED_THRESHOLD_REPLAY_V1\n"
            << "TARGET_ROWS " << target.size() << " ROW_FNV " << std::hex
            << target_ledger.state << std::dec << " C596 " << c596.size()
            << " C596_FNV " << std::hex << masks_fnv_agent(c596) << std::dec
            << '\n'
            << "SUPPORT_SIGN_CELLS " << c596.size() * target.size()
            << " SIGN_FNV " << std::hex << sign_ledger.state
            << " SUPPORT_LEDGER_FNV " << support_ledger.state << std::dec
            << " EQUALITIES " << equality_cells << '\n'
            << "D350 " << d350.size() << " FNV " << std::hex
            << masks_fnv_agent(d350) << std::dec << " THRESHOLD " << threshold
            << " JOINT_EXCLUDED " << joint.size() << '\n'
            << "HOSTILE_CARRIER " << hostile_carrier.size() << " FNV "
            << std::hex << masks_fnv_agent(hostile_carrier) << std::dec
            << " RANK8 " << hostile_rank8 << " RANK9 "
            << hostile_carrier.size() - hostile_rank8 << '\n'
            << "HOSTILE_REPLAY AFFECTED_ROWS " << hostile_audit.affected_rows
            << " SUPPORT_INCIDENCES " << hostile_audit.support_incidences
            << " TESTED_RELEVANT_BODIES " << hostile_audit.tested
            << " FAILURES " << hostile_audit.failures << " FAILURE_FNV "
            << std::hex << hostile_failure_ledger.state << " ROW_FNV "
            << hostile_audit.row_fnv << std::dec << '\n'
            << "RESIDUAL_RESPONDERS " << responders << " FNV " << std::hex
            << responder_ledger.state << std::dec << " TYPES " << types.size()
            << " TYPE_FNV " << std::hex << type_ledger.state << std::dec << '\n'
            << "RETENTION_COVER " << retain37.size() << " MASK_FNV "
            << std::hex << masks_fnv_agent(retain37)
            << " MASK_RESPONSE_FNV " << retain_response_ledger.state
            << std::dec << " PACKING " << packing.size() << " PACKING_FNV "
            << std::hex << packing_ledger.state << std::dec
            << " EXACT_MIN 37\n"
            << "DFINAL " << dfinal.size() << " FNV " << std::hex
            << masks_fnv_agent(dfinal) << std::dec
            << " IDENTITY D350_MINUS_H37\n"
            << "FINAL_CARRIER " << final_carrier.size() << " FNV "
            << std::hex << masks_fnv_agent(final_carrier) << std::dec
            << " RANK8 " << final_rank8 << " RANK9 "
            << final_carrier.size() - final_rank8 << " JOINT_RETAINED "
            << joint.size() << '\n'
            << "FINAL_REPLAY AFFECTED_ROWS " << final_audit.affected_rows
            << " SUPPORT_INCIDENCES " << final_audit.support_incidences
            << " TESTED_RELEVANT_BODIES " << final_audit.tested
            << " FAILURES " << final_audit.failures << " ROW_FNV "
            << std::hex << final_audit.row_fnv << std::dec << '\n'
            << "SCOPE INDEPENDENT_THRESHOLD_RECONSTRUCTION_RESIDUAL_"
               "HYPERGRAPH_COVER_PACKING_AND_RELEVANT_BODY_REPLAY_"
               "FIXED_391_ROWS_NO_PHYSICAL_ENTRY_NO_LRC14\n"
            << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "INDEPENDENT_REPAIRED_THRESHOLD_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}

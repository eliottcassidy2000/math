// Exact protected-nonjoint-mask quotient for one deletion after adjoining the
// selected endpoint-593 response, on the full 432-row preservation target.

#define ENDPOINT593_RESPONSE_CAPACITY_MAIN endpoint593_capacity_hidden_main
#include "endpoint593_response_capacity.cpp"
#undef ENDPOINT593_RESPONSE_CAPACITY_MAIN

#include <atomic>
#include <mutex>
#include <thread>

namespace {
constexpr u64 kSelectedAdditionFnv593 = UINT64_C(0x60873ef7a2b4ab90);

struct SingletonRow593 {
    PairAgent pair;
    u64 body_tests = 0;
    u64 exposed = 0;
    std::vector<std::pair<u32, u32>> singleton;  // body, sole witness
};

SingletonRow593 audit_singletons593(
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

    SingletonRow593 result;
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
        require(hits != 0, "augmented carrier failed target");
        if (hits == 1) result.singleton.push_back({body, sole});
    }
    require(result.body_tests == kBodyCount617, "body universe changed");
    return result;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 20,
                "usage: quotient JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE "
                "PREFIX391_CSV TOP594_CSV TOP593 ADDITION1 SINGLETON_CSV "
                "PROTECTED_TXT SAFE_DELETE1 WORKERS");
        const unsigned worker_count = std::stoul(argv[19]);
        require(worker_count >= 1 && worker_count <= 64,
                "workers escaped 1..64");
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> carrier = reconstruct_c3925_593(
            argv[2], argv[3], argv[4], argv[7], argv[8], argv[9], argv[10],
            argv[11]);
        (void)read_pairs_repaired(argv[5], 22647, kUniverseFnv);
        (void)read_pairs_repaired(argv[6], 1624, kOldUnionFnv);
        const auto prefix391 = read_pair_audit_rows593(argv[12], 391);
        const auto top594 = read_pair_audit_rows593(argv[13], 25);
        const auto top593_set = read_pairs_repaired(argv[14], 16, kTop593Fnv593);
        const std::vector<u32> addition = read_flexible_masks593(
            argv[15], 1, kSelectedAdditionFnv593);
        require(!joint_set.contains(addition.front()),
                "selected addition entered protected joint deck");

        std::set<std::pair<int, int>> target_set;
        for (PairAgent pair : prefix391)
            require(target_set.insert({pair.q, pair.r}).second,
                    "prefix row overlap");
        for (PairAgent pair : top594)
            require(target_set.insert({pair.q, pair.r}).second,
                    "top594 row overlap");
        for (const auto& pair : top593_set)
            require(target_set.insert(pair).second, "top593 row overlap");
        require(target_set.size() == 432, "full target count changed");
        std::vector<PairAgent> rows;
        Fnv target_ledger;
        for (const auto& [q, r] : target_set) {
            rows.push_back(PairAgent{q, r});
            target_ledger.add(q);
            target_ledger.add(r);
        }
        require(target_ledger.state == kFull432Fnv593,
                "full target identity changed");

        std::vector<u32> augmented = carrier;
        require(std::find(augmented.begin(), augmented.end(), addition.front()) ==
                    augmented.end(),
                "selected addition already in carrier");
        augmented.push_back(addition.front());
        require(augmented.size() == 3926, "augmented size changed");
        const u64 nonjoint_count = std::count_if(
            augmented.begin(), augmented.end(),
            [&](u32 mask) { return !joint_set.contains(mask); });
        require(nonjoint_count == 3505, "augmented nonjoint census changed");

        std::vector<SingletonRow593> audits(rows.size());
        std::atomic<std::size_t> cursor{0};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        auto worker = [&]() {
            try {
                while (true) {
                    const std::size_t index = cursor.fetch_add(1);
                    if (index >= rows.size()) return;
                    audits[index] = audit_singletons593(
                        rows[index], joint, joint_set, augmented);
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

        std::ofstream singleton_out(argv[16]);
        std::ofstream protected_out(argv[17]);
        std::ofstream safe_out(argv[18]);
        require(singleton_out && protected_out && safe_out,
                "cannot create singleton/deletion ledgers");
        singleton_out << "q,r,body_hex,witness_hex,witness_rank\n";
        std::map<u32, std::pair<u64, std::set<std::size_t>>> private_by_mask;
        Fnv singleton_ledger;
        u64 body_tests = 0, exposed = 0, private_obligations = 0;
        for (std::size_t row_index = 0; row_index < audits.size(); ++row_index) {
            const auto& audit = audits[row_index];
            body_tests += audit.body_tests;
            exposed += audit.exposed;
            for (const auto& [body, witness] : audit.singleton) {
                ++private_obligations;
                auto& item = private_by_mask[witness];
                ++item.first;
                item.second.insert(row_index);
                singleton_ledger.add(audit.pair.q);
                singleton_ledger.add(audit.pair.r);
                singleton_ledger.add(body);
                singleton_ledger.add(witness);
                singleton_out << audit.pair.q << ',' << audit.pair.r << ','
                              << hex8(body) << ',' << hex8(witness) << ','
                              << std::popcount(witness) << '\n';
            }
        }
        Fnv protected_ledger;
        for (const auto& [mask, item] : private_by_mask) {
            protected_out << hex8(mask) << '\n';
            protected_ledger.add(mask);
        }
        const std::set<u32> protected_set = [&]() {
            std::set<u32> result;
            for (const auto& [mask, item] : private_by_mask) result.insert(mask);
            return result;
        }();
        u32 safe_delete = UINT32_MAX;
        for (u32 mask : carrier)
            if (!joint_set.contains(mask) && !protected_set.contains(mask))
                safe_delete = std::min(safe_delete, mask);
        require(safe_delete != UINT32_MAX, "no safe original nonjoint deletion");
        safe_out << hex8(safe_delete) << '\n';
        require(singleton_out.good() && protected_out.good() && safe_out.good(),
                "singleton/deletion ledger write failed");

        std::cout << "LRC14_ENDPOINT593_FULL432_SINGLETON_QUOTIENT_V1\n"
                  << "BASE_C3925_FNV " << std::hex
                  << masks_fnv_agent(carrier) << " ADDITION "
                  << hex8(addition.front()) << " AUGMENTED_FNV "
                  << masks_fnv_agent(augmented) << std::dec << '\n'
                  << "TARGET_ROWS " << rows.size() << " TARGET_FNV "
                  << std::hex << target_ledger.state << std::dec
                  << " BODY_TESTS " << body_tests << " EXPOSED " << exposed
                  << '\n'
                  << "PRIVATE_OBLIGATIONS " << private_obligations
                  << " SINGLETON_FNV " << std::hex << singleton_ledger.state
                  << std::dec << " PROTECTED_NONJOINT "
                  << private_by_mask.size() << " PROTECTED_FNV " << std::hex
                  << protected_ledger.state << std::dec << '\n'
                  << "NONJOINT_TOTAL " << nonjoint_count << " SAFE "
                  << nonjoint_count - private_by_mask.size() << " UNSAFE "
                  << private_by_mask.size() << " JOINT_PROTECTED "
                  << joint.size() << '\n'
                  << "SELECTED_SAFE_DELETE " << hex8(safe_delete) << " RANK "
                  << std::popcount(safe_delete) << " DELETE_FNV " << std::hex
                  << masks_fnv_agent(std::vector<u32>{safe_delete}) << std::dec
                  << '\n';
        std::cout << "SCOPE EXACT_NONJOINT_SINGLE_DELETION_BOUNDARY_"
                     "AUGMENTED_CARRIER_FIXED_432_ROWS_JOINT_DECK_RETAINED_"
                     "NO_SIMULTANEOUS_DELETION_INFERENCE_NO_PHYSICAL_ENTRY_"
                     "NO_LRC14\nVERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT593_SINGLETON_QUOTIENT_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}

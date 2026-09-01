// Complete multiplicity-through-two audit of THM-4314's endpoint-590 layer.
// The all-carrier and protected-joint notions are computed separately in the
// same exhaustive row/body traversal. Scratch only; no physical-entry claim.

#define ENDPOINT591_TWOHIT_MAIN endpoint591_twohit_hidden_for_endpoint590
#include "04-computation/lrc14_endpoint591_complete_layer_closure_deletion_boundary_thm4314/endpoint591_twohit_census.cpp"
#undef ENDPOINT591_TWOHIT_MAIN

namespace {

constexpr u64 kTop590FnvStructure = UINT64_C(0x44aa8a793d162cf9);

struct LowWitness590 {
    u32 body = 0;
    unsigned hits = 0;
    u32 first = 0;
    u32 second = 0;
};

struct RowAudit590 {
    u64 bodies = 0;
    u64 all_zero = 0;
    u64 all_one = 0;
    u64 all_two = 0;
    u64 joint_exposed = 0;
    u64 protected_zero = 0;
    u64 protected_one = 0;
    u64 protected_two = 0;
    std::size_t active_joint = 0;
    std::size_t active_nonjoint = 0;
    std::vector<LowWitness590> all_low;
    std::vector<LowWitness590> protected_low;
};

std::vector<u32> reconstruct_final_carrier590(int argc, char** argv,
                                               std::vector<u32>& joint,
                                               std::vector<PairAgent>& rows) {
    require(argc == 20,
            "usage: endpoint590-low JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
            "OLD_UNION REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE "
            "ADDITION593 DELETE593 TOP590 COVER43 DELETE43 ALL_LOW_CSV "
            "PROTECTED_LOW_CSV WORKERS");
    joint = read_masks_agent(argv[1], kJointCountAgent, kJointFnvAgent, 8);
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
            "THM4311 exchanged carrier changed");

    const auto additions = read_cover591(argv[15]);
    const auto deletions = read_deletions591(argv[16]);
    require(masks_fnv_agent(additions) == kCover43Fnv &&
                masks_fnv_agent(deletions) == kDeletion43Fnv,
            "THM4313 exchange ledgers changed");
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
                "protected joint mask lost");

    const auto row_set = read_pairs_repaired(argv[14], 13,
                                             kTop590FnvStructure);
    for (const auto& [q, r] : row_set) {
        require(r == 590, "non-590 row escaped frontier");
        rows.push_back({q, r});
    }
    return carrier;
}

RowAudit590 census_row590(PairAgent row, const std::vector<u32>& joint,
                          const std::vector<u32>& carrier,
                          const std::unordered_set<u32>& joint_set) {
    const Geometry geometry = build_geometry(row.q, row.r);
    std::vector<u32> active_all, active_joint, active_nonjoint;
    for (u32 mask : carrier) {
        if (margin(geometry, mask).ticks < 0) continue;
        active_all.push_back(mask);
        if (joint_set.contains(mask)) active_joint.push_back(mask);
        else active_nonjoint.push_back(mask);
    }

    RowAudit590 result;
    result.active_joint = active_joint.size();
    result.active_nonjoint = active_nonjoint.size();
    const u32 limit = u32{1} << 30;
    for (u32 body = (u32{1} << 9) - 1; body < limit;
         body = next_combination(body)) {
        ++result.bodies;
        std::array<u32, 3> witnesses{};
        unsigned hits = 0;
        for (u32 mask : active_all) {
            if ((mask & body) != 0) continue;
            witnesses[hits++] = mask;
            if (hits == 3) break;
        }
        if (hits <= 2) {
            if (hits == 2 && witnesses[1] < witnesses[0])
                std::swap(witnesses[0], witnesses[1]);
            result.all_low.push_back({body, hits, witnesses[0], witnesses[1]});
            result.all_zero += hits == 0;
            result.all_one += hits == 1;
            result.all_two += hits == 2;
        }

        bool joint_hit = false;
        for (u32 mask : active_joint)
            if ((mask & body) == 0) {
                joint_hit = true;
                break;
            }
        if (joint_hit) {
            require(hits != 0, "joint hit absent from all-carrier census");
            continue;
        }

        ++result.joint_exposed;
        std::array<u32, 3> nonjoint_witnesses{};
        unsigned nonjoint_hits = 0;
        for (u32 mask : active_nonjoint) {
            if ((mask & body) != 0) continue;
            nonjoint_witnesses[nonjoint_hits++] = mask;
            if (nonjoint_hits == 3) break;
        }
        if (nonjoint_hits <= 2) {
            if (nonjoint_hits == 2 &&
                nonjoint_witnesses[1] < nonjoint_witnesses[0])
                std::swap(nonjoint_witnesses[0], nonjoint_witnesses[1]);
            result.protected_low.push_back(
                {body, nonjoint_hits, nonjoint_witnesses[0],
                 nonjoint_witnesses[1]});
            result.protected_zero += nonjoint_hits == 0;
            result.protected_one += nonjoint_hits == 1;
            result.protected_two += nonjoint_hits == 2;
        }
        if (hits <= 2) {
            require(hits == nonjoint_hits && witnesses[0] == nonjoint_witnesses[0]
                        && witnesses[1] == nonjoint_witnesses[1],
                    "all/protected low census diverged without joint hit");
        }
    }
    require(result.bodies == kBodyCount617, "body universe changed");
    return result;
}

void write_low590(std::ofstream& output, const std::vector<PairAgent>& rows,
                  const std::vector<RowAudit590>& audits, bool protected_scope,
                  Fnv& ledger) {
    output << "q,r,body_hex,hits,first_mask_hex,second_mask_hex\n";
    for (std::size_t index = 0; index < rows.size(); ++index) {
        const auto& low = protected_scope ? audits[index].protected_low
                                          : audits[index].all_low;
        for (const auto& item : low) {
            output << rows[index].q << ',' << rows[index].r << ','
                   << hex8(item.body) << ',' << item.hits << ',';
            if (item.hits >= 1) output << hex8(item.first);
            output << ',';
            if (item.hits >= 2) output << hex8(item.second);
            output << '\n';
            ledger.add(rows[index].q); ledger.add(rows[index].r);
            ledger.add(item.body); ledger.add(item.hits);
            ledger.add(item.first); ledger.add(item.second);
        }
    }
    require(output.good(), "low-witness ledger write failed");
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const unsigned worker_count = argc == 20 ? std::stoul(argv[19]) : 0;
        require(worker_count >= 1 && worker_count <= 64,
                "workers escaped 1..64");
        std::vector<u32> joint;
        std::vector<PairAgent> rows;
        const std::vector<u32> carrier =
            reconstruct_final_carrier590(argc, argv, joint, rows);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());

        std::vector<RowAudit590> audits(rows.size());
        std::atomic<std::size_t> cursor{0};
        std::exception_ptr worker_error;
        std::mutex error_mutex;
        auto worker = [&]() {
            try {
                while (true) {
                    const std::size_t index = cursor.fetch_add(1);
                    if (index >= rows.size()) return;
                    audits[index] = census_row590(rows[index], joint, carrier,
                                                  joint_set);
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

        std::ofstream all_output(argv[17]);
        std::ofstream protected_output(argv[18]);
        require(all_output && protected_output, "cannot create output ledgers");
        Fnv all_ledger, protected_ledger;
        write_low590(all_output, rows, audits, false, all_ledger);
        write_low590(protected_output, rows, audits, true, protected_ledger);

        u64 all_zero = 0, all_one = 0, all_two = 0;
        u64 exposed = 0, protected_zero = 0, protected_one = 0,
            protected_two = 0;
        std::cout << "LRC14_ENDPOINT590_LOW_WITNESS_CENSUS_V1\n"
                  << "CARRIER 3925 FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec
                  << " RANK8 3818 RANK9 107 JOINT 421 ROWS 13 ROW_FNV "
                  << std::hex << kTop590FnvStructure << std::dec
                  << " BODY_TESTS " << 13 * kBodyCount617 << '\n';
        for (std::size_t index = 0; index < rows.size(); ++index) {
            const auto& audit = audits[index];
            all_zero += audit.all_zero; all_one += audit.all_one;
            all_two += audit.all_two; exposed += audit.joint_exposed;
            protected_zero += audit.protected_zero;
            protected_one += audit.protected_one;
            protected_two += audit.protected_two;
            std::cout << "ROW " << rows[index].q << ',' << rows[index].r
                      << " ACTIVE_J " << audit.active_joint
                      << " ACTIVE_N " << audit.active_nonjoint
                      << " ALL_ZERO " << audit.all_zero
                      << " ALL_ONE " << audit.all_one
                      << " ALL_TWO " << audit.all_two
                      << " J_EXPOSED " << audit.joint_exposed
                      << " N_ZERO " << audit.protected_zero
                      << " N_ONE " << audit.protected_one
                      << " N_TWO " << audit.protected_two << '\n';
        }
        std::cout << "ALL_SUMMARY ZERO " << all_zero << " ONE " << all_one
                  << " TWO " << all_two << " LOW_FNV " << std::hex
                  << all_ledger.state << std::dec << '\n'
                  << "PROTECTED_SUMMARY EXPOSED " << exposed << " ZERO "
                  << protected_zero << " ONE " << protected_one << " TWO "
                  << protected_two << " LOW_FNV " << std::hex
                  << protected_ledger.state << std::dec << '\n'
                  << "SCOPE FIXED_THM4313_CARRIER_COMPLETE_ENDPOINT590_LAYER_"
                     "ALL_VERSUS_PROTECTED_JOINT_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT590_LOW_WITNESS_ERROR " << error.what() << '\n';
        return 1;
    }
}

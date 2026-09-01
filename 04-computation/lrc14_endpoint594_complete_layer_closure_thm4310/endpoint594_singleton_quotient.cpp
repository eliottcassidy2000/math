// Exact one-deletion quotient of the THM-4309 3,925-mask carrier on the
// endpoint-594 layer, promoted for THM-4310.

#define ENDPOINT595_REPAIRED_RAW_MAIN endpoint594_quotient_hidden_main
#include "04-computation/lrc14_endpoint595_support_threshold_residual_hypergraph_compression_thm4309/strict_repair_raw_audit.cpp"
#undef ENDPOINT595_REPAIRED_RAW_MAIN

#include <map>

namespace {
constexpr u64 kTargetFnvQ = UINT64_C(0xcce015c81f7121d9);
constexpr u64 kDeleteFnvQ = UINT64_C(0xff4c932f9a7adac8);
constexpr u64 kCarrierFnvQ = UINT64_C(0x6fbd0bffcf0ed78b);

std::vector<u32> read_flexible_q(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open deletion ledger");
    std::vector<u32> masks;
    std::set<u32> distinct;
    std::string token;
    while (input >> token) {
        const u32 mask = parse_mask_agent(token);
        require((std::popcount(mask) == 8 || std::popcount(mask) == 9) &&
                    distinct.insert(mask).second,
                "deletion rank/distinctness changed");
        masks.push_back(mask);
    }
    require(input.eof(), "deletion read failed");
    return masks;
}

struct PrivateSummary {
    u64 obligations = 0;
    u32 rows = 0;
};
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 15,
                "usage: quotient JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION REPAIRS76 OLD_ADD4 OLD_DELETE73 ADD10 "
                "FINAL_DELETE TARGET25 SINGLETON_CSV SUMMARY_CSV");
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const auto repairs = read_mixed617(argv[7], 76, kRepairFnv);
        const auto old_additions = read_mixed617(argv[8], 4, kOldAdditionFnv);
        const auto old_deletions = read_mixed617(argv[9], 73, kOldDeleteFnv);
        const auto additions = read_mixed617(argv[10], 10, kNewAdditionFnv);
        const auto deletions = read_flexible_q(argv[11]);
        Fnv deletion_fnv;
        for (u32 mask : deletions) deletion_fnv.add(mask);
        require(deletions.size() == 5104 && deletion_fnv.state == kDeleteFnvQ,
                "final deletion identity changed");

        std::vector<u32> augmented = build_mixed_carrier(argv[2], argv[3], argv[4]);
        std::set<u32> distinct(augmented.begin(), augmented.end());
        for (u32 mask : repairs) {
            require(distinct.insert(mask).second, "repair overlap");
            augmented.push_back(mask);
        }
        require(augmented.size() == 9088 && masks_fnv_agent(augmented) == kAugmentedFnv,
                "augmented carrier changed");
        const std::set<u32> old_deleted(old_deletions.begin(), old_deletions.end());
        std::vector<u32> c596;
        for (u32 mask : augmented)
            if (!old_deleted.contains(mask)) c596.push_back(mask);
        for (u32 mask : old_additions) c596.push_back(mask);
        require(c596.size() == 9019 && masks_fnv_agent(c596) == kC596Fnv,
                "C596 changed");
        const std::set<u32> deleted(deletions.begin(), deletions.end());
        std::vector<u32> carrier;
        for (u32 mask : c596)
            if (!deleted.contains(mask)) carrier.push_back(mask);
        distinct = std::set<u32>(carrier.begin(), carrier.end());
        for (u32 mask : additions) {
            require(distinct.insert(mask).second, "addition overlap");
            carrier.push_back(mask);
        }
        require(carrier.size() == 3925 && masks_fnv_agent(carrier) == kCarrierFnvQ,
                "C3925 changed");
        const u64 nonjoint_count = std::count_if(
            carrier.begin(), carrier.end(),
            [&](u32 mask) { return !joint_set.contains(mask); });
        require(nonjoint_count == 3504, "nonjoint carrier census changed");

        (void)read_pairs_repaired(argv[5], 22647, kUniverseFnv);
        (void)read_pairs_repaired(argv[6], 1624, kOldUnionFnv);
        const auto target_set = read_pairs_repaired(argv[12], 25, kTargetFnvQ);
        std::vector<PairAgent> rows;
        for (const auto& [q, r] : target_set) {
            require(r == 594, "target escaped endpoint 594");
            rows.push_back(PairAgent{q, r});
        }

        std::ofstream singleton_out(argv[13]);
        std::ofstream summary_out(argv[14]);
        require(singleton_out && summary_out, "cannot create quotient ledgers");
        singleton_out << "q,r,body_hex,witness_hex,witness_rank\n";
        summary_out << "witness_hex,witness_rank,private_obligations,row_count,row_bits_hex\n";

        std::map<u64, u64> hit_histogram;
        std::map<u32, PrivateSummary> private_by_mask;
        Fnv exposed_fnv, singleton_fnv;
        u64 exposed_total = 0, incidence_total = 0, body_total = 0;
        for (std::size_t row_index = 0; row_index < rows.size(); ++row_index) {
            const auto pair = rows[row_index];
            const Geometry geometry = build_geometry(pair.q, pair.r);
            std::vector<u32> active_joint;
            std::vector<u32> active_nonjoint;
            for (u32 mask : joint)
                if (margin(geometry, mask).ticks >= 0) active_joint.push_back(mask);
            for (u32 mask : carrier)
                if (margin(geometry, mask).ticks >= 0 && !joint_set.contains(mask))
                    active_nonjoint.push_back(mask);

            u64 row_bodies = 0;
            const u32 limit = u32{1} << 30;
            for (u32 body = (u32{1} << 9) - 1; body < limit;
                 body = next_combination(body)) {
                ++row_bodies;
                bool joint_hit = false;
                for (u32 mask : active_joint)
                    if ((mask & body) == 0) {
                        joint_hit = true;
                        break;
                    }
                if (joint_hit) continue;

                ++exposed_total;
                exposed_fnv.add(pair.q);
                exposed_fnv.add(pair.r);
                exposed_fnv.add(body);
                u64 hits = 0;
                u32 sole = 0;
                for (u32 mask : active_nonjoint)
                    if ((mask & body) == 0) {
                        ++hits;
                        sole = mask;
                    }
                require(hits != 0, "base carrier failed target");
                incidence_total += hits;
                ++hit_histogram[hits];
                if (hits == 1) {
                    auto& item = private_by_mask[sole];
                    ++item.obligations;
                    item.rows |= u32{1} << row_index;
                    singleton_fnv.add(pair.q);
                    singleton_fnv.add(pair.r);
                    singleton_fnv.add(body);
                    singleton_fnv.add(sole);
                    singleton_out << pair.q << ',' << pair.r << ',' << hex8(body)
                                  << ',' << hex8(sole) << ','
                                  << std::popcount(sole) << '\n';
                }
            }
            require(row_bodies == kBodyCount617, "body universe changed");
            body_total += row_bodies;
        }

        Fnv protected_fnv, summary_fnv;
        u64 private_obligations = 0;
        for (const auto& [mask, item] : private_by_mask) {
            private_obligations += item.obligations;
            protected_fnv.add(mask);
            summary_fnv.add(mask);
            summary_fnv.add(item.obligations);
            summary_fnv.add(item.rows);
            summary_out << hex8(mask) << ',' << std::popcount(mask) << ','
                        << item.obligations << ',' << std::popcount(item.rows)
                        << ',' << std::hex << item.rows << std::dec << '\n';
        }
        require(singleton_out.good() && summary_out.good(),
                "quotient ledger write failed");

        std::cout << "LRC14_ENDPOINT594_C3925_SINGLETON_DELETION_QUOTIENT_V1\n"
                  << "CARRIER 3925 FNV " << std::hex << kCarrierFnvQ << std::dec
                  << " TARGET_ROWS " << rows.size() << " FNV " << std::hex
                  << kTargetFnvQ
                  << std::dec << " BODY_TESTS " << body_total << '\n'
                  << "EXPOSED " << exposed_total << " EXPOSED_TRIPLE_FNV "
                  << std::hex << exposed_fnv.state << std::dec
                  << " HIT_INCIDENCES " << incidence_total << '\n'
                  << "PRIVATE_OBLIGATIONS " << private_obligations
                  << " SINGLETON_QUADRUPLE_FNV " << std::hex
                  << singleton_fnv.state << std::dec << " PROTECTED_MASKS "
                  << private_by_mask.size() << " PROTECTED_MASK_FNV "
                  << std::hex << protected_fnv.state << " SUMMARY_FNV "
                  << summary_fnv.state << std::dec << '\n'
                  << "HIT_HISTOGRAM";
        for (const auto& [hits, count] : hit_histogram)
            std::cout << ' ' << hits << ':' << count;
        std::cout << '\n';
        for (const auto& [mask, item] : private_by_mask)
            std::cout << "PROTECTED " << hex8(mask) << " RANK "
                      << std::popcount(mask) << " PRIVATE " << item.obligations
                      << " ROWS " << std::popcount(item.rows) << " ROW_BITS "
                      << std::hex << item.rows << std::dec << '\n';
        std::cout << "EXACT_NONJOINT_ONE_DELETION_QUOTIENT "
                  << "SAFE_NONJOINT_SINGLE_DELETIONS "
                  << nonjoint_count - private_by_mask.size()
                  << " UNSAFE_NONJOINT_SINGLE_DELETIONS "
                  << private_by_mask.size()
                  << " JOINT_SINGLE_DELETIONS_UNAUDITED "
                  << carrier.size() - nonjoint_count << '\n'
                  << "SCOPE EXACT_FOR_C3925_ON_FIXED_25_ROW_ENDPOINT594_LAYER_"
                     "NONJOINT_SINGLE_DELETIONS_ONLY_NO_OLD391_PRESERVATION_"
                     "INFERENCE_NO_PHYSICAL_ENTRY_NO_LRC14\nVERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT594_SINGLETON_QUOTIENT_ERROR " << error.what()
                  << '\n';
        return 1;
    }
}

// Fourth learned-carrier CEGAR round at the unique post-THM-4267 top edge.
//
// This program is deliberately a /tmp research artifact.  It rebuilds the
// canonical THM-4266 carrier from the frozen transcripts, appends only masks
// novel in the order-minimal active prefix for (520,688), checks the seed by
// exact atoms and an independent component formula, and descends through an
// explicitly supplied semantic residual.

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif
#define CARRIER_CEGAR_LIBRARY_ONLY
#include "04-computation/lrc14_three_round_learned_carrier_thm4266/carrier_cegar_descent.cpp"
#undef CARRIER_CEGAR_LIBRARY_ONLY
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

namespace {

void append_unseen_round4(const std::vector<u32>& source,
                          std::set<u32>& seen,
                          std::vector<u32>& target,
                          std::vector<u32>* novel = nullptr) {
    for (u32 mask : source) {
        if (seen.insert(mask).second) {
            target.push_back(mask);
            if (novel != nullptr) novel->push_back(mask);
        }
    }
}

u64 fnv_masks_round4(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

std::vector<u32> failed_bodies_round4(const std::vector<u32>& active) {
    std::vector<u32> bodies;
    bodies.reserve(EXPECTED_BODIES);
    u32 body = (UINT32_C(1) << 9) - 1;
    const u32 limit = UINT32_C(1) << 30;
    while (body != 0 && body < limit) {
        bodies.push_back(body);
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    require(bodies.size() == EXPECTED_BODIES,
            "failed-body universe changed");
    const unsigned hardware = std::thread::hardware_concurrency();
    const unsigned threadCount =
        std::max(1u, std::min(8u, hardware ? hardware : 1u));
    std::vector<std::vector<u32>> local(threadCount);
    std::vector<std::thread> workers;
    for (unsigned thread = 0; thread < threadCount; ++thread) {
        workers.emplace_back([&, thread]() {
            const std::size_t begin = bodies.size() * thread / threadCount;
            const std::size_t end = bodies.size() * (thread + 1) / threadCount;
            for (std::size_t index = begin; index < end; ++index) {
                const u32 candidate = bodies[index];
                bool covered = false;
                for (u32 repair : active) {
                    if ((candidate & repair) == 0) {
                        covered = true;
                        break;
                    }
                }
                if (!covered) local[thread].push_back(candidate);
            }
        });
    }
    for (auto& worker : workers) worker.join();
    std::vector<u32> failures;
    for (auto& part : local) {
        failures.insert(failures.end(), part.begin(), part.end());
    }
    std::sort(failures.begin(), failures.end());
    return failures;
}

struct BuiltCarrierRound4 {
    std::vector<u32> round3;
    std::set<u32> seen;
};

BuiltCarrierRound4 build_round3_carrier(
    const std::filesystem::path& packet,
    const PrefixCegar& prefix704,
    const PrefixCegar& prefix416,
    const PrefixCegar& prefix520,
    const PrefixCegar& prefix384) {
    require(prefix704.q == 416 && prefix704.r == 704 &&
                prefix704.masks.size() == 2608 &&
                prefix704.stated_fnv == UINT64_C(0x18ff663a123e684e),
            "round-one prefix changed");
    require(prefix416.q == 416 && prefix416.r == 700 &&
                prefix416.masks.size() == 5894 &&
                prefix416.stated_fnv == UINT64_C(0x701c0233c8f8abeb),
            "416,700 prefix changed");
    require(prefix520.q == 520 && prefix520.r == 700 &&
                prefix520.masks.size() == 4557 &&
                prefix520.stated_fnv == UINT64_C(0xd8466558bcd8ef9d),
            "520,700 prefix changed");
    require(prefix384.q == 384 && prefix384.r == 694 &&
                prefix384.masks.size() == 5307 &&
                prefix384.stated_fnv == UINT64_C(0x7b9a46c60e8514f8),
            "384,694 prefix changed");

    std::vector<u32> carrier;
    std::set<u32> seen;
    for (const EdgeCegar& expected : FROZEN_BAND_CEGAR) {
        const PrefixCegar prefix = read_prefix_cegar(
            packet / "replay_band" /
            (std::to_string(expected.q) + "_" +
             std::to_string(expected.r) + ".out"));
        require(prefix.q == expected.q && prefix.r == expected.r,
                "band transcript mismatch");
        append_unseen_round4(prefix.masks, seen, carrier);
    }
    require(carrier.size() == 4675 &&
                fnv_masks_round4(carrier) ==
                    UINT64_C(0xce4e76ec11df057c),
            "base carrier changed");

    append_unseen_round4(prefix704.masks, seen, carrier);
    require(carrier.size() == 4733 &&
                fnv_masks_round4(carrier) ==
                    UINT64_C(0xa7b046289655c733),
            "round-one carrier changed");

    const std::set<u32> round1_seen = seen;
    std::vector<u32> novel416;
    std::vector<u32> novel520;
    for (u32 mask : prefix416.masks) {
        if (!round1_seen.count(mask)) novel416.push_back(mask);
    }
    for (u32 mask : prefix520.masks) {
        if (!round1_seen.count(mask)) novel520.push_back(mask);
    }
    require(novel416.size() == 2970 && novel520.size() == 1194,
            "round-two source novelty changed");
    append_unseen_round4(novel416, seen, carrier);
    append_unseen_round4(novel520, seen, carrier);
    require(carrier.size() == 7986 &&
                fnv_masks_round4(carrier) ==
                    UINT64_C(0xbaef1d2f49444638),
            "round-two carrier changed");

    append_unseen_round4(prefix384.masks, seen, carrier);
    require(carrier.size() == 8319 &&
                fnv_masks_round4(carrier) ==
                    UINT64_C(0xe08b227730f6793c),
            "round-three carrier changed");
    return {std::move(carrier), std::move(seen)};
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 10,
                "usage: round4 REPLAY_BAND PREFIX_416_704 PREFIX_416_700 "
                "PREFIX_520_700 PREFIX_384_694 PREFIX_520_688 RESIDUAL_CSV "
                "START_ENDPOINT MIN_ENDPOINT");
        const PrefixCegar prefix704 = read_prefix_cegar(argv[2]);
        const PrefixCegar prefix416 = read_prefix_cegar(argv[3]);
        const PrefixCegar prefix520 = read_prefix_cegar(argv[4]);
        const PrefixCegar prefix384 = read_prefix_cegar(argv[5]);
        const PrefixCegar prefix688 = read_prefix_cegar(argv[6]);
        BuiltCarrierRound4 built = build_round3_carrier(
            argv[1], prefix704, prefix416, prefix520, prefix384);
        require(prefix688.q == 520 && prefix688.r == 688 &&
                    prefix688.masks.size() == 5398 &&
                    prefix688.stated_fnv == UINT64_C(0x6ab471da88c8e1d1),
                "520,688 discovery prefix changed");
        const int startEndpoint = std::stoi(argv[8]);
        const int minEndpoint = std::stoi(argv[9]);
        require(startEndpoint >= minEndpoint, "bad endpoint range");

        const ScanResult fullPrefixScan = scan_bodies(prefix688.masks);
        require(fullPrefixScan.failures == 0 &&
                    fullPrefixScan.max_checks == prefix688.masks.size(),
                "520,688 prefix lost order-minimal closure");

        std::vector<u32> novel;
        std::vector<u32> round4 = built.round3;
        append_unseen_round4(prefix688.masks, built.seen, round4, &novel);
        require(novel.size() == 199 && round4.size() == 8518,
                "round-four novelty count changed");
        const u64 novelFnv = fnv_masks_round4(novel);
        const u64 round4Fnv = fnv_masks_round4(round4);
        require(novelFnv == UINT64_C(0x772162cabf9167f9) &&
                    round4Fnv == UINT64_C(0x1603e3fe970f8428),
                "round-four carrier fingerprint changed");

        const std::vector<Cell> poolCells = build_pool_cells();
        require(poolCells.size() == 7133, "pool cell count changed");
        const CarrierScanCegar oldSeed = scan_carrier_cegar(
            poolCells, 520, 688, built.round3);
        const CarrierScanCegar newSeed = scan_carrier_cegar(
            poolCells, 520, 688, round4);
        require(oldSeed.active == 5934 && oldSeed.body.failures == 2 &&
                    oldSeed.body.first_failure == UINT32_C(0x07187008),
                "round-three seed resistance changed");
        require(newSeed.body.failures == 0,
                "round-four carrier does not close seed");

        const i64 seedG = std::gcd(520, 688);
        const PrimitivePair seedPrimitive = build_primitive(65, 86);
        const AtomData seedAtoms = build_cocycle_atoms(
            poolCells, seedPrimitive, seedG);
        const i128 seedDenominator =
            static_cast<i128>(seedPrimitive.grid) * seedG * COMMON;
        FnvLocal atomComponentLedger;
        u64 activeNovel = 0;
        u64 disjointFirst = 0;
        u64 disjointSecond = 0;
        u32 firstWitness = 0;
        u32 secondFailure = 0;
        {
            std::vector<u32> activeOld;
            for (u32 mask : built.round3) {
                const i128 margin = static_cast<i128>(63) *
                                        atom_mass_cegar(seedAtoms, mask) -
                                    static_cast<i128>(4) * seedDenominator;
                if (margin >= 0) activeOld.push_back(mask);
            }
            const std::vector<u32> failures =
                failed_bodies_round4(activeOld);
            require(failures.size() == 2 &&
                        failures.front() == oldSeed.body.first_failure,
                    "old active-deck body census changed");
            secondFailure = failures.back();
            require(secondFailure != 0,
                    "failed-body reconstruction changed");
        }
        for (u32 mask : novel) {
            const i128 mass = atom_mass_cegar(seedAtoms, mask);
            const i128 margin = static_cast<i128>(63) * mass -
                                static_cast<i128>(4) * seedDenominator;
            require(margin >= 0,
                    "source-prefix novel mask inactive at seed");
            const ComponentReport component = component_report(
                poolCells, seedPrimitive, seedG, mask);
            require(component.exact_mass_num == mass,
                    "round-four atom/component mismatch");
            ++activeNovel;
            if ((mask & oldSeed.body.first_failure) == 0) {
                ++disjointFirst;
                if (firstWitness == 0) firstWitness = mask;
            }
            if ((mask & secondFailure) == 0) ++disjointSecond;
            atomComponentLedger.add(mask);
            fnv_add_i128_cegar(atomComponentLedger, mass);
            fnv_add_i128_cegar(atomComponentLedger, margin);
        }
        require(firstWitness != 0 && disjointFirst != 0 &&
                    disjointSecond != 0,
                "round-four novelty misses old hostile bodies");

        std::cout << "CARRIER_CEGAR_ROUND4_V1\n";
        std::cout << "ROUND3 UNION " << built.round3.size() << " FNV "
                  << std::hex << fnv_masks_round4(built.round3) << std::dec
                  << "\n";
        std::cout << "SEED_FULL_DISCOVERY PAIR 520,688 PREFIX "
                  << prefix688.masks.size() << " PREFIX_FNV " << std::hex
                  << prefix688.stated_fnv << std::dec << " OVERLAP_ROUND3 "
                  << prefix688.masks.size() - novel.size() << " NOVEL "
                  << novel.size() << " NOVEL_FNV " << std::hex << novelFnv
                  << std::dec << " PREFIX_MAX_CHECKS "
                  << fullPrefixScan.max_checks << "\n";
        std::cout << "OLD_SEED ACTIVE " << oldSeed.active << " FAILURES "
                  << oldSeed.body.failures << " FIRST_FAILURE " << std::hex
                  << oldSeed.body.first_failure << " SECOND_FAILURE "
                  << secondFailure << std::dec << " FIRST_LABELS {"
                  << labels(oldSeed.body.first_failure) << "} SECOND_LABELS {"
                  << labels(secondFailure) << "}\n";
        std::cout << "ROUND4 UNION " << round4.size() << " FNV " << std::hex
                  << round4Fnv << " ATOM_COMPONENT_LEDGER "
                  << atomComponentLedger.state << std::dec
                  << " ACTIVE_NOVEL " << activeNovel
                  << " DISJOINT_FIRST " << disjointFirst
                  << " DISJOINT_SECOND " << disjointSecond
                  << " FIRST_WITNESS " << std::hex << firstWitness << std::dec
                  << " FIRST_WITNESS_LABELS {" << labels(firstWitness) << "}"
                  << " NEW_SEED_ACTIVE " << newSeed.active
                  << " NEW_SEED_FAILURES " << newSeed.body.failures
                  << " NEW_SEED_CHECKS " << newSeed.body.checks
                  << " NEW_SEED_MAX_CHECKS " << newSeed.body.max_checks
                  << "\n";

        const std::vector<EdgeCegar> edges = read_edges_cegar(argv[7]);
        int currentLayer = std::numeric_limits<int>::max();
        u64 rows = 0;
        u64 resistant = 0;
        u64 activeSum = 0;
        u64 checks = 0;
        u64 maxChecks = 0;
        FnvLocal layerLedger;
        std::vector<std::string> hostileRows;
        bool found = false;
        auto flush = [&]() {
            if (currentLayer == std::numeric_limits<int>::max()) return false;
            std::cout << "LAYER " << currentLayer << " ROWS " << rows
                      << " RESISTANT " << resistant << " ACTIVE_SUM "
                      << activeSum << " CHECKS " << checks << " MAX_CHECKS "
                      << maxChecks << " ROW_LEDGER " << std::hex
                      << layerLedger.state << std::dec << "\n";
            for (const std::string& row : hostileRows) std::cout << row << '\n';
            return resistant != 0;
        };
        for (const EdgeCegar& edge : edges) {
            if (edge.r > startEndpoint || edge.r < minEndpoint) continue;
            if (currentLayer != edge.r) {
                if (flush()) {
                    found = true;
                    break;
                }
                currentLayer = edge.r;
                rows = resistant = activeSum = checks = maxChecks = 0;
                layerLedger = FnvLocal{};
                hostileRows.clear();
            }
            const CarrierScanCegar scan = scan_carrier_cegar(
                poolCells, edge.q, edge.r, round4);
            ++rows;
            if (scan.body.failures != 0) ++resistant;
            activeSum += scan.active;
            checks += scan.body.checks;
            maxChecks = std::max(maxChecks, scan.body.max_checks);
            layerLedger.add(static_cast<u64>(edge.q));
            layerLedger.add(static_cast<u64>(edge.r));
            layerLedger.add(scan.active);
            layerLedger.add(scan.body.failures);
            layerLedger.add(scan.body.checks);
            layerLedger.add(scan.body.max_checks);
            layerLedger.add(scan.body.first_failure);
            layerLedger.add(scan.best_disjoint_repair);
            fnv_add_i128_cegar(layerLedger, scan.best_disjoint_margin);
            fnv_add_i128_cegar(layerLedger, scan.denominator);
            if (scan.body.failures != 0) {
                std::ostringstream row;
                row << "RESISTANT PAIR " << edge.q << ',' << edge.r
                    << " ACTIVE " << scan.active << " FAILURES "
                    << scan.body.failures << " FIRST_FAILURE " << std::hex
                    << scan.body.first_failure << std::dec
                    << " FIRST_FAILURE_LABELS {"
                    << labels(scan.body.first_failure) << "}"
                    << " BEST_DISJOINT_REPAIR " << std::hex
                    << scan.best_disjoint_repair << std::dec
                    << " BEST_DISJOINT_LABELS {"
                    << labels(scan.best_disjoint_repair) << "}"
                    << " BEST_DISJOINT_MARGIN_NUM "
                    << decimal(scan.best_disjoint_margin) << " DEN "
                    << decimal(scan.denominator);
                hostileRows.push_back(row.str());
            }
        }
        if (!found) found = flush();
        std::cout << "STOP FIRST_RESISTANT_LAYER "
                  << (found ? currentLayer : -1) << " VERDICT "
                  << (found ? "BOUNDARY_FOUND" : "RANGE_CLOSED") << "\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "CARRIER_CEGAR_ROUND4_ERROR " << error.what() << '\n';
        return 1;
    }
}

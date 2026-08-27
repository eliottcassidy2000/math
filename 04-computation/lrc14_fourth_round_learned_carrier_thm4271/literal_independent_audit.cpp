// Detached literal-joint-wall audit for the THM-4271 round-four carrier.
// No endpoint atoms or cocycle antiderivatives are used for activation.

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif
#define CARRIER_CLEANROOM_CONTROL_LIBRARY_ONLY
#include "04-computation/lrc14_three_round_learned_carrier_thm4266/cleanroom_resistance_controls.cpp"
#undef CARRIER_CLEANROOM_CONTROL_LIBRARY_ONLY
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

namespace {

void append_unseen_literal(const std::vector<u32>& source,
                           std::set<u32>& seen,
                           std::vector<u32>& target,
                           std::vector<u32>* novel = nullptr) {
    for (u32 repair : source) {
        if (seen.insert(repair).second) {
            target.push_back(repair);
            if (novel != nullptr) novel->push_back(repair);
        }
    }
}

u64 fnv_literal(const std::vector<u32>& masks) {
    Fnv ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

std::vector<u32> active_literal(int q, int r,
                                const std::vector<u32>& carrier) {
    const Geometry geometry = build_joint_geometry(q, r);
    std::vector<u32> active;
    active.reserve(carrier.size());
    for (u32 repair : carrier) {
        i64 mass = 0;
        for (const Cell& cell : geometry.cells) {
            if (cell.pair_safe && (cell.failed_pool & ~repair) == 0) {
                mass += cell.width;
            }
        }
        if (static_cast<i128>(63) * mass >=
            static_cast<i128>(4) * geometry.grid) {
            active.push_back(repair);
        }
    }
    return active;
}

std::vector<u32> failures_literal(const std::vector<u32>& bodies,
                                  const std::vector<u32>& active) {
    const unsigned hardware = std::thread::hardware_concurrency();
    const unsigned threadCount =
        std::max(1u, std::min(8u, hardware == 0 ? 1u : hardware));
    std::vector<std::vector<u32>> local(threadCount);
    std::vector<std::thread> workers;
    for (unsigned thread = 0; thread < threadCount; ++thread) {
        workers.emplace_back([&, thread]() {
            const std::size_t begin = bodies.size() * thread / threadCount;
            const std::size_t end = bodies.size() * (thread + 1) / threadCount;
            for (std::size_t index = begin; index < end; ++index) {
                const u32 body = bodies[index];
                bool covered = false;
                for (u32 repair : active) {
                    if ((body & repair) == 0) {
                        covered = true;
                        break;
                    }
                }
                if (!covered) local[thread].push_back(body);
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

struct ClosedLiteral {
    u64 cells = 0;
    u64 active = 0;
    BodyAudit body;
};

ClosedLiteral closed_literal(int q, int r,
                             const std::vector<u32>& carrier,
                             const std::vector<u32>& bodies) {
    const Geometry geometry = build_joint_geometry(q, r);
    const std::vector<u32> active = active_literal(q, r, carrier);
    const BodyAudit body = audit_bodies(bodies, active);
    require(body.bodies == EXPECTED_BODIES && body.failures == 0,
            "literal carrier unexpectedly fails");
    return {static_cast<u64>(geometry.cells.size()), active.size(), body};
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 7,
                "usage: literal REPLAY_BAND PREFIX_416_704 PREFIX_416_700 "
                "PREFIX_520_700 PREFIX_384_694 PREFIX_520_688");
        const std::filesystem::path band(argv[1]);
        const Parsed prefix704 = parse_probe(argv[2]);
        const Parsed prefix416 = parse_probe(argv[3]);
        const Parsed prefix520 = parse_probe(argv[4]);
        const Parsed prefix384 = parse_probe(argv[5]);
        const Parsed prefix688 = parse_probe(argv[6]);
        require(prefix704.q == 416 && prefix704.r == 704 &&
                    prefix704.deck.size() == 2608 &&
                    prefix416.q == 416 && prefix416.r == 700 &&
                    prefix416.deck.size() == 5894 &&
                    prefix520.q == 520 && prefix520.r == 700 &&
                    prefix520.deck.size() == 4557 &&
                    prefix384.q == 384 && prefix384.r == 694 &&
                    prefix384.deck.size() == 5307 &&
                    prefix688.q == 520 && prefix688.r == 688 &&
                    prefix688.deck.size() == 5398 &&
                    prefix688.declared_prefix_fnv ==
                        UINT64_C(0x6ab471da88c8e1d1),
                "literal source prefix identity changed");

        std::vector<u32> carrier;
        std::set<u32> seen;
        for (const auto& [q, r] : EXPECTED_PAIRS) {
            const Parsed parsed = parse_probe(
                band / (std::to_string(q) + "_" +
                        std::to_string(r) + ".out"));
            require(parsed.q == q && parsed.r == r,
                    "literal band transcript mismatch");
            append_unseen_literal(parsed.deck, seen, carrier);
        }
        require(carrier.size() == 4675 &&
                    fnv_literal(carrier) ==
                        UINT64_C(0xce4e76ec11df057c),
                "literal base carrier changed");
        append_unseen_literal(prefix704.deck, seen, carrier);
        require(carrier.size() == 4733 &&
                    fnv_literal(carrier) ==
                        UINT64_C(0xa7b046289655c733),
                "literal round-one carrier changed");
        const std::set<u32> round1Seen = seen;
        std::vector<u32> novel416;
        std::vector<u32> novel520;
        for (u32 repair : prefix416.deck) {
            if (!round1Seen.count(repair)) novel416.push_back(repair);
        }
        for (u32 repair : prefix520.deck) {
            if (!round1Seen.count(repair)) novel520.push_back(repair);
        }
        append_unseen_literal(novel416, seen, carrier);
        append_unseen_literal(novel520, seen, carrier);
        require(carrier.size() == 7986 &&
                    fnv_literal(carrier) ==
                        UINT64_C(0xbaef1d2f49444638),
                "literal round-two carrier changed");
        append_unseen_literal(prefix384.deck, seen, carrier);
        require(carrier.size() == 8319 &&
                    fnv_literal(carrier) ==
                        UINT64_C(0xe08b227730f6793c),
                "literal round-three carrier changed");
        const std::vector<u32> round3 = carrier;

        std::vector<u32> novel688;
        append_unseen_literal(prefix688.deck, seen, carrier, &novel688);
        require(novel688.size() == 199 && carrier.size() == 8518 &&
                    fnv_literal(novel688) ==
                        UINT64_C(0x772162cabf9167f9) &&
                    fnv_literal(carrier) ==
                        UINT64_C(0x1603e3fe970f8428),
                "literal round-four carrier changed");

        std::vector<u32> bodies;
        bodies.reserve(EXPECTED_BODIES);
        enumerate_bodies(0, 9, 0, bodies);
        require(bodies.size() == EXPECTED_BODIES,
                "literal body universe changed");

        std::cout << "THM4271_ROUND4_LITERAL_INDEPENDENT_AUDIT_V1\n";
        std::cout << "ROUND3 UNION " << round3.size() << " FNV " << std::hex
                  << fnv_literal(round3) << std::dec << " ROUND4_NOVEL "
                  << novel688.size() << " NOVEL_FNV " << std::hex
                  << fnv_literal(novel688) << std::dec << " ROUND4_UNION "
                  << carrier.size() << " UNION_FNV " << std::hex
                  << fnv_literal(carrier) << std::dec << "\n";

        Fnv labelledPrefix;
        labelledPrefix.add(520);
        labelledPrefix.add(688);
        labelledPrefix.add(prefix688.deck.size());
        const PairAudit fullPrefix = audit_pair(
            prefix688, bodies, labelledPrefix);

        const DirectScanControl oldSeed = direct_scan_control(
            520, 688, round3, bodies);
        require(oldSeed.active == 5934 && oldSeed.failures == 2 &&
                    oldSeed.first_failure == UINT32_C(0x07187008) &&
                    oldSeed.best_repair == UINT32_C(0x10458a01),
                "literal old seed resistance changed");
        const std::vector<u32> activeOld = active_literal(520, 688, round3);
        const std::vector<u32> failedOld = failures_literal(bodies, activeOld);
        require(failedOld == std::vector<u32>{UINT32_C(0x07187008),
                                              UINT32_C(0x27503008)},
                "literal old failed bodies changed");

        const std::vector<u32> activeNovel = active_literal(
            520, 688, novel688);
        require(activeNovel.size() == novel688.size(),
                "literal novel source mask inactive");
        u64 disjointFirst = 0;
        u64 disjointSecond = 0;
        for (u32 repair : activeNovel) {
            if ((repair & failedOld[0]) == 0) ++disjointFirst;
            if ((repair & failedOld[1]) == 0) ++disjointSecond;
        }
        require(disjointFirst == 1 && disjointSecond == 1,
                "literal novel disjointness changed");

        const ClosedLiteral newSeed = closed_literal(
            520, 688, carrier, bodies);
        require(newSeed.active == 6133,
                "literal new seed active count changed");
        std::cout << "OLD_SEED PAIR 520,688 JOINT_CELLS "
                  << build_joint_geometry(520,688).cells.size()
                  << " ACTIVE " << oldSeed.active << " FAILURES "
                  << oldSeed.failures << " FIRST_FAILURE " << std::hex
                  << failedOld[0] << " SECOND_FAILURE " << failedOld[1]
                  << std::dec << "\n";
        std::cout << "NOVEL_TRANSFER ACTIVE " << activeNovel.size()
                  << " DISJOINT_FIRST " << disjointFirst
                  << " DISJOINT_SECOND " << disjointSecond << "\n";
        std::cout << "NEW_SEED PAIR 520,688 JOINT_CELLS " << newSeed.cells
                  << " ACTIVE " << newSeed.active << " BODIES "
                  << newSeed.body.bodies << " FAILURES "
                  << newSeed.body.failures << " CHECKS "
                  << newSeed.body.checks << " MAX_CHECKS "
                  << newSeed.body.max_checks << "\n";

        check_and_print_control(
            256, 671, carrier, bodies, 4290, 26,
            UINT32_C(0x06166401), UINT32_C(0x39201280),
            "89508336872746176", "43867507598953758720");
        check_and_print_control(
            384, 671, carrier, bodies, 5988, 1,
            UINT32_C(0x0d186401), UINT32_C(0x10479200),
            "721295701293847680", "65801261398430638080");

        std::cout << "FULL_PREFIX JOINT_CELLS " << fullPrefix.joint_cells
                  << " REPAIRS " << fullPrefix.prefix_size
                  << " BODY_CHECKS " << fullPrefix.body_checks << "\n";
        std::cout << "VERDICT PASS LITERAL_PREFIX_SEED_TRANSFER_BOUNDARY\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "THM4271_LITERAL_ERROR " << error.what() << '\n';
        return 1;
    }
}

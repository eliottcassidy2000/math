// Literal-wall and recursive-body controls for the second and third CEGAR
// rounds.

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#endif
#define CARRIER_CLEANROOM_CONTROL_LIBRARY_ONLY
#include "cleanroom_resistance_controls.cpp"
#undef CARRIER_CLEANROOM_CONTROL_LIBRARY_ONLY
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

namespace {

struct ClosedControl {
    u64 active = 0;
    BodyAudit body;
    u64 joint_cells = 0;
};

ClosedControl direct_closed_control(int q, int r,
                                    const std::vector<u32>& carrier,
                                    const std::vector<u32>& bodies) {
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
    const BodyAudit body = audit_bodies(bodies, active);
    require(body.bodies == EXPECTED_BODIES && body.failures == 0,
            "round-two literal carrier does not close seed");
    return {active.size(), body,
            static_cast<u64>(geometry.cells.size())};
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 5 || argc == 6,
                "usage: cleanroom-round2 REPLAY_BAND PREFIX_416_704 "
                "PREFIX_416_700 PREFIX_520_700 [PREFIX_384_694]");
        const std::filesystem::path band(argv[1]);
        std::vector<u32> base;
        std::set<u32> seen;
        Fnv base_ledger;
        for (const auto& [q, r] : EXPECTED_PAIRS) {
            const Parsed parsed = parse_probe(
                band / (std::to_string(q) + "_" +
                        std::to_string(r) + ".out"));
            require(parsed.q == q && parsed.r == r,
                    "round-two band transcript mismatch");
            for (u32 repair : parsed.deck) {
                if (seen.insert(repair).second) {
                    base.push_back(repair);
                    base_ledger.add(repair);
                }
            }
        }
        require(base.size() == 4675 &&
                    base_ledger.state == UINT64_C(0xce4e76ec11df057c),
                "round-two literal base carrier changed");

        const Parsed prefix_704 = parse_probe(argv[2]);
        const Parsed prefix_416 = parse_probe(argv[3]);
        const Parsed prefix_520 = parse_probe(argv[4]);
        require(prefix_704.q == 416 && prefix_704.r == 704 &&
                    prefix_704.deck.size() == 2608 &&
                    prefix_416.q == 416 && prefix_416.r == 700 &&
                    prefix_416.deck.size() == 5894 &&
                    prefix_520.q == 520 && prefix_520.r == 700 &&
                    prefix_520.deck.size() == 4557,
                "round-two literal prefix identity changed");

        std::vector<u32> round1 = base;
        for (u32 repair : prefix_704.deck) {
            if (seen.insert(repair).second) round1.push_back(repair);
        }
        require(round1.size() == 4733,
                "round-two literal round-one count changed");

        std::set<u32> round1_seen = seen;
        std::vector<u32> novel_416;
        std::vector<u32> novel_520;
        for (u32 repair : prefix_416.deck) {
            if (!round1_seen.count(repair)) novel_416.push_back(repair);
        }
        for (u32 repair : prefix_520.deck) {
            if (!round1_seen.count(repair)) novel_520.push_back(repair);
        }
        require(novel_416.size() == 2970 && novel_520.size() == 1194,
                "round-two literal novel counts changed");

        std::vector<u32> round2 = round1;
        for (u32 repair : novel_416) {
            require(seen.insert(repair).second,
                    "round-two literal 416 novelty changed");
            round2.push_back(repair);
        }
        u64 only_520 = 0;
        for (u32 repair : novel_520) {
            if (seen.insert(repair).second) {
                round2.push_back(repair);
                ++only_520;
            }
        }
        require(only_520 == 283 && round2.size() == 7986,
                "round-two literal union changed");
        Fnv round2_ledger;
        for (u32 repair : round2) round2_ledger.add(repair);
        require(round2_ledger.state == UINT64_C(0xbaef1d2f49444638),
                "round-two literal union FNV changed");

        const bool has_round3 = argc == 6;
        Parsed prefix_384;
        std::vector<u32> round3 = round2;
        u64 round3_novel = 0;
        Fnv round3_novel_ledger;
        u64 round3_fnv = 0;
        if (has_round3) {
            prefix_384 = parse_probe(argv[5]);
            require(prefix_384.q == 384 && prefix_384.r == 694 &&
                        prefix_384.deck.size() == 5307 &&
                        prefix_384.declared_prefix_fnv ==
                            UINT64_C(0x7b9a46c60e8514f8),
                    "round-three literal prefix changed");
            for (u32 repair : prefix_384.deck) {
                if (seen.insert(repair).second) {
                    round3.push_back(repair);
                    ++round3_novel;
                    round3_novel_ledger.add(repair);
                }
            }
            require(round3_novel == 333 && round3.size() == 8319,
                    "round-three literal union count changed");
            require(round3_novel_ledger.state ==
                        UINT64_C(0x1d1422eaa770226d),
                    "round-three literal novel FNV changed");
            Fnv round3_ledger;
            for (u32 repair : round3) round3_ledger.add(repair);
            require(round3_ledger.state == UINT64_C(0xe08b227730f6793c),
                    "round-three literal union FNV changed");
            round3_fnv = round3_ledger.state;
        }

        std::vector<u32> bodies;
        bodies.reserve(EXPECTED_BODIES);
        enumerate_bodies(0, 9, 0, bodies);
        require(bodies.size() == EXPECTED_BODIES,
                "round-two literal body universe changed");

        std::cout << (has_round3
                          ? "CARRIER_CEGAR_ROUND23_CLEANROOM_V1\n"
                          : "CARRIER_CEGAR_ROUND2_CLEANROOM_V1\n");
        std::cout << "ROUND2 UNION " << round2.size() << " FNV " << std::hex
                  << round2_ledger.state << std::dec
                  << " NOVEL_416 " << novel_416.size()
                  << " NOVEL_520 " << novel_520.size()
                  << " APPENDED_520_ONLY " << only_520 << "\n";

        Fnv labelled_416;
        labelled_416.add(416);
        labelled_416.add(700);
        labelled_416.add(prefix_416.deck.size());
        const PairAudit full_416 = audit_pair(
            prefix_416, bodies, labelled_416);
        Fnv labelled_520;
        labelled_520.add(520);
        labelled_520.add(700);
        labelled_520.add(prefix_520.deck.size());
        const PairAudit full_520 = audit_pair(
            prefix_520, bodies, labelled_520);
        PairAudit full_384;
        if (has_round3) {
            Fnv labelled_384;
            labelled_384.add(384);
            labelled_384.add(694);
            labelled_384.add(prefix_384.deck.size());
            full_384 = audit_pair(prefix_384, bodies, labelled_384);
        }

        const ClosedControl closed_416 = direct_closed_control(
            416, 700, round2, bodies);
        const ClosedControl closed_520 = direct_closed_control(
            520, 700, round2, bodies);
        require(closed_416.active == 5894 && closed_520.active == 6514,
                "round-two literal active counts disagree with atoms");
        std::cout << "ROUND2_CARRIER_PAIR 416,700 JOINT_CELLS "
                  << closed_416.joint_cells << " ACTIVE " << closed_416.active
                  << " BODIES " << closed_416.body.bodies
                  << " FAILURES " << closed_416.body.failures
                  << " CHECKS " << closed_416.body.checks
                  << " MAX_CHECKS " << closed_416.body.max_checks << "\n";
        std::cout << "ROUND2_CARRIER_PAIR 520,700 JOINT_CELLS "
                  << closed_520.joint_cells << " ACTIVE " << closed_520.active
                  << " BODIES " << closed_520.body.bodies
                  << " FAILURES " << closed_520.body.failures
                  << " CHECKS " << closed_520.body.checks
                  << " MAX_CHECKS " << closed_520.body.max_checks << "\n";

        check_and_print_control(
            384, 694, round2, bodies, 5966, 1, UINT32_C(0x0d186401),
            UINT32_C(0x30209248), "764085760585922880",
            "34028372138979778560");
        if (has_round3) {
            const ClosedControl closed_384 = direct_closed_control(
                384, 694, round3, bodies);
            require(closed_384.active == 6299,
                    "round-three literal seed active count changed");
            std::cout << "ROUND3 UNION " << round3.size() << " FNV "
                      << std::hex << round3_fnv << std::dec
                      << " NOVEL " << round3_novel << "\n";
            std::cout << "ROUND3_CARRIER_PAIR 384,694 JOINT_CELLS "
                      << closed_384.joint_cells << " ACTIVE "
                      << closed_384.active << " BODIES "
                      << closed_384.body.bodies << " FAILURES "
                      << closed_384.body.failures << " CHECKS "
                      << closed_384.body.checks << " MAX_CHECKS "
                      << closed_384.body.max_checks << "\n";
            check_and_print_control(
                520, 688, round3, bodies, 5934, 2,
                UINT32_C(0x07187008), UINT32_C(0x10458a01),
                "71286709268960640", "11420425087469798400");
        }
        std::cout << "FULL_PREFIX_TOTAL_CELLS "
                  << full_416.joint_cells + full_520.joint_cells +
                         full_384.joint_cells
                  << " FULL_PREFIX_TOTAL_REPAIRS "
                  << full_416.prefix_size + full_520.prefix_size +
                         full_384.prefix_size
                  << "\n";
        std::cout << (has_round3
                          ? "VERDICT PASS ROUND23_LITERAL_WALL_FULL_PREFIXES_"
                            "CARRIER_CLOSURES_NEXT_RESISTANCE\n"
                          : "VERDICT PASS ROUND2_LITERAL_WALL_FULL_PREFIXES_"
                            "CARRIER_CLOSURES_NEXT_RESISTANCE\n");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "CARRIER_CEGAR_ROUND2_CLEANROOM_ERROR "
                  << error.what() << '\n';
        return 1;
    }
}

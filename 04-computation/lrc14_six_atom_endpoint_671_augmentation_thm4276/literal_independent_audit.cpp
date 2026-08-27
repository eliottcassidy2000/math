// Detached literal-joint-wall audit for the compact round-five augmentation.
// Activation is computed from fresh joint walls, never endpoint atoms.

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

constexpr std::array<u32, 26> ROUND4_FAIL_256 = {{
    0x06166401,0x07067001,0x07106409,0x07126088,0x07126401,0x07162401,
    0x07163400,0x07166008,0x0d107401,0x0d10e401,0x0d146401,0x0d186401,
    0x0d246401,0x0d506401,0x0d906401,0x0f106401,0x0f142401,0x0f142408,
    0x0f143400,0x0f146008,0x15923400,0x17162400,0x17922008,0x1d106401,
    0x1d902401,0x1f142400
}};
constexpr u32 ROUND4_FAIL_384 = UINT32_C(0x0d186401);
constexpr std::array<u32, 6> AUGMENTATION = {{
    0x00289285,0x0260812c,0x18689040,
    0x20c0c124,0x302c1006,0x30580888
}};
constexpr std::array<u32, 6> EXPECTED_PATTERNS = {{
    0x026a0080,0x05904d00,0x000000bd,
    0x02270060,0x0000e21c,0x00001002
}};

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

u32 literal_pattern(const std::vector<u32>& active256,
                    const std::vector<u32>& active384,
                    u32 repair) {
    u32 pattern = 0;
    if (std::find(active256.begin(), active256.end(), repair) !=
        active256.end()) {
        for (unsigned index = 0; index < ROUND4_FAIL_256.size(); ++index) {
            if ((repair & ROUND4_FAIL_256[index]) == 0) {
                pattern |= UINT32_C(1) << index;
            }
        }
    }
    if (std::find(active384.begin(), active384.end(), repair) !=
            active384.end() &&
        (repair & ROUND4_FAIL_384) == 0) {
        pattern |= UINT32_C(1) << 26;
    }
    return pattern;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 9,
                "usage: literal REPLAY_BAND PREFIX_416_704 PREFIX_416_700 "
                "PREFIX_520_700 PREFIX_384_694 PREFIX_520_688 "
                "PREFIX_256_671 PREFIX_384_671");
        const std::filesystem::path band(argv[1]);
        const Parsed prefix704 = parse_probe(argv[2]);
        const Parsed prefix416 = parse_probe(argv[3]);
        const Parsed prefix520 = parse_probe(argv[4]);
        const Parsed prefix384 = parse_probe(argv[5]);
        const Parsed prefix688 = parse_probe(argv[6]);
        const Parsed prefix256b = parse_probe(argv[7]);
        const Parsed prefix384b = parse_probe(argv[8]);
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
                    prefix256b.q == 256 && prefix256b.r == 671 &&
                    prefix256b.deck.size() == 13086 &&
                    prefix256b.declared_prefix_fnv ==
                        UINT64_C(0x06c1f44ec0661d8c) &&
                    prefix384b.q == 384 && prefix384b.r == 671 &&
                    prefix384b.deck.size() == 6986 &&
                    prefix384b.declared_prefix_fnv ==
                        UINT64_C(0x3a02a5774d3641ab),
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
        append_unseen_literal(prefix688.deck, seen, carrier);
        require(carrier.size() == 8518 &&
                    fnv_literal(carrier) ==
                        UINT64_C(0x1603e3fe970f8428),
                "literal round-four carrier changed");
        const std::vector<u32> round4 = carrier;

        std::vector<u32> compact = round4;
        std::vector<u32> novel;
        append_unseen_literal(
            std::vector<u32>(AUGMENTATION.begin(), AUGMENTATION.end()),
            seen, compact, &novel);
        require(novel == std::vector<u32>(AUGMENTATION.begin(),
                                          AUGMENTATION.end()) &&
                    compact.size() == 8524 &&
                    fnv_literal(compact) ==
                        UINT64_C(0x5ddb84a44f5d2ad7),
                "literal compact carrier changed");

        std::vector<u32> bodies;
        bodies.reserve(EXPECTED_BODIES);
        enumerate_bodies(0, 9, 0, bodies);
        require(bodies.size() == EXPECTED_BODIES,
                "literal body universe changed");

        const std::vector<u32> activeOld256 = active_literal(
            256, 671, round4);
        const std::vector<u32> activeOld384 = active_literal(
            384, 671, round4);
        require(activeOld256.size() == 4290 &&
                    activeOld384.size() == 5988,
                "literal old seed active census changed");
        const std::vector<u32> failed256 = failures_literal(
            bodies, activeOld256);
        const std::vector<u32> failed384 = failures_literal(
            bodies, activeOld384);
        require(failed256 == std::vector<u32>(ROUND4_FAIL_256.begin(),
                                              ROUND4_FAIL_256.end()) &&
                    failed384 == std::vector<u32>{ROUND4_FAIL_384},
                "literal obligation universe changed");

        u32 unionPattern = 0;
        Fnv patternLedger;
        for (std::size_t index = 0; index < AUGMENTATION.size(); ++index) {
            const u32 pattern = literal_pattern(
                active_literal(256, 671,
                               std::vector<u32>{AUGMENTATION[index]}),
                active_literal(384, 671,
                               std::vector<u32>{AUGMENTATION[index]}),
                AUGMENTATION[index]);
            require(pattern == EXPECTED_PATTERNS[index],
                    "literal compact pattern changed");
            unionPattern |= pattern;
            patternLedger.add(AUGMENTATION[index]);
            patternLedger.add(pattern);
        }
        require(unionPattern == UINT32_C(0x07ffffff) &&
                    patternLedger.state ==
                        UINT64_C(0x5e7cc39cf6cee4be),
                "literal selected-pattern certificate changed");

        const ClosedLiteral closed256 = closed_literal(
            256, 671, compact, bodies);
        const ClosedLiteral closed384 = closed_literal(
            384, 671, compact, bodies);
        require(closed256.active == 4296 && closed384.active == 5994,
                "literal compact seed active counts changed");

        std::cout << "ROUND5_COMPACT_LITERAL_INDEPENDENT_AUDIT_V1\n";
        std::cout << "ROUND4 UNION " << round4.size() << " FNV " << std::hex
                  << fnv_literal(round4) << std::dec << " COMPACT_NOVEL "
                  << novel.size() << " COMPACT_UNION " << compact.size()
                  << " COMPACT_FNV " << std::hex << fnv_literal(compact)
                  << std::dec << "\n";
        std::cout << "OBLIGATIONS PAIR256 " << failed256.size()
                  << " PAIR384 " << failed384.size() << "\n";
        std::cout << "SELECTED_PATTERN_LEDGER " << std::hex
                  << patternLedger.state << " UNION " << unionPattern
                  << std::dec << "\n";
        std::cout << "SEED 256,671 ACTIVE " << closed256.active
                  << " BODIES " << closed256.body.bodies << " FAILURES "
                  << closed256.body.failures << " CHECKS "
                  << closed256.body.checks << " MAX_CHECKS "
                  << closed256.body.max_checks << "\n";
        std::cout << "SEED 384,671 ACTIVE " << closed384.active
                  << " BODIES " << closed384.body.bodies << " FAILURES "
                  << closed384.body.failures << " CHECKS "
                  << closed384.body.checks << " MAX_CHECKS "
                  << closed384.body.max_checks << "\n";

        Fnv labelled256;
        labelled256.add(256);
        labelled256.add(671);
        labelled256.add(prefix256b.deck.size());
        const PairAudit full256 = audit_pair(
            prefix256b, bodies, labelled256);
        Fnv labelled384;
        labelled384.add(384);
        labelled384.add(671);
        labelled384.add(prefix384b.deck.size());
        const PairAudit full384 = audit_pair(
            prefix384b, bodies, labelled384);

        check_and_print_control(
            256, 670, compact, bodies, 2172, 1659,
            UINT32_C(0x01147409), UINT32_C(0x2ee00040),
            "472438974077602560", "21901065641802547200");
        check_and_print_control(
            384, 670, compact, bodies, 3851, 46,
            UINT32_C(0x043e6400), UINT32_C(0x3a001809),
            "1258637401934962560", "32851598462703820800");

        std::cout << "FULL_PREFIXES JOINT_CELLS "
                  << full256.joint_cells + full384.joint_cells
                  << " REPAIRS " << full256.prefix_size + full384.prefix_size
                  << " BODY_CHECKS "
                  << full256.body_checks + full384.body_checks << "\n";
        std::cout << "VERDICT PASS LITERAL_ATLAS_SEEDS_PREFIXES_BOUNDARY\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ROUND5_COMPACT_LITERAL_ERROR " << error.what() << '\n';
        return 1;
    }
}

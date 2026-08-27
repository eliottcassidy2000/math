// Literal-wall controls for the carrier resistance witnesses.
//
// This program reconstructs the original 4,675-mask carrier and the 58-mask
// CEGAR enrichment, but computes all activation margins from a fresh joint
// wall arrangement.  It does not call the endpoint-cocycle or atom routines.

#define main thm4254_cleanroom_unused_main
#include "lrc14_endpoint_cascade_direct_wall_body_audit.cpp"
#undef main

#include <set>

namespace {

struct DirectScanControl {
    u64 active = 0;
    u64 bodies = 0;
    u64 failures = 0;
    u64 checks = 0;
    u64 max_checks = 0;
    u32 first_failure = std::numeric_limits<u32>::max();
    u32 best_repair = 0;
    i128 best_margin = 0;
    i64 grid = 0;
};

DirectScanControl direct_scan_control(int q, int r,
                                      const std::vector<u32>& carrier,
                                      const std::vector<u32>& bodies) {
    const Geometry geometry = build_joint_geometry(q, r);
    std::vector<u32> active;
    std::vector<i128> margins;
    active.reserve(carrier.size());
    margins.reserve(carrier.size());
    for (u32 repair : carrier) {
        i64 mass = 0;
        for (const Cell& cell : geometry.cells) {
            if (cell.pair_safe && (cell.failed_pool & ~repair) == 0) {
                mass += cell.width;
            }
        }
        const i128 margin = static_cast<i128>(63) * mass -
                            static_cast<i128>(4) * geometry.grid;
        margins.push_back(margin);
        if (margin >= 0) active.push_back(repair);
    }

    const unsigned hardware = std::thread::hardware_concurrency();
    const unsigned thread_count =
        std::max(1u, std::min(8u, hardware == 0 ? 1u : hardware));
    std::vector<DirectScanControl> locals(thread_count);
    std::vector<std::thread> threads;
    for (unsigned thread = 0; thread < thread_count; ++thread) {
        threads.emplace_back([&, thread]() {
            DirectScanControl local;
            const std::size_t begin = bodies.size() * thread / thread_count;
            const std::size_t end = bodies.size() * (thread + 1) / thread_count;
            for (std::size_t index = begin; index < end; ++index) {
                const u32 body = bodies[index];
                u64 used = 0;
                bool covered = false;
                for (u32 repair : active) {
                    ++used;
                    if ((body & repair) == 0) {
                        covered = true;
                        break;
                    }
                }
                ++local.bodies;
                local.checks += used;
                if (!covered) {
                    ++local.failures;
                    local.first_failure =
                        std::min(local.first_failure, body);
                }
                local.max_checks = std::max(local.max_checks, used);
            }
            locals[thread] = local;
        });
    }
    for (auto& thread : threads) thread.join();

    DirectScanControl answer;
    answer.active = active.size();
    answer.grid = geometry.grid;
    for (const DirectScanControl& local : locals) {
        answer.bodies += local.bodies;
        answer.failures += local.failures;
        answer.checks += local.checks;
        answer.max_checks = std::max(answer.max_checks, local.max_checks);
        answer.first_failure =
            std::min(answer.first_failure, local.first_failure);
    }
    require(answer.failures != 0, "control row unexpectedly closes");

    bool first = true;
    for (std::size_t index = 0; index < carrier.size(); ++index) {
        if ((carrier[index] & answer.first_failure) != 0) continue;
        if (first || margins[index] > answer.best_margin ||
            (margins[index] == answer.best_margin &&
             carrier[index] < answer.best_repair)) {
            first = false;
            answer.best_margin = margins[index];
            answer.best_repair = carrier[index];
        }
    }
    require(!first && answer.best_margin < 0,
            "direct resistance witness is inconsistent");
    return answer;
}

void check_and_print_control(int q, int r,
                             const std::vector<u32>& carrier,
                             const std::vector<u32>& bodies,
                             u64 expected_active,
                             u64 expected_failures,
                             u32 expected_first,
                             u32 expected_best,
                             const std::string& atom_margin,
                             const std::string& atom_denominator) {
    const DirectScanControl scan =
        direct_scan_control(q, r, carrier, bodies);
    require(scan.active == expected_active &&
                scan.bodies == EXPECTED_BODIES &&
                scan.failures == expected_failures &&
                scan.first_failure == expected_first &&
                scan.best_repair == expected_best,
            "literal control disagrees with endpoint-atom scan");
    const i128 atom_num = -parse_i128(atom_margin);
    const i128 atom_den = parse_i128(atom_denominator);
    require(scan.best_margin * atom_den ==
                atom_num * static_cast<i128>(scan.grid),
            "literal best margin disagrees with endpoint-atom ratio");
    std::cout << "PAIR " << q << ',' << r
              << " CARRIER " << carrier.size()
              << " ACTIVE " << scan.active
              << " BODIES " << scan.bodies
              << " FAILURES " << scan.failures
              << " CHECKS " << scan.checks
              << " MAX_CHECKS " << scan.max_checks
              << " FIRST_FAILURE " << std::hex << scan.first_failure
              << std::dec << " LABELS {" << labels(scan.first_failure) << "}"
              << " BEST_REPAIR " << std::hex << scan.best_repair << std::dec
              << " LABELS {" << labels(scan.best_repair) << "}"
              << " LITERAL_BEST_MARGIN " << decimal(scan.best_margin)
              << '/' << scan.grid
              << " REDUCED " << fraction(scan.best_margin, scan.grid)
              << " CROSS_ATOM PASS\n";
}

void print_novel_transfer_control(int q, int r,
                                  const std::vector<u32>& novel,
                                  u32 hostile_body) {
    const Geometry geometry = build_joint_geometry(q, r);
    u64 active = 0;
    u64 disjoint = 0;
    u64 active_disjoint = 0;
    bool first = true;
    i128 best_margin = 0;
    u32 best_repair = 0;
    for (u32 repair : novel) {
        i64 mass = 0;
        for (const Cell& cell : geometry.cells) {
            if (cell.pair_safe && (cell.failed_pool & ~repair) == 0) {
                mass += cell.width;
            }
        }
        const i128 margin = static_cast<i128>(63) * mass -
                            static_cast<i128>(4) * geometry.grid;
        if (margin >= 0) ++active;
        if ((repair & hostile_body) == 0) {
            ++disjoint;
            if (margin >= 0) ++active_disjoint;
            if (first || margin > best_margin ||
                (margin == best_margin && repair < best_repair)) {
                first = false;
                best_margin = margin;
                best_repair = repair;
            }
        }
    }
    require(!first, "no novel repair is disjoint from hostile body");
    std::cout << "NOVEL_TRANSFER PAIR " << q << ',' << r
              << " NOVEL " << novel.size()
              << " ACTIVE " << active
              << " DISJOINT_FROM_FIRST " << disjoint
              << " ACTIVE_DISJOINT " << active_disjoint
              << " BEST_NOVEL_DISJOINT_REPAIR " << std::hex << best_repair
              << std::dec << " LABELS {" << labels(best_repair) << "}"
              << " BEST_NOVEL_DISJOINT_MARGIN " << decimal(best_margin)
              << '/' << geometry.grid
              << " REDUCED " << fraction(best_margin, geometry.grid) << "\n";
}

}  // namespace

#ifndef CARRIER_CLEANROOM_CONTROL_LIBRARY_ONLY
int main(int argc, char** argv) {
    try {
        require(argc == 3,
                "usage: resistance-controls REPLAY_BAND FULL_416_704");
        const std::filesystem::path band(argv[1]);
        std::vector<u32> base;
        std::set<u32> seen;
        Fnv base_ledger;
        for (const auto& [q, r] : EXPECTED_PAIRS) {
            const Parsed parsed = parse_probe(
                band / (std::to_string(q) + "_" +
                        std::to_string(r) + ".out"));
            require(parsed.q == q && parsed.r == r,
                    "band transcript mismatch");
            for (u32 repair : parsed.deck) {
                if (seen.insert(repair).second) {
                    base.push_back(repair);
                    base_ledger.add(repair);
                }
            }
        }
        require(base.size() == 4675 &&
                    base_ledger.state == UINT64_C(0xce4e76ec11df057c),
                "literal-control base carrier changed");

        const Parsed seed = parse_probe(argv[2]);
        require(seed.q == 416 && seed.r == 704 && seed.deck.size() == 2608,
                "literal-control seed prefix changed");
        std::vector<u32> enriched = base;
        std::vector<u32> novel_masks;
        Fnv novel_ledger;
        u64 novel = 0;
        for (u32 repair : seed.deck) {
            if (seen.insert(repair).second) {
                enriched.push_back(repair);
                novel_masks.push_back(repair);
                novel_ledger.add(repair);
                ++novel;
            }
        }
        require(novel == 58 && enriched.size() == 4733 &&
                    novel_ledger.state == UINT64_C(0xf25cf1f803d09f8a),
                "literal-control enrichment changed");

        std::vector<u32> bodies;
        bodies.reserve(EXPECTED_BODIES);
        enumerate_bodies(0, 9, 0, bodies);
        require(bodies.size() == EXPECTED_BODIES,
                "literal-control body universe changed");

        std::cout << "CARRIER_CEGAR_LITERAL_RESISTANCE_CONTROLS_V1\n";
        std::cout << "BASE " << base.size() << " FNV " << std::hex
                  << base_ledger.state << " NOVEL " << std::dec << novel
                  << " NOVEL_FNV " << std::hex << novel_ledger.state
                  << std::dec << " ENRICHED " << enriched.size() << "\n";
        check_and_print_control(
            416, 704, base, bodies, 3330, 1, UINT32_C(0x2f903000),
            UINT32_C(0x10458a01), "20520733194561600",
            "2337203273714749440");
        check_and_print_control(
            416, 700, enriched, bodies, 2924, 38, UINT32_C(0x051a7400),
            UINT32_C(0x10458a01), "150665738470752960",
            "18591389677276416000");
        check_and_print_control(
            520, 700, enriched, bodies, 3363, 1, UINT32_C(0x2d183001),
            UINT32_C(0x12244604), "62470194522372000",
            "4647847419319104000");
        print_novel_transfer_control(
            416, 704, novel_masks, UINT32_C(0x2f903000));
        print_novel_transfer_control(
            416, 700, novel_masks, UINT32_C(0x051a7400));
        print_novel_transfer_control(
            520, 700, novel_masks, UINT32_C(0x2d183001));
        std::cout << "VERDICT PASS THREE_LITERAL_WALL_CONTROLS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "CARRIER_CEGAR_LITERAL_CONTROL_ERROR "
                  << error.what() << '\n';
        return 1;
    }
}
#endif

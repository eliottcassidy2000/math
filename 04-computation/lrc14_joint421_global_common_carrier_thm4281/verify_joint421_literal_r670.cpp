// Detached literal-joint-wall audit of the locked 421-mask joint deck.
// This translation unit intentionally imports the canonical direct-wall
// checker, not the endpoint-cocycle/atom implementation.

#define CARRIER_CLEANROOM_CONTROL_LIBRARY_ONLY
#include "04-computation/lrc14_three_round_learned_carrier_thm4266/cleanroom_resistance_controls.cpp"
#undef CARRIER_CLEANROOM_CONTROL_LIBRARY_ONLY

#include <fstream>

namespace {

constexpr std::size_t EXPECTED_JOINT = 421;
constexpr u64 EXPECTED_JOINT_FNV = UINT64_C(0x20d63dd42fe8150e);

constexpr std::array<u32, 6> THM4276_AUGMENTATION = {{
    0x00289285, 0x0260812c, 0x18689040,
    0x20c0c124, 0x302c1006, 0x30580888}};

void append_unseen_literal421(const std::vector<u32>& source,
                              std::set<u32>& seen,
                              std::vector<u32>& target) {
    for (u32 repair : source) {
        if (seen.insert(repair).second) target.push_back(repair);
    }
}

u64 fnv_literal421(const std::vector<u32>& masks) {
    Fnv ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

std::vector<u32> active_literal421(int q, int r,
                                   const std::vector<u32>& carrier) {
    const Geometry geometry = build_joint_geometry(q, r);
    std::vector<u32> active;
    active.reserve(carrier.size());
    for (u32 repair : carrier) {
        i64 mass = 0;
        for (const Cell& cell : geometry.cells) {
            if (cell.pair_safe && (cell.failed_pool & ~repair) == 0)
                mass += cell.width;
        }
        if (static_cast<i128>(63) * mass >=
            static_cast<i128>(4) * geometry.grid) {
            active.push_back(repair);
        }
    }
    return active;
}

std::vector<u32> failures_literal421(const std::vector<u32>& bodies,
                                     const std::vector<u32>& active) {
    const unsigned hardware = std::thread::hardware_concurrency();
    const unsigned threads =
        std::max(1u, std::min(8u, hardware == 0 ? 1u : hardware));
    std::vector<std::vector<u32>> local(threads);
    std::vector<std::thread> workers;
    for (unsigned lane = 0; lane < threads; ++lane) {
        workers.emplace_back([&, lane]() {
            const std::size_t begin = bodies.size() * lane / threads;
            const std::size_t end = bodies.size() * (lane + 1) / threads;
            for (std::size_t index = begin; index < end; ++index) {
                bool covered = false;
                for (u32 repair : active) {
                    if ((repair & bodies[index]) == 0) {
                        covered = true;
                        break;
                    }
                }
                if (!covered) local[lane].push_back(bodies[index]);
            }
        });
    }
    for (auto& worker : workers) worker.join();
    std::vector<u32> failures;
    for (auto& part : local)
        failures.insert(failures.end(), part.begin(), part.end());
    std::sort(failures.begin(), failures.end());
    return failures;
}

std::vector<u32> read_joint421(const std::filesystem::path& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open joint deck");
    std::vector<u32> masks;
    std::set<u32> distinct;
    std::string word;
    while (in >> word) {
        const u64 wide = std::stoull(word, nullptr, 16);
        require(wide < (UINT64_C(1) << 30), "joint mask leaves pool");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8, "joint mask is not rank eight");
        require(distinct.insert(mask).second, "joint deck repeats a mask");
        masks.push_back(mask);
    }
    require(masks.size() == EXPECTED_JOINT &&
                fnv_literal421(masks) == EXPECTED_JOINT_FNV,
            "joint-deck ledger changed");
    return masks;
}

u64 body_fnv_literal(const std::vector<u32>& bodies) {
    Fnv ledger;
    for (u32 body : bodies) ledger.add(body);
    return ledger.state;
}

void enumerate_allowed_rank8(const std::vector<unsigned>& allowed,
                             std::size_t begin, unsigned need, u32 mask,
                             std::vector<u32>& output) {
    if (need == 0) {
        output.push_back(mask);
        return;
    }
    if (allowed.size() - begin < need) return;
    for (std::size_t index = begin;
         index + need <= allowed.size(); ++index) {
        enumerate_allowed_rank8(allowed, index + 1, need - 1,
                                mask | (u32{1} << allowed[index]), output);
    }
}

std::vector<u32> repairs_disjoint_from(u32 forbidden) {
    std::vector<unsigned> allowed;
    for (unsigned bit = 0; bit < 30; ++bit)
        if ((forbidden & (u32{1} << bit)) == 0) allowed.push_back(bit);
    std::vector<u32> repairs;
    enumerate_allowed_rank8(allowed, 0, 8, 0, repairs);
    return repairs;
}

struct ResponseAudit {
    u64 covered = 0;
    u64 failures = 0;
    u64 incidences = 0;
    Fnv ledger;
};

void add_response_row(ResponseAudit& audit, u64 row,
                      const std::vector<u32>& obligations,
                      const std::vector<u32>& active) {
    for (std::size_t index = 0; index < obligations.size(); ++index) {
        u64 hits = 0;
        for (u32 repair : active) hits += (repair & obligations[index]) == 0;
        audit.covered += hits != 0;
        audit.failures += hits == 0;
        audit.incidences += hits;
        audit.ledger.add(row);
        audit.ledger.add(index);
        audit.ledger.add(obligations[index]);
        audit.ledger.add(hits);
    }
}

}  // namespace

#ifndef JOINT421_LITERAL_LIBRARY_ONLY
int main(int argc, char** argv) {
    try {
        require(argc == 8,
                "usage: joint421-literal REPLAY_BAND PREFIX_416_704 PREFIX_416_700 PREFIX_520_700 PREFIX_384_694 PREFIX_520_688 JOINT_DECK");
        const std::filesystem::path band(argv[1]);
        const Parsed prefix704 = parse_probe(argv[2]);
        const Parsed prefix416 = parse_probe(argv[3]);
        const Parsed prefix520 = parse_probe(argv[4]);
        const Parsed prefix384 = parse_probe(argv[5]);
        const Parsed prefix688 = parse_probe(argv[6]);
        require(prefix704.q == 416 && prefix704.r == 704 &&
                    prefix704.deck.size() == 2608 &&
                    prefix704.declared_prefix_fnv == UINT64_C(0x18ff663a123e684e) &&
                    prefix416.q == 416 && prefix416.r == 700 &&
                    prefix416.deck.size() == 5894 &&
                    prefix416.declared_prefix_fnv == UINT64_C(0x701c0233c8f8abeb) &&
                    prefix520.q == 520 && prefix520.r == 700 &&
                    prefix520.deck.size() == 4557 &&
                    prefix520.declared_prefix_fnv == UINT64_C(0xd8466558bcd8ef9d) &&
                    prefix384.q == 384 && prefix384.r == 694 &&
                    prefix384.deck.size() == 5307 &&
                    prefix384.declared_prefix_fnv == UINT64_C(0x7b9a46c60e8514f8) &&
                    prefix688.q == 520 && prefix688.r == 688 &&
                    prefix688.deck.size() == 5398 &&
                    prefix688.declared_prefix_fnv == UINT64_C(0x6ab471da88c8e1d1),
                "source-prefix identity changed");

        std::vector<u32> carrier;
        std::set<u32> seen;
        for (const auto& [q, r] : EXPECTED_PAIRS) {
            const Parsed parsed = parse_probe(
                band / (std::to_string(q) + "_" + std::to_string(r) + ".out"));
            require(parsed.q == q && parsed.r == r,
                    "literal band transcript mismatch");
            append_unseen_literal421(parsed.deck, seen, carrier);
        }
        require(carrier.size() == 4675 &&
                    fnv_literal421(carrier) == UINT64_C(0xce4e76ec11df057c),
                "base carrier changed");
        append_unseen_literal421(prefix704.deck, seen, carrier);
        require(carrier.size() == 4733 &&
                    fnv_literal421(carrier) == UINT64_C(0xa7b046289655c733),
                "round-one carrier changed");
        const std::set<u32> round1_seen = seen;
        std::vector<u32> novel416;
        std::vector<u32> novel520;
        for (u32 repair : prefix416.deck)
            if (!round1_seen.count(repair)) novel416.push_back(repair);
        for (u32 repair : prefix520.deck)
            if (!round1_seen.count(repair)) novel520.push_back(repair);
        append_unseen_literal421(novel416, seen, carrier);
        append_unseen_literal421(novel520, seen, carrier);
        require(carrier.size() == 7986 &&
                    fnv_literal421(carrier) == UINT64_C(0xbaef1d2f49444638),
                "round-two carrier changed");
        append_unseen_literal421(prefix384.deck, seen, carrier);
        require(carrier.size() == 8319 &&
                    fnv_literal421(carrier) == UINT64_C(0xe08b227730f6793c),
                "round-three carrier changed");
        append_unseen_literal421(prefix688.deck, seen, carrier);
        require(carrier.size() == 8518 &&
                    fnv_literal421(carrier) == UINT64_C(0x1603e3fe970f8428),
                "THM-4271 carrier changed");
        append_unseen_literal421(
            std::vector<u32>(THM4276_AUGMENTATION.begin(), THM4276_AUGMENTATION.end()),
            seen, carrier);
        require(carrier.size() == 8524 &&
                    fnv_literal421(carrier) == UINT64_C(0x5ddb84a44f5d2ad7),
                "THM-4276 carrier changed");

        const std::vector<u32> joint = read_joint421(argv[7]);
        for (u32 mask : joint)
            require(!seen.count(mask), "joint mask already in THM-4276 carrier");
        std::vector<u32> augmented = carrier;
        append_unseen_literal421(joint, seen, augmented);
        require(augmented.size() == 8945, "augmented carrier count changed");

        std::vector<u32> bodies;
        bodies.reserve(EXPECTED_BODIES);
        enumerate_bodies(0, 9, 0, bodies);
        require(bodies.size() == EXPECTED_BODIES,
                "literal body universe changed");

        const std::vector<u32> active_base_a = active_literal421(256, 670, carrier);
        const std::vector<u32> active_base_b = active_literal421(384, 670, carrier);
        require(active_base_a.size() == 2172 && active_base_b.size() == 3851,
                "base r670 active census changed");
        const std::vector<u32> failures_a = failures_literal421(bodies, active_base_a);
        const std::vector<u32> failures_b = failures_literal421(bodies, active_base_b);
        require(failures_a.size() == 1659 && failures_b.size() == 46 &&
                    body_fnv_literal(failures_a) == UINT64_C(0x970f004b0f9e2edb) &&
                    body_fnv_literal(failures_b) == UINT64_C(0x6c4fd2dab94cf6a9),
                "base r670 obligation ledger changed");

        const std::vector<u32> active_joint_a = active_literal421(256, 670, joint);
        const std::vector<u32> active_joint_b = active_literal421(384, 670, joint);
        ResponseAudit response;
        add_response_row(response, 0, failures_a, active_joint_a);
        add_response_row(response, 1, failures_b, active_joint_b);
        require(active_joint_a.size() == 397 &&
                    active_joint_b.size() == 367 &&
                    response.covered == 1705 && response.failures == 0 &&
                    response.incidences == 7030 &&
                    response.ledger.state == UINT64_C(0x4ab5ca017f17efeb),
                "joint deck misses an r670 obligation");

        const std::vector<u32> active_augmented_a = active_literal421(256, 670, augmented);
        const std::vector<u32> active_augmented_b = active_literal421(384, 670, augmented);
        const BodyAudit closed_a = audit_bodies(bodies, active_augmented_a);
        const BodyAudit closed_b = audit_bodies(bodies, active_augmented_b);
        require(closed_a.bodies == EXPECTED_BODIES && closed_a.failures == 0 &&
                    closed_b.bodies == EXPECTED_BODIES && closed_b.failures == 0,
                "augmented literal r670 carrier does not close");
        const BodyAudit joint_body = audit_bodies(bodies, joint);
        require(joint_body.bodies == EXPECTED_BODIES && joint_body.failures == 0,
                "joint deck misses a nine-body combinatorially");
        require(active_augmented_a.size() == 2569 &&
                    active_augmented_b.size() == 4218 &&
                    closed_a.checks == UINT64_C(613731713) &&
                    closed_a.max_checks == 2497 &&
                    closed_b.checks == UINT64_C(541846700) &&
                    closed_b.max_checks == 4205 &&
                    joint_body.checks == UINT64_C(405170384) &&
                    joint_body.max_checks == 421,
                "frozen r670/body scan ledger changed");

        // First endpoint-669 hostile row and the ordered two-mask repair
        // supplied by the continuation search.  Re-evaluate both masks and
        // the final closure solely on the fresh literal joint grid.
        const std::vector<u32> active_augmented_669 =
            active_literal421(256, 669, augmented);
        const std::vector<u32> failures_669 =
            failures_literal421(bodies, active_augmented_669);
        const std::vector<u32> expected_failures_669 = {
            UINT32_C(0x0f047001), UINT32_C(0x1f302001),
            UINT32_C(0x1f803001), UINT32_C(0x37003088)};
        require(active_augmented_669.size() == 3359 &&
                    failures_669 == expected_failures_669 &&
                    body_fnv_literal(failures_669) ==
                        UINT64_C(0x6a4e5c7582d44240),
                "endpoint-669 literal failure count changed");
        const std::vector<u32> witnesses_669 = {
            UINT32_C(0x003884c8), UINT32_C(0x00c4c124)};
        const std::array<u64, 2> expected_patterns_669 = {
            UINT64_C(0x5), UINT64_C(0xa)};
        for (u32 witness : witnesses_669) {
            require(std::popcount(witness) == 8 && !seen.count(witness),
                    "endpoint-669 witness invalid or already present");
        }
        u64 union_669 = 0;
        Fnv witness_pattern_ledger;
        for (std::size_t index = 0; index < witnesses_669.size(); ++index) {
            const u32 witness = witnesses_669[index];
            const std::vector<u32> active_single =
                active_literal421(256, 669, std::vector<u32>{witness});
            require(active_single == std::vector<u32>{witness},
                    "endpoint-669 witness is literally inactive");
            u64 pattern = 0;
            for (std::size_t bit = 0; bit < failures_669.size(); ++bit)
                if ((witness & failures_669[bit]) == 0)
                    pattern |= u64{1} << bit;
            require(pattern == expected_patterns_669[index],
                    "endpoint-669 witness response changed");
            union_669 |= pattern;
            witness_pattern_ledger.add(witness);
            witness_pattern_ledger.add(pattern);
        }
        require(union_669 == UINT64_C(0xf),
                "ordered endpoint-669 witnesses miss a failure");
        std::vector<u32> augmented_669 = augmented;
        std::set<u32> seen_669 = seen;
        append_unseen_literal421(witnesses_669, seen_669, augmented_669);
        const std::vector<u32> active_closed_669 =
            active_literal421(256, 669, augmented_669);
        const BodyAudit closed_669 = audit_bodies(bodies, active_closed_669);
        require(closed_669.bodies == EXPECTED_BODIES &&
                    closed_669.failures == 0,
                "two-mask endpoint-669 literal augmentation does not close");
        require(witness_pattern_ledger.state ==
                    UINT64_C(0xa4359ec88357f64d) &&
                    active_closed_669.size() == 3361 &&
                    closed_669.checks == UINT64_C(569666955) &&
                    closed_669.max_checks == 3361,
                "frozen endpoint-669 literal ledger changed");

        // Descend two further endpoint units.  The inherited carrier here is
        // exactly the 8,945-mask THM-4276+joint deck plus the ordered r669
        // minimum-two augmentation just verified above.
        const std::vector<u32> active_augmented_667 =
            active_literal421(256, 667, augmented_669);
        const std::vector<u32> failures_667 =
            failures_literal421(bodies, active_augmented_667);
        require(failures_667.size() == 18 &&
                    body_fnv_literal(failures_667) ==
                        UINT64_C(0xce3d9a9b4b429af8),
                "endpoint-667 literal failure ledger changed");
        const std::vector<u32> witnesses_667 = {
            UINT32_C(0x02a05206), UINT32_C(0x10e05240)};
        const std::array<u64, 2> expected_patterns_667 = {
            UINT64_C(0x11c07), UINT64_C(0x2efff)};
        u64 union_667 = 0;
        Fnv witness_667_pattern_ledger;
        for (std::size_t index = 0; index < witnesses_667.size(); ++index) {
            const u32 witness = witnesses_667[index];
            require(std::popcount(witness) == 8 && !seen_669.count(witness),
                    "endpoint-667 witness invalid or already present");
            const std::vector<u32> active_single =
                active_literal421(256, 667, std::vector<u32>{witness});
            require(active_single == std::vector<u32>{witness},
                    "endpoint-667 witness is literally inactive");
            u64 pattern = 0;
            for (std::size_t bit = 0; bit < failures_667.size(); ++bit)
                if ((witness & failures_667[bit]) == 0)
                    pattern |= u64{1} << bit;
            require(pattern == expected_patterns_667[index],
                    "endpoint-667 witness response changed");
            union_667 |= pattern;
            witness_667_pattern_ledger.add(witness);
            witness_667_pattern_ledger.add(pattern);
        }
        require(union_667 == UINT64_C(0x3ffff),
                "ordered endpoint-667 witnesses miss a failure");

        // Detached lower bound for the exact minimum-two assertion.  Bodies
        // 13 and 16 form the frozen two-obligation packing certificate.  A
        // one-mask augmentation would have to be disjoint from their union;
        // enumerate that entire C(17,8) fibre and test literal activation.
        const u32 packing_union_667 =
            failures_667[13] | failures_667[16];
        require(failures_667[13] == UINT32_C(0x0f142401) &&
                    failures_667[16] == UINT32_C(0x151aa400) &&
                    packing_union_667 == UINT32_C(0x1f1ea401) &&
                    std::popcount(packing_union_667) == 13,
                "endpoint-667 packing identity changed");
        const std::vector<u32> packing_candidates_667 =
            repairs_disjoint_from(packing_union_667);
        require(packing_candidates_667.size() == 24310,
                "endpoint-667 packing fibre changed");
        require(fnv_literal421(packing_candidates_667) ==
                    UINT64_C(0x0962a60e00e97c69),
                "endpoint-667 packing candidate ledger changed");
        const std::vector<u32> packing_active_667 =
            active_literal421(256, 667, packing_candidates_667);
        require(packing_active_667.empty(),
                "endpoint-667 packing fibre admits a one-mask repair");
        const u64 packing_candidate_fnv_667 =
            fnv_literal421(packing_candidates_667);
        std::vector<u32> augmented_667 = augmented_669;
        std::set<u32> seen_667 = seen_669;
        append_unseen_literal421(witnesses_667, seen_667, augmented_667);
        const std::vector<u32> active_closed_667 =
            active_literal421(256, 667, augmented_667);
        const BodyAudit closed_667 = audit_bodies(bodies, active_closed_667);
        require(closed_667.bodies == EXPECTED_BODIES &&
                    closed_667.failures == 0,
                "two-mask endpoint-667 literal augmentation does not close");
        require(active_augmented_667.size() == 3982 &&
                    witness_667_pattern_ledger.state ==
                        UINT64_C(0x1d015ec917d3849b) &&
                    active_closed_667.size() == 3984 &&
                    closed_667.checks == UINT64_C(572210490) &&
                    closed_667.max_checks == 3984 &&
                    augmented_667.size() == 8949 &&
                    fnv_literal421(augmented_667) ==
                        UINT64_C(0x1496ebecca72684b),
                "frozen endpoint-667 literal ledger changed");

        // Final positive descent step.  The single failed body at (520,665)
        // is evaluated before the witness is appended; the appended carrier
        // is then rescanned over the complete nine-body universe.
        const std::vector<u32> active_augmented_665 =
            active_literal421(520, 665, augmented_667);
        const std::vector<u32> failures_665 =
            failures_literal421(bodies, active_augmented_665);
        const std::vector<u32> expected_failures_665 = {
            UINT32_C(0x151a3400)};
        require(active_augmented_665.size() == 5156 &&
                    failures_665 == expected_failures_665 &&
                    body_fnv_literal(failures_665) ==
                        UINT64_C(0x183de0540d550fac),
                "endpoint-665 literal hostile ledger changed");
        const u32 witness_665 = UINT32_C(0x0004409f);
        require(std::popcount(witness_665) == 8 &&
                    !seen_667.count(witness_665) &&
                    (witness_665 & failures_665.front()) == 0,
                "endpoint-665 witness identity/disjointness changed");
        const std::vector<u32> active_single_665 =
            active_literal421(520, 665, std::vector<u32>{witness_665});
        require(active_single_665 == std::vector<u32>{witness_665},
                "endpoint-665 witness is literally inactive");
        std::vector<u32> augmented_665 = augmented_667;
        std::set<u32> seen_665 = seen_667;
        append_unseen_literal421(std::vector<u32>{witness_665},
                                 seen_665, augmented_665);
        require(augmented_665.size() == 8950 &&
                    fnv_literal421(augmented_665) ==
                        UINT64_C(0xf07022300e266930),
                "endpoint-665 augmented carrier ledger changed");
        const std::vector<u32> active_closed_665 =
            active_literal421(520, 665, augmented_665);
        const BodyAudit closed_665 = audit_bodies(bodies, active_closed_665);
        require(active_closed_665.size() == 5157 &&
                    closed_665.bodies == EXPECTED_BODIES &&
                    closed_665.failures == 0 &&
                    closed_665.checks == UINT64_C(531031657) &&
                    closed_665.max_checks == 5157,
                "endpoint-665 literal closure ledger changed");

        // Final one-mask continuation.  First retain the exact two-body
        // hostile boundary, then check the least response-atlas witness on
        // the detached literal grid and close the full body universe.
        const std::vector<u32> active_hostile_664 =
            active_literal421(256, 664, augmented_665);
        const std::vector<u32> failures_664 =
            failures_literal421(bodies, active_hostile_664);
        const std::vector<u32> expected_failures_664 = {
            UINT32_C(0x0d143408), UINT32_C(0x2d143400)};
        require(active_hostile_664.size() == 4527 &&
                    failures_664 == expected_failures_664 &&
                    body_fnv_literal(failures_664) ==
                        UINT64_C(0xe9872402af5e2795),
                "endpoint-664 final hostile ledger changed");
        const u32 witness_664 = UINT32_C(0x00c0c125);
        require(std::popcount(witness_664) == 8 &&
                    !seen_665.count(witness_664),
                "endpoint-664 witness identity changed");
        const std::vector<u32> active_single_664 =
            active_literal421(256, 664, std::vector<u32>{witness_664});
        require(active_single_664 == std::vector<u32>{witness_664} &&
                    (witness_664 & failures_664[0]) == 0 &&
                    (witness_664 & failures_664[1]) == 0,
                "endpoint-664 witness is inactive or misses an obligation");
        std::vector<u32> augmented_664 = augmented_665;
        std::set<u32> seen_664 = seen_665;
        append_unseen_literal421(std::vector<u32>{witness_664},
                                 seen_664, augmented_664);
        require(augmented_664.size() == 8951 &&
                    fnv_literal421(augmented_664) ==
                        UINT64_C(0x188f82ab9dd1695a),
                "endpoint-664 augmented carrier ledger changed");
        const std::vector<u32> active_closed_664 =
            active_literal421(256, 664, augmented_664);
        const BodyAudit closed_664 = audit_bodies(bodies, active_closed_664);
        require(active_closed_664.size() == 4528 &&
                    closed_664.bodies == EXPECTED_BODIES &&
                    closed_664.failures == 0 &&
                    closed_664.checks == UINT64_C(536018179) &&
                    closed_664.max_checks == 4528,
                "endpoint-664 literal closure ledger changed");

        Fnv active_a_ledger;
        for (u32 mask : active_joint_a) active_a_ledger.add(mask);
        Fnv active_b_ledger;
        for (u32 mask : active_joint_b) active_b_ledger.add(mask);
        std::cout << "JOINT421_LITERAL_DESCENT_V4\n";
        std::cout << "THM4276_CARRIER " << carrier.size() << " FNV " << std::hex
                  << fnv_literal421(carrier) << std::dec << '\n';
        std::cout << "JOINT_DECK " << joint.size() << " FNV " << std::hex
                  << fnv_literal421(joint) << std::dec << " OVERLAP 0\n";
        std::cout << "AUGMENTED_CARRIER " << augmented.size() << " FNV "
                  << std::hex << fnv_literal421(augmented) << std::dec << '\n';
        std::cout << "BASE_OBLIGATIONS A " << failures_a.size() << " FNV "
                  << std::hex << body_fnv_literal(failures_a) << std::dec
                  << " B " << failures_b.size() << " FNV " << std::hex
                  << body_fnv_literal(failures_b)
                  << std::dec << '\n';
        std::cout << "JOINT_ACTIVITY A " << active_joint_a.size() << " FNV "
                  << std::hex << active_a_ledger.state << std::dec << " B "
                  << active_joint_b.size() << " FNV " << std::hex
                  << active_b_ledger.state
                  << std::dec << '\n';
        std::cout << "R670_RESPONSE OBLIGATIONS " << response.covered
                  << " FAILURES " << response.failures << " INCIDENCES "
                  << response.incidences << " LEDGER " << std::hex
                  << response.ledger.state << std::dec << '\n';
        std::cout << "PAIR 256,670 AUGMENTED_ACTIVE " << active_augmented_a.size()
                  << " BODIES " << closed_a.bodies << " FAILURES "
                  << closed_a.failures << " CHECKS " << closed_a.checks
                  << " MAX_CHECKS " << closed_a.max_checks << '\n';
        std::cout << "PAIR 384,670 AUGMENTED_ACTIVE " << active_augmented_b.size()
                  << " BODIES " << closed_b.bodies << " FAILURES "
                  << closed_b.failures << " CHECKS " << closed_b.checks
                  << " MAX_CHECKS " << closed_b.max_checks << '\n';
        std::cout << "JOINT_BODY_SCAN BODIES " << joint_body.bodies
                  << " FAILURES " << joint_body.failures << " CHECKS "
                  << joint_body.checks << " MAX_CHECKS "
                  << joint_body.max_checks << '\n';
        std::cout << "ENDPOINT669_HOSTILE PAIR 256,669 ACTIVE "
                  << active_augmented_669.size() << " FAILURES "
                  << failures_669.size() << " FAILURE_FNV " << std::hex
                  << body_fnv_literal(failures_669) << " MASKS";
        for (u32 body : failures_669) std::cout << ' ' << body;
        std::cout << std::dec << '\n';
        std::cout << "ENDPOINT669_ORDERED_MIN2 " << std::hex
                  << witnesses_669[0] << ',' << witnesses_669[1]
                  << " RESPONSE_UNION " << union_669
                  << " PATTERN_LEDGER " << witness_pattern_ledger.state
                  << std::dec << " FINAL_ACTIVE " << active_closed_669.size()
                  << " BODIES " << closed_669.bodies << " FAILURES "
                  << closed_669.failures << " CHECKS " << closed_669.checks
                  << " MAX_CHECKS " << closed_669.max_checks << '\n';
        std::cout << "ENDPOINT667_HOSTILE PAIR 256,667 ACTIVE "
                  << active_augmented_667.size() << " FAILURES "
                  << failures_667.size() << " FAILURE_FNV " << std::hex
                  << body_fnv_literal(failures_667) << " MASKS";
        for (u32 body : failures_667) std::cout << ' ' << body;
        std::cout << std::dec << '\n';
        std::cout << "ENDPOINT667_ORDERED_MIN2 " << std::hex
                  << witnesses_667[0] << ',' << witnesses_667[1]
                  << " PATTERNS " << expected_patterns_667[0] << ','
                  << expected_patterns_667[1] << " RESPONSE_UNION "
                  << union_667 << " PATTERN_LEDGER "
                  << witness_667_pattern_ledger.state << std::dec
                  << " FINAL_ACTIVE " << active_closed_667.size()
                  << " BODIES " << closed_667.bodies << " FAILURES "
                  << closed_667.failures << " CHECKS " << closed_667.checks
                  << " MAX_CHECKS " << closed_667.max_checks << '\n';
        std::cout << "ENDPOINT667_LITERAL_MIN2_LOWER_BOUND PACKING 13,16"
                  << " UNION " << std::hex << packing_union_667 << std::dec
                  << " CANDIDATES " << packing_candidates_667.size()
                  << " CANDIDATE_FNV " << std::hex
                  << packing_candidate_fnv_667 << std::dec
                  << " LITERAL_ACTIVE " << packing_active_667.size() << '\n';
        std::cout << "ENDPOINT667_FINAL_CARRIER " << augmented_667.size()
                  << " FNV " << std::hex << fnv_literal421(augmented_667)
                  << std::dec << '\n';
        std::cout << "ENDPOINT665_HOSTILE PAIR 520,665 ACTIVE "
                  << active_augmented_665.size() << " FAILURES "
                  << failures_665.size() << " FAILURE_FNV " << std::hex
                  << body_fnv_literal(failures_665) << " MASK "
                  << failures_665.front() << std::dec << '\n';
        std::cout << "ENDPOINT665_WITNESS " << std::hex << witness_665
                  << std::dec << " ACTIVE 1 DISJOINT 1 FINAL_CARRIER "
                  << augmented_665.size() << " FINAL_FNV " << std::hex
                  << fnv_literal421(augmented_665) << std::dec
                  << " FINAL_ACTIVE " << active_closed_665.size()
                  << " BODIES " << closed_665.bodies << " FAILURES "
                  << closed_665.failures << " CHECKS " << closed_665.checks
                  << " MAX_CHECKS " << closed_665.max_checks << '\n';
        std::cout << "ENDPOINT664_HOSTILE PAIR 256,664 ACTIVE "
                  << active_hostile_664.size() << " FAILURES "
                  << failures_664.size() << " FAILURE_FNV " << std::hex
                  << body_fnv_literal(failures_664) << " MASKS";
        for (u32 body : failures_664) std::cout << ' ' << body;
        std::cout << std::dec << '\n';
        std::cout << "ENDPOINT664_WITNESS " << std::hex << witness_664
                  << std::dec << " ACTIVE 1 DISJOINT 2 FINAL_CARRIER "
                  << augmented_664.size() << " FINAL_FNV " << std::hex
                  << fnv_literal421(augmented_664) << std::dec
                  << " FINAL_ACTIVE " << active_closed_664.size()
                  << " BODIES " << closed_664.bodies << " FAILURES "
                  << closed_664.failures << " CHECKS " << closed_664.checks
                  << " MAX_CHECKS " << closed_664.max_checks << '\n';
        std::cout << "VERDICT PASS DETACHED_LITERAL_WALL_DESCENT_THROUGH_664\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "JOINT421_LITERAL_R670_ERROR " << error.what() << '\n';
        return 1;
    }
}
#endif

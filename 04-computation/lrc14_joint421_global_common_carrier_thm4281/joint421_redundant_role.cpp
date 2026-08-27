#define JOINT421_LITERAL_LIBRARY_ONLY
#include "verify_joint421_literal_r670.cpp"
#undef JOINT421_LITERAL_LIBRARY_ONLY

int main(int argc, char** argv) {
    try {
        require(argc == 8, "usage: role REPLAY_BAND P704 P416 P520 P384 P688 JOINT");
        const std::filesystem::path band(argv[1]);
        const Parsed prefix704 = parse_probe(argv[2]);
        const Parsed prefix416 = parse_probe(argv[3]);
        const Parsed prefix520 = parse_probe(argv[4]);
        const Parsed prefix384 = parse_probe(argv[5]);
        const Parsed prefix688 = parse_probe(argv[6]);

        std::vector<u32> carrier;
        std::set<u32> seen;
        for (const auto& [q, r] : EXPECTED_PAIRS) {
            const Parsed parsed = parse_probe(
                band / (std::to_string(q) + "_" + std::to_string(r) + ".out"));
            require(parsed.q == q && parsed.r == r,
                    "literal band transcript mismatch");
            append_unseen_literal421(parsed.deck, seen, carrier);
        }
        append_unseen_literal421(prefix704.deck, seen, carrier);
        const std::set<u32> round1_seen = seen;
        std::vector<u32> novel416;
        std::vector<u32> novel520;
        for (u32 repair : prefix416.deck)
            if (!round1_seen.count(repair)) novel416.push_back(repair);
        for (u32 repair : prefix520.deck)
            if (!round1_seen.count(repair)) novel520.push_back(repair);
        append_unseen_literal421(novel416, seen, carrier);
        append_unseen_literal421(novel520, seen, carrier);
        append_unseen_literal421(prefix384.deck, seen, carrier);
        append_unseen_literal421(prefix688.deck, seen, carrier);
        append_unseen_literal421(
            std::vector<u32>(THM4276_AUGMENTATION.begin(),
                             THM4276_AUGMENTATION.end()),
            seen, carrier);
        require(carrier.size() == 8524 &&
                    fnv_literal421(carrier) == UINT64_C(0x5ddb84a44f5d2ad7),
                "THM-4276 carrier changed");

        const std::vector<u32> joint = read_joint421(argv[7]);
        constexpr u32 target = UINT32_C(0x003c900c);
        const auto target_position =
            std::find(joint.begin(), joint.end(), target);
        require(target_position != joint.end() &&
                    std::distance(joint.begin(), target_position) == 318,
                "target deck position changed");
        std::vector<u32> core420 = joint;
        core420.erase(core420.begin() + 318);
        require(core420.size() == 420 &&
                    fnv_literal421(core420) ==
                        UINT64_C(0xe72c3c8b50ec6c6e),
                "420-core ledger changed");
        std::vector<u32> bodies;
        bodies.reserve(EXPECTED_BODIES);
        enumerate_bodies(0, 9, 0, bodies);
        require(bodies.size() == EXPECTED_BODIES,
                "literal body universe changed");
        const BodyAudit core_body = audit_bodies(bodies, core420);
        require(core_body.bodies == EXPECTED_BODIES &&
                    core_body.failures == 0 &&
                    core_body.checks == UINT64_C(405168772) &&
                    core_body.max_checks == 420,
                "420 core does not cover every nine-body");

        const std::vector<u32> active_base_a =
            active_literal421(256, 670, carrier);
        const std::vector<u32> active_base_b =
            active_literal421(384, 670, carrier);
        const std::vector<u32> failures_a =
            failures_literal421(bodies, active_base_a);
        const std::vector<u32> failures_b =
            failures_literal421(bodies, active_base_b);
        require(failures_a.size() == 1659 && failures_b.size() == 46 &&
                    body_fnv_literal(failures_a) ==
                        UINT64_C(0x970f004b0f9e2edb) &&
                    body_fnv_literal(failures_b) ==
                        UINT64_C(0x6c4fd2dab94cf6a9),
                "base obligations changed");

        const std::vector<u32> active_a = active_literal421(256, 670, joint);
        const std::vector<u32> active_b = active_literal421(384, 670, joint);
        const bool target_active_a =
            std::find(active_a.begin(), active_a.end(), target) != active_a.end();
        const bool target_active_b =
            std::find(active_b.begin(), active_b.end(), target) != active_b.end();

        u64 target_incidence_a = 0;
        u64 target_incidence_b = 0;
        u64 target_private_a = 0;
        u64 target_private_b = 0;
        u64 uncovered_without_target_a = 0;
        u64 uncovered_without_target_b = 0;
        Fnv private_ledger;
        auto audit = [&](unsigned row, const std::vector<u32>& obligations,
                         const std::vector<u32>& active, u64& incidence,
                         u64& private_count, u64& uncovered) {
            for (u32 body : obligations) {
                u64 hits = 0;
                bool target_hit = false;
                for (u32 repair : active) {
                    if ((repair & body) != 0) continue;
                    ++hits;
                    target_hit = target_hit || repair == target;
                }
                incidence += target_hit;
                if (target_hit && hits == 1) {
                    ++private_count;
                    private_ledger.add(row);
                    private_ledger.add(body);
                }
                uncovered += target_hit && hits == 1;
            }
        };
        audit(0, failures_a, active_a, target_incidence_a,
              target_private_a, uncovered_without_target_a);
        audit(1, failures_b, active_b, target_incidence_b,
              target_private_b, uncovered_without_target_b);
        require(target_active_a && target_active_b &&
                    target_incidence_a == 28 && target_incidence_b == 0 &&
                    target_private_a == 5 && target_private_b == 0 &&
                    uncovered_without_target_a == 5 &&
                    uncovered_without_target_b == 0 &&
                    private_ledger.state == UINT64_C(0x3a11d5bef4e6318d),
                "target endpoint-role ledger changed");

        std::cout << "JOINT421_REDUNDANT_ROLE_V1\n"
                  << "TARGET " << std::hex << target << std::dec
                  << " DECK_INDEX 318 ACTIVE_A " << target_active_a
                  << " ACTIVE_B " << target_active_b << '\n'
                  << "CORE420 " << core420.size() << " FNV " << std::hex
                  << fnv_literal421(core420) << std::dec
                  << " BODY_SCAN " << core_body.bodies
                  << " FAILURES " << core_body.failures
                  << " CHECKS " << core_body.checks
                  << " MAX_CHECKS " << core_body.max_checks << '\n'
                  << "INCIDENCE_A " << target_incidence_a
                  << " INCIDENCE_B " << target_incidence_b
                  << " TOTAL " << target_incidence_a + target_incidence_b
                  << '\n'
                  << "PRIVATE_A " << target_private_a
                  << " PRIVATE_B " << target_private_b
                  << " TOTAL " << target_private_a + target_private_b
                  << " PRIVATE_FNV " << std::hex << private_ledger.state
                  << std::dec << '\n'
                  << "DELETE_UNCOVERED_A " << uncovered_without_target_a
                  << " DELETE_UNCOVERED_B " << uncovered_without_target_b
                  << " TOTAL "
                  << uncovered_without_target_a + uncovered_without_target_b
                  << '\n'
                  << "VERDICT PASS EXACT_LITERAL_ENDPOINT_ROLE\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "JOINT421_REDUNDANT_ROLE_ERROR " << error.what() << '\n';
        return 1;
    }
}

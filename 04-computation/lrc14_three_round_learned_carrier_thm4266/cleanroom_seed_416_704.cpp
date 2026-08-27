// One-row wrapper around the THM-4254 literal-wall clean-room auditor.
// It verifies the newly discovered (416,704) full-prefix transcript without
// using endpoint cocycles, atoms, colex ranking, or a superset transform.

#define main thm4254_cleanroom_unused_main
#include "lrc14_endpoint_cascade_direct_wall_body_audit.cpp"
#undef main

int main(int argc, char** argv) {
    try {
        require(argc == 2, "usage: cleanroom-seed FULL_PREFIX_TRANSCRIPT");
        const Parsed parsed = parse_probe(argv[1]);
        require(parsed.q == 416 && parsed.r == 704,
                "unexpected clean-room seed pair");
        require(parsed.deck.size() == 2608 &&
                parsed.declared_prefix_fnv == UINT64_C(0x18ff663a123e684e),
                "seed prefix fingerprint changed");

        std::vector<u32> bodies;
        bodies.reserve(EXPECTED_BODIES);
        enumerate_bodies(0, 9, 0, bodies);
        require(bodies.size() == EXPECTED_BODIES,
                "recursive body universe changed");

        std::cout << "CARRIER_CEGAR_CLEANROOM_SEED_V1\n";
        Fnv labelled;
        labelled.add(static_cast<u64>(parsed.q));
        labelled.add(static_cast<u64>(parsed.r));
        labelled.add(parsed.deck.size());
        const PairAudit audit = audit_pair(parsed, bodies, labelled);
        std::cout << "SUMMARY PAIR 416,704 JOINT_CELLS "
                  << audit.joint_cells << " PREFIX " << audit.prefix_size
                  << " BODY_CHECKS " << audit.body_checks
                  << " LABELLED_FNV " << std::hex << labelled.state
                  << std::dec << "\n";
        std::cout << "VERDICT PASS LITERAL_WALL_DIRECT_MASS "
                     "RECURSIVE_BODY_ORDER_MINIMAL\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "CARRIER_CEGAR_CLEANROOM_ERROR " << error.what() << '\n';
        return 1;
    }
}

// Scratch-only dump of the exact active endpoint-589 q=50 carrier hypergraph.

#define ENDPOINT589_EXCHANGED_CARRIER_AUDIT_MAIN endpoint589_hidden_for_q50_structure
#include "endpoint589_exchanged_carrier_audit.cpp"
#undef ENDPOINT589_EXCHANGED_CARRIER_AUDIT_MAIN

int main(int argc, char** argv) {
    try {
        require(argc == 19,
                "usage: dump JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE OLD_UNION "
                "REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE ADDITION593 "
                "DELETE593 COVER43 DELETE43 COVER9 DELETE9 OUTPUT");
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
        const std::vector<u32> base =
            reconstruct_final_carrier_quotient(argv, joint);
        const std::unordered_set<u32> base_set(base.begin(), base.end());
        const std::vector<u32> additions = read_cover_quotient(argv[16]);
        const std::vector<u32> deletions = read_delete9_589(argv[17]);
        require(masks_fnv_agent(additions) == kAddition9Fnv &&
                    masks_fnv_agent(deletions) == kDeletion9Fnv,
                "exchange ledgers changed");
        const std::unordered_set<u32> deletion_set(deletions.begin(),
                                                   deletions.end());
        std::vector<u32> carrier;
        for (u32 mask : base)
            if (!deletion_set.contains(mask)) carrier.push_back(mask);
        for (u32 mask : additions) {
            require(!base_set.contains(mask), "addition overlaps base");
            carrier.push_back(mask);
        }
        require(carrier.size() == 3925 &&
                    masks_fnv_agent(carrier) == kExchangedCarrier589Fnv,
                "exchanged carrier changed");
        const Geometry geometry = build_geometry(50, 589);
        std::ofstream output(argv[18]);
        require(static_cast<bool>(output), "cannot create active ledger");
        output << "mask_hex,rank,joint,margin_ticks\n";
        u64 active = 0, active_joint = 0;
        Fnv ledger;
        for (u32 mask : carrier) {
            const auto value = margin(geometry, mask);
            if (value.ticks < 0) continue;
            ++active;
            active_joint += joint_set.contains(mask);
            ledger.add(mask);
            output << hex8(mask) << ',' << std::popcount(mask) << ','
                   << joint_set.contains(mask) << ',' << decimal(value.ticks)
                   << '\n';
        }
        require(output.good() && active == 1398 && active_joint == 118 &&
                    ledger.state == UINT64_C(0xc075113890c7f5e1),
                "active q50 identity changed");
        std::cout << "ENDPOINT589_Q50_ACTIVE_CARRIER_V1\n"
                  << "ACTIVE " << active << " JOINT " << active_joint
                  << " NONJOINT " << active - active_joint << " FNV "
                  << std::hex << ledger.state << std::dec << '\n'
                  << "SCOPE FINITE_EXACT_FIXED_EXCHANGED_CARRIER_Q50\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Q50_ACTIVE_DUMP_ERROR " << error.what() << '\n';
        return 1;
    }
}

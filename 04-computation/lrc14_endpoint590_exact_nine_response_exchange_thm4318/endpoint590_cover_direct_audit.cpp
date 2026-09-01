// Direct arithmetic verification of a frozen nine-mask response cover for
// all 100 endpoint-590 failures. No solver state is trusted by this replay.

#define ENDPOINT590_RESPONSE_STRUCTURE_MAIN endpoint590_response_hidden_for_cover_audit
#include "endpoint590_response_structure.cpp"
#undef ENDPOINT590_RESPONSE_STRUCTURE_MAIN

namespace {
struct CoverMask590 {
    u32 mask = 0;
    unsigned declared_rank = 0;
    unsigned declared_weight = 0;
};

std::vector<CoverMask590> read_cover_candidate590(
    const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open endpoint590 cover");
    std::string line;
    require(std::getline(input, line) &&
                line == "mask_hex,rank,weight,w0,w1",
            "endpoint590 cover header changed");
    std::vector<CoverMask590> result;
    std::set<u32> distinct;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        std::string token, ignored0, ignored1;
        CoverMask590 entry;
        fields >> token >> entry.declared_rank >> entry.declared_weight
               >> ignored0 >> ignored1;
        entry.mask = parse_mask_agent(token);
        require(fields && entry.declared_rank == 9 &&
                    std::popcount(entry.mask) == entry.declared_rank &&
                    distinct.insert(entry.mask).second,
                "endpoint590 cover rank/distinctness changed");
        result.push_back(entry);
    }
    require(input.eof() && result.size() == 9,
            "endpoint590 cover count changed");
    return result;
}
}

int main(int argc, char** argv) {
    try {
        require(argc == 19,
                "usage: direct JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE "
                "OLD_UNION REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE "
                "ADDITION593 DELETE593 COVER43 DELETE43 FAILURES COVER "
                "INCIDENCE_CSV");
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::vector<u32> carrier =
            reconstruct_final_carrier590(argc, argv, joint);
        const std::unordered_set<u32> carrier_set(carrier.begin(),
                                                  carrier.end());
        const std::vector<u32> bodies = read_failure_bodies590(argv[16]);
        const std::vector<CoverMask590> cover =
            read_cover_candidate590(argv[17]);
        const Geometry geometry = build_geometry(210, 590);

        Fnv cover_fnv;
        for (const auto& entry : cover) {
            require(!carrier_set.contains(entry.mask),
                    "cover candidate already belongs to carrier");
            require(margin(geometry, entry.mask).ticks >= 0,
                    "cover candidate is inactive at (210,590)");
            cover_fnv.add(entry.mask);
        }
        require(cover_fnv.state == UINT64_C(0xd1cf49e4b811b958),
                "endpoint590 nine-mask cover identity changed");

        std::vector<unsigned> body_hits(bodies.size(), 0);
        std::vector<unsigned> mask_incidences(cover.size(), 0);
        std::vector<Fnv> mask_response_fnv(cover.size());
        for (std::size_t index = 0; index < cover.size(); ++index) {
            for (std::size_t body_index = 0; body_index < bodies.size();
                 ++body_index) {
                if ((cover[index].mask & bodies[body_index]) != 0) continue;
                ++body_hits[body_index];
                ++mask_incidences[index];
                mask_response_fnv[index].add(210);
                mask_response_fnv[index].add(590);
                mask_response_fnv[index].add(bodies[body_index]);
            }
            require(mask_incidences[index] == cover[index].declared_weight,
                    "direct incidence differs from cover declaration");
        }

        std::array<u64, 10> hit_distribution{};
        Fnv hit_ledger;
        unsigned minimum_hits = 10;
        unsigned maximum_hits = 0;
        u64 total_incidences = 0;
        for (std::size_t index = 0; index < bodies.size(); ++index) {
            const unsigned hits = body_hits[index];
            require(hits >= 1 && hits <= cover.size(),
                    "nine-mask cover leaves a failure uncovered");
            ++hit_distribution[hits];
            minimum_hits = std::min(minimum_hits, hits);
            maximum_hits = std::max(maximum_hits, hits);
            total_incidences += hits;
            hit_ledger.add(210); hit_ledger.add(590);
            hit_ledger.add(bodies[index]); hit_ledger.add(hits);
        }

        std::ofstream incidence_out(argv[18]);
        require(static_cast<bool>(incidence_out),
                "cannot create endpoint590 incidence ledger");
        incidence_out << "mask_hex,rank,responses,response_fnv\n";
        Fnv incidence_ledger;
        for (std::size_t index = 0; index < cover.size(); ++index) {
            incidence_out << hex8(cover[index].mask) << ','
                          << cover[index].declared_rank << ','
                          << mask_incidences[index] << ',' << std::hex
                          << mask_response_fnv[index].state << std::dec << '\n';
            incidence_ledger.add(cover[index].mask);
            incidence_ledger.add(cover[index].declared_rank);
            incidence_ledger.add(mask_incidences[index]);
            incidence_ledger.add(mask_response_fnv[index].state);
        }
        require(incidence_out.good(), "incidence ledger write failed");
        require(mask_incidences ==
                    std::vector<unsigned>{11, 6, 14, 4, 18, 16, 24, 15, 28} &&
                    total_incidences == 136 && minimum_hits == 1 &&
                    maximum_hits == 3 && hit_distribution[1] == 69 &&
                    hit_distribution[2] == 26 && hit_distribution[3] == 5 &&
                    hit_ledger.state == UINT64_C(0x57b3865558b0767b) &&
                    incidence_ledger.state == UINT64_C(0x921483aaf38df968),
                "endpoint590 nine-mask direct census changed");

        std::cout << "LRC14_ENDPOINT590_NINE_RESPONSE_DIRECT_AUDIT_V1\n"
                  << "CARRIER " << carrier.size() << " FNV " << std::hex
                  << masks_fnv_agent(carrier) << std::dec << '\n'
                  << "FAILURES " << bodies.size() << " FAILURE_FNV "
                  << std::hex << kFailureFnv590 << std::dec << '\n'
                  << "COVER " << cover.size() << " RANK9 " << cover.size()
                  << " COVER_FNV " << std::hex << cover_fnv.state << std::dec
                  << '\n';
        for (std::size_t index = 0; index < cover.size(); ++index) {
            std::cout << "MASK " << hex8(cover[index].mask) << " RANK "
                      << cover[index].declared_rank << " RESPONSES "
                      << mask_incidences[index] << " RESPONSE_FNV " << std::hex
                      << mask_response_fnv[index].state << std::dec << '\n';
        }
        std::cout << "TOTAL_INCIDENCES " << total_incidences << " HIT_RANGE "
                  << minimum_hits << ".." << maximum_hits << " HIT_DISTRIBUTION";
        for (unsigned hits = 1; hits <= cover.size(); ++hits)
            if (hit_distribution[hits])
                std::cout << ' ' << hits << ':' << hit_distribution[hits];
        std::cout << " HIT_LEDGER_FNV " << std::hex << hit_ledger.state
                  << " INCIDENCE_LEDGER_FNV " << incidence_ledger.state
                  << std::dec << '\n'
                  << "UNCOVERED 0\n"
                  << "SCOPE DIRECT_ARITHMETIC_FIXED_FAILURE_COVER_"
                     "ADDITIONS_ONLY_NO_DELETIONS_NO_EXCHANGE_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT590_DIRECT_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}

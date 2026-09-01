// Complete rank-8/rank-9 response census for the 100 failures left by the
// unchanged THM-4313/THM-4314 carrier at (210,590). This only describes the
// fixed labelled response universe; its greedy cover is an upper bound.

#define ENDPOINT590_FINAL_CARRIER_AUDIT_MAIN endpoint590_capacity_hidden_for_response_structure
#include "endpoint590_final_carrier_audit.cpp"
#undef ENDPOINT590_FINAL_CARRIER_AUDIT_MAIN

#include <array>
#include <map>

namespace {
constexpr std::size_t kFailureCount590 = 100;
constexpr u64 kFailureFnv590 = UINT64_C(0x8d19cba1e86e53b5);

struct Bits590 {
    u64 low = 0;
    u64 high = 0;

    auto operator<=>(const Bits590&) const = default;
};

struct Signature590 {
    u64 count8 = 0;
    u64 count9 = 0;
    u32 least8 = 0;
    u32 least9 = 0;
};

unsigned weight590(Bits590 bits) {
    return std::popcount(bits.low) + std::popcount(bits.high);
}

unsigned intersection_weight590(Bits590 left, Bits590 right) {
    return std::popcount(left.low & right.low) +
           std::popcount(left.high & right.high);
}

Bits590 intersect590(Bits590 left, Bits590 right) {
    return {left.low & right.low, left.high & right.high};
}

void subtract590(Bits590& left, Bits590 right) {
    left.low &= ~right.low;
    left.high &= ~right.high;
}

std::vector<u32> read_failure_bodies590(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open endpoint590 failures");
    std::string line;
    require(std::getline(input, line) && line == "q,r,body_hex",
            "endpoint590 failure header changed");
    std::vector<u32> bodies;
    std::set<u32> distinct;
    Fnv ledger;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream fields(line);
        int q = 0, r = 0;
        std::string token;
        fields >> q >> r >> token;
        require(fields && q == 210 && r == 590,
                "endpoint590 failure escaped (210,590)");
        const u32 body = parse_mask_agent(token);
        require(std::popcount(body) == 9 && distinct.insert(body).second,
                "endpoint590 failure body rank/distinctness changed");
        bodies.push_back(body);
        ledger.add(q); ledger.add(r); ledger.add(body);
    }
    require(input.eof() && bodies.size() == kFailureCount590 &&
                ledger.state == kFailureFnv590,
            "endpoint590 failure identity changed");
    return bodies;
}

std::vector<u32> reconstruct_final_carrier590(
    int argc, char** argv, const std::vector<u32>& joint) {
    require(argc == 19,
            "usage: response JOINT BASE8951 ADD45 SUFFIX9 UNIVERSE OLD_UNION "
            "REPAIRS76 ADD4 DELETE73 ADD10 FINAL_DELETE ADDITION593 DELETE593 "
            "COVER43 DELETE43 FAILURES RESPONSE_CSV GREEDY_CSV");
    const std::unordered_set<u32> joint_set(joint.begin(), joint.end());
    const std::vector<u32> base = reconstruct_c3925_593(
        argv[2], argv[3], argv[4], argv[7], argv[8], argv[9], argv[10],
        argv[11]);
    (void)read_pairs_repaired(argv[5], 22647, kUniverseFnv);
    (void)read_pairs_repaired(argv[6], 1624, kOldUnionFnv);
    const auto add593 = read_flexible_masks593(argv[12], 1, kAddition593Fnv);
    const auto del593 = read_flexible_masks593(argv[13], 1, kDeletion593Fnv);
    std::vector<u32> carrier593;
    for (u32 mask : base)
        if (mask != del593.front()) carrier593.push_back(mask);
    require(carrier593.size() + 1 == base.size(), "THM4311 deletion absent");
    carrier593.push_back(add593.front());
    require(carrier593.size() == 3925 &&
                masks_fnv_agent(carrier593) == kExchangedCarrierFnv,
            "THM4311 carrier changed");
    const std::unordered_set<u32> carrier593_set(carrier593.begin(),
                                                 carrier593.end());
    const auto additions = read_cover590(argv[14]);
    const auto deletions = read_deletions590(argv[15]);
    require(masks_fnv_agent(additions) == kCover43Fnv &&
                masks_fnv_agent(deletions) == kDeletion43Fnv,
            "THM4313 exchange ledger changed");
    const std::unordered_set<u32> deletion_set(deletions.begin(),
                                               deletions.end());
    std::vector<u32> carrier;
    for (u32 mask : carrier593) {
        if (deletion_set.contains(mask)) {
            require(!joint_set.contains(mask), "THM4313 deleted joint mask");
        } else {
            carrier.push_back(mask);
        }
    }
    require(carrier.size() + deletions.size() == carrier593.size(),
            "THM4313 deletion absent");
    for (u32 mask : additions) {
        require(!joint_set.contains(mask) && !carrier593_set.contains(mask),
                "THM4313 addition overlaps inherited carrier/joint deck");
        carrier.push_back(mask);
    }
    require(carrier.size() == 3925 &&
                masks_fnv_agent(carrier) == kFinalCarrier590Fnv &&
                std::set<u32>(carrier.begin(), carrier.end()).size() ==
                    carrier.size(),
            "THM4313/THM4314 final carrier changed");
    return carrier;
}
}

#ifndef ENDPOINT590_RESPONSE_STRUCTURE_MAIN
#define ENDPOINT590_RESPONSE_STRUCTURE_MAIN main
#endif

int ENDPOINT590_RESPONSE_STRUCTURE_MAIN(int argc, char** argv) {
    try {
        const auto joint = read_masks_agent(argv[1], kJointCountAgent,
                                            kJointFnvAgent, 8);
        const std::vector<u32> carrier =
            reconstruct_final_carrier590(argc, argv, joint);
        const std::unordered_set<u32> carrier_set(carrier.begin(),
                                                  carrier.end());
        const std::vector<u32> bodies = read_failure_bodies590(argv[16]);
        const Geometry geometry = build_geometry(210, 590);

        std::map<Bits590, Signature590> signatures;
        std::array<u64, 2> active_count{};
        std::array<u64, 2> responder_count{};
        std::array<Fnv, 2> responder_fnv;
        std::array<unsigned, 2> maximum_weight{};
        std::array<u32, 2> least_maximum{};
        std::array<u64, kFailureCount590> multiplicity{};
        for (unsigned rank : {8U, 9U}) {
            const std::size_t slot = rank - 8;
            u64 universe_count = 0;
            const u32 limit = UINT32_C(1) << 30;
            for (u32 mask = (UINT32_C(1) << rank) - 1; mask < limit;
                 mask = next_combination(mask)) {
                ++universe_count;
                if (margin(geometry, mask).ticks < 0) continue;
                ++active_count[slot];
                Bits590 response;
                for (std::size_t index = 0; index < bodies.size(); ++index) {
                    if ((mask & bodies[index]) != 0) continue;
                    if (index < 64)
                        response.low |= UINT64_C(1) << index;
                    else
                        response.high |= UINT64_C(1) << (index - 64);
                    ++multiplicity[index];
                }
                const unsigned weight = weight590(response);
                if (weight == 0) continue;
                require(!carrier_set.contains(mask),
                        "existing carrier mask responds to frozen failure");
                ++responder_count[slot];
                responder_fnv[slot].add(mask);
                if (weight > maximum_weight[slot] ||
                    (weight == maximum_weight[slot] &&
                     (least_maximum[slot] == 0 || mask < least_maximum[slot]))) {
                    maximum_weight[slot] = weight;
                    least_maximum[slot] = mask;
                }
                auto& record = signatures[response];
                u64& count = rank == 8 ? record.count8 : record.count9;
                u32& least = rank == 8 ? record.least8 : record.least9;
                ++count;
                if (least == 0 || mask < least) least = mask;
            }
            require(universe_count == (rank == 8 ? UINT64_C(5852925)
                                                  : UINT64_C(14307150)),
                    "rank universe count changed");
        }

        Fnv multiplicity_fnv;
        u64 minimum_multiplicity = std::numeric_limits<u64>::max();
        u64 maximum_multiplicity = 0;
        for (std::size_t index = 0; index < bodies.size(); ++index) {
            multiplicity_fnv.add(210); multiplicity_fnv.add(590);
            multiplicity_fnv.add(bodies[index]);
            multiplicity_fnv.add(multiplicity[index]);
            minimum_multiplicity = std::min(minimum_multiplicity,
                                            multiplicity[index]);
            maximum_multiplicity = std::max(maximum_multiplicity,
                                            multiplicity[index]);
        }
        require(minimum_multiplicity > 0,
                "some endpoint590 failure has no rank8/9 response");
        require(active_count == std::array<u64, 2>{1073813, 6446114} &&
                    responder_count == std::array<u64, 2>{36285, 568812} &&
                    responder_fnv[0].state == UINT64_C(0x12fda2562de3e9a6) &&
                    responder_fnv[1].state == UINT64_C(0xb531e8d2059a473c) &&
                    maximum_weight == std::array<unsigned, 2>{19, 29} &&
                    least_maximum == std::array<u32, 2>{UINT32_C(0x0ba04140),
                                                        UINT32_C(0x0a0a4416)} &&
                    signatures.size() == 14368 &&
                    minimum_multiplicity == 5570 &&
                    maximum_multiplicity == 32922 &&
                    multiplicity_fnv.state == UINT64_C(0xa8e79ecbc16f3c1c),
                "endpoint590 response census changed");

        std::ofstream response_out(argv[17]);
        std::ofstream greedy_out(argv[18]);
        require(response_out && greedy_out, "cannot create response ledgers");
        response_out << "w0,w1,weight,count8,least8,count9,least9\n";
        for (const auto& [bits, record] : signatures) {
            response_out << std::hex << std::setw(16) << std::setfill('0')
                         << bits.low << ',' << std::setw(16) << bits.high
                         << std::dec << ',' << weight590(bits) << ','
                         << record.count8 << ',';
            if (record.least8) response_out << hex8(record.least8);
            response_out << ',' << record.count9 << ',';
            if (record.least9) response_out << hex8(record.least9);
            response_out << '\n';
        }

        Bits590 uncovered{UINT64_MAX, (UINT64_C(1) << 36) - 1};
        unsigned covered = 0;
        unsigned step = 0;
        greedy_out << "step,mask_hex,rank,newly_covered,total_covered,w0,w1\n";
        while (weight590(uncovered) != 0) {
            unsigned best_gain = 0;
            u32 best_mask = 0;
            unsigned best_rank = 0;
            Bits590 best_response;
            for (const auto& [bits, record] : signatures) {
                const unsigned gain = intersection_weight590(bits, uncovered);
                if (gain == 0) continue;
                u32 representative = 0;
                unsigned rank = 0;
                for (const auto [candidate, candidate_rank] :
                     {std::pair{record.least8, 8U},
                      std::pair{record.least9, 9U}}) {
                    if (candidate != 0 &&
                        (representative == 0 || candidate < representative)) {
                        representative = candidate;
                        rank = candidate_rank;
                    }
                }
                if (gain > best_gain ||
                    (gain == best_gain && representative < best_mask)) {
                    best_gain = gain;
                    best_mask = representative;
                    best_rank = rank;
                    best_response = bits;
                }
            }
            require(best_gain != 0, "greedy response cover stalled");
            subtract590(uncovered, best_response);
            covered += best_gain;
            greedy_out << ++step << ',' << hex8(best_mask) << ',' << best_rank
                       << ',' << best_gain << ',' << covered << ',' << std::hex
                       << std::setw(16) << std::setfill('0') << best_response.low
                       << ',' << std::setw(16) << best_response.high << std::dec
                       << '\n';
        }
        require(step == 10, "endpoint590 greedy cover size changed");
        require(response_out.good() && greedy_out.good(),
                "response ledger write failed");

        std::cout << "LRC14_ENDPOINT590_RESPONSE_STRUCTURE_V1\n"
                  << "FAILURES " << bodies.size() << " FAILURE_FNV "
                  << std::hex << kFailureFnv590 << std::dec << '\n'
                  << "RANK8 ACTIVE " << active_count[0] << " RESPONDERS "
                  << responder_count[0] << " RESPONDER_FNV " << std::hex
                  << responder_fnv[0].state << std::dec << " MAX_WEIGHT "
                  << maximum_weight[0] << " LEAST_MAX " << hex8(least_maximum[0])
                  << '\n'
                  << "RANK9 ACTIVE " << active_count[1] << " RESPONDERS "
                  << responder_count[1] << " RESPONDER_FNV " << std::hex
                  << responder_fnv[1].state << std::dec << " MAX_WEIGHT "
                  << maximum_weight[1] << " LEAST_MAX " << hex8(least_maximum[1])
                  << '\n'
                  << "SIGNATURES " << signatures.size()
                  << " MULTIPLICITY_RANGE " << minimum_multiplicity << ".."
                  << maximum_multiplicity << " MULTIPLICITY_FNV " << std::hex
                  << multiplicity_fnv.state << std::dec << '\n'
                  << "GREEDY_COVER_SIZE " << step << '\n'
                  << "SCOPE COMPLETE_RANK8_RANK9_RESPONSE_CENSUS_GREEDY_"
                     "UPPER_BOUND_NOT_MINIMUM_NO_EXCHANGE_NO_PHYSICAL_ENTRY_"
                     "NO_LRC14\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ENDPOINT590_RESPONSE_ERROR " << error.what() << '\n';
        return 1;
    }
}

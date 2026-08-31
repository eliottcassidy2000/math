// Scratch-only detached audit of the exact four-mask (100,629) repair.

#define main r632_detached_hostile_main
#include "r632_detached_hostile_survivor.cpp"
#undef main

namespace {

constexpr std::array<u32, 28> kBodies629 = {
    0x031c5208, 0x035c5008, 0x037c4008, 0x039c3008,
    0x039e1008, 0x064c5401, 0x070c7008, 0x07143408,
    0x071e1400, 0x072c6008, 0x07485401, 0x074c5001,
    0x078c1088, 0x078c3008, 0x07983008, 0x17087008,
    0x17107008, 0x17286008, 0x17306008, 0x17883008,
    0x23346008, 0x27841088, 0x27903008, 0x27c41008,
    0x33105088, 0x33107008, 0x33306001, 0x33306008};

constexpr std::array<u32, 4> kWitness629 = {
    0x002ac4c0, 0x3882a082, 0x0041c325, 0x08c28e40};

constexpr std::array<u32, 4> kResponses629 = {
    0x0c00000, 0x0000c27, 0x0687198, 0xf178240};

u32 response629(u32 mask) {
    u32 response = 0;
    for (std::size_t i = 0; i < kBodies629.size(); ++i)
        if ((mask & kBodies629[i]) == 0) response |= u32{1} << i;
    return response;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 2, "usage: r629_witness RESPONSES_CSV");
        const Geometry geometry = build_geometry(100, 629);
        Fnv body_ledger;
        for (u32 body : kBodies629) {
            require(std::popcount(body) == 9, "body rank changed");
            body_ledger.add(body);
        }
        require(body_ledger.state == UINT64_C(0x61ee3b61e7de8594),
                "body identity changed");

        std::ofstream csv(argv[1]);
        require(static_cast<bool>(csv), "cannot create response CSV");
        csv << "body_ordinal,body_hex,body_labels,hit_count,consumer_masks\n";
        Fnv witness_ledger;
        Fnv response_ledger;
        u32 joined = 0;
        std::array<u64, 28> hit_counts{};
        std::array<std::string, 28> consumers{};
        std::cout << "R629_MIXED_WITNESS_DETACHED_V1\n"
                  << "PAIR 100,629 GRID " << geometry.grid << " CELLS "
                  << geometry.cells << " FAILURE_CLASSES_RANK_LE9 "
                  << geometry.classes.size() << '\n';
        for (std::size_t i = 0; i < kWitness629.size(); ++i) {
            const u32 mask = kWitness629[i];
            const unsigned rank = std::popcount(mask);
            require(rank == (i == 0 ? 8 : 9), "witness rank changed");
            const Margin value = margin(geometry, mask);
            const u32 response = response629(mask);
            require(value.ticks > 0 && response == kResponses629[i],
                    "witness activity/response changed");
            witness_ledger.add(mask);
            response_ledger.add(mask);
            response_ledger.add(rank);
            response_ledger.add(response);
            response_ledger.add(static_cast<u64>(value.mass));
            add_i128(response_ledger, value.ticks);
            joined |= response;
            for (std::size_t body = 0; body < kBodies629.size(); ++body)
                if ((response >> body) & 1u) {
                    ++hit_counts[body];
                    if (!consumers[body].empty()) consumers[body] += ';';
                    consumers[body] += hex8(mask);
                }
            std::cout << "WITNESS " << i << " MASK " << hex8(mask)
                      << " RANK " << rank << " LABELS {" << labels(mask)
                      << "} RESPONSE " << std::hex << std::setw(7)
                      << std::setfill('0') << response << std::dec
                      << std::setfill(' ') << " COVER "
                      << std::popcount(response) << " MASS " << value.mass
                      << " MARGIN_TICKS63 " << decimal(value.ticks) << '\n';
        }
        require(joined == (u32{1} << 28) - 1,
                "witness does not cover all bodies");
        u64 incidences = 0;
        for (std::size_t i = 0; i < kBodies629.size(); ++i) {
            require(hit_counts[i] > 0, "body lacks consumer");
            incidences += hit_counts[i];
            csv << i << ',' << hex8(kBodies629[i]) << ",\"{"
                << labels(kBodies629[i]) << "}\"," << hit_counts[i]
                << ",\"" << consumers[i] << "\"\n";
            std::cout << "BODY " << i << ' ' << hex8(kBodies629[i])
                      << " LABELS {" << labels(kBodies629[i]) << "} HITS "
                      << hit_counts[i] << " CONSUMERS " << consumers[i]
                      << '\n';
        }
        require(csv.good(), "response CSV write failed");
        std::cout << "WITNESS_FNV " << std::hex << witness_ledger.state
                  << " RESPONSE_FNV " << response_ledger.state << std::dec
                  << " INCIDENCES " << incidences << '\n'
                  << "LOWER_CERTIFICATE separate exact rational dual 7/2 "
                     "over all 419 response classes\n"
                  << "SCOPE IMPORT_FREE_LITERAL_WALL_FIXED_PAIR_LABELLED_"
                     "FOUR_MASK_UPPER_CERTIFICATE_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_MIXED_FOUR_MASK_WITNESS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "R629_WITNESS_ERROR " << error.what() << '\n';
        return 1;
    }
}

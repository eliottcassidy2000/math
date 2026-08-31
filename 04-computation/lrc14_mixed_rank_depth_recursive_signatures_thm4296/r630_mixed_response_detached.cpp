#define main r632_detached_hostile_main
#include "r632_detached_hostile_survivor.cpp"
#undef main

namespace {

constexpr std::array<u32, 2> kBodies630 = {0x27841208, 0x27c81008};

template <class Callback>
void choose_rank_rec(const std::vector<unsigned char>& positions,
                     std::size_t start, unsigned chosen, unsigned target,
                     u32 mask, Callback& callback) {
    if (chosen == target) {
        callback(mask);
        return;
    }
    const std::size_t needed = target - chosen;
    for (std::size_t index = start;
         index + needed <= positions.size(); ++index)
        choose_rank_rec(positions, index + 1, chosen + 1, target,
                        mask | (u32{1} << positions[index]), callback);
}

template <class Callback>
u64 enumerate_rank(u32 forbidden, unsigned rank, Callback callback) {
    std::vector<unsigned char> positions;
    for (unsigned bit = 0; bit < 30; ++bit)
        if ((forbidden & (u32{1} << bit)) == 0)
            positions.push_back(static_cast<unsigned char>(bit));
    if (positions.size() < rank) return 0;
    u64 count = 0;
    auto counted = [&](u32 mask) {
        require(std::popcount(mask) == rank && (mask & forbidden) == 0,
                "rank response enumeration escaped");
        callback(mask);
        ++count;
    };
    choose_rank_rec(positions, 0, 0, rank, 0, counted);
    return count;
}

struct ResponseAudit630 {
    unsigned rank = 0;
    std::array<u64, 4> counts{};
    std::array<u32, 4> least{};
    std::array<bool, 4> have{};
    std::array<Fnv, 4> response_ledgers{};
    u64 raw_candidates = 0;
    u64 active_candidates = 0;
    Fnv ledger;
};

ResponseAudit630 audit_rank(const Geometry& geometry, unsigned rank) {
    std::map<u32, unsigned> raw;
    for (std::size_t index = 0; index < kBodies630.size(); ++index)
        enumerate_rank(kBodies630[index], rank,
                       [&](u32 mask) { raw[mask] |= 1u << index; });
    ResponseAudit630 out;
    out.rank = rank;
    out.raw_candidates = raw.size();
    for (const auto& [mask, response] : raw) {
        const Margin value = margin(geometry, mask);
        if (value.ticks < 0) continue;
        ++out.active_candidates;
        ++out.counts[response];
        out.response_ledgers[response].add(mask);
        if (!out.have[response] || mask < out.least[response]) {
            out.have[response] = true;
            out.least[response] = mask;
        }
        out.ledger.add(mask);
        out.ledger.add(response);
        out.ledger.add(static_cast<u64>(value.mass));
        add_i128(out.ledger, value.ticks);
    }
    return out;
}

}  // namespace


int main() {
    try {
        const Geometry geometry = build_geometry(100, 630);
        std::cout << "R630_MIXED_RESPONSE_DETACHED_V1\n"
                  << "PAIR 100,630 GRID " << geometry.grid << " CELLS "
                  << geometry.cells << " FAILURE_CLASSES_RANK_LE9 "
                  << geometry.classes.size() << " BODIES "
                  << hex8(kBodies630[0]) << ',' << hex8(kBodies630[1])
                  << '\n';
        u64 total_full = 0;
        u32 least_full = 0;
        unsigned least_rank = 0;
        Fnv aggregate;
        for (unsigned rank : {8u, 9u}) {
            const ResponseAudit630 audit = audit_rank(geometry, rank);
            require(audit.counts[1] > 0 && audit.counts[2] > 0,
                    "one singleton response absent");
            aggregate.add(rank);
            aggregate.add(audit.raw_candidates);
            aggregate.add(audit.active_candidates);
            aggregate.add(audit.ledger.state);
            for (unsigned response = 1; response <= 3; ++response) {
                aggregate.add(response);
                aggregate.add(audit.counts[response]);
                aggregate.add(audit.least[response]);
            }
            std::cout << "RANK " << rank << " RAW_CANDIDATES "
                      << audit.raw_candidates << " ACTIVE_CANDIDATES "
                      << audit.active_candidates << " RESPONSE1 "
                      << audit.counts[1] << " LEAST1 " << hex8(audit.least[1])
                      << " RESPONSE2 " << audit.counts[2] << " LEAST2 "
                      << hex8(audit.least[2]) << " FULL_RESPONDERS "
                      << audit.counts[3] << " LEAST_FULL "
                      << hex8(audit.least[3]) << " FULL_FNV " << std::hex
                      << audit.response_ledgers[3].state << " LEDGER_FNV "
                      << audit.ledger.state << std::dec << '\n';
            total_full += audit.counts[3];
            if (audit.counts[3] > 0 &&
                (least_rank == 0 || audit.least[3] < least_full)) {
                least_full = audit.least[3];
                least_rank = rank;
            }
        }
        require(total_full > 0, "no mixed full responder");
        const Margin witness = margin(geometry, least_full);
        require(witness.ticks > 0 &&
                    (least_full & kBodies630[0]) == 0 &&
                    (least_full & kBodies630[1]) == 0,
                "least full witness invalid");
        std::cout << "MIXED_FULL_RESPONDERS " << total_full
                  << " LEAST_FULL " << hex8(least_full) << " RANK "
                  << least_rank << " LABELS {" << labels(least_full)
                  << "} MASS " << witness.mass << " MARGIN_TICKS63 "
                  << decimal(witness.ticks) << " AGGREGATE_FNV " << std::hex
                  << aggregate.state << std::dec << '\n'
                  << "MINIMUM_ADDITIONS 1 LOWER nonempty-obligation-set "
                     "UPPER least-full-responder\n"
                  << "SCOPE IMPORT_FREE_LITERAL_WALL_FIXED_PAIR_LABELLED_"
                     "RANK8_OR_RANK9_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_MIXED_MINIMUM_ONE\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "R630_RESPONSE_ERROR " << error.what() << '\n';
        return 1;
    }
}

// Scratch-only detached rank-8/rank-9 response quotient and exact cover for
// the four labelled failures at (100,628).

#define main r632_detached_hostile_main
#include "r632_detached_hostile_survivor.cpp"
#undef main

namespace {

constexpr std::array<u32, 4> kBodies628 = {
    0x05346408, 0x15306408, 0x17581400, 0x27d01008};

template <class Callback>
void choose_rank628(const std::vector<unsigned char>& positions,
                    std::size_t start, unsigned chosen, unsigned target,
                    u32 mask, Callback& callback) {
    if (chosen == target) {
        callback(mask);
        return;
    }
    const std::size_t needed = target - chosen;
    for (std::size_t index = start;
         index + needed <= positions.size(); ++index)
        choose_rank628(positions, index + 1, chosen + 1, target,
                       mask | (u32{1} << positions[index]), callback);
}

template <class Callback>
void enumerate_disjoint628(u32 body, unsigned rank, Callback callback) {
    std::vector<unsigned char> positions;
    for (unsigned bit = 0; bit < 30; ++bit)
        if ((body & (u32{1} << bit)) == 0)
            positions.push_back(static_cast<unsigned char>(bit));
    require(positions.size() == 21, "body complement changed");
    choose_rank628(positions, 0, 0, rank, 0, callback);
}

struct Class628 {
    u64 count8 = 0;
    u64 count9 = 0;
    u32 least8 = 0;
    u32 least9 = 0;
};

struct Rank628 {
    unsigned rank = 0;
    u64 raw = 0;
    u64 active = 0;
    Fnv ledger;
};

Rank628 audit_rank628(const Geometry& geometry, unsigned rank,
                      std::array<Class628, 16>& classes) {
    std::map<u32, unsigned> raw;
    for (std::size_t i = 0; i < kBodies628.size(); ++i)
        enumerate_disjoint628(kBodies628[i], rank,
            [&](u32 mask) { raw[mask] |= 1u << i; });
    Rank628 audit{rank, raw.size(), 0, {}};
    for (const auto& [mask, response] : raw) {
        const Margin value = margin(geometry, mask);
        if (value.ticks < 0) continue;
        ++audit.active;
        audit.ledger.add(mask);
        audit.ledger.add(response);
        audit.ledger.add(static_cast<u64>(value.mass));
        add_i128(audit.ledger, value.ticks);
        Class628& c = classes[response];
        if (rank == 8) {
            ++c.count8;
            if (c.count8 == 1 || mask < c.least8) c.least8 = mask;
        } else {
            ++c.count9;
            if (c.count9 == 1 || mask < c.least9) c.least9 = mask;
        }
    }
    return audit;
}

struct Solution628 {
    unsigned size = 100;
    unsigned rank8 = 0;
    std::vector<unsigned> responses;
    std::vector<u32> masks;
};

bool better628(const Solution628& a, const Solution628& b) {
    if (a.size != b.size) return a.size < b.size;
    if (a.rank8 != b.rank8) return a.rank8 > b.rank8;
    return a.masks < b.masks;
}

}  // namespace

int main() {
    try {
        const Geometry geometry = build_geometry(100, 628);
        Fnv body_ledger;
        for (u32 body : kBodies628) {
            require(std::popcount(body) == 9, "body rank changed");
            body_ledger.add(body);
        }
        require(body_ledger.state == UINT64_C(0x2e7f9ddccc403b9d),
                "body ledger changed");
        std::array<Class628, 16> classes{};
        const Rank628 rank8 = audit_rank628(geometry, 8, classes);
        const Rank628 rank9 = audit_rank628(geometry, 9, classes);

        Solution628 best;
        const unsigned subset_limit = 1u << 15;
        for (unsigned subset = 1; subset < subset_limit; ++subset) {
            Solution628 candidate;
            candidate.size = std::popcount(subset);
            if (candidate.size > best.size) continue;
            unsigned joined = 0;
            bool valid = true;
            for (unsigned response = 1; response < 16; ++response) {
                if (((subset >> (response - 1)) & 1u) == 0) continue;
                const Class628& c = classes[response];
                if (c.count8 + c.count9 == 0) {
                    valid = false;
                    break;
                }
                const bool choose8 = c.count8 > 0;
                candidate.rank8 += choose8;
                candidate.responses.push_back(response);
                candidate.masks.push_back(choose8 ? c.least8 : c.least9);
                joined |= response;
            }
            if (!valid || joined != 15) continue;
            std::sort(candidate.masks.begin(), candidate.masks.end());
            if (better628(candidate, best)) best = candidate;
        }
        require(best.size < 100, "no response cover");

        std::cout << "R628_MIXED_RESPONSE_DETACHED_V1\n"
                  << "PAIR 100,628 GRID " << geometry.grid << " CELLS "
                  << geometry.cells << " FAILURE_CLASSES_RANK_LE9 "
                  << geometry.classes.size() << " BODIES_FNV " << std::hex
                  << body_ledger.state << std::dec << '\n';
        for (std::size_t i = 0; i < kBodies628.size(); ++i)
            std::cout << "BODY " << i << ' ' << hex8(kBodies628[i])
                      << " LABELS {" << labels(kBodies628[i]) << "}\n";
        for (const Rank628* audit : {&rank8, &rank9})
            std::cout << "RANK " << audit->rank << " RAW_CANDIDATES "
                      << audit->raw << " ACTIVE_CANDIDATES " << audit->active
                      << " RESPONSE_FNV " << std::hex << audit->ledger.state
                      << std::dec << '\n';
        Fnv class_ledger;
        unsigned class_count = 0;
        for (unsigned response = 1; response < 16; ++response) {
            const Class628& c = classes[response];
            if (c.count8 + c.count9 == 0) continue;
            ++class_count;
            class_ledger.add(response);
            class_ledger.add(c.count8);
            class_ledger.add(c.count9);
            class_ledger.add(c.least8);
            class_ledger.add(c.least9);
            std::cout << "CLASS RESPONSE " << std::hex << response << std::dec
                      << " COVER " << std::popcount(response)
                      << " COUNT8 " << c.count8 << " LEAST8 "
                      << hex8(c.least8) << " COUNT9 " << c.count9
                      << " LEAST9 " << hex8(c.least9) << '\n';
        }
        std::cout << "CLASSES " << class_count << " CLASS_FNV " << std::hex
                  << class_ledger.state << std::dec << '\n';
        Fnv witness_ledger;
        for (std::size_t i = 0; i < best.responses.size(); ++i) {
            const unsigned response = best.responses[i];
            const Class628& c = classes[response];
            const bool choose8 = c.count8 > 0;
            const u32 mask = choose8 ? c.least8 : c.least9;
            const Margin value = margin(geometry, mask);
            require(value.ticks >= 0, "selected witness inactive");
            witness_ledger.add(mask);
            witness_ledger.add(response);
            witness_ledger.add(static_cast<u64>(value.mass));
            add_i128(witness_ledger, value.ticks);
            std::cout << "WITNESS " << hex8(mask) << " RANK "
                      << std::popcount(mask) << " RESPONSE " << std::hex
                      << response << std::dec << " COVER "
                      << std::popcount(response) << " MASS " << value.mass
                      << " MARGIN_TICKS63 " << decimal(value.ticks) << '\n';
        }
        std::cout << "MINIMUM_ADDITIONS " << best.size
                  << " LOWER exhaustive-all-32767-response-class-subsets"
                  << " PREFERRED_RANK8 " << best.rank8 << " WITNESS_FNV "
                  << std::hex << witness_ledger.state << std::dec << '\n'
                  << "SCOPE IMPORT_FREE_LITERAL_WALL_FIXED_PAIR_LABELLED_"
                     "RANK8_OR_RANK9_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_MIXED_RESPONSE_MINIMUM\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "R628_RESPONSE_ERROR " << error.what() << '\n';
        return 1;
    }
}

// Detached completeness audit for the response-class universe used by the
// (520,663) signature surgery.  Endpoint activity is reconstructed from the
// literal safe arcs of the two actual speeds on their common grid.  It does
// not call the primitive-pair/cocycle activity builder used by the primary
// atlas.  The rank-eight superset accumulation and body-complement enumerator
// are also local implementations.

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-function"
#pragma clang diagnostic ignored "-Wunneeded-internal-declaration"
#endif
#define CARRIER_CEGAR_LIBRARY_ONLY
#include "04-computation/lrc14_three_round_learned_carrier_thm4266/carrier_cegar_descent.cpp"
#undef CARRIER_CEGAR_LIBRARY_ONLY
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include <array>
#include <bit>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace {

constexpr int kQ = 520;
constexpr int kR = 663;
constexpr u64 kRank8Count = UINT64_C(5852925);
constexpr u64 kBodyComplements = UINT64_C(203490);
constexpr std::size_t kObligationCount = 53;
constexpr u64 kObligationFnv = UINT64_C(0x1e5d6dfe0b676151);
constexpr u64 kExpectedZetaOperations = UINT64_C(152170690);
constexpr u64 kExpectedActive = UINT64_C(2879147);
constexpr u64 kExpectedActiveFnv = UINT64_C(0x6dbc8f16ef5c3ff7);
constexpr u64 kExpectedCandidateCount = UINT64_C(3490314);
constexpr u64 kExpectedActiveNonempty = UINT64_C(1253628);
constexpr u64 kExpectedActiveIncidences = UINT64_C(2720981);
constexpr u64 kExpectedClassCount = UINT64_C(7124);
constexpr u64 kExpectedNonemptyClassCount = UINT64_C(7123);
constexpr u64 kExpectedMaximalCount = UINT64_C(810);
constexpr u64 kExpectedClassFnv = UINT64_C(0x2158656972d58de7);
constexpr u64 kExpectedActiveResponseFnv = UINT64_C(0xa8eb70069447e610);

std::array<std::array<u64, 9>, 31> detached_choose{};

void init_detached_choose() {
    for (int n = 0; n <= 30; ++n) {
        detached_choose[n][0] = 1;
        for (int k = 1; k <= 8; ++k) {
            detached_choose[n][k] =
                n == 0 ? 0 : detached_choose[n - 1][k] +
                                   detached_choose[n - 1][k - 1];
        }
    }
    require(detached_choose[30][8] == kRank8Count,
            "detached rank-eight universe changed");
}

u64 detached_rank8(u32 mask) {
    require(std::popcount(mask) == 8, "detached rank called off rank eight");
    u64 rank = 0;
    int ordinal = 1;
    for (int bit = 0; bit < 30; ++bit) {
        if ((mask & (u32{1} << bit)) == 0) continue;
        rank += detached_choose[bit][ordinal++];
    }
    require(ordinal == 9 && rank < kRank8Count,
            "detached colex rank escaped universe");
    return rank;
}

u32 detached_unrank8(u64 rank) {
    require(rank < kRank8Count, "detached unrank escaped universe");
    const u64 original = rank;
    u32 mask = 0;
    int ceiling = 29;
    for (int ordinal = 8; ordinal >= 1; --ordinal) {
        while (ceiling >= 0 && detached_choose[ceiling][ordinal] > rank)
            --ceiling;
        require(ceiling >= ordinal - 1, "detached colex unrank failed");
        mask |= u32{1} << ceiling;
        rank -= detached_choose[ceiling][ordinal];
        --ceiling;
    }
    require(rank == 0 && detached_rank8(mask) == original,
            "detached colex round trip failed");
    return mask;
}

i128 detached_lcm(i128 left, i128 right) {
    require(left > 0 && right > 0, "nonpositive detached LCM input");
    return left / gcd_i128(left, right) * right;
}

struct LiteralArc {
    i128 left = 0;
    i128 right = 0;
    i128 prefix = 0;
};

struct LiteralPair {
    i128 grid = 0;
    i128 pool_scale = 0;
    i128 total_safe = 0;
    std::vector<LiteralArc> arcs;
};

std::vector<std::pair<i128, i128>> literal_speed_arcs(int speed,
                                                       i128 grid) {
    const i128 denominator = static_cast<i128>(14) * speed;
    require(speed > 0 && grid % denominator == 0,
            "literal grid is not speed-divisible");
    const i128 unit = grid / denominator;
    std::vector<std::pair<i128, i128>> arcs;
    arcs.reserve(speed);
    for (int tooth = 0; tooth < speed; ++tooth) {
        arcs.push_back({static_cast<i128>(14 * tooth + 1) * unit,
                        static_cast<i128>(14 * tooth + 13) * unit});
    }
    return arcs;
}

LiteralPair build_literal_pair_detached(int q, int r) {
    require(q > 0 && q < r, "invalid detached literal pair");
    LiteralPair out;
    out.grid = detached_lcm(
        detached_lcm(static_cast<i128>(COMMON),
                     static_cast<i128>(14) * q),
        static_cast<i128>(14) * r);
    require(out.grid % COMMON == 0,
            "detached literal grid lost pool scale");
    out.pool_scale = out.grid / COMMON;
    const auto q_arcs = literal_speed_arcs(q, out.grid);
    const auto r_arcs = literal_speed_arcs(r, out.grid);
    std::size_t qi = 0;
    std::size_t ri = 0;
    i128 running = 0;
    while (qi < q_arcs.size() && ri < r_arcs.size()) {
        const i128 left = std::max(q_arcs[qi].first, r_arcs[ri].first);
        const i128 right = std::min(q_arcs[qi].second, r_arcs[ri].second);
        if (left < right) {
            if (!out.arcs.empty() && out.arcs.back().right == left) {
                out.arcs.back().right = right;
            } else {
                out.arcs.push_back({left, right, running});
            }
            running += right - left;
        }
        if (q_arcs[qi].second < r_arcs[ri].second) {
            ++qi;
        } else if (r_arcs[ri].second < q_arcs[qi].second) {
            ++ri;
        } else {
            ++qi;
            ++ri;
        }
    }
    require(!out.arcs.empty() && running > 0,
            "detached literal pair has no safe intersection");
    out.total_safe = running;
    return out;
}

i128 literal_prefix_detached(const LiteralPair& pair, i128 tick) {
    require(0 <= tick && tick <= pair.grid,
            "detached literal prefix outside circle");
    std::size_t low = 0;
    std::size_t high = pair.arcs.size();
    while (low < high) {
        const std::size_t middle = low + (high - low) / 2;
        if (pair.arcs[middle].right <= tick)
            low = middle + 1;
        else
            high = middle;
    }
    if (low == pair.arcs.size()) return pair.total_safe;
    const LiteralArc& arc = pair.arcs[low];
    return arc.prefix +
           (tick > arc.left ? std::min(tick, arc.right) - arc.left : 0);
}

struct LostBody {
    u32 body = 0;
    unsigned deleted_response = 0;
};

std::vector<LostBody> read_lost_bodies(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open lost-body CSV");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "ordinal,body_hex,deleted_response_hex",
            "lost-body header changed");
    std::vector<LostBody> bodies;
    FnvLocal ledger;
    while (std::getline(input, line)) {
        require(!line.empty(), "empty lost-body row");
        const std::size_t comma1 = line.find(',');
        const std::size_t comma2 = line.find(',', comma1 + 1);
        require(comma1 != std::string::npos && comma2 != std::string::npos &&
                    line.find(',', comma2 + 1) == std::string::npos,
                "malformed lost-body row");
        const u64 ordinal = std::stoull(line.substr(0, comma1));
        const u32 body = static_cast<u32>(
            std::stoul(line.substr(comma1 + 1, comma2 - comma1 - 1),
                       nullptr, 16));
        const unsigned response = static_cast<unsigned>(
            std::stoul(line.substr(comma2 + 1), nullptr, 16));
        require(ordinal == bodies.size() && std::popcount(body) == 9 &&
                    body < (u32{1} << 30) && response > 0 && response < 32 &&
                    (bodies.empty() || bodies.back().body < body),
                "lost-body identity/order changed");
        ledger.add(body);
        ledger.add(response);
        bodies.push_back({body, response});
    }
    require(bodies.size() == kObligationCount &&
                ledger.state == kObligationFnv,
            "lost-body count/FNV changed");
    return bodies;
}

template <class Callback>
void choose_positions(const std::vector<unsigned char>& positions,
                      int start, int remaining, u32 mask,
                      Callback& callback) {
    if (remaining == 0) {
        callback(mask);
        return;
    }
    for (int index = start;
         index <= static_cast<int>(positions.size()) - remaining; ++index) {
        choose_positions(positions, index + 1, remaining - 1,
                         mask | (u32{1} << positions[index]), callback);
    }
}

template <class Callback>
u64 enumerate_rank8_supersets(u32 atom, Callback callback) {
    const int arity = std::popcount(atom);
    require(arity <= 8, "detached atom exceeds rank eight");
    std::vector<unsigned char> complement;
    complement.reserve(30 - arity);
    for (int bit = 0; bit < 30; ++bit)
        if ((atom & (u32{1} << bit)) == 0)
            complement.push_back(static_cast<unsigned char>(bit));
    u64 emitted = 0;
    auto checked = [&](u32 extra) {
        const u32 repair = atom | extra;
        require(std::popcount(repair) == 8 && (repair & atom) == atom,
                "detached superset enumeration failed");
        callback(repair, detached_rank8(repair));
        ++emitted;
    };
    choose_positions(complement, 0, 8 - arity, 0, checked);
    return emitted;
}

template <class Callback>
u64 enumerate_disjoint_rank8(u32 body, Callback callback) {
    require(std::popcount(body) == 9, "detached body is not rank nine");
    std::vector<unsigned char> complement;
    complement.reserve(21);
    for (int bit = 0; bit < 30; ++bit)
        if ((body & (u32{1} << bit)) == 0)
            complement.push_back(static_cast<unsigned char>(bit));
    require(complement.size() == 21, "detached body complement changed");
    u64 emitted = 0;
    auto checked = [&](u32 repair) {
        require(std::popcount(repair) == 8 && (repair & body) == 0,
                "detached disjoint enumeration failed");
        callback(repair, detached_rank8(repair));
        ++emitted;
    };
    choose_positions(complement, 0, 8, 0, checked);
    require(emitted == kBodyComplements,
            "detached disjoint-repair count changed");
    return emitted;
}

struct ResponseClass {
    u64 multiplicity = 0;
    u32 least_mask = 0;
    unsigned cover = 0;
    bool maximal = false;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 4,
                "usage: literal-response-audit LOST_CSV CLASS_TSV ACTIVE_TSV");
        init_detached_choose();
        const std::vector<LostBody> obligations = read_lost_bodies(argv[1]);
        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "fixed-pool cell count changed");
        const LiteralPair pair = build_literal_pair_detached(kQ, kR);
        require(pair.grid == COMMON && pair.pool_scale == 1,
                "target literal grid unexpectedly refines fixed-pool grid");

        // Literal cell integration is the independent endpoint-activity path.
        std::map<u32, i128> atom_mass;
        i128 all_cell_mass = 0;
        u64 eligible_cells = 0;
        for (const Cell& cell : cells) {
            const i128 left = static_cast<i128>(cell.left) * pair.pool_scale;
            const i128 right = static_cast<i128>(cell.right) * pair.pool_scale;
            const i128 mass = literal_prefix_detached(pair, right) -
                              literal_prefix_detached(pair, left);
            require(mass >= 0, "negative detached literal cell mass");
            all_cell_mass += mass;
            if (std::popcount(cell.failed_pool) <= 8) {
                atom_mass[cell.failed_pool] += mass;
                ++eligible_cells;
            }
        }
        require(all_cell_mass == pair.total_safe,
                "literal cell masses do not telescope");

        // Accumulate every literal atom into all rank-eight supersets using a
        // single, local colex array (the primary uses an eight-way cocycle
        // accumulator).  Zero-mass atom keys are deliberately retained.
        std::vector<i128> literal_mass(kRank8Count, 0);
        u64 zeta_operations = 0;
        for (const auto& [atom, mass] : atom_mass) {
            zeta_operations += enumerate_rank8_supersets(
                atom, [&](u32, u64 rank) { literal_mass[rank] += mass; });
        }
        require(zeta_operations == kExpectedZetaOperations,
                "detached literal zeta operation count changed");

        std::vector<unsigned char> active(kRank8Count, 0);
        u64 active_count = 0;
        FnvLocal active_ledger;
        for (u64 rank = 0; rank < kRank8Count; ++rank) {
            const bool is_active = static_cast<i128>(63) * literal_mass[rank] >=
                                   static_cast<i128>(4) * pair.grid;
            active[rank] = static_cast<unsigned char>(is_active);
            if (is_active) {
                ++active_count;
                active_ledger.add(detached_unrank8(rank));
            }
        }
        require(active_count == kExpectedActive &&
                    active_ledger.state == kExpectedActiveFnv,
                "detached literal active universe disagrees with primary");

        // Build the union of all rank-eight masks disjoint from at least one
        // obligation, and independently attach the complete 53-bit response.
        std::vector<u64> response_by_rank(kRank8Count, 0);
        u64 complement_checks = 0;
        for (std::size_t obligation = 0; obligation < obligations.size();
             ++obligation) {
            complement_checks += enumerate_disjoint_rank8(
                obligations[obligation].body, [&](u32, u64 rank) {
                    response_by_rank[rank] |= UINT64_C(1) << obligation;
                });
        }
        require(complement_checks == obligations.size() * kBodyComplements,
                "detached complement accounting changed");
        u64 candidate_count = 0;
        for (u64 pattern : response_by_rank) candidate_count += pattern != 0;
        require(candidate_count == kExpectedCandidateCount,
                "detached unique candidate-union count changed");

        std::map<u64, ResponseClass> classes;
        u64 active_incidences = 0;
        u64 active_nonempty = 0;
        u64 response_union = 0;
        FnvLocal active_response_ledger;
        for (u64 rank = 0; rank < kRank8Count; ++rank) {
            if (!active[rank]) continue;
            const u32 repair = detached_unrank8(rank);
            const u64 pattern = response_by_rank[rank];
            response_union |= pattern;
            active_incidences += std::popcount(pattern);
            active_nonempty += pattern != 0;
            active_response_ledger.add(repair);
            active_response_ledger.add(pattern);
            ResponseClass& response_class = classes[pattern];
            ++response_class.multiplicity;
            if (response_class.multiplicity == 1 ||
                repair < response_class.least_mask)
                response_class.least_mask = repair;
        }
        const u64 full_response = (UINT64_C(1) << obligations.size()) - 1;
        require(response_union == full_response &&
                    active_nonempty == kExpectedActiveNonempty &&
                    active_incidences == kExpectedActiveIncidences &&
                    active_response_ledger.state == kExpectedActiveResponseFnv &&
                    classes.size() == kExpectedClassCount &&
                    classes.contains(0) &&
                    classes.size() - 1 == kExpectedNonemptyClassCount,
                "detached active-response quotient disagrees with primary");

        std::vector<u64> descending;
        descending.reserve(classes.size() - 1);
        for (auto& [pattern, response_class] : classes) {
            response_class.cover = std::popcount(pattern);
            if (pattern != 0) descending.push_back(pattern);
        }
        std::sort(descending.begin(), descending.end(),
                  [](u64 left, u64 right) {
                      const unsigned left_cover = std::popcount(left);
                      const unsigned right_cover = std::popcount(right);
                      if (left_cover != right_cover)
                          return left_cover > right_cover;
                      return left < right;
                  });
        std::vector<u64> maximal;
        for (u64 pattern : descending) {
            bool dominated = false;
            for (u64 prior : maximal) {
                if ((pattern & ~prior) == 0) {
                    dominated = true;
                    break;
                }
            }
            if (!dominated) maximal.push_back(pattern);
        }
        require(maximal.size() == kExpectedMaximalCount,
                "detached maximal response count disagrees with primary");
        for (u64 pattern : maximal) classes.at(pattern).maximal = true;

        std::ofstream class_output(argv[2]);
        require(static_cast<bool>(class_output),
                "cannot create detached class TSV");
        class_output << "class\tpattern_hex\tleast_mask_hex\tmultiplicity"
                        "\tcover\tmaximal\n";
        std::map<u64, u64> class_id;
        FnvLocal class_ledger;
        u64 class_index = 0;
        for (const auto& [pattern, response_class] : classes) {
            class_id.emplace(pattern, class_index);
            class_output << class_index << '\t' << std::hex << std::setw(14)
                         << std::setfill('0') << pattern << '\t'
                         << std::setw(8) << response_class.least_mask << std::dec
                         << std::setfill(' ') << '\t'
                         << response_class.multiplicity << '\t'
                         << response_class.cover << '\t'
                         << response_class.maximal << '\n';
            class_ledger.add(pattern);
            class_ledger.add(response_class.least_mask);
            class_ledger.add(response_class.multiplicity);
            class_ledger.add(response_class.cover);
            class_ledger.add(response_class.maximal);
            ++class_index;
        }
        require(class_output.good() && class_ledger.state == kExpectedClassFnv,
                "detached class TSV/FNV disagrees with primary");

        std::ofstream active_output(argv[3]);
        require(static_cast<bool>(active_output),
                "cannot create detached active-mask TSV");
        active_output << "colex_rank\tmask_hex\tclass\tpattern_hex\n";
        for (u64 rank = 0; rank < kRank8Count; ++rank) {
            if (!active[rank]) continue;
            const u32 repair = detached_unrank8(rank);
            const u64 pattern = response_by_rank[rank];
            active_output << rank << '\t' << std::hex << std::setw(8)
                          << std::setfill('0') << repair << std::dec
                          << std::setfill(' ') << '\t' << class_id.at(pattern)
                          << '\t' << std::hex << std::setw(14)
                          << std::setfill('0') << pattern << std::dec
                          << std::setfill(' ') << '\n';
        }
        require(active_output.good(), "failed writing detached active TSV");

        std::cout << "THM4281_LITERAL_RESPONSE_UNIVERSE_AUDIT_V1\n"
                  << "TARGET " << kQ << ',' << kR << " GRID "
                  << decimal(pair.grid) << " POOL_SCALE "
                  << decimal(pair.pool_scale) << " SAFE_TICKS "
                  << decimal(pair.total_safe) << " LITERAL_ARCS "
                  << pair.arcs.size() << '\n'
                  << "POOL_CELLS " << cells.size() << " ELIGIBLE_CELLS "
                  << eligible_cells << " ATOM_CLASSES " << atom_mass.size()
                  << " ZETA_OPERATIONS " << zeta_operations << '\n'
                  << "RANK8_UNIVERSE " << kRank8Count << " ACTIVE "
                  << active_count << " ACTIVE_FNV " << std::hex
                  << active_ledger.state << std::dec << '\n'
                  << "OBLIGATIONS " << obligations.size()
                  << " COMPLEMENT_CHECKS " << complement_checks
                  << " UNIQUE_RESPONSE_CANDIDATES " << candidate_count
                  << " ACTIVE_NONEMPTY " << active_nonempty
                  << " ACTIVE_INCIDENCES " << active_incidences << '\n'
                  << "RESPONSE_CLASSES_ALL " << classes.size()
                  << " NONEMPTY " << classes.size() - 1 << " MAXIMAL "
                  << maximal.size() << " CLASS_FNV " << std::hex
                  << class_ledger.state << " ACTIVE_RESPONSE_FNV "
                  << active_response_ledger.state << std::dec << '\n'
                  << "RESPONSE_UNION " << std::hex << response_union
                  << std::dec << " COVERED " << std::popcount(response_union)
                  << " OF " << obligations.size() << '\n'
                  << "VERDICT PASS DETACHED_LITERAL_COMPLETE_RESPONSE_UNIVERSE\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "LITERAL_RESPONSE_AUDIT_ERROR " << error.what() << '\n';
        return 1;
    }
}

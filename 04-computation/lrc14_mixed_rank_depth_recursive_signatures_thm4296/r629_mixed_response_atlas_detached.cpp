// Scratch-only detached rank-8/rank-9 response quotient for the first
// failure layer after the (100,630) repair.  The included engine reconstructs
// the fixed-pair walls literally and imports no maintained project source.

#define main r632_detached_hostile_main
#include "r632_detached_hostile_survivor.cpp"
#undef main

#include <unordered_map>

namespace {

constexpr u64 kRank9Count629 = UINT64_C(14307150);

struct Class629 {
    u64 count8 = 0;
    u64 count9 = 0;
    u32 least8 = 0;
    u32 least9 = 0;
};

struct RankAudit629 {
    unsigned rank = 0;
    u64 raw_incidences = 0;
    u64 raw_candidates = 0;
    u64 active_candidates = 0;
    u64 active_incidences = 0;
    u64 maximum_cover = 0;
    u32 least_maximum = 0;
    u64 full_responders = 0;
    u32 least_full = 0;
    Fnv response_ledger;
};

std::array<std::array<u64, 10>, 31> choose629{};

void init_choose629() {
    for (unsigned n = 0; n <= 30; ++n) {
        choose629[n][0] = 1;
        for (unsigned k = 1; k <= 9; ++k) {
            if (k > n) choose629[n][k] = 0;
            else if (k == n) choose629[n][k] = 1;
            else choose629[n][k] =
                choose629[n - 1][k - 1] + choose629[n - 1][k];
        }
    }
    require(choose629[30][8] == kRank8Count &&
                choose629[30][9] == kRank9Count629,
            "combination universe changed");
}

u64 colex629(u32 mask, unsigned rank) {
    require(std::popcount(mask) == rank, "colex rank changed");
    u64 answer = 0;
    unsigned ordinal = 1;
    for (unsigned bit = 0; bit < 30; ++bit)
        if ((mask >> bit) & 1u) answer += choose629[bit][ordinal++];
    require(ordinal == rank + 1 && answer < choose629[30][rank],
            "colex escaped");
    return answer;
}

template <class Callback>
void choose_rank629(const std::vector<unsigned char>& positions,
                    std::size_t start, unsigned chosen, unsigned target,
                    u32 mask, Callback& callback) {
    if (chosen == target) {
        callback(mask);
        return;
    }
    const std::size_t needed = target - chosen;
    for (std::size_t index = start;
         index + needed <= positions.size(); ++index)
        choose_rank629(positions, index + 1, chosen + 1, target,
                       mask | (u32{1} << positions[index]), callback);
}

template <class Callback>
u64 enumerate_disjoint629(u32 body, unsigned rank, Callback callback) {
    std::vector<unsigned char> positions;
    for (unsigned bit = 0; bit < 30; ++bit)
        if ((body & (u32{1} << bit)) == 0)
            positions.push_back(static_cast<unsigned char>(bit));
    require(positions.size() == 21, "body complement size changed");
    u64 count = 0;
    auto counted = [&](u32 mask) {
        require(std::popcount(mask) == rank && (mask & body) == 0,
                "disjoint enumeration escaped");
        callback(mask);
        ++count;
    };
    choose_rank629(positions, 0, 0, rank, 0, counted);
    require(count == choose629[21][rank], "disjoint count changed");
    return count;
}

std::vector<u32> read_r629_failures(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open r629 failure ledger");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,body_hex",
            "r629 failure header changed");
    std::vector<u32> bodies;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        int q = 0;
        int r = 0;
        std::string token;
        row >> q >> r >> token;
        require(row && q == 100 && r == 629,
                "malformed/non-r629 failure row");
        const u32 body = static_cast<u32>(std::stoul(token, nullptr, 16));
        require(std::popcount(body) == 9, "failure body rank changed");
        bodies.push_back(body);
    }
    require(bodies.size() == 28 &&
                std::is_sorted(bodies.begin(), bodies.end()) &&
                std::adjacent_find(bodies.begin(), bodies.end()) == bodies.end(),
            "r629 failure set changed");
    return bodies;
}

RankAudit629 audit_rank629(const Geometry& geometry,
                           const std::vector<u32>& failures,
                           unsigned rank,
                           std::map<u32, Class629>& classes) {
    std::vector<u32> raw(choose629[30][rank]);
    RankAudit629 audit;
    audit.rank = rank;
    for (std::size_t index = 0; index < failures.size(); ++index)
        audit.raw_incidences += enumerate_disjoint629(
            failures[index], rank, [&](u32 mask) {
                raw[colex629(mask, rank)] |= u32{1} << index;
            });

    const u32 limit = u32{1} << 30;
    u64 ordinal = 0;
    for (u32 mask = (u32{1} << rank) - 1; mask < limit;
         mask = next_combination(mask), ++ordinal) {
        require(ordinal == colex629(mask, rank),
                "numeric/colex order changed");
        const u32 response = raw[ordinal];
        if (response == 0) continue;
        ++audit.raw_candidates;
        const Margin value = margin(geometry, mask);
        if (value.ticks < 0) continue;
        ++audit.active_candidates;
        const unsigned cover = std::popcount(response);
        audit.active_incidences += cover;
        audit.response_ledger.add(mask);
        audit.response_ledger.add(response);
        audit.response_ledger.add(static_cast<u64>(value.mass));
        add_i128(audit.response_ledger, value.ticks);
        Class629& c = classes[response];
        if (rank == 8) {
            ++c.count8;
            if (c.count8 == 1 || mask < c.least8) c.least8 = mask;
        } else {
            ++c.count9;
            if (c.count9 == 1 || mask < c.least9) c.least9 = mask;
        }
        if (cover > audit.maximum_cover ||
            (cover == audit.maximum_cover &&
             (audit.least_maximum == 0 || mask < audit.least_maximum))) {
            audit.maximum_cover = cover;
            audit.least_maximum = mask;
        }
        if (cover == failures.size()) {
            ++audit.full_responders;
            if (audit.full_responders == 1 || mask < audit.least_full)
                audit.least_full = mask;
        }
    }
    require(ordinal == choose629[30][rank], "rank universe count changed");
    return audit;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 3, "usage: r629_atlas FAILURES MIXED_TSV");
        init_choose629();
        const std::vector<u32> failures = read_r629_failures(argv[1]);
        const Geometry geometry = build_geometry(100, 629);
        Fnv failure_ledger;
        for (u32 body : failures) failure_ledger.add(body);

        std::map<u32, Class629> classes;
        const RankAudit629 rank8 =
            audit_rank629(geometry, failures, 8, classes);
        const RankAudit629 rank9 =
            audit_rank629(geometry, failures, 9, classes);

        std::ofstream output(argv[2]);
        require(static_cast<bool>(output), "cannot create mixed atlas");
        output << "a_hex\tb_hex\tcover\tcount\tleast_mask\tleast_rank"
                  "\tcount8\tcount9\tmass\tmargin_ticks63\n";
        Fnv class_ledger;
        u64 classes_with8 = 0;
        u64 classes_with9 = 0;
        u64 prefer8 = 0;
        u32 covered = 0;
        for (const auto& [response, c] : classes) {
            require(response != 0 && c.count8 + c.count9 > 0,
                    "empty response class");
            const bool choose8 = c.count8 > 0;
            const u32 least = choose8 ? c.least8 : c.least9;
            const unsigned least_rank = choose8 ? 8 : 9;
            const Margin value = margin(geometry, least);
            require(value.ticks >= 0, "representative inactive");
            classes_with8 += c.count8 > 0;
            classes_with9 += c.count9 > 0;
            prefer8 += choose8;
            covered |= response;
            class_ledger.add(response);
            class_ledger.add(c.count8);
            class_ledger.add(c.count9);
            class_ledger.add(least);
            class_ledger.add(least_rank);
            class_ledger.add(static_cast<u64>(value.mass));
            add_i128(class_ledger, value.ticks);
            output << std::hex << std::setw(16) << std::setfill('0')
                   << static_cast<u64>(response) << '\t'
                   << std::setw(16) << UINT64_C(0) << std::dec
                   << std::setfill(' ') << '\t' << std::popcount(response)
                   << '\t' << c.count8 + c.count9 << '\t' << hex8(least)
                   << '\t' << least_rank << '\t' << c.count8 << '\t'
                   << c.count9 << '\t' << value.mass << '\t'
                   << decimal(value.ticks) << '\n';
        }
        require(output.good(), "mixed atlas write failed");
        require(covered == ((u32{1} << failures.size()) - 1),
                "some failure lacks any mixed response");

        std::cout << "R629_MIXED_RESPONSE_ATLAS_DETACHED_V1\n"
                  << "PAIR 100,629 GRID " << geometry.grid << " CELLS "
                  << geometry.cells << " PAIR_TICKS " << geometry.pair_ticks
                  << " FAILURE_CLASSES_RANK_LE9 " << geometry.classes.size()
                  << '\n'
                  << "FAILURES " << failures.size() << " FNV " << std::hex
                  << failure_ledger.state << std::dec << '\n';
        for (std::size_t i = 0; i < failures.size(); ++i)
            std::cout << "BODY " << i << ' ' << hex8(failures[i])
                      << " LABELS {" << labels(failures[i]) << "}\n";
        for (const RankAudit629* audit : {&rank8, &rank9})
            std::cout << "RANK " << audit->rank << " RAW_INCIDENCES "
                      << audit->raw_incidences << " RAW_CANDIDATES "
                      << audit->raw_candidates << " ACTIVE_CANDIDATES "
                      << audit->active_candidates << " ACTIVE_INCIDENCES "
                      << audit->active_incidences << " RESPONSE_FNV "
                      << std::hex << audit->response_ledger.state << std::dec
                      << " MAX_COVER " << audit->maximum_cover
                      << " LEAST_MAX " << hex8(audit->least_maximum)
                      << " FULL_RESPONDERS " << audit->full_responders
                      << " LEAST_FULL " << hex8(audit->least_full) << '\n';
        std::cout << "MIXED_CLASSES " << classes.size()
                  << " CLASSES_WITH_RANK8 " << classes_with8
                  << " CLASSES_WITH_RANK9 " << classes_with9
                  << " REPRESENTATIVES_PREFER_RANK8 " << prefer8
                  << " CLASS_FNV " << std::hex << class_ledger.state
                  << std::dec << '\n'
                  << "SCOPE IMPORT_FREE_LITERAL_WALL_FIXED_PAIR_LABELLED_"
                     "RANK8_OR_RANK9_RESPONSE_QUOTIENT_NO_PHYSICAL_ENTRY_"
                     "NO_LRC14\n"
                  << "VERDICT PASS EXACT_MIXED_RESPONSE_QUOTIENT\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "R629_ATLAS_ERROR " << error.what() << '\n';
        return 1;
    }
}

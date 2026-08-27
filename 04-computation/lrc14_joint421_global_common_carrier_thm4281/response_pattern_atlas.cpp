// Private THM-4276 continuation: quotient every rank-eight repair by its
// labelled response to all failed bodies at (256,670) and (384,670).
//
// The expensive Cartesian product is avoided exactly: a failed nine-body has
// a 21-point complement, hence precisely C(21,8)=203490 disjoint repairs.
// We enumerate those complements, retain the endpoint-activation coordinate,
// and only materialize response vectors for repairs with nonempty response.

#define CASCADE_LIBRARY_ONLY
#include "04-computation/lrc14_three_round_learned_carrier_thm4266/cascade_pair_exhaustive_primary.cpp"
#undef CASCADE_LIBRARY_ONLY

#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>

namespace {

constexpr std::size_t EXPECTED_OBLIGATIONS_A = 1659;
constexpr std::size_t EXPECTED_OBLIGATIONS_B = 46;
constexpr std::size_t EXPECTED_OBLIGATIONS =
    EXPECTED_OBLIGATIONS_A + EXPECTED_OBLIGATIONS_B;
constexpr std::size_t RESPONSE_WORDS =
    (EXPECTED_OBLIGATIONS + 63) / 64;
constexpr u64 DISJOINT_REPAIRS_PER_BODY = UINT64_C(203490);

struct ObligationRow {
    int q = 0;
    int r = 0;
    u64 active_carrier = 0;
    u64 stated_failure_fnv = 0;
    std::vector<u32> bodies;
};

std::vector<std::string> split(const std::string& value, char delimiter) {
    std::vector<std::string> pieces;
    std::stringstream stream(value);
    std::string piece;
    while (std::getline(stream, piece, delimiter)) pieces.push_back(piece);
    return pieces;
}

ObligationRow parse_obligation_line(const std::string& line) {
    std::istringstream in(line);
    std::string tag;
    std::string pair_tag;
    std::string pair;
    std::string active_tag;
    std::string failures_tag;
    std::string fnv_tag;
    std::string fnv_hex;
    std::string masks_tag;
    std::string masks_hex;
    u64 stated_failures = 0;
    ObligationRow row;
    in >> tag >> pair_tag >> pair >> active_tag >> row.active_carrier
       >> failures_tag >> stated_failures >> fnv_tag >> fnv_hex
       >> masks_tag >> masks_hex;
    require(in && tag == "OBLIGATION_ROW" && pair_tag == "PAIR" &&
                active_tag == "ACTIVE" && failures_tag == "FAILURES" &&
                fnv_tag == "FAILURE_FNV" && masks_tag == "MASKS_HEX",
            "malformed obligation row");
    const std::vector<std::string> qr = split(pair, ',');
    require(qr.size() == 2, "malformed obligation pair");
    row.q = std::stoi(qr[0]);
    row.r = std::stoi(qr[1]);
    row.stated_failure_fnv = std::stoull(fnv_hex, nullptr, 16);
    for (const std::string& token : split(masks_hex, ',')) {
        require(!token.empty(), "empty failed-body token");
        row.bodies.push_back(static_cast<u32>(std::stoul(token, nullptr, 16)));
    }
    require(row.bodies.size() == stated_failures,
            "failed-body count disagrees with transcript");
    require(std::is_sorted(row.bodies.begin(), row.bodies.end()),
            "failed bodies are not ordered");
    FnvLocal ledger;
    for (u32 body : row.bodies) {
        require(std::popcount(body) == 9 && body < (u32{1} << 30),
                "invalid failed body");
        ledger.add(body);
    }
    require(ledger.state == row.stated_failure_fnv,
            "failed-body FNV disagrees with transcript");
    return row;
}

std::vector<ObligationRow> read_obligations(const std::string& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open obligation transcript");
    std::vector<ObligationRow> rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("OBLIGATION_ROW ", 0) == 0) {
            rows.push_back(parse_obligation_line(line));
        }
    }
    require(rows.size() == 2, "obligation row count changed");
    require(rows[0].q == 256 && rows[0].r == 670 &&
                rows[0].active_carrier == 2172 &&
                rows[0].bodies.size() == EXPECTED_OBLIGATIONS_A &&
                rows[0].stated_failure_fnv == UINT64_C(0x970f004b0f9e2edb),
            "(256,670) obligation ledger changed");
    require(rows[1].q == 384 && rows[1].r == 670 &&
                rows[1].active_carrier == 3851 &&
                rows[1].bodies.size() == EXPECTED_OBLIGATIONS_B &&
                rows[1].stated_failure_fnv == UINT64_C(0x6c4fd2dab94cf6a9),
            "(384,670) obligation ledger changed");
    return rows;
}

struct ActiveUniverse {
    std::vector<unsigned char> active;
    u64 count = 0;
    u64 fnv = 0;
    i128 denominator = 0;
    u64 zeta_operations = 0;
};

ActiveUniverse build_active_universe(const std::vector<Cell>& cells,
                                     int q, int r) {
    const i64 g = std::gcd(q, r);
    const PrimitivePair primitive = build_primitive(q / g, r / g);
    const AtomData atoms = build_cocycle_atoms(cells, primitive, g);
    std::vector<std::pair<u32, i128>> atom_list(atoms.mass.begin(),
                                                atoms.mass.end());
    const unsigned hardware = std::thread::hardware_concurrency();
    const unsigned thread_count =
        std::max(1u, std::min(8u, hardware ? hardware : 1u));
    std::vector<std::vector<i128>> local(
        thread_count, std::vector<i128>(EXPECTED_REPAIRS, 0));
    std::vector<u64> local_operations(thread_count, 0);
    std::vector<std::thread> workers;
    for (unsigned thread = 0; thread < thread_count; ++thread) {
        workers.emplace_back([&, thread]() {
            const std::size_t begin = atom_list.size() * thread / thread_count;
            const std::size_t end = atom_list.size() * (thread + 1) /
                                    thread_count;
            for (std::size_t index = begin; index < end; ++index) {
                const auto [mask, value] = atom_list[index];
                add_supersets_pair(mask, 8 - std::popcount(mask), 0, 0,
                                   value, local[thread],
                                   local_operations[thread]);
            }
        });
    }
    for (std::thread& worker : workers) worker.join();
    std::vector<i128> masses(EXPECTED_REPAIRS, 0);
    for (unsigned thread = 0; thread < thread_count; ++thread) {
        for (std::size_t rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            masses[rank] += local[thread][rank];
        }
    }
    local.clear();
    local.shrink_to_fit();
    ActiveUniverse out;
    out.denominator = static_cast<i128>(primitive.grid) * g * COMMON;
    out.zeta_operations = std::accumulate(local_operations.begin(),
                                          local_operations.end(), u64{0});
    require(out.zeta_operations == UINT64_C(152170690),
            "pair zeta operation count changed");
    out.active.assign(EXPECTED_REPAIRS, 0);
    FnvLocal ledger;
    u32 repair = (u32{1} << 8) - 1;
    const u32 limit = u32{1} << 30;
    u64 rank = 0;
    while (repair != 0 && repair < limit) {
        const bool active = static_cast<i128>(63) * masses[rank] >=
                            static_cast<i128>(4) * out.denominator;
        out.active[rank] = active;
        if (active) {
            ++out.count;
            ledger.add(repair);
        }
        ++rank;
        const u32 next = next_combination(repair);
        if (next <= repair) break;
        repair = next;
    }
    require(rank == EXPECTED_REPAIRS, "repair enumeration changed");
    out.fnv = ledger.state;
    return out;
}

template <class Callback>
void enumerate_complement_rec(const std::array<unsigned char, 21>& positions,
                              int start, int chosen, u32 mask, u64 rank,
                              Callback& callback) {
    if (chosen == 8) {
        callback(mask, rank);
        return;
    }
    const int needed = 8 - chosen;
    for (int index = start; index <= 21 - needed; ++index) {
        const int bit = positions[index];
        enumerate_complement_rec(
            positions, index + 1, chosen + 1, mask | (u32{1} << bit),
            rank + choose8_local[bit][chosen + 1], callback);
    }
}

template <class Callback>
void enumerate_disjoint_repairs(u32 body, Callback callback) {
    std::array<unsigned char, 21> positions{};
    std::size_t count = 0;
    for (int bit = 0; bit < 30; ++bit) {
        if ((body & (u32{1} << bit)) == 0) positions[count++] = bit;
    }
    require(count == positions.size(), "body complement size changed");
    u64 emitted = 0;
    auto checked = [&](u32 repair, u64 rank) {
        require(rank < EXPECTED_REPAIRS && std::popcount(repair) == 8 &&
                    (repair & body) == 0,
                "complement enumeration failure");
        callback(repair, rank);
        ++emitted;
    };
    enumerate_complement_rec(positions, 0, 0, 0, 0, checked);
    require(emitted == DISJOINT_REPAIRS_PER_BODY,
            "disjoint repair count changed");
}

u32 unrank_colex8(u64 rank) {
    require(rank < EXPECTED_REPAIRS, "unrank outside repair universe");
    const u64 original = rank;
    u32 mask = 0;
    int upper = 29;
    for (int ordinal = 8; ordinal >= 1; --ordinal) {
        while (upper >= 0 && choose8_local[upper][ordinal] > rank) --upper;
        require(upper >= ordinal - 1, "colex unrank failed");
        mask |= u32{1} << upper;
        rank -= choose8_local[upper][ordinal];
        --upper;
    }
    require(rank == 0 && colex_rank8_local(mask) == original,
            "colex unrank round trip failed");
    return mask;
}

bool pattern_less(const u64* left, const u64* right) {
    for (std::size_t word = RESPONSE_WORDS; word-- > 0;) {
        if (left[word] != right[word]) return left[word] < right[word];
    }
    return false;
}

bool pattern_equal(const u64* left, const u64* right) {
    for (std::size_t word = 0; word < RESPONSE_WORDS; ++word) {
        if (left[word] != right[word]) return false;
    }
    return true;
}

u64 pattern_popcount(const u64* pattern) {
    u64 count = 0;
    for (std::size_t word = 0; word < RESPONSE_WORDS; ++word) {
        count += std::popcount(pattern[word]);
    }
    return count;
}

}  // namespace

#ifndef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc >= 3 && argc <= 6,
                "usage: response-pattern-atlas OBLIGATIONS_OUT ATLAS_TSV "
                "[NONEMPTY_MASKS] [ALLOWED_MASKS] [MASK_CLASS_TSV]");
        init_choose8_local();
        const std::vector<ObligationRow> rows = read_obligations(argv[1]);
        std::vector<u32> allowed_masks;
        if (argc >= 5) {
            std::ifstream allowed_input(argv[4]);
            require(static_cast<bool>(allowed_input),
                    "cannot open allowed-mask input");
            std::string word;
            FnvLocal allowed_ledger;
            while (allowed_input >> word) {
                const u32 mask =
                    static_cast<u32>(std::stoul(word, nullptr, 16));
                require(std::popcount(mask) == 8 &&
                            mask < (u32{1} << 30),
                        "invalid allowed mask");
                allowed_masks.push_back(mask);
                allowed_ledger.add(mask);
            }
            require(allowed_masks.size() == 212718 &&
                        allowed_ledger.state ==
                            UINT64_C(0x1952773ca107f9ea),
                    "rectangle-common allowed-mask ledger changed");
            require(std::is_sorted(allowed_masks.begin(),
                                   allowed_masks.end()) &&
                        std::adjacent_find(allowed_masks.begin(),
                                           allowed_masks.end()) ==
                            allowed_masks.end(),
                    "allowed masks not strictly ordered");
            std::cout << "ALLOWED_MASKS " << allowed_masks.size()
                      << " FNV " << std::hex << allowed_ledger.state
                      << std::dec << '\n';
        }
        std::vector<u32> bodies;
        std::vector<unsigned char> row_of_obligation;
        for (std::size_t row = 0; row < rows.size(); ++row) {
            for (u32 body : rows[row].bodies) {
                bodies.push_back(body);
                row_of_obligation.push_back(static_cast<unsigned char>(row));
            }
        }
        require(bodies.size() == EXPECTED_OBLIGATIONS,
                "total obligation count changed");

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cell count changed");
        ActiveUniverse active_a = build_active_universe(cells, 256, 670);
        std::cout << "ACTIVE_ROW 256,670 COUNT " << active_a.count
                  << " FNV " << std::hex << active_a.fnv << std::dec
                  << " ZETA_OPERATIONS " << active_a.zeta_operations << '\n';
        ActiveUniverse active_b = build_active_universe(cells, 384, 670);
        std::cout << "ACTIVE_ROW 384,670 COUNT " << active_b.count
                  << " FNV " << std::hex << active_b.fnv << std::dec
                  << " ZETA_OPERATIONS " << active_b.zeta_operations << '\n';

        const unsigned hardware = std::thread::hardware_concurrency();
        const unsigned thread_count =
            std::max(1u, std::min(8u, hardware ? hardware : 1u));
        std::vector<std::vector<u32>> covering_ranks(EXPECTED_OBLIGATIONS);
        const std::size_t mark_words = (EXPECTED_REPAIRS + 63) / 64;
        std::vector<std::vector<u64>> local_marks(
            thread_count, std::vector<u64>(mark_words, 0));
        std::vector<u64> local_incidences(thread_count, 0);
        std::vector<std::thread> workers;
        for (unsigned thread = 0; thread < thread_count; ++thread) {
            workers.emplace_back([&, thread]() {
                const std::size_t begin = bodies.size() * thread / thread_count;
                const std::size_t end = bodies.size() * (thread + 1) /
                                        thread_count;
                for (std::size_t obligation = begin; obligation < end;
                     ++obligation) {
                    std::vector<u32>& hits = covering_ranks[obligation];
                    hits.reserve(65536);
                    const std::vector<unsigned char>& active =
                        row_of_obligation[obligation] == 0
                            ? active_a.active : active_b.active;
                    enumerate_disjoint_repairs(
                        bodies[obligation], [&](u32, u64 rank) {
                            if (!active[rank]) return;
                            hits.push_back(static_cast<u32>(rank));
                            local_marks[thread][rank >> 6] |=
                                UINT64_C(1) << (rank & 63);
                        });
                    local_incidences[thread] += hits.size();
                }
            });
        }
        for (std::thread& worker : workers) worker.join();
        std::vector<u64> candidate_marks(mark_words, 0);
        for (unsigned thread = 0; thread < thread_count; ++thread) {
            for (std::size_t word = 0; word < mark_words; ++word) {
                candidate_marks[word] |= local_marks[thread][word];
            }
        }
        local_marks.clear();
        local_marks.shrink_to_fit();
        const u64 incidences = std::accumulate(local_incidences.begin(),
                                               local_incidences.end(), u64{0});
        u64 candidate_count = 0;
        for (u64 word : candidate_marks) candidate_count += std::popcount(word);
        std::cout << "OBLIGATIONS " << bodies.size()
                  << " COMPLEMENT_CHECKS "
                  << bodies.size() * DISJOINT_REPAIRS_PER_BODY
                  << " ACTIVE_INCIDENCES " << incidences
                  << " NONEMPTY_REPAIRS " << candidate_count
                  << " EMPTY_REPAIRS " << EXPECTED_REPAIRS - candidate_count
                  << '\n';

        std::vector<int32_t> candidate_by_rank(EXPECTED_REPAIRS, -1);
        std::vector<u32> rank_by_candidate;
        rank_by_candidate.reserve(candidate_count);
        std::ofstream candidate_masks;
        if (argc >= 4) {
            candidate_masks.open(argv[3]);
            require(static_cast<bool>(candidate_masks),
                    "cannot create nonempty-mask output");
        }
        FnvLocal candidate_ledger;
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            if ((candidate_marks[rank >> 6] &
                 (UINT64_C(1) << (rank & 63))) == 0) continue;
            candidate_by_rank[rank] =
                static_cast<int32_t>(rank_by_candidate.size());
            rank_by_candidate.push_back(static_cast<u32>(rank));
            const u32 mask = unrank_colex8(rank);
            candidate_ledger.add(mask);
            if (argc >= 4) {
                candidate_masks << std::hex << std::setw(8)
                                << std::setfill('0') << mask << std::dec
                                << '\n';
            }
        }
        require(rank_by_candidate.size() == candidate_count,
                "candidate indexing changed");
        require(argc < 4 || candidate_masks.good(),
                "failed writing nonempty masks");
        std::cout << "NONEMPTY_MASK_FNV " << std::hex
                  << candidate_ledger.state << std::dec << '\n';
        candidate_marks.clear();
        candidate_marks.shrink_to_fit();

        std::vector<u64> patterns(candidate_count * RESPONSE_WORDS, 0);
        workers.clear();
        for (unsigned thread = 0; thread < thread_count; ++thread) {
            workers.emplace_back([&, thread]() {
                for (std::size_t word = thread; word < RESPONSE_WORDS;
                     word += thread_count) {
                    const std::size_t begin = word * 64;
                    const std::size_t end =
                        std::min(begin + 64, EXPECTED_OBLIGATIONS);
                    for (std::size_t obligation = begin; obligation < end;
                         ++obligation) {
                        const u64 bit = UINT64_C(1) << (obligation & 63);
                        for (u32 rank : covering_ranks[obligation]) {
                            const int32_t candidate = candidate_by_rank[rank];
                            require(candidate >= 0,
                                    "cover incidence lost its candidate");
                            patterns[static_cast<std::size_t>(candidate) *
                                         RESPONSE_WORDS + word] |= bit;
                        }
                    }
                }
            });
        }
        for (std::thread& worker : workers) worker.join();
        covering_ranks.clear();
        covering_ranks.shrink_to_fit();
        candidate_by_rank.clear();
        candidate_by_rank.shrink_to_fit();

        std::vector<u32> order;
        order.reserve(allowed_masks.empty() ? candidate_count
                                            : allowed_masks.size());
        for (u32 candidate = 0; candidate < candidate_count; ++candidate) {
            const u32 mask = unrank_colex8(rank_by_candidate[candidate]);
            if (allowed_masks.empty() ||
                std::binary_search(allowed_masks.begin(), allowed_masks.end(),
                                   mask)) {
                order.push_back(candidate);
            }
        }
        require(allowed_masks.empty() || order.size() == allowed_masks.size(),
                "allowed mask missing from nonempty response universe");
        std::sort(order.begin(), order.end(), [&](u32 left, u32 right) {
            return pattern_less(&patterns[static_cast<std::size_t>(left) *
                                          RESPONSE_WORDS],
                                &patterns[static_cast<std::size_t>(right) *
                                          RESPONSE_WORDS]);
        });

        std::ofstream atlas(argv[2]);
        require(static_cast<bool>(atlas), "cannot create atlas TSV");
        atlas << "class\tleast_mask_hex\tmultiplicity\tcover\twords_hi_to_lo\n";
        FnvLocal class_ledger;
        u64 class_count = 0;
        u64 max_cover = 0;
        u64 singleton_classes = 0;
        std::array<u64, RESPONSE_WORDS> union_pattern{};
        std::vector<int32_t> class_by_candidate(candidate_count, -1);
        for (std::size_t begin = 0; begin < order.size();) {
            std::size_t end = begin + 1;
            const u64* pattern =
                &patterns[static_cast<std::size_t>(order[begin]) *
                          RESPONSE_WORDS];
            while (end < order.size() &&
                   pattern_equal(
                       pattern,
                       &patterns[static_cast<std::size_t>(order[end]) *
                                 RESPONSE_WORDS])) {
                ++end;
            }
            u32 least_rank = rank_by_candidate[order[begin]];
            for (std::size_t index = begin + 1; index < end; ++index) {
                least_rank = std::min(least_rank,
                                      rank_by_candidate[order[index]]);
            }
            const u32 least_mask = unrank_colex8(least_rank);
            const u64 multiplicity = end - begin;
            const u64 cover = pattern_popcount(pattern);
            max_cover = std::max(max_cover, cover);
            singleton_classes += cover == 1;
            class_ledger.add(least_mask);
            class_ledger.add(multiplicity);
            class_ledger.add(cover);
            for (std::size_t word = 0; word < RESPONSE_WORDS; ++word) {
                class_ledger.add(pattern[word]);
                union_pattern[word] |= pattern[word];
            }
            atlas << class_count << '\t' << std::hex << std::setw(8)
                  << std::setfill('0') << least_mask << std::dec << '\t'
                  << multiplicity << '\t' << cover << '\t';
            for (std::size_t word = RESPONSE_WORDS; word-- > 0;) {
                if (word + 1 != RESPONSE_WORDS) atlas << ',';
                atlas << std::hex << std::setw(16) << std::setfill('0')
                      << pattern[word] << std::dec;
            }
            atlas << '\n';
            for (std::size_t index = begin; index < end; ++index) {
                class_by_candidate[order[index]] =
                    static_cast<int32_t>(class_count);
            }
            ++class_count;
            begin = end;
        }
        require(atlas.good(), "failed while writing atlas TSV");
        if (argc == 6) {
            require(!allowed_masks.empty(),
                    "mask-class output requires allowed-mask input");
            std::ofstream mapping(argv[5]);
            require(static_cast<bool>(mapping),
                    "cannot create mask-class TSV");
            mapping << "mask_hex\tclass\n";
            FnvLocal mapping_ledger;
            u64 mapped = 0;
            for (u32 candidate = 0; candidate < candidate_count;
                 ++candidate) {
                const int32_t response_class =
                    class_by_candidate[candidate];
                if (response_class < 0) continue;
                const u32 mask =
                    unrank_colex8(rank_by_candidate[candidate]);
                mapping << std::hex << std::setw(8) << std::setfill('0')
                        << mask << std::dec << '\t' << response_class
                        << '\n';
                mapping_ledger.add(mask);
                mapping_ledger.add(response_class);
                ++mapped;
            }
            require(mapping.good() && mapped == allowed_masks.size(),
                    "failed writing complete mask-class map");
            std::cout << "MASK_CLASS_ROWS " << mapped << " FNV "
                      << std::hex << mapping_ledger.state << std::dec
                      << '\n';
        }
        u64 covered_obligations = 0;
        for (u64 word : union_pattern)
            covered_obligations += std::popcount(word);
        std::cout << "SELECTED_REPAIRS " << order.size()
                  << " OBLIGATIONS_COVERED " << covered_obligations
                  << " FAILURES "
                  << EXPECTED_OBLIGATIONS - covered_obligations << '\n';
        std::cout << "RESPONSE_CLASSES_NONEMPTY " << class_count
                  << " RESPONSE_CLASSES_ALL " << class_count + 1
                  << " CLASS_FNV " << std::hex << class_ledger.state
                  << std::dec << " MAX_COVER " << max_cover
                  << " SINGLETON_CLASSES " << singleton_classes << '\n';
        std::cout << "EMPTY_CLASS MULTIPLICITY "
                  << EXPECTED_REPAIRS - candidate_count << '\n';
        std::cout << "VERDICT PASS EXACT_LABELLED_PAIR_RESPONSE_QUOTIENT\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "RESPONSE_PATTERN_ATLAS_ERROR " << error.what() << '\n';
        return 1;
    }
}
#endif

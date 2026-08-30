#define CASCADE_LIBRARY_ONLY
#include "04-computation/lrc14_three_round_learned_carrier_thm4266/cascade_pair_exhaustive_primary.cpp"
#undef CASCADE_LIBRARY_ONLY

#include <bit>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <unordered_map>

namespace {

struct SigRow {
    int q = 0;
    int r = 0;
    std::array<u64, 7> words{};
};

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> out;
    std::istringstream in(line);
    std::string field;
    while (std::getline(in, field, ',')) out.push_back(field);
    return out;
}

std::vector<u32> read_deck(const std::string& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open deck");
    std::vector<u32> deck;
    std::set<u32> seen;
    std::string token;
    FnvLocal ledger;
    while (in >> token) {
        const u32 mask = static_cast<u32>(std::stoul(token, nullptr, 16));
        require(std::popcount(mask) == 8 && mask < (u32{1} << 30) &&
                    seen.insert(mask).second,
                "invalid deck mask");
        deck.push_back(mask);
        ledger.add(mask);
    }
    require(deck.size() == 421 &&
                ledger.state == UINT64_C(0x20d63dd42fe8150e),
            "deck identity changed");
    return deck;
}

std::vector<SigRow> read_signatures(const std::string& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open signature csv");
    std::string line;
    require(static_cast<bool>(std::getline(in, line)) &&
                line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "signature header changed");
    std::vector<SigRow> rows;
    int previous_q = -1;
    int previous_r = -1;
    while (std::getline(in, line)) {
        const auto f = split_csv(line);
        require(f.size() == 10, "bad signature row");
        SigRow row;
        row.q = std::stoi(f[0]);
        row.r = std::stoi(f[1]);
        unsigned count = 0;
        for (int i = 0; i < 7; ++i) {
            row.words[i] = std::stoull(f[3 + i], nullptr, 16);
            count += std::popcount(row.words[i]);
        }
        require(count == static_cast<unsigned>(std::stoul(f[2])),
                "inactive count mismatch");
        require(rows.empty() || std::pair{previous_q, previous_r} <
                                    std::pair{row.q, row.r},
                "signature rows not ordered");
        rows.push_back(row);
        previous_q = row.q;
        previous_r = row.r;
    }
    require(rows.size() == 24223, "signature universe changed");
    return rows;
}

std::array<u64, 7> singleton_signature(std::size_t index) {
    require(index < 421, "bad singleton index");
    std::array<u64, 7> out{};
    out[index / 64] |= UINT64_C(1) << (index % 64);
    return out;
}

std::vector<u32> private_bodies(const std::vector<u32>& deck,
                                std::size_t target_index,
                                u64& target_disjoint,
                                u64& retained_checks) {
    const u32 target = deck[target_index];
    std::vector<u32> out;
    target_disjoint = 0;
    retained_checks = 0;
    u64 bodies = 0;
    u32 body = (u32{1} << 9) - 1;
    while (body < (u32{1} << 30)) {
        ++bodies;
        if ((body & target) == 0) {
            ++target_disjoint;
            bool covered = false;
            for (std::size_t i = 0; i < deck.size(); ++i) {
                if (i == target_index) continue;
                ++retained_checks;
                if ((body & deck[i]) == 0) {
                    covered = true;
                    break;
                }
            }
            if (!covered) out.push_back(body);
        }
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    require(bodies == UINT64_C(14307150) &&
                target_disjoint == UINT64_C(497420),
            "body enumeration changed");
    return out;
}

struct PairAtoms {
    int q = 0;
    int r = 0;
    AtomData atoms;
    i128 denominator = 0;
};

PairAtoms build_pair_atoms(const std::vector<Cell>& cells, int q, int r) {
    PairAtoms out;
    out.q = q;
    out.r = r;
    const i64 g = std::gcd(q, r);
    const PrimitivePair primitive = build_primitive(q / g, r / g);
    out.atoms = build_cocycle_atoms(cells, primitive, g);
    out.denominator = static_cast<i128>(primitive.grid) * g * COMMON;
    return out;
}

i128 margin(const PairAtoms& pair, u32 repair) {
    i128 mass = 0;
    for (const auto& [atom, value] : pair.atoms.mass) {
        if ((atom & ~repair) == 0) mass += value;
    }
    return static_cast<i128>(63) * mass -
           static_cast<i128>(4) * pair.denominator;
}

u64 mask_fnv(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 6,
                "usage: singleton-fibre-cegar DECK SIGNATURES INDEX Q R");
        init_choose8_local();
        const std::vector<u32> deck = read_deck(argv[1]);
        const std::vector<SigRow> rows = read_signatures(argv[2]);
        const std::size_t target_index = std::stoul(argv[3]);
        const int target_q = std::stoi(argv[4]);
        const int target_r = std::stoi(argv[5]);
        require(target_index < deck.size(), "target index out of range");
        require(target_index == 396 && target_q == 366 && target_r == 644,
                "target surgery identity changed");
        const auto wanted = singleton_signature(target_index);
        std::vector<std::pair<int, int>> fibre;
        bool target_found = false;
        for (const SigRow& row : rows) {
            if (row.words != wanted) continue;
            fibre.emplace_back(row.q, row.r);
            target_found |= row.q == target_q && row.r == target_r;
        }
        require(target_found && !fibre.empty(), "target not in singleton fibre");

        u64 target_disjoint = 0;
        u64 retained_checks = 0;
        const std::vector<u32> obligations = private_bodies(
            deck, target_index, target_disjoint, retained_checks);
        require(!obligations.empty(), "deletion exposes no private body");
        u32 body_union = 0;
        FnvLocal obligation_ledger;
        for (u32 body : obligations) {
            body_union |= body;
            obligation_ledger.add(body);
        }

        std::vector<u32> one_mask_responders;
        u32 repair = (u32{1} << 8) - 1;
        while (repair < (u32{1} << 30)) {
            if ((repair & body_union) == 0) one_mask_responders.push_back(repair);
            const u32 next = next_combination(repair);
            if (next <= repair) break;
            repair = next;
        }
        require(!one_mask_responders.empty(), "no one-mask body responder");

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        std::vector<PairAtoms> fibre_atoms;
        fibre_atoms.reserve(fibre.size());
        for (const auto& [q, r] : fibre) {
            fibre_atoms.push_back(build_pair_atoms(cells, q, r));
        }

        const std::size_t words = (fibre.size() + 63) / 64;
        struct Pattern {
            std::vector<u64> bits;
            u32 least_mask = 0;
            u64 multiplicity = 0;
        };
        std::map<std::vector<u64>, Pattern> patterns;
        u64 active_incidences = 0;
        u64 equalities = 0;
        std::vector<u32> fibre_wide;
        for (u32 candidate : one_mask_responders) {
            std::vector<u64> bits(words, 0);
            std::size_t active = 0;
            for (std::size_t i = 0; i < fibre_atoms.size(); ++i) {
                const i128 m = margin(fibre_atoms[i], candidate);
                equalities += m == 0;
                if (m >= 0) {
                    bits[i / 64] |= UINT64_C(1) << (i % 64);
                    ++active;
                    ++active_incidences;
                }
            }
            if (active == fibre.size()) fibre_wide.push_back(candidate);
            Pattern& p = patterns[bits];
            ++p.multiplicity;
            if (p.least_mask == 0 || candidate < p.least_mask)
                p.least_mask = candidate;
            p.bits = std::move(bits);
        }

        std::vector<unsigned char> covered(fibre.size(), 0);
        std::vector<std::pair<u32, std::vector<u64>>> greedy;
        while (true) {
            std::size_t remaining = 0;
            for (unsigned char bit : covered) remaining += !bit;
            if (remaining == 0) break;
            std::size_t best_gain = 0;
            u32 best_mask = 0;
            std::vector<u64> best_bits;
            for (const auto& [key, p] : patterns) {
                std::size_t gain = 0;
                for (std::size_t i = 0; i < fibre.size(); ++i) {
                    gain += !covered[i] &&
                            ((p.bits[i / 64] >> (i % 64)) & 1u);
                }
                if (gain > best_gain ||
                    (gain == best_gain && gain > 0 && p.least_mask < best_mask)) {
                    best_gain = gain;
                    best_mask = p.least_mask;
                    best_bits = p.bits;
                }
            }
            if (best_gain == 0) break;
            for (std::size_t i = 0; i < fibre.size(); ++i) {
                if ((best_bits[i / 64] >> (i % 64)) & 1u) covered[i] = 1;
            }
            greedy.emplace_back(best_mask, std::move(best_bits));
        }
        std::size_t covered_count = 0;
        for (unsigned char bit : covered) covered_count += bit;

        FnvLocal fibre_ledger;
        for (const auto& [q, r] : fibre) {
            fibre_ledger.add(q);
            fibre_ledger.add(r);
        }
        FnvLocal pattern_ledger;
        for (const auto& [bits, p] : patterns) {
            pattern_ledger.add(p.least_mask);
            pattern_ledger.add(p.multiplicity);
            for (u64 word : bits) pattern_ledger.add(word);
        }
        require(fibre.size() == 36 &&
                    fibre_ledger.state == UINT64_C(0x3d92ab45b46a72c0) &&
                    deck[target_index] == UINT32_C(0x082022c9) &&
                    target_disjoint == UINT64_C(497420) &&
                    retained_checks == UINT64_C(13677671) &&
                    obligations.size() == 8 &&
                    obligation_ledger.state == UINT64_C(0xe3fd616a4aa21839) &&
                    body_union == UINT32_C(0x33dfdc06) &&
                    one_mask_responders.size() == 495 &&
                    mask_fnv(one_mask_responders) ==
                        UINT64_C(0x3c76be2ab12086ed) &&
                    active_incidences == 1193 && equalities == 0 &&
                    patterns.size() == 83 &&
                    pattern_ledger.state == UINT64_C(0xbb249b742c391810) &&
                    fibre_wide == std::vector<u32>{UINT32_C(0x042022c9),
                                                   UINT32_C(0x0c202289)} &&
                    mask_fnv(fibre_wide) == UINT64_C(0x72f8883feafa2f41) &&
                    greedy.size() == 1 &&
                    greedy.front().first == UINT32_C(0x042022c9) &&
                    covered_count == 36,
                "singleton response audit changed");
        std::cout << "THM4286_SINGLETON_FIBRE_CEGAR_V1\n"
                  << "TARGET " << target_q << ',' << target_r
                  << " INDEX " << target_index << " MASK " << std::hex
                  << std::setw(8) << std::setfill('0') << deck[target_index]
                  << std::dec << std::setfill(' ') << '\n'
                  << "FIBRE " << fibre.size() << " FNV " << std::hex
                  << fibre_ledger.state << std::dec << " ROWS";
        for (const auto& [q, r] : fibre) std::cout << ' ' << q << ',' << r;
        std::cout << '\n'
                  << "BODY_SCAN TARGET_DISJOINT " << target_disjoint
                  << " RETAINED_CHECKS " << retained_checks
                  << " OBLIGATIONS " << obligations.size() << " FNV "
                  << std::hex << obligation_ledger.state << std::dec
                  << " UNION " << std::hex << std::setw(8)
                  << std::setfill('0') << body_union << std::dec
                  << std::setfill(' ') << " UNION_SIZE "
                  << std::popcount(body_union) << '\n'
                  << "ONE_MASK_RESPONDERS " << one_mask_responders.size()
                  << " FNV " << std::hex << mask_fnv(one_mask_responders)
                  << std::dec << " ACTIVITY_TESTS "
                  << one_mask_responders.size() * fibre.size()
                  << " ACTIVE_INCIDENCES " << active_incidences
                  << " EQUALITIES " << equalities << '\n'
                  << "ACTIVITY_CLASSES " << patterns.size() << " FNV "
                  << std::hex << pattern_ledger.state << std::dec
                  << " FIBRE_WIDE " << fibre_wide.size() << " FNV "
                  << std::hex << mask_fnv(fibre_wide) << std::dec;
        if (!fibre_wide.empty()) {
            std::cout << " LEAST " << std::hex << std::setw(8)
                      << std::setfill('0') << fibre_wide.front() << std::dec
                      << std::setfill(' ');
        }
        std::cout << " MASKS";
        for (u32 mask : fibre_wide) {
            std::cout << ' ' << std::hex << std::setw(8)
                      << std::setfill('0') << mask << std::dec
                      << std::setfill(' ');
        }
        std::cout << '\n';
        for (const auto& [bits, p] : patterns) {
            std::size_t active = 0;
            for (u64 word : bits) active += std::popcount(word);
            if (active == 0) continue;
            std::cout << "CLASS LEAST " << std::hex << std::setw(8)
                      << std::setfill('0') << p.least_mask << std::dec
                      << std::setfill(' ') << " MULT " << p.multiplicity
                      << " ACTIVE " << active << " ROWS";
            for (std::size_t i = 0; i < fibre.size(); ++i) {
                if ((bits[i / 64] >> (i % 64)) & 1u)
                    std::cout << ' ' << fibre[i].first << ',' << fibre[i].second;
            }
            std::cout << '\n';
        }
        std::cout
                  << "GREEDY_CLASSES " << greedy.size() << " COVERED "
                  << covered_count << " OF " << fibre.size() << " MASKS";
        for (const auto& [mask, bits] : greedy) {
            std::cout << ' ' << std::hex << std::setw(8)
                      << std::setfill('0') << mask << std::dec
                      << std::setfill(' ');
        }
        std::cout << "\nVERDICT PASS DERIVED_SINGLETON_FIBRE_RESPONSE_AUDIT\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SINGLETON_FIBRE_CEGAR_ERROR " << error.what() << '\n';
        return 1;
    }
}

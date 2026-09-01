// Independent index-265 audit using the atom-by-atom margin path rather than
// the colex/zeta activity-intersection path used by index265_exact.cpp.

#define CASCADE_LIBRARY_ONLY
#include "04-computation/lrc14_three_round_learned_carrier_thm4266/cascade_pair_exhaustive_primary.cpp"
#undef CASCADE_LIBRARY_ONLY

#include <filesystem>
#include <fstream>
#include <set>

namespace {

constexpr std::size_t kIndex = 265;
constexpr u32 kOld = UINT32_C(0x20820649);
constexpr std::array<u32, 8> kExpectedBodies = {
    UINT32_C(0x09392104), UINT32_C(0x0d30e080),
    UINT32_C(0x0d382104), UINT32_C(0x15386080),
    UINT32_C(0x186c9080), UINT32_C(0x19786000),
    UINT32_C(0x1d489080), UINT32_C(0x1f087000)};
constexpr std::array<u32, 2> kWitnesses = {
    UINT32_C(0x22020e09), UINT32_C(0x00868489)};

using PairD = std::pair<int, int>;

struct SigRowD {
    int q = 0;
    int r = 0;
    std::array<u64, 7> words{};
};

std::vector<std::string> split_csvD(const std::string& line) {
    std::vector<std::string> out;
    std::istringstream input(line);
    std::string field;
    while (std::getline(input, field, ',')) out.push_back(field);
    return out;
}

std::vector<u32> read_deckD(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open deck");
    std::vector<u32> deck;
    std::set<u32> distinct;
    FnvLocal ledger;
    std::string token;
    while (input >> token) {
        const u32 mask = static_cast<u32>(std::stoul(token, nullptr, 16));
        require(mask < (u32{1} << 30) && std::popcount(mask) == 8 &&
                    distinct.insert(mask).second,
                "invalid deck mask");
        deck.push_back(mask);
        ledger.add(mask);
    }
    require(input.eof() && deck.size() == 421 &&
                ledger.state == UINT64_C(0x20d63dd42fe8150e),
            "deck identity changed");
    return deck;
}

std::vector<SigRowD> read_signaturesD(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open signatures");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "signature header changed");
    std::vector<SigRowD> out;
    PairD previous{-1, -1};
    while (std::getline(input, line)) {
        const auto fields = split_csvD(line);
        require(fields.size() == 10, "bad signature row");
        SigRowD row;
        row.q = std::stoi(fields[0]);
        row.r = std::stoi(fields[1]);
        unsigned weight = 0;
        for (unsigned word = 0; word < 7; ++word) {
            row.words[word] = std::stoull(fields[3 + word], nullptr, 16);
            weight += std::popcount(row.words[word]);
        }
        require(weight == std::stoul(fields[2]) &&
                    (out.empty() || previous < PairD{row.q, row.r}),
                "signature count/order changed");
        out.push_back(row);
        previous = {row.q, row.r};
    }
    require(input.eof() && out.size() == 24223,
            "signature atlas size changed");
    return out;
}

std::set<PairD> read_liveD(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open live residual");
    std::set<PairD> out;
    std::string line;
    while (std::getline(input, line)) {
        const auto fields = split_csvD(line);
        require(fields.size() == 2, "bad live row");
        require(out.insert({std::stoi(fields[0]), std::stoi(fields[1])}).second,
                "duplicate live row");
    }
    require(input.eof() && out.size() == 22647,
            "live residual size changed");
    return out;
}

std::array<u64, 7> singleton_signatureD(std::size_t index) {
    require(index < 421, "bad signature index");
    std::array<u64, 7> out{};
    out[index / 64] |= UINT64_C(1) << (index % 64);
    return out;
}

std::vector<u32> private_bodiesD(const std::vector<u32>& deck,
                                 std::size_t target_index,
                                 u64& target_disjoint,
                                 u64& retained_checks) {
    target_disjoint = retained_checks = 0;
    std::vector<u32> out;
    u32 body = (u32{1} << 9) - 1;
    u64 body_count = 0;
    while (body < (u32{1} << 30)) {
        ++body_count;
        if ((body & deck[target_index]) == 0) {
            ++target_disjoint;
            bool covered = false;
            for (std::size_t index = 0; index < deck.size(); ++index) {
                if (index == target_index) continue;
                ++retained_checks;
                if ((body & deck[index]) == 0) {
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
    require(body_count == UINT64_C(14307150) &&
                target_disjoint == UINT64_C(497420),
            "body universe changed");
    return out;
}

struct PairAtomsD {
    int q = 0;
    int r = 0;
    AtomData atoms;
    i128 denominator = 0;
};

PairAtomsD build_pair_atomsD(const std::vector<Cell>& cells, int q, int r) {
    PairAtomsD out;
    out.q = q;
    out.r = r;
    const i64 g = std::gcd(q, r);
    const PrimitivePair primitive = build_primitive(q / g, r / g);
    out.atoms = build_cocycle_atoms(cells, primitive, g);
    out.denominator = static_cast<i128>(primitive.grid) * g * COMMON;
    return out;
}

i128 marginD(const PairAtomsD& pair, u32 repair) {
    i128 mass = 0;
    for (const auto& [atom, value] : pair.atoms.mass)
        if ((atom & ~repair) == 0) mass += value;
    return static_cast<i128>(63) * mass -
           static_cast<i128>(4) * pair.denominator;
}

void add_nonnegative_i128D(FnvLocal& ledger, i128 value) {
    require(value >= 0, "negative margin in positive ledger");
    const __uint128_t wide = static_cast<__uint128_t>(value);
    ledger.add(static_cast<u64>(wide));
    ledger.add(static_cast<u64>(wide >> 64));
}

struct BodyAuditD {
    u64 bodies = 0;
    u64 checks = 0;
    u64 failures = 0;
    u64 fnv = 0;
};

BodyAuditD rescanD(const std::vector<u32>& deck) {
    BodyAuditD out;
    FnvLocal ledger;
    u32 body = (u32{1} << 9) - 1;
    while (body < (u32{1} << 30)) {
        ++out.bodies;
        bool hit = false;
        for (u32 mask : deck) {
            ++out.checks;
            if ((body & mask) == 0) {
                hit = true;
                break;
            }
        }
        out.failures += !hit;
        ledger.add(body);
        ledger.add(hit ? 1 : 0);
        const u32 next = next_combination(body);
        if (next <= body) break;
        body = next;
    }
    require(out.bodies == UINT64_C(14307150), "body count changed");
    out.fnv = ledger.state;
    return out;
}

void choose_from_bits(const std::vector<unsigned>& bits, std::size_t start,
                      unsigned need, u32 mask, std::vector<u32>& out) {
    if (need == 0) {
        out.push_back(mask);
        return;
    }
    for (std::size_t index = start; index + need <= bits.size(); ++index)
        choose_from_bits(bits, index + 1, need - 1,
                         mask | (u32{1} << bits[index]), out);
}

unsigned response(u32 mask, const std::vector<u32>& bodies) {
    unsigned out = 0;
    for (unsigned index = 0; index < bodies.size(); ++index)
        if ((mask & bodies[index]) == 0) out |= 1u << index;
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 5,
                "usage: direct DECK SIGNATURES LIVE THREADS");
        const std::vector<u32> deck = read_deckD(argv[1]);
        require(deck[kIndex] == kOld, "index-265 old mask changed");
        const std::vector<SigRowD> signature_rows = read_signaturesD(argv[2]);
        const std::set<PairD> live = read_liveD(argv[3]);
        const unsigned thread_count = static_cast<unsigned>(std::stoul(argv[4]));
        require(thread_count >= 1 && thread_count <= 8,
                "thread count outside 1..8");

        const auto wanted = singleton_signatureD(kIndex);
        std::vector<PairD> ideal;
        FnvLocal ideal_ledger;
        for (const SigRowD& row : signature_rows) {
            const PairD pair{row.q, row.r};
            if (!live.contains(pair) || row.words != wanted) continue;
            ideal.push_back(pair);
            ideal_ledger.add(row.q);
            ideal_ledger.add(row.r);
        }
        require(ideal.size() == 367 &&
                    ideal_ledger.state == UINT64_C(0xd422b161d94ebae4),
                "index-265 ideal changed");

        u64 target_disjoint = 0;
        u64 retained_checks = 0;
        const std::vector<u32> bodies = private_bodiesD(
            deck, kIndex, target_disjoint, retained_checks);
        require(bodies == std::vector<u32>(kExpectedBodies.begin(),
                                           kExpectedBodies.end()),
                "independent private-body census changed");
        FnvLocal body_ledger;
        for (u32 body : bodies) body_ledger.add(body);

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool-cell count changed");
        std::vector<PairAtomsD> atoms;
        atoms.reserve(ideal.size());
        for (const auto& [q, r] : ideal)
            atoms.push_back(build_pair_atomsD(cells, q, r));

        // Literal verification of the 420 retained joint masks.  This does
        // not trust the signature bit other than for selecting the row set.
        u64 retained_activity_tests = 0;
        u64 retained_equalities = 0;
        FnvLocal retained_margin_ledger;
        for (std::size_t row = 0; row < atoms.size(); ++row) {
            require(marginD(atoms[row], kOld) < 0,
                    "old mask is not strictly inactive on ideal row");
            for (std::size_t index = 0; index < deck.size(); ++index) {
                if (index == kIndex) continue;
                const i128 value = marginD(atoms[row], deck[index]);
                ++retained_activity_tests;
                retained_equalities += value == 0;
                require(value >= 0,
                        "retained joint mask is inactive on ideal row");
                retained_margin_ledger.add(ideal[row].first);
                retained_margin_ledger.add(ideal[row].second);
                retained_margin_ledger.add(deck[index]);
                add_nonnegative_i128D(retained_margin_ledger, value);
            }
        }

        // The integral lower packing from the complete quotient is {B1,B7}.
        // A one-mask cover would have to be disjoint from their union.  Scan
        // the complete choose(17,8)=24,310 universe of such masks directly and
        // exhibit a negative ideal-row margin for every candidate.
        const u32 packed_union = bodies[1] | bodies[7];
        std::vector<unsigned> available;
        for (unsigned bit = 0; bit < 30; ++bit)
            if ((packed_union & (u32{1} << bit)) == 0)
                available.push_back(bit);
        require(available.size() == 17,
                "lower-packing complement size changed");
        std::vector<u32> candidates;
        choose_from_bits(available, 0, 8, 0, candidates);
        require(candidates.size() == 24310 &&
                    std::is_sorted(candidates.begin(), candidates.end()) == false,
                "candidate generator changed");
        std::sort(candidates.begin(), candidates.end());
        require(std::adjacent_find(candidates.begin(), candidates.end()) ==
                    candidates.end(),
                "duplicate lower candidate");

        std::vector<int> first_blocker(candidates.size(), -1);
        std::vector<u64> test_counts(candidates.size(), 0);
        std::atomic<std::size_t> next{0};
        std::vector<std::thread> workers;
        for (unsigned thread = 0; thread < thread_count; ++thread) {
            workers.emplace_back([&]() {
                while (true) {
                    const std::size_t candidate_index = next.fetch_add(1);
                    if (candidate_index >= candidates.size()) break;
                    for (std::size_t row = 0; row < atoms.size(); ++row) {
                        ++test_counts[candidate_index];
                        if (marginD(atoms[row], candidates[candidate_index]) < 0) {
                            first_blocker[candidate_index] = static_cast<int>(row);
                            break;
                        }
                    }
                }
            });
        }
        for (std::thread& worker : workers) worker.join();
        FnvLocal candidate_ledger;
        FnvLocal blocker_ledger;
        u64 lower_tests = 0;
        std::set<PairD> blocker_rows;
        for (std::size_t index = 0; index < candidates.size(); ++index) {
            require(first_blocker[index] >= 0,
                    "one-mask lower-certificate survivor found");
            const PairD blocker = ideal[first_blocker[index]];
            candidate_ledger.add(candidates[index]);
            blocker_ledger.add(candidates[index]);
            blocker_ledger.add(blocker.first);
            blocker_ledger.add(blocker.second);
            lower_tests += test_counts[index];
            blocker_rows.insert(blocker);
        }

        FnvLocal witness_margin_ledger;
        u64 witness_tests = 0;
        u64 witness_equalities = 0;
        unsigned witness_union = 0;
        for (u32 witness : kWitnesses) {
            const unsigned pattern = response(witness, bodies);
            witness_union |= pattern;
            for (std::size_t row = 0; row < atoms.size(); ++row) {
                const i128 value = marginD(atoms[row], witness);
                ++witness_tests;
                witness_equalities += value == 0;
                require(value >= 0, "cover witness is not common-active");
                witness_margin_ledger.add(witness);
                witness_margin_ledger.add(ideal[row].first);
                witness_margin_ledger.add(ideal[row].second);
                add_nonnegative_i128D(witness_margin_ledger, value);
                add_nonnegative_i128D(witness_margin_ledger,
                                     63 * atoms[row].denominator);
            }
            std::cout << "WITNESS " << std::hex << std::setw(8)
                      << std::setfill('0') << witness << " RESPONSE "
                      << std::setw(2) << pattern << std::dec
                      << std::setfill(' ') << '\n';
        }
        require(witness_union == 0xff,
                "two witnesses do not cover all obligations");

        std::vector<u32> rebuilt;
        for (std::size_t index = 0; index < deck.size(); ++index)
            if (index != kIndex) rebuilt.push_back(deck[index]);
        rebuilt.insert(rebuilt.end(), kWitnesses.begin(), kWitnesses.end());
        require(rebuilt.size() == 422 &&
                    std::set<u32>(rebuilt.begin(), rebuilt.end()).size() == 422,
                "rebuilt deck identity malformed");
        FnvLocal rebuilt_ledger;
        for (u32 mask : rebuilt) rebuilt_ledger.add(mask);
        const BodyAuditD body_audit = rescanD(rebuilt);
        require(body_audit.failures == 0,
                "independent rebuilt-deck body scan failed");

        std::cout << "INDEX265_DIRECT_ATOM_AUDIT_V1\n"
                  << "IDEAL ROWS " << ideal.size() << " FNV " << std::hex
                  << ideal_ledger.state << std::dec << '\n'
                  << "PRIVATE_BODIES " << bodies.size() << " FNV " << std::hex
                  << body_ledger.state << std::dec << " TARGET_DISJOINT "
                  << target_disjoint << " RETAINED_CHECKS " << retained_checks
                  << '\n'
                  << "RETAINED_LITERAL_CELLS " << retained_activity_tests
                  << " EQUALITIES " << retained_equalities << " MARGIN_FNV "
                  << std::hex << retained_margin_ledger.state << std::dec << '\n'
                  << "LOWER_PAIR B1,B7 SUPPORT 82 UNION " << std::hex
                  << packed_union << std::dec << " COMPLEMENT_LABELS "
                  << available.size() << " CANDIDATES " << candidates.size()
                  << " CANDIDATE_FNV " << std::hex << candidate_ledger.state
                  << std::dec << " ACTIVITY_TESTS " << lower_tests
                  << " SURVIVORS 0 BLOCKER_ROWS " << blocker_rows.size()
                  << " BLOCKER_FNV " << std::hex << blocker_ledger.state
                  << std::dec << '\n'
                  << "WITNESS_CELLS " << witness_tests << " EQUALITIES "
                  << witness_equalities << " MARGIN_FNV " << std::hex
                  << witness_margin_ledger.state << std::dec << '\n'
                  << "EXACT_MINIMUM 2 LOWER_NO_ONE_MASK UPPER_TWO_MASKS\n"
                  << "REBUILT_DECK " << rebuilt.size() << " FNV " << std::hex
                  << rebuilt_ledger.state << std::dec << " BODY_SCAN "
                  << body_audit.bodies << " CHECKS " << body_audit.checks
                  << " FAILURES " << body_audit.failures << " BODY_FNV "
                  << std::hex << body_audit.fnv << std::dec << '\n'
                  << "SCOPE INDEPENDENT_DIRECT_ATOMS_FIXED_POOL_COMMON_DECK_"
                     "NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS INDEX265_EXACT_MINIMUM_TWO\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "INDEX265_DIRECT_ERROR " << error.what() << '\n';
        return 1;
    }
}

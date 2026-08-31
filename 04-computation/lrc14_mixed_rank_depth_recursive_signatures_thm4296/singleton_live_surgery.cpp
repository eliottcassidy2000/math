// Exploratory exact singleton inactive-signature-ideal surgery inside the
// post-THM4287 residual.  This imports only already-audited arithmetic helpers;
// all candidate, activity, margin, and rebuilt-body conclusions are recomputed.

#define CASCADE_LIBRARY_ONLY
#include "04-computation/lrc14_three_round_learned_carrier_thm4266/cascade_pair_exhaustive_primary.cpp"
#undef CASCADE_LIBRARY_ONLY

#include <fstream>
#include <map>
#include <set>
#include <sstream>

namespace {

using Pair = std::pair<int, int>;

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
    Pair previous{-1, -1};
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
        require(count == static_cast<unsigned>(std::stoul(f[2])) &&
                    (rows.empty() || previous < Pair{row.q, row.r}),
                "signature count/order changed");
        rows.push_back(row);
        previous = {row.q, row.r};
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
    target_disjoint = retained_checks = 0;
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
    for (const auto& [atom, value] : pair.atoms.mass)
        if ((atom & ~repair) == 0) mass += value;
    return static_cast<i128>(63) * mass -
           static_cast<i128>(4) * pair.denominator;
}

u64 mask_fnv(const std::vector<u32>& masks) {
    FnvLocal ledger;
    for (u32 mask : masks) ledger.add(mask);
    return ledger.state;
}

std::set<Pair> read_live(const std::string& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open live residual");
    std::set<Pair> out;
    std::string line;
    while (std::getline(in, line)) {
        const auto fields = split_csv(line);
        require(fields.size() == 2, "bad live row");
        const Pair pair{std::stoi(fields[0]), std::stoi(fields[1])};
        require(pair.first > 0 && pair.first < pair.second &&
                    out.insert(pair).second,
                "bad/duplicate live pair");
    }
    require(out.size() == 22647, "post-THM4287 residual count changed");
    return out;
}

struct BodyAudit {
    u64 bodies = 0;
    u64 checks = 0;
    u64 failures = 0;
    u64 fnv = 0;
};

BodyAudit rescan(const std::vector<u32>& deck) {
    BodyAudit out;
    FnvLocal ledger;
    u32 body = (u32{1} << 9) - 1;
    while (body < (u32{1} << 30)) {
        ++out.bodies;
        bool hit = false;
        for (u32 mask : deck) {
            ++out.checks;
            if ((mask & body) == 0) {
                hit = true;
                break;
            }
        }
        if (!hit) ++out.failures;
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

void add_nonnegative_i128(FnvLocal& ledger, i128 value) {
    require(value >= 0, "negative margin in positive ledger");
    const __uint128_t wide = static_cast<__uint128_t>(value);
    ledger.add(static_cast<u64>(wide));
    ledger.add(static_cast<u64>(wide >> 64));
}

struct FullWitness {
    u32 mask = 0;
    i128 min_num = 0;
    i128 min_den = 1;
    Pair weakest{};
    u64 margin_fnv = 0;
};

} // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        if (argc == 2 && std::string(argv[1]) == "--help") {
            std::cout
                << "usage: singleton-live DECK SIGNATURES LIVE INDEX TARGET_Q,R\n"
                << "writes one deterministic singleton-surgery report to stdout\n";
            return 0;
        }
        require(argc == 6,
                "usage: singleton-live DECK SIGNATURES LIVE INDEX TARGET_Q,R");
        init_choose8_local();
        const std::vector<u32> deck = read_deck(argv[1]);
        const std::vector<SigRow> rows = read_signatures(argv[2]);
        const std::set<Pair> live = read_live(argv[3]);
        const std::size_t target_index = std::stoul(argv[4]);
        const auto target_fields = split_csv(argv[5]);
        require(target_fields.size() == 2, "bad target pair");
        const Pair target{std::stoi(target_fields[0]),
                          std::stoi(target_fields[1])};
        require(target_index < deck.size() && live.contains(target),
                "bad target index/pair");

        const auto wanted = singleton_signature(target_index);
        std::vector<Pair> ideal;
        bool target_found = false;
        FnvLocal ideal_ledger;
        for (const SigRow& row : rows) {
            const Pair pair{row.q, row.r};
            if (row.words != wanted || !live.contains(pair)) continue;
            ideal.push_back(pair);
            ideal_ledger.add(row.q);
            ideal_ledger.add(row.r);
            target_found |= pair == target;
        }
        require(target_found && !ideal.empty(), "target absent from live ideal");

        u64 target_disjoint = 0;
        u64 retained_checks = 0;
        const std::vector<u32> obligations = private_bodies(
            deck, target_index, target_disjoint, retained_checks);
        u32 body_union = 0;
        FnvLocal obligation_ledger;
        for (u32 body : obligations) {
            body_union |= body;
            obligation_ledger.add(body);
        }

        std::vector<u32> responders;
        if (obligations.empty()) {
            responders.push_back(0);
        } else if (std::popcount(body_union) <= 22) {
            u32 repair = (u32{1} << 8) - 1;
            while (repair < (u32{1} << 30)) {
                if ((repair & body_union) == 0) responders.push_back(repair);
                const u32 next = next_combination(repair);
                if (next <= repair) break;
                repair = next;
            }
        }

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool-cell count changed");
        std::vector<PairAtoms> atoms;
        atoms.reserve(ideal.size());
        for (const auto& [q, r] : ideal)
            atoms.push_back(build_pair_atoms(cells, q, r));

        std::vector<FullWitness> full;
        u64 activity_tests = 0;
        u64 activity_equalities = 0;
        for (u32 candidate : responders) {
            if (candidate == 0) continue;
            bool all = true;
            i128 minimum_num = 0;
            i128 minimum_den = 1;
            Pair weakest{};
            bool have_minimum = false;
            FnvLocal margin_ledger;
            for (const PairAtoms& pair : atoms) {
                const i128 value = margin(pair, candidate);
                ++activity_tests;
                activity_equalities += value == 0;
                if (value < 0) {
                    all = false;
                    break;
                }
                if (value >= 0 && (!have_minimum ||
                        value * minimum_den < minimum_num *
                                                   (63 * pair.denominator))) {
                    minimum_num = value;
                    minimum_den = 63 * pair.denominator;
                    weakest = {pair.q, pair.r};
                    have_minimum = true;
                }
                margin_ledger.add(pair.q);
                margin_ledger.add(pair.r);
                margin_ledger.add(candidate);
                const i128 shifted = value >= 0 ? value : -value;
                margin_ledger.add(value >= 0 ? 1 : 0);
                add_nonnegative_i128(margin_ledger, shifted);
                add_nonnegative_i128(margin_ledger, 63 * pair.denominator);
            }
            if (all) {
                const i128 divisor = gcd_i128(minimum_num, minimum_den);
                full.push_back({candidate, minimum_num / divisor,
                                minimum_den / divisor, weakest,
                                margin_ledger.state});
            }
        }

        std::vector<u32> rebuilt;
        BodyAudit body_audit;
        FnvLocal rebuilt_ledger;
        if (obligations.empty() || !full.empty()) {
            for (std::size_t index = 0; index < deck.size(); ++index)
                if (index != target_index) rebuilt.push_back(deck[index]);
            if (!obligations.empty()) rebuilt.push_back(full.front().mask);
            std::set<u32> distinct(rebuilt.begin(), rebuilt.end());
            require(distinct.size() == rebuilt.size(),
                    "rebuilt deck has duplicate mask");
            for (u32 mask : rebuilt) rebuilt_ledger.add(mask);
            body_audit = rescan(rebuilt);
            require(body_audit.failures == 0, "rebuilt deck fails body cover");
        }

        std::cout << "SINGLETON_LIVE_SIGNATURE_IDEAL_SURGERY_V1\n"
                  << "LIVE_ROWS " << live.size() << " TARGET "
                  << target.first << ',' << target.second << " J "
                  << target_index << " OLD_MASK " << std::hex
                  << std::setw(8) << std::setfill('0') << deck[target_index]
                  << std::dec << std::setfill(' ') << '\n'
                  << "IDEAL_ROWS " << ideal.size() << " FNV " << std::hex
                  << ideal_ledger.state << std::dec << " ROWS";
        for (const auto& [q, r] : ideal) std::cout << ' ' << q << ',' << r;
        std::cout << "\nOBLIGATIONS " << obligations.size() << " FNV "
                  << std::hex << obligation_ledger.state << " UNION "
                  << std::setw(8) << std::setfill('0') << body_union
                  << std::dec << std::setfill(' ') << " UNION_SIZE "
                  << std::popcount(body_union) << " TARGET_DISJOINT "
                  << target_disjoint << " RETAINED_CHECKS "
                  << retained_checks << " BODIES";
        for (u32 body : obligations)
            std::cout << ' ' << std::hex << std::setw(8)
                      << std::setfill('0') << body << std::dec
                      << std::setfill(' ');
        std::cout << "\nFULL_BODY_RESPONDERS " << responders.size()
                  << " FNV " << std::hex << mask_fnv(responders) << std::dec
                  << " ACTIVITY_TESTS " << activity_tests
                  << " EQUALITIES " << activity_equalities
                  << " COMMON_ACTIVE " << full.size() << " MASKS";
        for (const FullWitness& witness : full)
            std::cout << ' ' << std::hex << std::setw(8)
                      << std::setfill('0') << witness.mask << std::dec
                      << std::setfill(' ');
        std::cout << '\n';
        for (const FullWitness& witness : full) {
            std::cout << "WITNESS " << std::hex << std::setw(8)
                      << std::setfill('0') << witness.mask << std::dec
                      << std::setfill(' ') << " WEAKEST_PAIR "
                      << witness.weakest.first << ',' << witness.weakest.second
                      << " NORMALIZED_GAP " << decimal(witness.min_num)
                      << '/' << decimal(witness.min_den) << " MARGIN_FNV "
                      << std::hex << witness.margin_fnv << std::dec << '\n';
        }
        if (!rebuilt.empty()) {
            std::cout << "REBUILT_DECK " << rebuilt.size() << " FNV "
                      << std::hex << rebuilt_ledger.state << std::dec
                      << " BODY_SCAN " << body_audit.bodies << " CHECKS "
                      << body_audit.checks << " FAILURES "
                      << body_audit.failures << " BODY_FNV " << std::hex
                      << body_audit.fnv << std::dec << '\n';
        }
        std::cout << "EXACT_MINIMUM ";
        if (obligations.empty()) std::cout << "0";
        else if (!full.empty()) std::cout << "1";
        else std::cout << ">1";
        std::cout << " CONSUMER COMMON_DECK_NODE NEW_LIVE_ROWS "
                  << (full.empty() && !obligations.empty() ? 0 : ideal.size())
                  << "\nVERDICT PASS FINITE_EXACT_EXPLORATORY\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "SINGLETON_LIVE_ERROR " << error.what() << '\n';
        return 1;
    }
}

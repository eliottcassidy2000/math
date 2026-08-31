// Exact full common-activity response quotient for the post-THM4287
// singleton inactive-signature ideal J={19}.

#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "04-computation/lrc14_joint421_global_common_carrier_thm4281/response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

#include <array>
#include <bit>
#include <fstream>
#include <iomanip>
#include <map>
#include <queue>
#include <set>
#include <sstream>

namespace {

using Pair19 = std::pair<int, int>;
constexpr std::size_t TARGET_INDEX19 = 19;

std::vector<std::string> split19(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream in(line);
    std::string field;
    while (std::getline(in, field, ',')) out.push_back(field);
    return out;
}

std::vector<u32> read_deck19(const std::string& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open deck");
    std::vector<u32> deck;
    std::set<u32> seen;
    FnvLocal ledger;
    std::string token;
    while (in >> token) {
        const u32 mask = static_cast<u32>(std::stoul(token, nullptr, 16));
        require(mask < (u32{1} << 30) && std::popcount(mask) == 8 &&
                    seen.insert(mask).second,
                "bad deck mask");
        deck.push_back(mask);
        ledger.add(mask);
    }
    require(deck.size() == 421 &&
                ledger.state == UINT64_C(0x20d63dd42fe8150e) &&
                deck[TARGET_INDEX19] == UINT32_C(0x1804aa01),
            "joint deck identity changed");
    return deck;
}

std::set<Pair19> read_live19(const std::string& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open live residual");
    std::set<Pair19> out;
    std::string line;
    while (std::getline(in, line)) {
        const auto f = split19(line);
        require(f.size() == 2, "bad live row");
        require(out.insert({std::stoi(f[0]), std::stoi(f[1])}).second,
                "duplicate live row");
    }
    require(out.size() == 22647, "live residual count changed");
    return out;
}

std::vector<Pair19> read_ideal19(const std::string& path,
                                 const std::set<Pair19>& live,
                                 u64& ledger_value) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open signatures");
    std::string line;
    require(static_cast<bool>(std::getline(in, line)) &&
                line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "signature header changed");
    std::vector<Pair19> ideal;
    FnvLocal ledger;
    std::size_t all_rows = 0;
    while (std::getline(in, line)) {
        ++all_rows;
        const auto f = split19(line);
        require(f.size() == 10, "bad signature row");
        const Pair19 pair{std::stoi(f[0]), std::stoi(f[1])};
        if (!live.contains(pair) || std::stoul(f[2]) != 1) continue;
        const std::size_t word = TARGET_INDEX19 / 64;
        const u64 bit = UINT64_C(1) << (TARGET_INDEX19 % 64);
        bool wanted = true;
        for (std::size_t i = 0; i < 7; ++i) {
            const u64 value = std::stoull(f[3 + i], nullptr, 16);
            if (value != (i == word ? bit : 0)) wanted = false;
        }
        if (!wanted) continue;
        ideal.push_back(pair);
        ledger.add(pair.first);
        ledger.add(pair.second);
    }
    require(all_rows == 24223 && ideal.size() == 36 &&
                std::find(ideal.begin(), ideal.end(), Pair19{338,636}) !=
                    ideal.end(),
            "J19 ideal changed");
    ledger_value = ledger.state;
    return ideal;
}

u32 next19(u32 value) {
    const u32 low = value & (0u - value);
    const u32 ripple = value + low;
    return ripple | (((value ^ ripple) >> 2) / low);
}

std::vector<u32> obligations19(const std::vector<u32>& deck,
                               u64& checks, u64& fnv_value,
                               u32& body_union) {
    std::vector<u32> out;
    checks = 0;
    body_union = 0;
    FnvLocal ledger;
    u32 body = (u32{1} << 9) - 1;
    for (std::size_t ordinal = 0; ordinal < 14307150; ++ordinal) {
        if ((body & deck[TARGET_INDEX19]) == 0) {
            bool retained = false;
            for (std::size_t index = 0; index < deck.size(); ++index) {
                if (index == TARGET_INDEX19) continue;
                ++checks;
                if ((body & deck[index]) == 0) {
                    retained = true;
                    break;
                }
            }
            if (!retained) {
                out.push_back(body);
                body_union |= body;
                ledger.add(body);
            }
        }
        if (ordinal + 1 < 14307150) body = next19(body);
    }
    require(body == UINT32_C(0x3fe00000) && out.size() == 8 &&
                body_union == UINT32_C(0x27fb548a),
            "J19 obligations changed");
    fnv_value = ledger.state;
    return out;
}

i128 direct_mass19(const AtomData& atoms, u32 mask) {
    i128 mass = 0;
    for (const auto& [atom, value] : atoms.mass)
        if ((atom & ~mask) == 0) mass += value;
    return mass;
}

struct BodyAudit19 {
    u64 checks = 0;
    u64 failures = 0;
    u64 fnv = 0;
};

BodyAudit19 rescan19(const std::vector<u32>& deck) {
    BodyAudit19 out;
    FnvLocal ledger;
    u32 body = (u32{1} << 9) - 1;
    for (std::size_t ordinal = 0; ordinal < 14307150; ++ordinal) {
        bool hit = false;
        for (u32 mask : deck) {
            ++out.checks;
            if ((body & mask) == 0) {
                hit = true;
                break;
            }
        }
        if (!hit) ++out.failures;
        ledger.add(body);
        ledger.add(hit ? 1 : 0);
        if (ordinal + 1 < 14307150) body = next19(body);
    }
    require(body == UINT32_C(0x3fe00000), "rebuilt body endpoint changed");
    out.fnv = ledger.state;
    return out;
}

void add_i128_19(FnvLocal& ledger, i128 value) {
    require(value >= 0, "negative i128 ledger value");
    const __uint128_t wide = static_cast<__uint128_t>(value);
    ledger.add(static_cast<u64>(wide));
    ledger.add(static_cast<u64>(wide >> 64));
}

} // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 4, "usage: j19-full DECK SIGNATURES LIVE");
        init_choose8_local();
        const std::vector<u32> deck = read_deck19(argv[1]);
        const std::set<Pair19> live = read_live19(argv[3]);
        u64 ideal_fnv = 0;
        const std::vector<Pair19> ideal = read_ideal19(argv[2], live, ideal_fnv);
        u64 obligation_checks = 0;
        u64 obligation_fnv = 0;
        u32 body_union = 0;
        const std::vector<u32> obligations = obligations19(
            deck, obligation_checks, obligation_fnv, body_union);

        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        std::vector<unsigned char> common_active(EXPECTED_REPAIRS, 1);
        u64 common_stage_fnv = 0;
        for (std::size_t row = 0; row < ideal.size(); ++row) {
            const ActiveUniverse active =
                build_active_universe(cells, ideal[row].first, ideal[row].second);
            for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank)
                common_active[rank] &= active.active[rank];
            FnvLocal stage;
            stage.add(row);
            stage.add(ideal[row].first);
            stage.add(ideal[row].second);
            for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank)
                if (common_active[rank]) stage.add(unrank_colex8(rank));
            common_stage_fnv = stage.state;
            std::cout << "INTERSECTION_ROW " << row + 1 << " OF "
                      << ideal.size() << " PAIR " << ideal[row].first << ','
                      << ideal[row].second << " ACTIVE "
                      << std::count(common_active.begin(), common_active.end(), 1)
                      << " STAGE_FNV " << std::hex << common_stage_fnv
                      << std::dec << '\n';
        }

        struct Class19 { u64 count = 0; u32 least = 0; };
        std::map<unsigned, Class19> classes;
        FnvLocal common_ledger;
        FnvLocal response_ledger;
        u64 common_count = 0;
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            if (!common_active[rank]) continue;
            const u32 mask = unrank_colex8(rank);
            ++common_count;
            common_ledger.add(mask);
            unsigned pattern = 0;
            for (unsigned index = 0; index < obligations.size(); ++index)
                if ((mask & obligations[index]) == 0)
                    pattern |= 1u << index;
            response_ledger.add(mask);
            response_ledger.add(pattern);
            Class19& cls = classes[pattern];
            ++cls.count;
            if (cls.count == 1 || mask < cls.least) cls.least = mask;
        }

        constexpr unsigned full_pattern = 0xff;
        std::array<int, 256> distance;
        std::array<unsigned, 256> predecessor{};
        std::array<unsigned, 256> used{};
        distance.fill(-1);
        std::queue<unsigned> queue;
        distance[0] = 0;
        queue.push(0);
        while (!queue.empty()) {
            const unsigned state = queue.front();
            queue.pop();
            for (const auto& [pattern, cls] : classes) {
                (void)cls;
                if (pattern == 0) continue;
                const unsigned joined = state | pattern;
                if (distance[joined] >= 0) continue;
                distance[joined] = distance[state] + 1;
                predecessor[joined] = state;
                used[joined] = pattern;
                queue.push(joined);
            }
        }
        require(distance[full_pattern] > 0, "common-active quotient not coverable");
        std::vector<u32> witness;
        for (unsigned state = full_pattern; state != 0;
             state = predecessor[state])
            witness.push_back(classes.at(used[state]).least);
        std::reverse(witness.begin(), witness.end());

        u64 equality_count = 0;
        FnvLocal margin_ledger;
        for (u32 mask : witness) {
            bool have_min = false;
            i128 min_num = 0;
            i128 min_den = 1;
            Pair19 weakest{};
            for (const Pair19 pair : ideal) {
                const i64 g = std::gcd(pair.first, pair.second);
                const PrimitivePair primitive =
                    build_primitive(pair.first / g, pair.second / g);
                const AtomData atoms = build_cocycle_atoms(cells, primitive, g);
                const i128 denominator =
                    static_cast<i128>(primitive.grid) * g * COMMON;
                const i128 raw = static_cast<i128>(63) *
                                     direct_mass19(atoms, mask) -
                                 static_cast<i128>(4) * denominator;
                require(raw >= 0, "witness lost common activity");
                equality_count += raw == 0;
                const i128 gap_den = static_cast<i128>(63) * denominator;
                if (!have_min || raw * min_den < min_num * gap_den) {
                    min_num = raw;
                    min_den = gap_den;
                    weakest = pair;
                    have_min = true;
                }
                margin_ledger.add(mask);
                margin_ledger.add(pair.first);
                margin_ledger.add(pair.second);
                add_i128_19(margin_ledger, raw);
                add_i128_19(margin_ledger, gap_den);
            }
            const i128 divisor = gcd_i128(min_num, min_den);
            std::cout << "WITNESS_MARGIN MASK " << std::hex << std::setw(8)
                      << std::setfill('0') << mask << std::dec
                      << std::setfill(' ') << " WEAKEST_PAIR "
                      << weakest.first << ',' << weakest.second
                      << " NORMALIZED_GAP " << decimal(min_num/divisor)
                      << '/' << decimal(min_den/divisor) << '\n';
        }

        std::vector<u32> rebuilt;
        for (std::size_t index = 0; index < deck.size(); ++index)
            if (index != TARGET_INDEX19) rebuilt.push_back(deck[index]);
        rebuilt.insert(rebuilt.end(), witness.begin(), witness.end());
        std::set<u32> distinct(rebuilt.begin(), rebuilt.end());
        require(distinct.size() == rebuilt.size(), "rebuilt masks not distinct");
        FnvLocal rebuilt_ledger;
        for (u32 mask : rebuilt) rebuilt_ledger.add(mask);
        const BodyAudit19 body_audit = rescan19(rebuilt);
        require(body_audit.failures == 0, "rebuilt deck fails body cover");

        std::cout << "J19_FULL_COMMON_RESPONSE_V1\n"
                  << "IDEAL 19 ROWS " << ideal.size() << " FNV " << std::hex
                  << ideal_fnv << std::dec << " MAX_ENDPOINT 636\n"
                  << "OBLIGATIONS " << obligations.size() << " FNV "
                  << std::hex << obligation_fnv << " UNION " << body_union
                  << std::dec << " RETAINED_CHECKS " << obligation_checks
                  << "\nCOMMON_ACTIVE " << common_count << " FNV "
                  << std::hex << common_ledger.state << " RESPONSE_FNV "
                  << response_ledger.state << std::dec << " CLASSES "
                  << classes.size() << " FULL_RESPONDER "
                  << (classes.contains(full_pattern) ? 1 : 0) << '\n';
        for (const auto& [pattern, cls] : classes)
            std::cout << "CLASS " << std::hex << pattern << " LEAST "
                      << std::setw(8) << std::setfill('0') << cls.least
                      << std::dec << std::setfill(' ') << " COVER "
                      << std::popcount(pattern) << " COUNT " << cls.count
                      << '\n';
        std::cout << "EXACT_MINIMUM " << distance[full_pattern]
                  << " WITNESS";
        for (u32 mask : witness)
            std::cout << ' ' << std::hex << std::setw(8)
                      << std::setfill('0') << mask << std::dec
                      << std::setfill(' ');
        std::cout << " EQUALITIES " << equality_count << " MARGIN_FNV "
                  << std::hex << margin_ledger.state << std::dec << '\n'
                  << "REBUILT_DECK " << rebuilt.size() << " FNV "
                  << std::hex << rebuilt_ledger.state << std::dec
                  << " BODY_SCAN 14307150 CHECKS " << body_audit.checks
                  << " FAILURES " << body_audit.failures << " BODY_FNV "
                  << std::hex << body_audit.fnv << std::dec << '\n'
                  << "CONSUMER COMMON_DECK_NODE NEW_LIVE_ROWS 36 "
                     "PROOF_GRAPH_OVERLAP_CURRENT 0\n"
                  << "SCOPE FIXED_POOL_FINITE_EXACT NO_PHYSICAL_ENTRY\n"
                  << "VERDICT PASS EXACT_FULL_COMMON_RESPONSE_QUOTIENT\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "J19_FULL_COMMON_RESPONSE_ERROR " << error.what() << '\n';
        return 1;
    }
}

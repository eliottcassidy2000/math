// Scratch-only structurally independent audit of the proposed index-297
// recursive signature ideal.  We reuse the maintained import-free literal
// wall primitives from the detached index-19 audit, but independently derive
// the ideal, private bodies, common response quotient, dual, normalized
// margins, rebuilt deck, body cover, and old-union overlap.

#define main thm4296_detached_j19_unused_main
#include "04-computation/lrc14_mixed_rank_depth_recursive_signatures_thm4296/detached_j19_audit.cpp"
#undef main

namespace {

constexpr std::size_t kTarget297 = 297;
constexpr std::size_t kIdealCount297 = 42;
constexpr u64 kIdealFnv297 = UINT64_C(0x211843ee21a19170);
constexpr u32 kOldMask297 = UINT32_C(0x08b204c0);
constexpr std::array<u32, 2> kWitness297 = {
    UINT32_C(0x089284c0), UINT32_C(0x08330481)};
constexpr std::array<u32, 4> kBodies297 = {
    UINT32_C(0x17087001), UINT32_C(0x2248d208),
    UINT32_C(0x27446008), UINT32_C(0x3548a008)};
constexpr u64 kObligationFnv297 = UINT64_C(0x81e5ca5f045bb2cf);
constexpr u32 kObligationUnion297 = UINT32_C(0x374cf209);
constexpr u64 kCommonCount297 = UINT64_C(2208823);
constexpr u64 kCommonFnv297 = UINT64_C(0xdfc36b5aa486b0ea);
constexpr u64 kResponseFnv297 = UINT64_C(0xb1665afb732ece83);
constexpr u64 kRebuiltFnv297 = UINT64_C(0x60d75261322593ac);

std::vector<u32> read_deck297(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open joint deck");
    std::vector<u32> deck;
    std::set<u32> distinct;
    Fnv ledger;
    std::string token;
    while (input >> token) {
        std::size_t used = 0;
        const u64 wide = std::stoull(token, &used, 16);
        require(used == token.size() && wide < (UINT64_C(1) << 30),
                "invalid deck token");
        const u32 mask = static_cast<u32>(wide);
        require(std::popcount(mask) == 8 && distinct.insert(mask).second,
                "invalid or repeated deck mask");
        deck.push_back(mask);
        ledger.add(mask);
    }
    require(input.eof() && deck.size() == kDeckCount &&
                ledger.state == kDeckFnv,
            "joint deck identity changed");
    require(deck[kTarget297] == kOldMask297,
            "index-297 deck mask differs from candidate");
    return deck;
}

struct Ideal297 {
    std::vector<Pair> rows;
    u64 fnv = 0;
    int maximum_endpoint = 0;
};

Ideal297 read_ideal297(const std::filesystem::path& path,
                       const std::set<Pair>& live,
                       const std::filesystem::path& output_path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open signature atlas");
    std::string line;
    require(static_cast<bool>(std::getline(input, line)) &&
                line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "signature header changed");
    Ideal297 out;
    Fnv ledger;
    std::size_t atlas_rows = 0;
    while (std::getline(input, line)) {
        ++atlas_rows;
        const std::vector<std::string> fields = split(line, ',');
        require(fields.size() == 10, "malformed signature row");
        const Pair pair{std::stoi(fields[0]), std::stoi(fields[1])};
        std::array<u64, 7> words{};
        unsigned signature_weight = 0;
        for (std::size_t word = 0; word < words.size(); ++word) {
            words[word] = std::stoull(fields[3 + word], nullptr, 16);
            signature_weight += std::popcount(words[word]);
        }
        require(signature_weight == std::stoul(fields[2]),
                "inactive count disagrees with signature words");
        if (!live.contains(pair) || signature_weight != 1) continue;
        const std::size_t target_word = kTarget297 / 64;
        const u64 target_bit = UINT64_C(1) << (kTarget297 % 64);
        bool exact = true;
        for (std::size_t word = 0; word < words.size(); ++word)
            exact &= words[word] == (word == target_word ? target_bit : 0);
        if (!exact) continue;
        out.rows.push_back(pair);
        out.maximum_endpoint = std::max(out.maximum_endpoint, pair.r);
        ledger.add(pair.q);
        ledger.add(pair.r);
    }
    require(input.eof() && atlas_rows == kAtlasRows,
            "signature atlas row count changed");
    out.fnv = ledger.state;
    require(out.rows.size() == kIdealCount297 && out.fnv == kIdealFnv297,
            "index-297 ideal identity differs from candidate");

    std::ofstream ideal_out(output_path, std::ios::binary);
    require(static_cast<bool>(ideal_out), "cannot create ideal ledger");
    for (Pair pair : out.rows) ideal_out << pair.q << ',' << pair.r << '\n';
    require(ideal_out.good(), "ideal-ledger write failed");
    return out;
}

struct Obligations297 {
    std::vector<u32> bodies;
    u64 target_disjoint = 0;
    u64 retained_checks = 0;
    u64 fnv = 0;
    u32 body_union = 0;
};

Obligations297 private_obligations297(const std::vector<u32>& deck) {
    Obligations297 out;
    Fnv ledger;
    u32 body = (u32{1} << 9) - 1;
    for (u64 ordinal = 0; ordinal < kBodyCount; ++ordinal) {
        if ((body & deck[kTarget297]) == 0) {
            ++out.target_disjoint;
            bool retained = false;
            for (std::size_t index = 0; index < deck.size(); ++index) {
                if (index == kTarget297) continue;
                ++out.retained_checks;
                if ((body & deck[index]) == 0) {
                    retained = true;
                    break;
                }
            }
            if (!retained) {
                out.bodies.push_back(body);
                out.body_union |= body;
                ledger.add(body);
            }
        }
        if (ordinal + 1 < kBodyCount) body = next_combination(body);
    }
    require(body == UINT32_C(0x3fe00000), "body enumeration endpoint changed");
    out.fnv = ledger.state;
    require(out.target_disjoint == UINT64_C(497420) &&
                out.bodies == std::vector<u32>(kBodies297.begin(), kBodies297.end()) &&
                out.fnv == kObligationFnv297 &&
                out.body_union == kObligationUnion297,
            "index-297 private bodies differ from candidate");
    return out;
}

struct RowActivity297 {
    Pair pair;
    i64 grid = 0;
    u64 cells = 0;
    u64 failure_classes = 0;
    u64 zeta_operations = 0;
    u64 active_count = 0;
    std::vector<u64> active;
    i128 old_margin = 0;
    std::vector<i128> rebuilt_margins;
};

RowActivity297 compute_activity297(Pair pair,
                                   const std::vector<u32>& rebuilt) {
    const Geometry geometry = build_geometry(pair);
    std::vector<i64> masses(kRepairCount, 0);
    u64 operations = 0;
    for (const auto& [failure, width] : geometry.widths) {
        const int arity = std::popcount(failure);
        if (arity <= 8)
            add_supersets(failure, 8 - arity, 0, 0, width,
                          masses, operations);
    }
    RowActivity297 out;
    out.pair = pair;
    out.grid = geometry.grid;
    out.cells = geometry.cells;
    out.failure_classes = geometry.widths.size();
    out.zeta_operations = operations;
    out.active.assign((kRepairCount + 63) / 64, 0);
    for (u64 rank = 0; rank < kRepairCount; ++rank) {
        const i128 margin = static_cast<i128>(63) * masses[rank] -
                            static_cast<i128>(4) * geometry.grid;
        if (margin >= 0) {
            out.active[rank / 64] |= UINT64_C(1) << (rank % 64);
            ++out.active_count;
        }
    }
    out.old_margin = static_cast<i128>(63) * masses[colex_rank(kOldMask297)] -
                     static_cast<i128>(4) * geometry.grid;
    out.rebuilt_margins.reserve(rebuilt.size());
    for (u32 mask : rebuilt)
        out.rebuilt_margins.push_back(
            static_cast<i128>(63) * masses[colex_rank(mask)] -
            static_cast<i128>(4) * geometry.grid);
    return out;
}

std::set<Pair> read_pair_set297(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open pair set");
    std::set<Pair> out;
    std::string line;
    while (std::getline(input, line)) {
        require(!line.empty(), "blank pair-set row");
        const std::vector<std::string> fields = split(line, ',');
        require(fields.size() == 2, "malformed pair-set row");
        const Pair pair{std::stoi(fields[0]), std::stoi(fields[1])};
        require(out.insert(pair).second, "duplicate pair-set row");
    }
    return out;
}

struct ResponseClass297 {
    u64 count = 0;
    u32 least = 0;
};

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 9,
                "usage: audit DECK SIGNATURES LIVE TYPED_UNION "
                "IDEAL_OUT MARGINS_OUT THREADS EXPECTED_UNION_COUNT");
        init_choose();
        const std::vector<u32> deck = read_deck297(argv[1]);
        const std::set<Pair> live = read_live(argv[3]);
        const Ideal297 ideal = read_ideal297(argv[2], live, argv[5]);
        const Obligations297 obligations = private_obligations297(deck);

        std::vector<u32> rebuilt;
        rebuilt.reserve(422);
        for (std::size_t index = 0; index < deck.size(); ++index)
            if (index != kTarget297) rebuilt.push_back(deck[index]);
        rebuilt.insert(rebuilt.end(), kWitness297.begin(), kWitness297.end());
        require(rebuilt.size() == 422 &&
                    std::set<u32>(rebuilt.begin(), rebuilt.end()).size() == 422,
                "rebuilt deck rank/distinctness changed");
        Fnv rebuilt_ledger;
        for (u32 mask : rebuilt) rebuilt_ledger.add(mask);
        require(rebuilt_ledger.state == kRebuiltFnv297,
                "rebuilt deck FNV differs from candidate");

        const unsigned thread_count = static_cast<unsigned>(std::stoul(argv[7]));
        require(thread_count >= 1 && thread_count <= 8,
                "thread count outside 1..8");
        std::vector<RowActivity297> rows(ideal.rows.size());
        std::atomic<std::size_t> next{0};
        std::vector<std::thread> workers;
        for (unsigned thread = 0; thread < thread_count; ++thread) {
            workers.emplace_back([&]() {
                while (true) {
                    const std::size_t index = next.fetch_add(1);
                    if (index >= ideal.rows.size()) break;
                    rows[index] = compute_activity297(ideal.rows[index], rebuilt);
                }
            });
        }
        for (std::thread& worker : workers) worker.join();

        std::vector<u64> common((kRepairCount + 63) / 64, ~UINT64_C(0));
        if (kRepairCount % 64 != 0)
            common.back() = (UINT64_C(1) << (kRepairCount % 64)) - 1;
        std::ofstream margins_out(argv[6], std::ios::binary);
        require(static_cast<bool>(margins_out), "cannot create margin ledger");
        margins_out << "q,r,grid,old_margin,witness0_margin,witness1_margin\n";

        Fnv activity_ledger;
        Fnv strict_ledger;
        u64 strict_cells = 0;
        u64 equalities = 0;
        u64 operation_total = 0;
        bool global_first = true;
        i128 global_min_num = 0;
        i128 global_min_den = 1;
        Pair global_weakest{};
        u32 global_weakest_mask = 0;
        for (std::size_t row = 0; row < rows.size(); ++row) {
            require(rows[row].pair == ideal.rows[row] && rows[row].old_margin < 0,
                    "old index-297 mask is not strictly inactive");
            require(rows[row].rebuilt_margins.size() == rebuilt.size(),
                    "rebuilt margin vector changed");
            for (std::size_t index = 0; index < rebuilt.size(); ++index) {
                const i128 margin = rows[row].rebuilt_margins[index];
                require(margin >= 0, "rebuilt mask inactive on ideal row");
                ++strict_cells;
                equalities += margin == 0;
                strict_ledger.add(ideal.rows[row].q);
                strict_ledger.add(ideal.rows[row].r);
                strict_ledger.add(rebuilt[index]);
                add_i128(strict_ledger, margin);
                const i128 denominator = static_cast<i128>(63) * rows[row].grid;
                if (global_first || fraction_less(margin, denominator,
                                                  global_min_num,
                                                  global_min_den)) {
                    global_first = false;
                    global_min_num = margin;
                    global_min_den = denominator;
                    global_weakest = ideal.rows[row];
                    global_weakest_mask = rebuilt[index];
                }
            }
            for (std::size_t word = 0; word < common.size(); ++word)
                common[word] &= rows[row].active[word];
            operation_total += rows[row].zeta_operations;
            activity_ledger.add(ideal.rows[row].q);
            activity_ledger.add(ideal.rows[row].r);
            activity_ledger.add(rows[row].grid);
            activity_ledger.add(rows[row].cells);
            activity_ledger.add(rows[row].failure_classes);
            activity_ledger.add(rows[row].zeta_operations);
            activity_ledger.add(rows[row].active_count);
            add_i128(activity_ledger, rows[row].old_margin);
            const std::size_t w0 = rebuilt.size() - 2;
            const std::size_t w1 = rebuilt.size() - 1;
            margins_out << ideal.rows[row].q << ',' << ideal.rows[row].r << ','
                        << rows[row].grid << ',' << decimal(rows[row].old_margin)
                        << ',' << decimal(rows[row].rebuilt_margins[w0])
                        << ',' << decimal(rows[row].rebuilt_margins[w1]) << '\n';
        }
        require(margins_out.good(), "margin-ledger write failed");
        require(equalities == 0, "rebuilt deck has an activity equality");

        std::map<unsigned, ResponseClass297> classes;
        Fnv common_ledger;
        Fnv response_ledger;
        u64 common_count = 0;
        u32 repair = (u32{1} << 8) - 1;
        for (u64 rank = 0; rank < kRepairCount; ++rank) {
            if (common[rank / 64] & (UINT64_C(1) << (rank % 64))) {
                ++common_count;
                common_ledger.add(repair);
                unsigned response = 0;
                for (unsigned index = 0; index < obligations.bodies.size(); ++index)
                    if ((repair & obligations.bodies[index]) == 0)
                        response |= 1u << index;
                response_ledger.add(repair);
                response_ledger.add(response);
                ResponseClass297& cls = classes[response];
                ++cls.count;
                if (cls.count == 1 || repair < cls.least) cls.least = repair;
            }
            if (rank + 1 < kRepairCount) repair = next_combination(repair);
        }
        require(repair == UINT32_C(0x3fc00000) &&
                    common_count == kCommonCount297 &&
                    common_ledger.state == kCommonFnv297 &&
                    response_ledger.state == kResponseFnv297 &&
                    classes.size() == 10 && !classes.contains(0xf),
                "common response quotient differs from candidate");

        std::vector<unsigned> maximal;
        for (const auto& [pattern, cls] : classes) {
            if (pattern == 0) continue;
            bool dominated = false;
            for (const auto& [other, other_cls] : classes)
                if (other != pattern && (pattern & ~other) == 0) {
                    dominated = true;
                    break;
                }
            if (!dominated) maximal.push_back(pattern);
        }
        require(maximal == std::vector<unsigned>({0x5, 0xe}),
                "maximal response antichain differs from candidate");

        // Exact integral dual y_0=y_3=1.  Every realized response has load <=1.
        unsigned maximum_dual_load = 0;
        for (const auto& [pattern, cls] : classes) {
            const unsigned load = ((pattern & 0x1) != 0) + ((pattern & 0x8) != 0);
            maximum_dual_load = std::max(maximum_dual_load, load);
            require(load <= 1, "candidate dual violates a response class");
        }
        require(maximum_dual_load == 1, "candidate dual is unexpectedly slack");

        unsigned witness_union = 0;
        Fnv witness_margin_ledger;
        for (std::size_t witness_index = 0;
             witness_index < kWitness297.size(); ++witness_index) {
            const u32 witness = kWitness297[witness_index];
            const u64 witness_rank = colex_rank(witness);
            require(common[witness_rank / 64] &
                        (UINT64_C(1) << (witness_rank % 64)),
                    "candidate witness is not common-active");
            unsigned response = 0;
            for (unsigned index = 0; index < obligations.bodies.size(); ++index)
                if ((witness & obligations.bodies[index]) == 0)
                    response |= 1u << index;
            witness_union |= response;
            bool first = true;
            i128 minimum_num = 0;
            i128 minimum_den = 1;
            Pair weakest{};
            const std::size_t rebuilt_index = rebuilt.size() - 2 + witness_index;
            for (std::size_t row = 0; row < rows.size(); ++row) {
                const i128 margin = rows[row].rebuilt_margins[rebuilt_index];
                const i128 denominator = static_cast<i128>(63) * rows[row].grid;
                if (first || fraction_less(margin, denominator,
                                           minimum_num, minimum_den)) {
                    first = false;
                    minimum_num = margin;
                    minimum_den = denominator;
                    weakest = ideal.rows[row];
                }
                witness_margin_ledger.add(witness);
                witness_margin_ledger.add(ideal.rows[row].q);
                witness_margin_ledger.add(ideal.rows[row].r);
                add_i128(witness_margin_ledger, margin);
                add_i128(witness_margin_ledger, denominator);
            }
            const i128 divisor = gcd128(minimum_num, minimum_den);
            std::cout << "WITNESS " << std::hex << std::setw(8)
                      << std::setfill('0') << witness << " RESPONSE "
                      << response << std::dec << std::setfill(' ')
                      << " WEAKEST " << weakest.q << ',' << weakest.r
                      << " GAP " << decimal(minimum_num / divisor) << '/'
                      << decimal(minimum_den / divisor) << '\n';
        }
        require(witness_union == 0xf,
                "candidate witnesses do not cover all obligations");

        const BodyAudit body_audit = scan_body_cover(rebuilt);
        require(body_audit.failures == 0 &&
                    body_audit.fnv == UINT64_C(0x10c4c3ed46d44bf1),
                "rebuilt deck body scan failed");

        const std::set<Pair> typed_union = read_pair_set297(argv[4]);
        const std::size_t expected_union_count = std::stoull(argv[8]);
        require(typed_union.size() == expected_union_count,
                "old typed-union size changed");
        std::vector<Pair> overlap;
        Fnv overlap_ledger;
        for (Pair pair : ideal.rows)
            if (typed_union.contains(pair)) {
                overlap.push_back(pair);
                overlap_ledger.add(pair.q);
                overlap_ledger.add(pair.r);
            }
        require(overlap.empty(), "index-297 overlaps old typed union");

        const i128 global_divisor = gcd128(global_min_num, global_min_den);
        std::cout << "DETACHED_INDEX297_COMMON_DECK_AUDIT_V1\n"
                  << "IDEAL INDEX 297 ROWS " << ideal.rows.size()
                  << " FNV " << std::hex << ideal.fnv << std::dec
                  << " MAX_ENDPOINT " << ideal.maximum_endpoint << " ROWS";
        for (Pair pair : ideal.rows) std::cout << ' ' << pair.q << ',' << pair.r;
        std::cout << '\n'
                  << "OLD_MASK " << std::hex << std::setw(8)
                  << std::setfill('0') << deck[kTarget297] << std::dec
                  << std::setfill(' ') << " OBLIGATIONS "
                  << obligations.bodies.size() << " FNV " << std::hex
                  << obligations.fnv << " UNION " << obligations.body_union
                  << std::dec << " TARGET_DISJOINT "
                  << obligations.target_disjoint << " RETAINED_CHECKS "
                  << obligations.retained_checks << " BODIES";
        for (u32 body : obligations.bodies)
            std::cout << ' ' << std::hex << std::setw(8)
                      << std::setfill('0') << body << std::dec
                      << std::setfill(' ');
        std::cout << '\n'
                  << "ACTIVITY ROWS " << rows.size() << " REPAIRS_PER_ROW "
                  << kRepairCount << " ZETA_OPERATIONS " << operation_total
                  << " ACTIVITY_FNV " << std::hex << activity_ledger.state
                  << std::dec << '\n'
                  << "COMMON_ACTIVE " << common_count << " FNV " << std::hex
                  << common_ledger.state << " RESPONSE_FNV "
                  << response_ledger.state << std::dec << " CLASSES "
                  << classes.size() << " FULL_RESPONDERS 0\n";
        for (const auto& [pattern, cls] : classes)
            std::cout << "CLASS " << std::hex << pattern << " LEAST "
                      << std::setw(8) << std::setfill('0') << cls.least
                      << std::dec << std::setfill(' ') << " COVER "
                      << std::popcount(pattern) << " COUNT " << cls.count
                      << '\n';
        std::cout << "MAXIMAL_ANTICHAIN 5,e DUAL y0=1,y3=1 VALUE 2 "
                     "MAX_LOAD 1 EXACT_MINIMUM 2\n"
                  << "STRICT_REBUILT_CELLS " << strict_cells
                  << " EQUALITIES " << equalities << " FNV " << std::hex
                  << strict_ledger.state << " WITNESS_MARGIN_FNV "
                  << witness_margin_ledger.state << std::dec << '\n'
                  << "GLOBAL_WEAKEST MASK " << std::hex << std::setw(8)
                  << std::setfill('0') << global_weakest_mask << std::dec
                  << std::setfill(' ') << " PAIR " << global_weakest.q << ','
                  << global_weakest.r << " GAP "
                  << decimal(global_min_num / global_divisor) << '/'
                  << decimal(global_min_den / global_divisor) << '\n'
                  << "REBUILT_DECK " << rebuilt.size() << " FNV " << std::hex
                  << rebuilt_ledger.state << std::dec << " BODY_SCAN "
                  << kBodyCount << " CHECKS " << body_audit.checks
                  << " FAILURES " << body_audit.failures << " BODY_FNV "
                  << std::hex << body_audit.fnv << std::dec << '\n'
                  << "OLD_TYPED_UNION " << typed_union.size()
                  << " OVERLAP " << overlap.size() << " OVERLAP_FNV "
                  << std::hex << overlap_ledger.state << std::dec << '\n'
                  << "SCOPE DETACHED_LITERAL_FIXED_POOL_SEPARATE_422_MASK_DECK_"
                     "OLD_TYPED_UNION_OVERLAP_ONLY_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS EXACT_INDEX297_MINIMUM_TWO_COMMON_DECK\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "DETACHED_INDEX297_ERROR " << error.what() << '\n';
        return 1;
    }
}

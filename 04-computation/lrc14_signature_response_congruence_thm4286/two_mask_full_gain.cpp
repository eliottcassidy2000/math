#define RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY
#include "04-computation/lrc14_joint421_global_common_carrier_thm4281/response_pattern_atlas.cpp"
#undef RESPONSE_PATTERN_ATLAS_LIBRARY_ONLY

#include <fstream>
#include <set>
#include <sstream>

namespace {

struct Row {
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

std::vector<Row> read_rows(const std::string& path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open signatures");
    std::string line;
    require(static_cast<bool>(std::getline(in, line)) &&
                line == "q,r,inactive_count,w0,w1,w2,w3,w4,w5,w6",
            "signature header changed");
    std::vector<Row> rows;
    FnvLocal ledger;
    while (std::getline(in, line)) {
        const auto f = split_csv(line);
        require(f.size() == 10, "bad signature row");
        Row row;
        row.q = std::stoi(f[0]);
        row.r = std::stoi(f[1]);
        unsigned count = 0;
        for (int i = 0; i < 7; ++i) {
            row.words[i] = std::stoull(f[3 + i], nullptr, 16);
            count += std::popcount(row.words[i]);
        }
        require(count == static_cast<unsigned>(std::stoul(f[2])),
                "signature count changed");
        ledger.add(row.q); ledger.add(row.r);
        for (u64 word : row.words) ledger.add(word);
        rows.push_back(row);
    }
    require(rows.size() == 24223 &&
                ledger.state == UINT64_C(0x96507e0f0046a67a),
            "signature universe changed");
    std::cout << "SIGNATURE_ROWS " << rows.size() << " FNV " << std::hex
              << ledger.state << std::dec << '\n';
    return rows;
}

i128 atom_mass(const AtomData& atoms, u32 mask) {
    i128 mass = 0;
    for (const auto& [failure, value] : atoms.mass)
        if ((failure & ~mask) == 0) mass += value;
    return mass;
}

void add_i128(FnvLocal& ledger, i128 value) {
    const __uint128_t bits = static_cast<__uint128_t>(value);
    ledger.add(static_cast<u64>(bits));
    ledger.add(static_cast<u64>(bits >> 64));
}

bool fraction_less(i128 a, i128 b, i128 c, i128 d) {
    require(a >= 0 && b > 0 && c >= 0 && d > 0,
            "fraction comparison domain error");
    bool reversed = false;
    while (true) {
        const i128 qa = a / b;
        const i128 qc = c / d;
        if (qa != qc) return reversed ? qa > qc : qa < qc;
        a %= b;
        c %= d;
        if (a == 0 || c == 0) {
            if (a == 0 && c == 0) return false;
            const bool current_less = a == 0;
            return reversed ? !current_less : current_less;
        }
        std::swap(a, b);
        std::swap(c, d);
        reversed = !reversed;
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        require(argc == 4,
                "usage: full-gain SIGNATURES PRIOR_RESIDUAL MODE");
        const std::vector<Row> rows = read_rows(argv[1]);
        std::ifstream residual_in(argv[2]);
        require(static_cast<bool>(residual_in), "cannot open prior residual");
        std::set<std::pair<int, int>> prior_residual;
        FnvLocal residual_ledger;
        std::pair<int, int> previous{-1, -1};
        std::string residual_line;
        while (std::getline(residual_in, residual_line)) {
            const auto f = split_csv(residual_line);
            require(f.size() == 2, "bad prior residual row");
            const std::pair<int, int> pair{std::stoi(f[0]), std::stoi(f[1])};
            require(previous < pair, "prior residual not unique and ordered");
            prior_residual.insert(pair);
            residual_ledger.add(pair.first); residual_ledger.add(pair.second);
            previous = pair;
        }
        require(prior_residual.size() == 23373 &&
                    residual_ledger.state == UINT64_C(0xc6ab0ae49ee32273),
                "prior residual census changed");
        const std::string mode = argv[3];
        std::vector<std::size_t> deleted;
        std::vector<u32> appended;
        u64 expected_candidates = 0;
        u64 expected_candidate_fnv = 0;
        u64 expected_tests = 0;
        u64 expected_active_incidences = 0;
        u64 expected_activity_fnv = 0;
        u64 expected_common = 0;
        u64 expected_common_fnv = 0;
        u64 expected_net = 0;
        u64 expected_net_fnv = 0;
        int expected_minimum_q = 0, expected_minimum_r = 0;
        u32 expected_minimum_mask = 0;
        i128 expected_minimum_numerator = 0;
        i128 expected_minimum_denominator = 1;
        if (mode == "two-mask") {
            deleted = {107, 318, 374};
            appended = {UINT32_C(0x32043014), UINT32_C(0x20807016),
                        UINT32_C(0x128c8012)};
            expected_candidates = 112;
            expected_candidate_fnv = UINT64_C(0xd2fcb3bf714762ef);
            expected_tests = 336;
            expected_active_incidences = 330;
            expected_activity_fnv = UINT64_C(0x122dd5b2bb75e118);
            expected_common = 106;
            expected_common_fnv = UINT64_C(0x2dab535d17359dfb);
            expected_net = 36;
            expected_net_fnv = UINT64_C(0x8da60395e47e11a3);
            expected_minimum_q = 332; expected_minimum_r = 496;
            expected_minimum_mask = UINT32_C(0x128c8012);
            expected_minimum_numerator =
                static_cast<i128>(UINT64_C(77081685186123072));
            expected_minimum_denominator =
                static_cast<i128>(UINT64_C(10513328712007080960));
        } else if (mode == "index396") {
            deleted = {396};
            appended = {UINT32_C(0x042022c9)};
            expected_candidates = 36;
            expected_candidate_fnv = UINT64_C(0x3d92ab45b46a72c0);
            expected_tests = 36;
            expected_active_incidences = 36;
            expected_activity_fnv = UINT64_C(0xb25d97184aacd3bc);
            expected_common = 36;
            expected_common_fnv = UINT64_C(0x3d92ab45b46a72c0);
            expected_net = 36;
            expected_net_fnv = UINT64_C(0x3d92ab45b46a72c0);
            expected_minimum_q = 238; expected_minimum_r = 366;
            expected_minimum_mask = UINT32_C(0x042022c9);
            expected_minimum_numerator =
                static_cast<i128>(UINT64_C(179590440538495728));
            expected_minimum_denominator =
                static_cast<i128>(UINT64_C(11122656401155178880));
        } else {
            require(false, "unknown mode");
        }
        std::array<u64, 7> allowed{};
        for (std::size_t index : deleted)
            allowed[index / 64] |= UINT64_C(1) << (index % 64);
        const std::vector<Cell> cells = build_pool_cells();
        require(cells.size() == 7133, "pool cells changed");
        u64 candidates = 0, common = 0, net = 0, equalities = 0;
        u64 tests = 0, active_incidences = 0;
        FnvLocal candidate_ledger, common_ledger, net_ledger, activity_ledger;
        i128 minimum_positive = 0, minimum_grid = 1;
        int minimum_q = 0, minimum_r = 0;
        u32 minimum_mask = 0;
        bool minimum_set = false;
        std::cout << "THM4286_FULL_GAIN_V1 MODE " << mode << '\n';
        for (const Row& row : rows) {
            bool subset = true;
            for (int i = 0; i < 7; ++i)
                subset &= (row.words[i] & ~allowed[i]) == 0;
            if (!subset) continue;
            ++candidates;
            candidate_ledger.add(row.q); candidate_ledger.add(row.r);
            const i64 g = std::gcd(row.q, row.r);
            const PrimitivePair primitive = build_primitive(row.q / g,
                                                             row.r / g);
            const AtomData atoms = build_cocycle_atoms(cells, primitive, g);
            const i128 denominator =
                static_cast<i128>(primitive.grid) * g * COMMON;
            bool all_active = true;
            u32 activity_word = 0;
            std::vector<i128> margins(appended.size());
            for (std::size_t i = 0; i < appended.size(); ++i) {
                margins[i] = static_cast<i128>(63) *
                                 atom_mass(atoms, appended[i]) -
                             static_cast<i128>(4) * denominator;
                ++tests;
                equalities += margins[i] == 0;
                if (margins[i] >= 0) {
                    activity_word |= u32{1} << i;
                    ++active_incidences;
                } else {
                    all_active = false;
                }
            }
            activity_ledger.add(row.q); activity_ledger.add(row.r);
            activity_ledger.add(activity_word);
            for (i128 margin : margins) add_i128(activity_ledger, margin);
            add_i128(activity_ledger, denominator);
            std::cout << "CANDIDATE " << row.q << ',' << row.r
                      << " ACTIVITY " << std::hex << activity_word
                      << std::dec << " MARGINS";
            for (i128 margin : margins)
                std::cout << ' ' << decimal(margin);
            std::cout << " DEN " << decimal(denominator) << '\n';
            if (!all_active) continue;
            ++common;
            common_ledger.add(row.q); common_ledger.add(row.r);
            std::cout << "COMMON " << row.q << ',' << row.r << '\n';
            if (prior_residual.contains({row.q, row.r})) {
                ++net;
                net_ledger.add(row.q); net_ledger.add(row.r);
                std::cout << "NET " << row.q << ',' << row.r << '\n';
            }
            for (std::size_t i = 0; i < appended.size(); ++i) {
                if (!minimum_set || fraction_less(margins[i], denominator,
                                                  minimum_positive,
                                                  minimum_grid)) {
                    minimum_set = true;
                    minimum_positive = margins[i];
                    minimum_grid = denominator;
                    minimum_q = row.q; minimum_r = row.r;
                    minimum_mask = appended[i];
                }
            }
        }
        require(candidates == expected_candidates &&
                    candidate_ledger.state == expected_candidate_fnv &&
                    tests == expected_tests &&
                    active_incidences == expected_active_incidences &&
                    equalities == 0 &&
                    activity_ledger.state == expected_activity_fnv &&
                    common == expected_common &&
                    common_ledger.state == expected_common_fnv &&
                    net == expected_net &&
                    net_ledger.state == expected_net_fnv &&
                    minimum_q == expected_minimum_q &&
                    minimum_r == expected_minimum_r &&
                    minimum_mask == expected_minimum_mask &&
                    minimum_positive == expected_minimum_numerator &&
                    minimum_grid == expected_minimum_denominator,
                "full-gain summary changed");
        std::cout << "SUMMARY CANDIDATES " << candidates
                  << " CANDIDATE_FNV " << std::hex << candidate_ledger.state
                  << std::dec << " TESTS " << tests
                  << " ACTIVE_INCIDENCES " << active_incidences
                  << " EQUALITIES " << equalities
                  << " ACTIVITY_FNV " << std::hex << activity_ledger.state
                  << std::dec << " COMMON " << common << " COMMON_FNV "
                  << std::hex << common_ledger.state << std::dec << '\n'
                  << "NET_FROM_PRIOR_RESIDUAL " << net << " NET_FNV "
                  << std::hex << net_ledger.state << std::dec
                  << " NEW_RESIDUAL " << prior_residual.size() - net << '\n'
                  << "MIN_COMMON_GAP " << minimum_q << ',' << minimum_r
                  << " MASK " << std::hex << minimum_mask << std::dec
                  << " NUM " << decimal(minimum_positive) << " DEN "
                  << decimal(minimum_grid) << '\n'
                  << "VERDICT PASS EXACT_FULL_SIGNATURE_SUBSET_ACTIVITY\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FULL_GAIN_ERROR " << error.what() << '\n';
        return 1;
    }
}

#include "04-computation/lrc14_endpoint_carrier_signature_surgery_thm4283/carrier_scan_support.cpp"

#include <array>
#include <fstream>
#include <iomanip>

namespace {

constexpr std::size_t N176 = 176;
constexpr std::size_t WORDS176 = 3;

struct Failure176 {
    int q = 0;
    int r = 0;
    u32 body = 0;
    unsigned row = 0;
};

struct Pattern176 {
    std::array<u64, WORDS176> w{};
    u32 least = 0;
    u64 multiplicity = 1;
};

unsigned row_index176(int q, int r) {
    if (q == 100 && r == 636) return 0;
    if (q == 256 && r == 636) return 1;
    if (q == 100 && r == 632) return 2;
    if (q == 256 && r == 632) return 3;
    require(false, "unexpected failure row");
    return 0;
}

std::vector<Failure176> read_failures176(const char* path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open failure ledger");
    std::string line;
    std::getline(in, line);
    require(line == "q,r,body_hex", "failure header changed");
    std::vector<Failure176> out;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const std::size_t c1 = line.find(',');
        const std::size_t c2 = line.find(',', c1 + 1);
        require(c1 != std::string::npos && c2 != std::string::npos &&
                    line.find(',', c2 + 1) == std::string::npos,
                "malformed failure row");
        Failure176 f;
        f.q = std::stoi(line.substr(0, c1));
        f.r = std::stoi(line.substr(c1 + 1, c2 - c1 - 1));
        f.body = static_cast<u32>(std::stoul(line.substr(c2 + 1), nullptr, 16));
        f.row = row_index176(f.q, f.r);
        require(std::popcount(f.body) == 9, "failure body rank changed");
        out.push_back(f);
    }
    require(out.size() == N176, "failure count changed");
    const std::array<u64, 4> expected{64, 37, 8, 67};
    std::array<u64, 4> counts{};
    for (const Failure176& f : out) ++counts[f.row];
    require(counts == expected, "failure row counts changed");
    return out;
}

bool pattern_order176(const Pattern176& a, const Pattern176& b) {
    for (std::size_t i = WORDS176; i-- > 0;)
        if (a.w[i] != b.w[i]) return a.w[i] < b.w[i];
    return a.least < b.least;
}

bool same_pattern176(const Pattern176& a, const Pattern176& b) {
    return a.w == b.w;
}

u64 cover176(const Pattern176& p) {
    return std::popcount(p.w[0]) + std::popcount(p.w[1]) +
           std::popcount(p.w[2]);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 3, "usage: PROG failures176.csv classes.tsv");
        init_choose8_local();
        const std::vector<Failure176> failures = read_failures176(argv[1]);
        FnvLocal failure_ledger;
        for (const Failure176& f : failures) {
            failure_ledger.add(f.q);
            failure_ledger.add(f.r);
            failure_ledger.add(f.body);
        }
        const std::vector<Cell> cells = build_pool_cells();
        std::array<ActiveUniverse, 4> active{
            build_active_universe(cells, 100, 636),
            build_active_universe(cells, 256, 636),
            build_active_universe(cells, 100, 632),
            build_active_universe(cells, 256, 632)};
        std::cout << "LRC14_176_ROWAWARE_RESPONSE_V1\n"
                  << "FAILURES " << failures.size() << " LABELLED_FNV "
                  << std::hex << failure_ledger.state << std::dec << '\n';
        constexpr std::array<std::pair<int,int>,4> pairs{{
            {100,636},{256,636},{100,632},{256,632}}};
        for (std::size_t i = 0; i < pairs.size(); ++i)
            std::cout << "ACTIVE_ROW " << pairs[i].first << ','
                      << pairs[i].second << " COUNT " << active[i].count
                      << " FNV " << std::hex << active[i].fnv << std::dec
                      << '\n';

        std::vector<Pattern176> raw;
        raw.reserve(100000);
        u64 empty = 0;
        u64 incidences = 0;
        u64 max_cover = 0;
        u32 least_max = 0;
        FnvLocal nonempty_ledger;
        for (u64 rank = 0; rank < EXPECTED_REPAIRS; ++rank) {
            const u32 mask = unrank_colex8(rank);
            Pattern176 p;
            p.least = mask;
            for (std::size_t i = 0; i < failures.size(); ++i) {
                const Failure176& f = failures[i];
                if (!active[f.row].active[rank] || (mask & f.body) != 0)
                    continue;
                p.w[i / 64] |= UINT64_C(1) << (i % 64);
            }
            const u64 c = cover176(p);
            if (c == 0) {
                ++empty;
                continue;
            }
            incidences += c;
            nonempty_ledger.add(mask);
            for (u64 word : p.w) nonempty_ledger.add(word);
            if (c > max_cover || (c == max_cover &&
                                  (least_max == 0 || mask < least_max))) {
                max_cover = c;
                least_max = mask;
            }
            raw.push_back(p);
        }
        std::sort(raw.begin(), raw.end(), pattern_order176);
        std::vector<Pattern176> classes;
        for (const Pattern176& p : raw) {
            if (!classes.empty() && same_pattern176(classes.back(), p)) {
                ++classes.back().multiplicity;
                classes.back().least = std::min(classes.back().least, p.least);
            } else {
                classes.push_back(p);
            }
        }
        FnvLocal class_ledger;
        std::ofstream out(argv[2]);
        require(static_cast<bool>(out), "cannot create class atlas");
        out << "class\tleast_mask_hex\tmultiplicity\tcover\tw2\tw1\tw0\n";
        for (std::size_t i = 0; i < classes.size(); ++i) {
            const Pattern176& p = classes[i];
            class_ledger.add(p.least);
            class_ledger.add(p.multiplicity);
            class_ledger.add(cover176(p));
            for (u64 word : p.w) class_ledger.add(word);
            out << i << '\t' << std::hex << std::setw(8) << std::setfill('0')
                << p.least << std::dec << '\t' << p.multiplicity << '\t'
                << cover176(p) << '\t' << std::hex << std::setw(16)
                << p.w[2] << '\t' << std::setw(16) << p.w[1] << '\t'
                << std::setw(16) << p.w[0] << std::dec << '\n';
        }
        require(out.good(), "atlas write failed");
        std::cout << "NONEMPTY_MASKS " << raw.size() << " EMPTY_MASKS "
                  << empty << " RESPONSE_INCIDENCES " << incidences
                  << " RESPONSE_FNV " << std::hex << nonempty_ledger.state
                  << std::dec << '\n'
                  << "NONEMPTY_CLASSES " << classes.size()
                  << " ALL_CLASSES " << classes.size() + 1 << " CLASS_FNV "
                  << std::hex << class_ledger.state << std::dec
                  << " MAX_COVER " << max_cover << " LEAST_MAX " << std::hex
                  << std::setw(8) << std::setfill('0') << least_max << std::dec
                  << std::setfill(' ') << '\n'
                  << "VERDICT PASS FINITE_EXACT_ROWAWARE_RESPONSE_QUOTIENT\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "LRC14_176_ROWAWARE_RESPONSE_ERROR " << e.what() << '\n';
        return 1;
    }
}

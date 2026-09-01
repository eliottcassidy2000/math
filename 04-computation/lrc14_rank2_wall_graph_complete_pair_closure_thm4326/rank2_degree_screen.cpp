// Scratch degree-bound screen for the rank-two literal graph.

#define main rank2_allbody_hidden_main
#include "rank2_allbody_audit.cpp"
#undef main

#include <filesystem>
#include <fstream>
#include <set>

namespace {
struct Fnv {
    u64 state = UINT64_C(0xcbf29ce484222325);
    void add(u64 word) {
        for (unsigned byte = 0; byte < 8; ++byte) {
            state ^= (word >> (8 * byte)) & UINT64_C(0xff);
            state *= UINT64_C(0x100000001b3);
        }
    }
};

void add_i128(Fnv& ledger, i128 value) {
    const u128 bits = static_cast<u128>(value);
    ledger.add(static_cast<u64>(bits));
    ledger.add(static_cast<u64>(bits >> 64));
}

std::vector<std::pair<int, int>> read_pairs(const std::filesystem::path& path) {
    std::ifstream input(path);
    require(static_cast<bool>(input), "cannot open pair ledger");
    std::vector<std::pair<int, int>> pairs;
    std::set<std::pair<int, int>> distinct;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::size_t comma = line.find(',');
        require(comma != std::string::npos &&
                    line.find(',', comma + 1) == std::string::npos,
                "pair line must be q,r without header");
        const int q = std::stoi(line.substr(0, comma));
        const int r = std::stoi(line.substr(comma + 1));
        require(q > 0 && q < r && distinct.insert({q, r}).second,
                "pair ordering/distinctness failed");
        pairs.push_back({q, r});
    }
    require(input.eof() && !pairs.empty(), "empty/incomplete pair ledger");
    std::sort(pairs.begin(), pairs.end(), [](auto left, auto right) {
        if (left.second != right.second) return left.second > right.second;
        return left.first < right.first;
    });
    return pairs;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        require(argc == 3, "usage: degree_screen PAIRS_QR_CSV OUTPUT_CSV");
        const auto pairs = read_pairs(argv[1]);
        std::ofstream output(argv[2], std::ios::binary);
        require(static_cast<bool>(output), "cannot create screen ledger");
        output << "q,r,grid,rank2_total,degree_bound_mass,degree_bound_ticks,positive,top9_hex\n";
        Fnv ledger;
        u64 positive = 0;
        int highest_bad_r = -1;
        u64 bad_at_highest = 0;
        for (const auto [q, r] : pairs) {
            const Graph graph = build_graph(q, r);
            std::array<std::pair<i64, unsigned>, 30> degrees{};
            for (unsigned vertex = 0; vertex < 30; ++vertex)
                degrees[vertex] = {graph.incident[vertex], vertex};
            std::sort(degrees.begin(), degrees.end(), [](auto left, auto right) {
                if (left.first != right.first) return left.first > right.first;
                return left.second < right.second;
            });
            i64 bound_mass = graph.total;
            u32 top9 = 0;
            for (unsigned index = 0; index < 9; ++index) {
                bound_mass -= degrees[index].first;
                top9 |= u32{1} << degrees[index].second;
            }
            const i128 ticks = static_cast<i128>(63) * bound_mass -
                               static_cast<i128>(4) * graph.grid;
            const bool ok = ticks > 0;
            positive += ok;
            if (!ok) {
                if (highest_bad_r == -1) highest_bad_r = r;
                if (r == highest_bad_r) ++bad_at_highest;
            }
            output << q << ',' << r << ',' << decimal(graph.grid) << ','
                   << decimal(graph.total) << ',' << decimal(bound_mass) << ','
                   << decimal(ticks) << ',' << ok << ',' << hex8(top9) << '\n';
            ledger.add(q); ledger.add(r); add_i128(ledger, graph.grid);
            add_i128(ledger, graph.total); add_i128(ledger, bound_mass);
            add_i128(ledger, ticks); ledger.add(ok); ledger.add(top9);
        }
        require(output.good(), "screen ledger write failed");
        std::cout << "LRC14_RANK2_DEGREE_SCREEN_V1\n"
                  << "PAIRS " << pairs.size() << " POSITIVE " << positive
                  << " NONPOSITIVE " << pairs.size() - positive
                  << " HIGHEST_NONPOSITIVE_ENDPOINT " << highest_bad_r
                  << " BAD_AT_HIGHEST " << bad_at_highest << " LEDGER_FNV "
                  << std::hex << ledger.state << std::dec << '\n'
                  << "BOUND L2(B)>=TOTAL-SUM_NINE_LARGEST_WEIGHTED_DEGREES\n"
                  << "SCOPE FINITE_EXACT_DEGREE_LOWER_BOUND_NO_CLAIM_FOR_NONPOSITIVE_ROWS_NO_PHYSICAL_ENTRY_NO_LRC14\n"
                  << "VERDICT PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "RANK2_DEGREE_SCREEN_ERROR " << error.what() << '\n';
        return 1;
    }
}

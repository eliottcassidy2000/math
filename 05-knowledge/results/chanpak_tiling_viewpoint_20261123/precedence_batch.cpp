// Read the default upper-triangle ASCII tournament stream of nauty gentourng.
// Build and run from this directory:
//   c++ -std=c++17 -O3 precedence_batch.cpp -o /tmp/precedence-batch
//   gentourng -q 9  | /tmp/precedence-batch 9
//   gentourng -q 10 | /tmp/precedence-batch 10
// Each line lists the pairs i<j in lexicographic order, with 1 meaning i->j.
// A machine word evaluates the same Boolean parity computation in 64 boards.
// Every retained class is checked independently by enumerating all Hamilton words.

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
using Batch = std::uint64_t;
using Histogram = std::array<int, 12>;
constexpr int max_order = 10;
constexpr int max_subsets = 1 << max_order;
constexpr int width = 64;

const Histogram all_nine = {48, 0, 38, 4, 48, 0, 35, 2, 31, 0, 43, 4};
const Histogram path_nine = {23, 0, 0, 0, 26, 0, 3, 0, 13, 0, 4, 0};
const Histogram all_ten = {186, 2, 123, 8, 170, 0, 131, 4, 134, 2, 141, 10};
const Histogram path_ten = {74, 0, 0, 0, 76, 0, 16, 0, 59, 0, 10, 0};

int order, subset_count;
std::array<std::string, width> records;
Batch arc[max_order][max_order], odd_before[max_order][max_order];
Batch ends[max_subsets][max_order], cycles[max_subsets];
Histogram all{}, is_path{};
int out[max_order], preceding[max_order];

int lowest_vertex(int nonempty) {
    return __builtin_ctz(static_cast<unsigned>(nonempty));
}

int hamilton_words(int used, int next) {
    int total = 0;
    while (next) {
        int v = lowest_vertex(next);
        next &= next - 1;
        int enlarged = used | (1 << v);
        int count = enlarged == subset_count - 1 ? 1 : hamilton_words(enlarged, out[v] & ~enlarged);
        // Exactly these 'count' completions place all of 'used' before v here.
        if (count & 1) preceding[v] ^= used;
        total += count;
    }
    return total;
}

void check_class(int lane) {
    std::fill_n(out, order, 0);
    std::fill_n(preceding, order, 0);
    int pos = 0;
    for (int i = 0; i < order; ++i) {
        for (int j = i+1; j < order; ++j) {
            if (records[lane][pos++] == '1') out[i] |= 1 << j;
            else out[j] |= 1 << i;
        }
    }
    int height = hamilton_words(0, subset_count - 1);
    int ranks = 0;
    std::array<int, max_order> sorted{};
    for (int i = 0; i < order; ++i) {
        for (int j = i+1; j < order; ++j) {
            if (((preceding[j] >> i) & 1) != ((odd_before[i][j] >> lane) & 1)) {
                throw std::runtime_error("cycle and word precedence disagree");
            }
        }
        int rank = __builtin_popcount(static_cast<unsigned>(preceding[i]));
        ranks |= 1 << rank;
        sorted[rank] = i;
    }
    if (!(height & 1) || ranks != subset_count - 1) throw std::runtime_error("invalid retained order");
    bool path = true;
    for (int k = 1; k < order; ++k) {
        if (!((out[sorted[k-1]] >> sorted[k]) & 1)) path = false;
    }
    int bin = (height % 24) / 2;
    ++all[bin];
    is_path[bin] += path;
}

void process(int size) {
    Batch active = size == width ? ~Batch(0) : (Batch(1) << size) - 1;
    int pos = 0;
    for (int i = 0; i < order; ++i) {
        for (int j = i+1; j < order; ++j) {
            Batch forward = 0;
            for (int k = 0; k < size; ++k) forward |= Batch(records[k][pos] == '1') << k;
            ++pos;
            arc[i][j] = forward;
            arc[j][i] = active ^ forward;
        }
    }

    // ends[S][j] counts Hamilton paths from min(S) to j in S, modulo two.
    // Closing such a path at the root counts each directed cycle exactly once.
    cycles[0] = 0;
    for (int set = 1; set < subset_count; ++set) {
        int root = lowest_vertex(set);
        int nonroots = set & (set - 1);
        if (!nonroots) ends[set][root] = active;
        Batch closed = 0;
        for (int todo = nonroots; todo; todo &= todo - 1) {
            int j = lowest_vertex(todo);
            int smaller = set ^ (1 << j);
            int previous = smaller & (smaller - 1);
            Batch paths = previous ? 0 : arc[root][j];
            for (; previous; previous &= previous - 1) {
                int v = lowest_vertex(previous);
                paths ^= ends[smaller][v] & arc[v][j];
            }
            ends[set][j] = paths;
            closed ^= paths & arc[j][root];
        }
        cycles[set] = closed;
    }

    // Each surviving cycle contributes to every pair in its vertex support.
    for (int bit = 1; bit < subset_count; bit *= 2) {
        for (int start = 0; start < subset_count; start += 2*bit) {
            for (int j = 0; j < bit; ++j) cycles[start+j] ^= cycles[start+j+bit];
        }
    }
    for (int i = 0; i < order; ++i) {
        for (int j = i+1; j < order; ++j) {
            odd_before[i][j] = arc[i][j] ^ cycles[(1 << i) | (1 << j)];
        }
    }

    // For i<j<k, a directed triangle has R_ij = R_jk different from R_ik.
    for (int i = 0; i < order && active; ++i) {
        for (int j = i+1; j < order && active; ++j) {
            for (int k = j+1; k < order; ++k) {
                Batch first = odd_before[i][j], second = odd_before[j][k], outer = odd_before[i][k];
                active &= ~((first ^ outer) & ~(first ^ second));
            }
        }
    }
    while (active) {
        int lane = __builtin_ctzll(active);
        active &= active - 1;
        check_class(lane);
    }
}

void print(const char* label, const Histogram& counts) {
    std::cout << label;
    for (int count : counts) std::cout << ' ' << count;
    std::cout << '\n';
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 2 || (std::string(argv[1]) != "9" && std::string(argv[1]) != "10")) {
        std::cerr << "Usage: gentourng -q N | precedence-batch N  (N=9 or 10)\n";
        return 1;
    }
    order = std::stoi(argv[1]);
    subset_count = 1 << order;
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    int size = 0, total = 0;
    while (std::getline(std::cin, records[size])) {
        const auto& row = records[size];
        if (row.size() != static_cast<std::size_t>(order*(order-1)/2)
            || row.find_first_not_of("01") != std::string::npos) {
            std::cerr << "Invalid tournament at record " << total+1 << '\n';
            return 2;
        }
        ++total;
        if (++size == width) {
            process(size);
            size = 0;
        }
    }
    if (size) process(size);
    std::cout << "ordinary classes " << total << '\n';
    print("transitive", all);
    print("transitive path", is_path);
    if (total != (order == 9 ? 191536 : 9733056)
        || all != (order == 9 ? all_nine : all_ten)
        || is_path != (order == 9 ? path_nine : path_ten)) {
        std::cerr << "Incorrect stream size or histogram\n";
        return 3;
    }
}

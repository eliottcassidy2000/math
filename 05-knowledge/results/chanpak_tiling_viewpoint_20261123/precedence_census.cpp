// Input: gentourng -q n, whose default records give the upper triangle by rows.
#include <array>
#include <cstdlib>
#include <iostream>
#include <string>

using Mask = unsigned;
constexpr int max_order = 10;

struct Precedence {
    int n;
    std::array<Mask, max_order> rows{};
    std::array<Mask, max_order> predecessors{};

    int extend(Mask used, int last) {
        Mask remaining = ((1u << n) - 1) ^ used;
        if (!remaining) return 1;
        Mask choices = last < 0 ? remaining : remaining & rows[last];
        int paths = 0;
        while (choices) {
            Mask bit = choices & -choices;
            choices ^= bit;
            int vertex = __builtin_ctz(bit);
            int completions = extend(used | bit, vertex);
            if (completions % 2) predecessors[vertex] ^= used;
            paths += completions;
        }
        return paths;
    }

    bool transitive() const {
        Mask ranks = 0;
        for (int v = 0; v < n; ++v) ranks |= 1u << __builtin_popcount(predecessors[v]);
        return ranks == (1u << n) - 1;
    }
};

int main(int argc, char** argv) {
    int n = argc == 2 ? std::atoi(argv[1]) : 0;
    if (n < 2 || n > max_order) return 2;
    std::array<int,24> residues{};
    int total = 0;
    std::string record;
    while (std::cin >> record) {
        if (record.size() != static_cast<size_t>(n*(n-1)/2) ||
            record.find_first_not_of("01") != std::string::npos) return 1;
        Precedence t{n};
        int k = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = i+1; j < n; ++j) {
                if (record[k++] == '1') t.rows[i] |= 1u << j;
                else t.rows[j] |= 1u << i;
            }
        }
        ++total;
        int height = t.extend(0, -1);
        if (t.transitive()) ++residues[height % 24];
    }
    std::cout << "order " << n << " ordinary classes " << total << " bins 1,3,...,23:";
    for (int r = 1; r < 24; r+=2) std::cout << ' ' << residues[r];
    std::cout << '\n';
    if (n == 9) {
        std::array<int,24> expected{};
        int counts[] = {48,0,38,4,48,0,35,2,31,0,43,4};
        for (int r = 1; r < 24; r+=2) expected[r] = counts[(r-1)/2];
        if (total != 191536 || residues != expected) return 1;
    }
}

// Read gentourng's default upper-triangle strings and test H | 32767, Q2/H dyadic.
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>

int main(int argc, char** argv) {
    if (argc != 2) return 2;
    int n = std::stoi(argv[1]);
    if (n < 2 || n > 9) return 2;
    int full = (1 << n) - 1;
    int total = 0, hits = 0;
    std::map<int,int> candidates;
    std::string encoding;
    while (std::cin >> encoding) {
        if (static_cast<int>(encoding.size()) != n*(n-1)/2) return 2;
        std::array<int,9> outgoing{};
        int position = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = i+1; j < n; ++j) {
                char bit = encoding[position++];
                if (bit == '1') outgoing[i] |= 1 << j;
                else if (bit == '0') outgoing[j] |= 1 << i;
                else return 2;
            }
        }
        ++total;
        std::array<std::array<int,9>,512> starts;
        std::array<int,512> paths;
        paths[0] = 1;
        for (int set = 1; set <= full; ++set) {
            paths[set] = 0;
            for (int mask = set; mask; mask &= mask-1) {
                int i = __builtin_ctz(mask), rest = set ^ (1 << i);
                int count = (rest == 0);
                for (int next = rest & outgoing[i]; next; next &= next-1)
                    count += starts[rest][__builtin_ctz(next)];
                starts[set][i] = count;
                paths[set] += count;
            }
        }
        int h = paths[full];
        if (h == 1 || 32767 % h != 0) continue;
        ++candidates[h];
        std::int64_t q = 0;
        for (int set = 0; set <= full; ++set)
            q += std::int64_t(paths[set]) * paths[full ^ set];
        if (q % h == 0 && ((q/h) & ((q/h)-1)) == 0) {
            ++hits;
            std::cout << "dyadic arm " << encoding << " H " << h << " Q2 " << q << '\n';
        }
    }
    const int census[] = {0,1,1,2,4,12,56,456,6880,191536};
    std::cout << "order " << n << ", isotypes " << total << ", hits " << hits << '\n';
    for (const auto& [h, count] : candidates) std::cout << h << ": " << count << '\n';
    return (total == census[n] && hits == 0) ? 0 : 1;
}

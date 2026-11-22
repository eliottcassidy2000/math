// Exact synchronized spacetime-profile census for elementary cellular automata.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using U = unsigned __int128;
using V = std::uint64_t;
U number(const std::string &s) {
    U n = 0;
    for (char c : s)
        n = 10 * n + (c - '0');
    return n;
}
std::string digits(U n) {
    if (!n)
        return "0";
    std::string s;
    while (n) {
        s.push_back('0' + n % 10);
        n /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}
U f(V x) {
    V z = 0;
    while (x) {
        V h = x / 2;
        if (!(x & 1))
            z += h + (h & 1);
        x = h;
    }
    return z;
}
struct Scan {
    int n, length;
    U h, q, c, t, e, r, k, N;
    V programs = 0, inside = 0, matches = 0;
    void test(V x, V y) {
        ++programs;
        U lower = k + q * x * y + c * y + t * x + e * (y & 1) + r * (x & 1);
        if (lower <= N && N - lower <= h * (x + y)) {
            ++inside;
            if (lower + h * (f(x) + f(y)) == N) {
                ++matches;
                if(matches==1) std::cout << "hit n " << n << " words " << x << " " << y << "\n";
            }
        }
    }
};
int main(int argc, char **argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "usage: farey_elementary max_ring_width [lower] < charts.txt\n";
        return 2;
    }
    int max_width = std::stoi(argv[1]);
    if (max_width < 1 || max_width > 16)
        return 2;
    U target = argc == 3 ? number(argv[2]) : (U(1) << 60) * ((U(1) << 61) - 1) - 1;
    std::vector<Scan> scans;
    int n, steps = 0;
    std::array<std::string, 7> data;
    while (std::cin >> n >> data[0] >> data[1] >> data[2] >> data[3] >> data[4] >> data[5] >>
           data[6]) {
        if (n > 16)
            continue;
        Scan s{n,
               0,
               number(data[0]),
               number(data[1]),
               number(data[2]),
               number(data[3]),
               number(data[4]),
               number(data[5]),
               number(data[6]),
               target};
        for (int l = 1; l < 62; ++l) {
            V low = V(1) << (l - 1), high = 2 * low - 1;
            if (s.k + s.q * low * low + (s.t + s.c) * low > target)
                break;
            if (s.k + s.q * high * high + (s.t + s.c + 2 * s.h) * high + s.e + s.r >= target) {
                s.length = l;
                scans.push_back(s);
                steps = std::max(steps, l - 1);
            }
        }
    }
    std::vector<std::vector<unsigned>> at(steps + 1);
    for (unsigned i = 0; i < scans.size(); ++i)
        at[scans[i].length - 1].push_back(i);
    for (int w = 1; w <= max_width; ++w) {
        unsigned total = 1u << w, mask = total - 1;
        std::vector<unsigned> next(total);
        for (unsigned rule = 0; rule < 256; ++rule) {
            for (unsigned s = 0; s < total; ++s) {
                unsigned a = (s << 1 & mask) | (s >> (w - 1)), b = s,
                         c = (s >> 1) | ((s & 1) << (w - 1));
                unsigned value = 0;
                for (unsigned j = 0; j < 8; ++j)
                    if (rule >> j & 1) {
                        value |= (j & 4 ? a : ~a) & (j & 2 ? b : ~b) & (j & 1 ? c : ~c);
                    }
                next[s] = value & mask;
            }
            for (unsigned seed = 0; seed < total; ++seed) {
                unsigned s = seed;
                std::array<V, 16> tape;
                tape.fill(1);
                for (int i = 0; i <= steps; ++i) {
                    for (unsigned index : at[i])
                        for (int d = 0; d < w; ++d)
                            scans[index].test(tape[0], tape[d]);
                    if (i == steps)
                        break;
                    for (int d = 0; d < w; ++d)
                        tape[d] = 2 * tape[d] + ((s >> d) & 1);
                    s = next[s];
                }
            }
        }
    }
    for (const auto &s : scans)
        std::cout << "n " << s.n << " length " << s.length << " programs " << s.programs
                  << " in envelope " << s.inside << " hits " << s.matches << "\n";
}

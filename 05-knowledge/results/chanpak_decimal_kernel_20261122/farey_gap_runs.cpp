// Exact diagonal-line inversion for address differences with exactly three one-runs.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
using U = unsigned __int128;
using V = uint64_t;

U decimal(const std::string &s) {
    U n = 0;
    for (char c : s)
        n = 10 * n + c - '0';
    return n;
}
std::string digits(U n) {
    std::string s;
    do {
        s.push_back('0' + n % 10);
        n /= 10;
    } while (n);
    std::reverse(s.begin(), s.end());
    return s;
}
U reading(U n) {
    U r = 0;
    while (n) {
        bool even = (n & 1) == 0;
        n >>= 1;
        if (even)
            r += n + (n & 1);
    }
    return r;
}

struct Diagonal {
    U h, q, c, t, e, r, k, target, max_gap;
    V tested = 0, inside = 0, local = 0, hits = 0;
    void check(U gap) {
        ++tested;
        U leftover = target - k - t * gap;
        U slope = q * gap + t + c;
        double s = double(slope), amount = double(leftover);
        V step = V(2 * amount / (s + std::sqrt(s * s + 4 * double(q) * amount)));
        auto value = [&](V u) -> U { return U(u) * (slope + q * u); };
        while (step && value(step) > leftover)
            --step;
        while (value(step + 1) <= leftover)
            ++step;
        U x = gap + step, y = step, used = value(step) + r * (x & 1) + e * (y & 1);
        if (used > leftover || leftover - used > h * (x + y))
            return;
        ++inside;
        U charge = leftover - used;
        if (charge % (2 * h))
            return;
        ++local;
        if (charge == h * (reading(x) + reading(y))) {
            ++hits;
            std::cout << "hit in scan coordinates " << digits(x) << " " << digits(y) << "\n";
        }
    }
    void grow(int first, int remaining, U gap) {
        for (int bottom = first; (U(1) << bottom) <= max_gap - gap; ++bottom) {
            U lower = U(1) << bottom;
            for (int top = bottom + 1;; ++top) {
                U next = gap + (U(1) << top) - lower;
                if (next > max_gap)
                    break;
                if (remaining == 1)
                    check(next);
                else
                    grow(top + 1, remaining - 1, next);
            }
        }
    }
};
int main(int argc, char **argv) {
    if (argc < 3 || argc > 4)
        return 2;
    int selected = std::stoi(argv[1]), sign = std::stoi(argv[2]), n;
    if ((selected != 4 && selected != 8) || (sign != 0 && sign != 1))
        return 2;
    std::string h, q, c, t, e, r, k;
    while (std::cin >> n >> h >> q >> c >> t >> e >> r >> k) {
        if (n != selected)
            continue;
        U target = argc == 4 ? decimal(argv[3]) : (U(1) << 60) * ((U(1) << 61) - 1) - 1;
        Diagonal s{decimal(h), decimal(q), decimal(c), decimal(t), decimal(e),
                   decimal(r), decimal(k), target,     0};
        if (target < s.k)
            return 0;
        if (sign) {
            std::swap(s.c, s.t);
            std::swap(s.e, s.r);
        }
        s.max_gap = (target - s.k) / s.t;
        if (s.max_gap >= (U(1) << 111) || target >= (U(1) << 122))
            throw std::runtime_error("exceeds the supported E61 range");
        s.grow(0, 3, 0);
        std::cout << n << " sign " << sign << " tested " << s.tested << " envelope " << s.inside
                  << " local " << s.local << " hits " << s.hits << std::endl;
        return 0;
    }
    return 2;
}

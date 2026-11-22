// Unrestricted two-address search on a Farey chart. Input rows: n h q c t e r k.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

static_assert(std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits == 53);
using U = unsigned __int128;
using V = uint64_t;
using S = __int128;
const U TARGET = (U(1) << 60) * ((U(1) << 61) - 1) - 1;

U parse(std::string s) {
    U v = 0;
    for (char d : s)
        v = 10 * v + (d - '0');
    return v;
}
std::string show(U n) {
    if (n == 0)
        return "0";
    std::string s;
    while (n) {
        s.push_back('0' + n % 10);
        n /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}
struct ByteTables {
    std::array<unsigned, 256> value{}, reverse{};
    ByteTables() {
        for (unsigned i = 1; i < 256; ++i) {
            value[i] = value[i / 2] + ((i & 1) ? 0 : i / 2 + ((i / 2) & 1));
            reverse[i] = (reverse[i / 2] >> 1) | ((i & 1) << 7);
        }
    }
} const bytes;
U reading(U n) {
    U sum = 0;
    while (n > 255) {
        unsigned low = n & 255;
        n >>= 8;
        sum += bytes.value[low] + (255 - bytes.reverse[low]) * n + ((low < 128) ? (n & 1) : 0);
    }
    return sum + bytes.value[V(n)];
}
U ceildiv(U a, U b) { return a / b + (a % b != 0); }
V inverse(V a, V m) {
    if (m == 1)
        return 0;
    S x = 1, y = 0;
    V b = m;
    while (b) {
        V d = a / b, t = a % b;
        a = b;
        b = t;
        S z = x - S(d) * y;
        x = y;
        y = z;
    }
    if (a != 1)
        throw std::runtime_error("nonunit");
    return (x % S(m) + m) % m;
}
struct Chart {
    V h;
    U q, c, t, e, r, k;
};
struct Scan {
    Chart original, a;
    U N;
    bool swapped;
    V sx, ox, sy, oy;
    Scan(Chart o, U n, bool flip, V step_x, V offset_x, V step_y, V offset_y)
        : original(o), a(o), N(n), swapped(flip), sx(step_x), ox(offset_x), sy(step_y),
          oy(offset_y) {
        a.k = o.k + o.t * ox + o.c * oy + o.q * ox * oy;
        if (sx % 2 == 0) {
            a.k += o.r * (ox & 1);
            a.r = 0;
        }
        if (sy % 2 == 0) {
            a.k += o.e * (oy & 1);
            a.e = 0;
        }
        a.c = sy * (o.c + o.q * ox);
        a.t = sx * (o.t + o.q * oy);
        a.q = o.q * sx * sy;
    }
    V exact_rows = 0, candidates = 0, hits = 0;
    void verify(V x, U y, U s, U rem) {
        if (y % sy != oy)
            return;
        ++candidates;
        if (s * y + a.h * reading(y) == rem) {
            ++hits;
            std::cout << "HIT " << show(swapped ? y : x) << " " << show(swapped ? x : y) << "\n";
        }
    }
    // Recover a long opposite flank from congruence and a reversed fixed low word.
    void solve(V coord, int requested_parity = -1) {
        V x = coord * sx + ox;
        ++exact_rows;
        const U s = original.c + original.q * x,
                base = original.k + original.t * x + original.r * (x & 1) + a.h * reading(x);
        if (base > N)
            return;
        for (V j = 0; j < 2; ++j) {
            if ((sy % 2 == 0 && j != (oy & 1)) ||
                (sy == 1 && requested_parity >= 0 && j != V(requested_parity)) ||
                base + j * original.e > N)
                continue;
            U rem = N - base - j * original.e;
            U lo = ceildiv(rem, s + a.h), hi = rem / s;
            lo = std::max(lo, U(sy) * (coord + swapped) + oy);
            if (lo > hi)
                continue;
            if (hi - lo <= 10) {
                for (U y = lo + ((lo & 1) != j); y <= hi; y += 2)
                    verify(x, y, s, rem);
                continue;
            }
            V g = std::gcd(V(s % a.h), a.h);
            if (rem % g)
                continue;
            V m = a.h / g;
            V residue = m == 1 ? 0 : U(V((rem / g) % m)) * inverse(V((s / g) % m), m) % m;
            V B = 1;
            while (U(B) * B * m * 6 < hi - lo) {
                if (B >= (V(1) << 40))
                    throw std::runtime_error("wide split");
                B *= 2;
            }
            V invB = inverse(B % m, m);
            V v = 0;
            for (V l = 0; l < B; ++l) {
                if ((B == 1 || (l & 1) == j) && rem >= s * l + a.h * reading(l)) {
                    U zlo = lo > l ? ceildiv(lo - l, B) : 0, zhi = hi >= l ? (hi - l) / B : 0;
                    U t = rem - s * l - a.h * reading(l), den = s * B + U(a.h) * (B - 1 - v);
                    zlo = std::max(zlo, ceildiv(t, den + a.h));
                    zhi = std::min(zhi, t / den);
                    if (zlo <= zhi) {
                        V res = U((residue + m - l % m) % m) * invB % m;
                        U first = zlo + (res + m - V(zlo % m)) % m, step = m;
                        if (B == 1) {
                            if ((first & 1) != j)
                                first += m;
                            step = 2 * m;
                        }
                        for (U z = first; z <= zhi; z += step)
                            verify(x, B * z + l, s, rem);
                    }
                }
                if (l + 1 < B) {
                    V bit = B / 2;
                    while (v & bit) {
                        v ^= bit;
                        bit /= 2;
                    }
                    v ^= bit;
                }
            }
        }
    }
    // List every lattice point in a dyadic strip of the shifted hyperbola.
    // Adjacent Farey fractions of order W supply the slanted integer lines.
    struct Fraction {
        V u, v;
    };
    template <V P> void tangent(V low, V high, V W) {
        const V aq = a.q % P, ac = a.c % P, at = a.t % P, ae = a.e % P, ar = a.r % P, ak = a.k % P,
                level = N % P;
        const double C = double(a.c) / double(a.q), T = double(a.t) / double(a.q);
        const double D = double(N - a.k) / double(a.q);
        const double g =
            2.0 * double(a.h) / double(a.q) * (sx + sy * D / double(low) / double(low)) +
            2.0 * (double(a.e) + double(a.r) + double(a.h) * (ox + oy)) / double(a.q) / double(low);
        const double bend = 2.0 * (D + 1) / double(low) / double(low) / double(low) * W * W;
        const double lam0 = (D + C * T) / ((low + C) * (low + C));
        Fraction L{V(lam0), 1}, R{V(lam0) + 1, 1};
        while (L.v + R.v <= W) {
            Fraction m{L.u + R.u, L.v + R.v};
            if (double(m.u) / m.v <= lam0)
                L = m;
            else
                R = m;
        }
        for (V X = low; X <= high; X += W) {
            V count = std::min(W, high - X + 1);
            U den = a.q * X + a.c;
            V yb = V(D / (X + C) - T);
            S diff = S(N - a.k) - S(a.t * X) - S(den * yb);
            double f0 = double(diff) / double(den);
            double lambda = (double(yb) + f0 + T) / (X + C);
            while (double(L.u) / L.v > lambda) {
                V k = (W + R.v) / L.v;
                Fraction before{k * L.u - R.u, k * L.v - R.v};
                R = L;
                L = before;
            }
            const bool lower = std::abs(L.u - lambda * L.v) < std::abs(R.u - lambda * R.v);
            Fraction sel = lower ? L : R;
            V aa = sel.u, b = sel.v, inv = lower ? (b - R.v % b) % b : L.v % b;
            double tilt = (double(aa) - b * lambda) * count;
            double error = std::ldexp(1., -46) * (1. + b * (1. + std::abs(f0) + bend + g) +
                                                  double(count) * (aa + b * lambda));
            int64_t bot = std::ceil(b * (f0 - g) + std::min(0., tilt) - error);
            int64_t top = std::floor(b * (f0 + bend) + std::max(0., tilt) + error);
            for (int64_t line = bot; line <= top; ++line) {
                int64_t rv = line % int64_t(b);
                if (rv < 0)
                    rv += b;
                V first = (V(rv) * inv) % b;
                for (V d = first; d < count; d += b) {
                    V x = X + d;

                    V y = yb + (line - int64_t(aa * d)) / int64_t(b);
                    if (y < x + swapped)
                        continue;
                    // Cheap exact prime-divisor test before the full address reading.
                    if ((aq * (x % P) * (y % P) + at * (x % P) + ac * (y % P) + ar * (x & 1) +
                         ae * (y & 1) + ak) %
                            P !=
                        level)
                        continue;
                    V raw_x = sx * x + ox, raw_y = sy * y + oy;
                    U used = original.k + original.t * raw_x + original.r * (raw_x & 1) +
                             U(a.h) * reading(raw_x) + original.e * (raw_y & 1);
                    if (used <= N)
                        verify(raw_x, raw_y, original.c + original.q * raw_x, N - used);
                }
            }
        }
    }
    void run(V upper, V start = 0) {
        const double D = double(N - a.k) / double(a.q), c = double(a.c) / double(a.q),
                     t = double(a.t) / double(a.q);
        if (a.c >= a.q || a.t >= a.q || upper >= (V(1) << 52))
            throw std::runtime_error("unexpected scale");
        if (start == 0)
            solve(0);
        // Use rational tangents where the curvature and nonlinear width are small.
        // Larger cyclic residue steps take over near either axis.
        V line_low = upper + 1;
        for (V low = 1; low <= upper; low *= 2) {
            V high = std::min(upper, low * 2 - 1);
            double idealW = double(low) / std::cbrt(D);
            double spread = double(a.h) / double(a.q) * (sx + sy * D / double(low) / double(low));
            if (idealW >= 80 && spread * idealW < 4 && D / double(low) < 3e15) {
                if (high < start)
                    continue;
                V W = std::max(V(1), V(idealW));
                if (!(a.h % 241))
                    tangent<241>(std::max(low, start), high, W);
                else if (!(a.h % 65537))
                    tangent<65537>(std::max(low, start), high, W);
                else if (!(a.h % 257))
                    tangent<257>(std::max(low, start), high, W);
                else
                    tangent<29>(std::max(low, start), high, W);
                line_low = std::min(line_low, low);
            } else if (line_low <= upper)
                throw std::runtime_error("noncontiguous tangent strip");
        }
        upper = line_low - 1;
        // A divisor of h eliminates both nonlinear readings. Keep y in one slope
        // progression (per pair of parities) before using the geometric test.
        V M = 1;
        for (V d = 2; d <= 100000; ++d)
            if (a.h % d == 0)
                M = d;
        V jumps = 2 * M;
        for (V xi = 0; xi < jumps; ++xi) {

            V sm = (a.c % M + (a.q % M) * xi) % M, g = std::gcd(sm, M), mod = M / g;
            V inv = inverse((sm / g) % mod, mod);
            for (V j = 0; j < 2; ++j) {
                if (j && a.e % M == 0)
                    continue;
                V right = (N % M + 4 * M - (a.k % M) - (a.t % M) * U(xi % M) % M -
                           (a.r % M) * (xi % 2) - (a.e % M) * j) %
                          M;
                if (right % g)
                    continue;
                V yres = U(right / g) * inv % mod;
                const int selected = a.e % M == 0 ? -1 : int(j);
                const double scale = D / double(mod), off = (t + double(yres)) / double(mod);
                for (V low = 1; low <= upper;) {
                    V high = std::min(upper, low * 2 - 1);
                    if (high < start) {
                        low = high + 1;
                        continue;
                    }
                    double top = D / double(low);
                    double margin =
                        (2.0 / double(low) +
                         2.0 * double(a.h) / double(a.q) * (sx + sy * top / double(low)) +
                         2.0 * (double(a.e) + double(a.r) + double(a.h) * (ox + oy)) / double(a.q) /
                             double(low)) /
                            double(mod) +
                        std::ldexp(1.0, -47) * (top / double(mod) + 2);
                    V first = std::max(low, start);
                    first += (xi + jumps - first % jumps) % jumps;
                    if (margin >= 0.45) {
                        for (V x = first; x <= high; x += jumps)
                            solve(x, selected);
                    } else {
                        constexpr int batch = 128;
                        std::array<unsigned char, batch> maybe;
                        while (first <= high && high - first >= jumps * (batch - 1)) {
                            for (int i = 0; i < batch; ++i) {
                                double u = scale / (double(first + i * jumps) + c) - off;
                                maybe[i] = std::abs(u - std::nearbyint(u)) <= margin;
                            }
                            for (int i = 0; i < batch; ++i)
                                if (maybe[i])
                                    solve(first + i * jumps, selected);
                            first += jumps * batch;
                        }
                        for (V x = first; x <= high; x += jumps) {
                            double u = scale / (double(x) + c) - off;
                            if (std::abs(u - std::nearbyint(u)) <= margin)
                                solve(x, selected);
                        }
                    }
                    low = high + 1;
                }
            }
        }
    }
};
int main(int argc, char **argv) {
    if (argc < 2 || argc > 7) {
        std::cerr << "usage: farey_hyperbola n [E61|lower [direction [first [last [left_residue]]]]]"
                     " < charts.txt\n";
        return 2;
    }
    int direction = argc > 3 ? std::stoi(argv[3]) : -1;
    V start = argc > 4 ? std::stoull(argv[4]) : 0, stop = argc > 5 ? std::stoull(argv[5]) : ~V(0);
    int wanted = std::stoi(argv[1]), n;
    std::array<std::string, 7> f;
    while (std::cin >> n >> f[0] >> f[1] >> f[2] >> f[3] >> f[4] >> f[5] >> f[6]) {
        if (n != wanted)
            continue;
        Chart a = {V(parse(f[0])), parse(f[1]), parse(f[2]), parse(f[3]),
                   parse(f[4]),    parse(f[5]), parse(f[6])};
        U N = argc > 2 && std::string(argv[2]) != "E61" ? parse(argv[2]) : TARGET;
        if (a.k > N)
            return 0;
        V g = std::gcd(std::gcd(V(a.c % a.h), a.h), std::gcd(V(a.q % a.h), V(a.e % a.h)));
        V step = 2 * g;
        for (V residue = 0; residue < step; ++residue) {
            if ((a.k % g + a.t * residue % g + (a.r % g) * (residue & 1)) % g != N % g)
                continue;
            if (argc > 6 && residue != std::stoull(argv[6]))
                continue;
            for (bool swapped : {false, true}) {
                if (direction >= 0 && direction != int(swapped))
                    continue;
                Chart o = a;
                if (swapped) {
                    std::swap(o.t, o.c);
                    std::swap(o.r, o.e);
                }
                Scan scan(o, N, swapped, swapped ? 1 : step, swapped ? 0 : residue,
                          swapped ? step : 1, swapped ? residue : 0);
                if (scan.a.k > N)
                    continue;
                U bound = (N - scan.a.k) / scan.a.q;
                V b = std::sqrt(double(bound));
                while (U(b + 1) * (b + 1) <= bound)
                    ++b;
                while (U(b) * b > bound)
                    --b;
                scan.run(std::min(b, stop), start);
                std::cout << "n " << n << " x mod " << step << " = " << residue << " swapped "
                          << swapped << " smaller bound " << b << " exact rows " << scan.exact_rows
                          << " addresses " << scan.candidates << " hits " << scan.hits << std::endl;
            }
        }
        return 0;
    }
    return 2;
}

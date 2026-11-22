// Exact two-ended inverse on a bounded binary-run language. All chart arithmetic is integral.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using U = unsigned __int128;
using L = std::uint64_t;

int width(U a) {
    L high = a >> 64, low = a;
    return high ? 128 - __builtin_clzll(high) : low ? 64 - __builtin_clzll(low) : 0;
}

int ones(U a) {
    return __builtin_popcountll(L(a)) + __builtin_popcountll(L(a >> 64));
}

int changes(U a, int bits) {
    if (bits <= 1) return 0;
    return ones((a ^ (a >> 1)) & ((U(1) << (bits - 1)) - 1));
}

U reading(U y) {
    U value = 0;
    while (y) {
        bool zero = !(y & 1);
        y >>= 1;
        if (zero) value += y + (y & 1);
    }
    return value;
}

std::string decimal(U n) {
    if (n == 0) return "0";
    std::string s;
    while (n) {
        s += char('0' + n % 10);
        n /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

struct Chart {
    unsigned c, h, even, odd;
    int bits;
};
constexpr std::array<Chart, 3> top = {{{81,25,43,68,115},{411,75,367,392,113},{1021,125,915,940,112}}};

struct Search {
    Chart chart;
    U target;
    int max_runs;
    U reference;
    L nodes = 0, leaves = 0;
    std::vector<U> hits;
    std::array<unsigned,128> divisor{};
    std::array<std::vector<unsigned>,128> quotient_residue;

    Search(Chart chart_, U target_, int max_runs_, U reference_):
        chart(chart_), target(target_), max_runs(max_runs_), reference(reference_) {
        // Solve B*c*z = rem modulo h at each suffix width B=2^b.
        for (int b=1; b<=chart.bits; ++b) {
            unsigned multiplier = unsigned((U(1) << b) % chart.h) * chart.c % chart.h;
            divisor[b] = std::gcd(multiplier,chart.h);
            unsigned mod = chart.h / divisor[b];
            quotient_residue[b].resize(mod);
            for (unsigned i=0; i<mod; ++i) quotient_residue[b][(i*(multiplier/divisor[b]))%mod] = i;
        }
    }

    bool feasible(U y) const {
        return changes(y ^ reference, chart.bits) < max_runs;
    }

    void visit(int b, U r, U reverse, U fr, int low_changes) {
        ++nodes;
        if (low_changes >= max_runs) return;
        U B = U(1) << b;
        unsigned k = (r & 1) ? chart.odd : chart.even;
        U offset = k + chart.c * r + chart.h * fr;
        if (offset > target) return;
        U rem = target - offset;
        U den = chart.c * B + chart.h * (B - 1 - reverse);
        U lo = (rem + den + chart.h - 1) / (den + chart.h), hi = rem / den;
        int rest = chart.bits - b;
        lo = std::max(lo, rest ? U(1) << (rest-1) : U(0));
        hi = std::min(hi, (U(1) << rest) - 1);
        if (lo > hi) return;
        unsigned rhs = rem % chart.h;
        if (rhs % divisor[b]) return;
        unsigned mod = chart.h/divisor[b];
        unsigned residue = quotient_residue[b][rhs/divisor[b]];
        lo += (residue + mod - unsigned(lo % mod)) % mod;
        if (lo > hi) return;
        hi -= (unsigned(hi % mod) + mod - residue) % mod;
        int loose = width(lo ^ hi);
        int high_changes = changes((lo ^ (reference >> b)) >> loose, rest-loose);
        if (high_changes + low_changes >= max_runs) return;
        if (hi-lo <= 2 * mod) {
            for (U z=lo; z<=hi; z+=mod) {
                ++leaves;
                U y = z*B + r;
                if (feasible(y) && k + chart.c*y + chart.h*reading(y) == target) hits.push_back(y);
            }
            return;
        }
        for (unsigned d=0; d<=1; ++d) {
            int toggle = int(d ^ ((r>>(b-1)) & 1) ^ ((reference>>b)&1) ^ ((reference>>(b-1))&1));
            U next_fr = fr + d * (B-1-reverse + (r < B/2));
            visit(b+1, r + d*B, 2*reverse + d, next_fr, low_changes + toggle);
        }
    }

    void run() {
        for (unsigned bit=0; bit<=1; ++bit) {
            unsigned k = bit ? chart.odd : chart.even;
            unsigned base = k+chart.c*bit;
            if (target < base || (target-base) % (2*std::gcd(chart.c,chart.h))) continue;
            visit(1, bit, bit, 0, 0);
        }
    }
};

void self_test() {
    for (int bits=4; bits<=13; bits+=3) {
        for (int which : {0,1,2}) {
            Chart chart = top[which];
            chart.bits=bits;
            for (unsigned j=0; j<35; ++j) {
                U first = U(1) << (bits-1), y = first + ((319*j) % first);
                U target = (y%2 ? chart.odd : chart.even) + chart.c*y + chart.h*reading(y) + 2*(j%2);
                U ref = j % 2 ? U(0xaaaaaaaaU) : 0;
                Search s(chart, target, j%bits+1, ref);
                s.run();
                std::vector<U> direct;
                for (U a=first; a<2*first; ++a)
                    if (s.feasible(a) && (a%2 ? chart.odd : chart.even)+chart.c*a+chart.h*reading(a)==target)
                        direct.push_back(a);
                std::sort(s.hits.begin(), s.hits.end());
                if (s.hits!=direct) throw std::runtime_error("bounded control failed");
            }
        }
    }
}

int main(int argc, char** argv) {
    if (argc<2 || argc>3) {
        std::cerr << "usage: suffix_runs MAX_RUNS [periodic-reference-from-low-bit]\n";
        return 2;
    }
    int runs=std::stoi(argv[1]);
    std::string period=argc==3 ? argv[2] : "0";
    if (runs<1 || runs>115 || period.empty() || period.find_first_not_of("01")!=std::string::npos) return 2;
    self_test();
    U center = (U(1)<<60)*((U(1)<<61)-1);
    for (Chart ch : top) {
        U ref=0;
        for (int i=0;i<ch.bits;++i) ref|=U(period[i%period.size()]=='1')<<i;
        Search s(ch,center-1,runs,ref);
        s.run();
        std::cout << "width " << ch.bits << " run cap " << runs << " reference " << period
                  << " nodes " << s.nodes << " leaf candidates " << s.leaves
                  << " hits " << s.hits.size() << std::endl;
        for (U y:s.hits) std::cout << "address " << decimal(y) << std::endl;
    }
}

// Exact inverse on signed differences between the two reflected halves of an address.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

using Wide = __int128_t;
using Small = std::int64_t;
using Unsigned = std::uint64_t;

struct Chart { int c, h, even, odd; };

Wide read(Wide y) {
    Wide value = 0;
    while (y) {
        Wide next = y >> 1;
        if (!(y & 1)) value += next + (next & 1);
        y = next;
    }
    return value;
}

std::string digits(Wide n) {
    if (n == 0) return "0";
    std::string answer;
    while (n > 0) {
        answer.push_back(char('0' + n % 10));
        n /= 10;
    }
    std::reverse(answer.begin(), answer.end());
    return answer;
}

Unsigned reverse_bits(Unsigned x, int width) {
    Unsigned result = 0;
    for (int i = 0; i < width; ++i) {
        result = 2 * result + (x & 1);
        x >>= 1;
    }
    return result;
}

Wide height(const Chart& chart, Wide word) {
    return (word & 1 ? chart.odd : chart.even) + chart.c * word + chart.h * read(word);
}

// A word is s, then g central bits m, then the b-bit reversal of v.
// Keep v = sign*s + offset and invert the monotone quadratic in s.
std::vector<Wide> invert(const Chart& chart, Wide target, int b, int g, int m,
                         int sign, Small offset) {
    const Small B = Small(1) << b, A = Small(1) << g;
    Small lo = B/2, hi = B-1;
    if (sign == 1) {
        lo = std::max(lo, -offset);
        hi = std::min(hi, B-1-offset);
    } else {
        lo = std::max(lo, offset-(B-1));
        hi = std::min(hi, offset);
    }
    if (lo > hi) return {};
    Wide constant = Wide(chart.c + chart.h)*B - Wide(chart.h)*(1 + Wide(offset));
    auto leading = [&](Small s) {
        return (Wide(A)*s+m) * (constant-Wide(chart.h)*sign*s);
    };
    const Wide room = Wide(chart.c + chart.h)*(B-1) + Wide(chart.h)*(Wide(A)*B)
                    + std::max(chart.even, chart.odd) - std::min(chart.even, chart.odd);
    const Wide level = target - std::min(chart.even, chart.odd);
    if (leading(lo) > level || leading(hi) + room < level) return {};
    Small top = hi;
    if (leading(hi) > level) {
        // Estimate via the stable quadratic formula. The integer comparisons
        // below are decisive even on platforms with binary64 long double.
        const long double linear = (long double)(Wide(A)*constant - Wide(chart.h)*sign*m);
        const long double rhs = (long double)(level - m*constant);
        const long double a = -static_cast<long double>(A) * chart.h * sign;
        const long double disc = std::max(0.L, linear*linear + 4*a*rhs);
        long double estimate = 2*rhs/(linear + std::sqrt(disc));
        top = static_cast<Small>(std::clamp(estimate, (long double)lo, (long double)hi));
        while (top > lo && leading(top) > level) --top;
        while (top < hi && leading(top+1) <= level) ++top;
    }
    std::vector<Wide> solutions;
    for (Small s = top; s >= lo && leading(s)+room >= level; --s) {
        Small v = sign*s + offset;
        Wide word = (Wide(A)*s+m)*B + reverse_bits(v, b);
        int k = word & 1 ? chart.odd : chart.even;
        if ((target-k-chart.c*word) % chart.h == 0 && height(chart, word) == target)
            solutions.push_back(word);
    }
    return solutions;
}

struct Search {
    Chart chart;
    Wide target;
    int b, g, weight;
    unsigned long long offsets = 0, hits = 0;

    void visit(int next, int used, Small delta) {
        Small B = Small(1) << b;
        if (delta > -B && delta < B) {
            ++offsets;
            for (int m = 0; m < (1 << g); ++m)
                for (int sign : {-1,1}) {
                    Small offset = delta + (sign == -1 ? B-1 : 0);
                    for (Wide word : invert(chart, target, b, g, m, sign, offset)) {
                        ++hits;
                        std::cout << "hit " << digits(word) << " sign " << sign << " offset " << delta << '\n';
                    }
                }
        }
        if (used == weight) return;
        for (int bit = next; bit <= b; ++bit) {
            Small step = Small(1) << bit;
            visit(bit+2, used+1, delta + step);
            visit(bit+2, used+1, delta - step);
        }
    }
};

void control() {
    Chart chart{81,25,43,68};
    for (int b = 1; b <= 5; ++b) {
        Small B = Small(1) << b;
        for (int g : {0,1,2})
            for (int m=0; m<(1<<g); ++m)
                for (int sign : {-1,1})
                    for (Small offset=-B; offset<=2*B; ++offset) {
                        std::map<Wide,std::vector<Wide>> exact;
                        for (Small s=B-1; s>=B/2; --s) {
                            Small v=sign*s+offset;
                            if (v < 0 || v >= B) continue;
                            Wide word=((Wide(1<<g)*s+m)*B)+reverse_bits(v,b);
                            exact[height(chart,word)].push_back(word);
                        }
                        for (auto& [n, words] : exact) {
                            assert(invert(chart,n,b,g,m,sign,offset) == words);
                            if (!exact.count(n+2)) assert(invert(chart,n+2,b,g,m,sign,offset).empty());
                        }
                    }
    }
}

int main(int argc, char** argv) {
    int weight = argc == 2 ? std::stoi(argv[1]) : 4;
    if (weight < 0 || weight > 58) throw std::invalid_argument("signed digit weight must be 0..58");
    control();
    Wide target = (Wide(1)<<60)*((Wide(1)<<61)-1)-1;
    for (auto [width,chart] : std::vector<std::pair<int,Chart>>{
             {115,{81,25,43,68}}, {112,{1021,125,915,940}}}) {
        Search search{chart,target,width/2,width%2,weight};
        search.visit(0,0,0);
        std::cout << "width " << width << " signed weight <= " << weight << " offsets "
                  << search.offsets << " hits " << search.hits << std::endl;
    }
}

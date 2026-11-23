// Coefficients for a transitive replacement of one new vertex over the seven Heegner primes.
#include <array>
#include <iostream>

using Poly = std::array<int, 6>;
using Rows = std::array<int, 8>;

void accumulate(Poly& out, const Poly& in, bool edge, bool gap) {
    if (edge) for (int k=0; k<=5; ++k) out[k] += in[k];
    if (gap) for (int k=0; k<5; ++k) out[k+1] += in[k];
}

// A gap contains one nonempty path through the module. The variable tracks its visits.
Poly skeleton(const Rows& rows, int outmask) {
    std::array<std::array<Poly,8>,256> paths{};
    for (int i=0; i<8; ++i) {
        paths[1<<i][i][0] = 1;
        paths[1<<i][i][1] = (outmask>>i)&1;
    }
    for (int s=1; s<256; ++s) for (int b=s; b; b&=b-1) {
        int tail=__builtin_ctz(b);
        for (int rest=255^s; rest; rest&=rest-1) {
            int next=__builtin_ctz(rest);
            bool edge = (rows[tail]>>next)&1;
            bool gap = !(outmask&(1<<tail)) && (outmask&(1<<next));
            accumulate(paths[s|(1<<next)][next],paths[s][tail],edge,gap);
        }
    }
    Poly total{};
    for (int i=0; i<8; ++i) accumulate(total,paths[255][i],true,!(outmask&(1<<i)));
    return total;
}

int main() {
    const std::array<int,7> frame{4,9,26,33,75,23,47}; // 3,7,11,19,43,67,163
    for (int fixed=0; fixed<128; ++fixed) {
        Rows rows{};
        rows[7]=fixed;
        for (int i=0; i<7; ++i) rows[i]=frame[i]|((1-((fixed>>i)&1))<<7);
        for (int module=0; module<256; ++module) {
            Poly p=skeleton(rows,module);
            std::cout<<fixed<<' '<<module;
            for (int value:p) std::cout<<' '<<value;
            std::cout<<'\n';
        }
    }
}

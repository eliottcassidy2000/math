// Count binary addresses in the necessary size window matching the exact decimal suffix.
#include <algorithm>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using U = unsigned __int128;
constexpr int modulus = 10000;

std::string decimal(U n) {
    if (!n) return "0";
    std::string out;
    while (n) { out += char('0'+n%10); n/=10; }
    std::reverse(out.begin(), out.end());
    return out;
}

struct Chart { unsigned c,h,k0,k1; int width,tail; };

U at_most(U bound, Chart ch, U target) {
    if (bound < U(1)<<(ch.width-1)) return 0;
    bound = std::min(bound, (U(1)<<ch.width)-1);
    // f is always even. Track f/2 modulo exactly the precision it contributes.
    int m=modulus/std::gcd(modulus,int(2*ch.h));
    std::vector<U> before(modulus*m), after(before.size());
    int tight_y=1,tight_f=0;
    bool tight=true;
    for (int i=ch.width-2;i>=0;--i) {
        bool forced_one=i<ch.tail;
        std::fill(after.begin(),after.end(),U(0));
        for (int y=0;y<modulus;++y) {
            int even=(2*y)%modulus, step=((y+1)/2)%m;
            for (int f=0;f<m;++f) {
                U ways=before[y*m+f];
                if (!ways) continue;
                if (!forced_one) after[even*m+(f+step)%m]+=ways;
                after[(even+1)*m+f]+=ways;
            }
        }
        if (tight) {
            int bit=(bound>>i)&1;
            int even=(2*tight_y)%modulus;
            if (bit && !forced_one) after[even*m+(tight_f+(tight_y+1)/2)%m]++;
            if (!bit && forced_one) tight=false;
            if (!bit) tight_f=(tight_f+(tight_y+1)/2)%m;
            tight_y=even+bit;
        }
        before.swap(after);
    }
    if (tight) ++before[tight_y*m+tight_f];
    U total=0;
    for (int y=0;y<modulus;++y)
        for (int f=0;f<m;++f)
            if (((y%2?ch.k1:ch.k0)+ch.c*y+2*ch.h*f)%modulus==target%modulus)
                total+=before[y*m+f];
    return total;
}

int main() {
    for (Chart c : {Chart{81,25,43,68,10,0},Chart{1021,125,915,940,10,2}}) {
        for (unsigned a : {587u,791u}) {
            unsigned f=0;
            for (int x=a;x;x>>=1) if (!(x&1)) f+=(x/2)+(x/2%2);
            U target=(a%2?c.k1:c.k0)+c.c*a+c.h*f;
            U direct=0;
            for (unsigned y=512;y<=881;++y) {
                if ((y+1)% (1u<<c.tail)) continue;
                unsigned fy=0;
                for (unsigned x=y;x;x>>=1) if (!(x&1)) fy+=(x/2)+(x/2%2);
                if (((y%2?c.k1:c.k0)+c.c*y+c.h*fy)%modulus==target%modulus) ++direct;
            }
            if (at_most(881,c,target)!=direct) return 1;
        }
    }
    U target=(U(1)<<60)*((U(1)<<61)-1)-1;
    for (Chart c: {Chart{81,25,43,68,115,0},Chart{1021,125,915,940,112,5}}) {
        U lo=std::max(U(1)<<(c.width-1),(target-std::max(c.k0,c.k1)+c.c+c.h-1)/(c.c+c.h));
        U hi=std::min((U(1)<<c.width)-1,(target-std::min(c.k0,c.k1))/c.c);
        U count=at_most(hi,c,target)-at_most(lo-1,c,target);
        std::cout << "width " << c.width << ", window [" << decimal(lo) << ", " << decimal(hi)
                  << "], forced terminal ones " << c.tail << ", matching lower suffix 2175: "
                  << decimal(count) << std::endl;
    }
}

// Independent literal subset dynamic program. No inherited count formulas.
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
using U = std::uint64_t;
void plus(U& a,U b) {
    if (b > std::numeric_limits<U>::max()-a) throw std::overflow_error("count addition");
    a += b;
}
U times(U a,U b) {
    if (a && b > std::numeric_limits<U>::max()/a) throw std::overflow_error("count multiplication");
    return a*b;
}
int main() {
    unsigned n,r;
    while(std::cin>>n>>r) {
        if(n<1 || n>23 || r>=n) throw std::runtime_error("bad graph order/root");
        unsigned size=1u<<n, full=size-1;
        std::vector<unsigned> rows(n);
        for(auto& a:rows) std::cin>>a;
        for(unsigned i=0;i<n;++i) {
            if(rows[i]&~full || (rows[i]>>i&1)) throw std::runtime_error("invalid row");
            for(unsigned j=i+1;j<n;++j)
                if(((rows[i]>>j&1)+(rows[j]>>i&1))!=1) throw std::runtime_error("not tournament");
        }
        std::vector<U> dp(static_cast<std::size_t>(size)*n),h(size);
        for(unsigned v=0;v<n;++v)dp[(std::size_t(1u<<v))*n+v]=1;
        for(unsigned mask=1;mask<size;++mask) {
            unsigned ends=mask;
            while(ends) {
                unsigned v=__builtin_ctz(ends);ends&=ends-1;
                U a=dp[std::size_t(mask)*n+v];
                if(!a)continue;
                plus(h[mask],a);
                unsigned av=rows[v]&~mask&full;
                while(av) {
                    unsigned u=__builtin_ctz(av);av&=av-1;
                    plus(dp[std::size_t(mask|(1u<<u))*n+u],a);
                }
            }
        }
        h[0]=1;
        U H=h[full],E=dp[std::size_t(full)*n+r],d=h[full^(1u<<r)];
        std::fill(dp.begin(),dp.end(),0);
        dp[std::size_t(1u<<r)*n+r]=1;
        U cover=0;
        for(unsigned mask=1;mask<size;++mask) {
            if(!(mask>>r&1))continue;
            unsigned ends=mask;
            while(ends) {
                unsigned v=__builtin_ctz(ends);ends&=ends-1;
                U a=dp[std::size_t(mask)*n+v];
                if(!a)continue;
                plus(cover,times(a,h[full^mask]));
                unsigned av=rows[v]&~mask&full;
                while(av) {
                    unsigned u=__builtin_ctz(av);av&=av-1;
                    plus(dp[std::size_t(mask|(1u<<u))*n+u],a);
                }
            }
        }
        std::cout<<H<<' '<<E<<' '<<cover<<' '<<d<<'\n';
    }
}

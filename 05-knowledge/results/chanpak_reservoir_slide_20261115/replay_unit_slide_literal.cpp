#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
using Word=unsigned;
using Count=std::uint64_t;

struct Graph {
    Word mark=0;
    std::vector<Word> rows;
    void append(const std::vector<Word>& core, Word m) {
        const unsigned n=rows.size(), k=core.size();
        const Word newer=((1u<<k)-1)<<n;
        for(auto& row:rows)row|=newer;
        for(Word row:core)rows.push_back(row<<n);
        mark|=m<<n;
    }
    void tape(const std::string& s) {
        for(char c:s)append({0},c=='1');
    }
    std::vector<Word> apex() {
        for(unsigned i=0;i<rows.size();++i)if(!((mark>>i)&1))rows[i]|=1u<<rows.size();
        rows.push_back(mark);
        return rows;
    }
};

Count count(const std::vector<Word>& G) {
    const unsigned n=G.size();
    const Word all=(1u<<n)-1, stride=1u<<(n-1);
    // compress a subset containing i by erasing bit i before indexing
    std::vector<Count> dp(std::size_t(n)*stride);
    auto pos=[&](unsigned i, Word s) {
        const Word lo=(1u<<i)-1;
        return std::size_t(i)*stride+((s>>1)&~lo)+(s&lo);
    };
    auto add=[](Count x, Count y) {
        if(y>std::numeric_limits<Count>::max()-x)throw std::runtime_error("path count overflow");
        return x+y;
    };
    for(Word s=1;s<=all;++s)for(Word t=s;t;t&=t-1) {
        const unsigned i=__builtin_ctz(t);
        const Word rem=s^(1u<<i);
        Count val=rem==0;
        for(Word u=rem&G[i];u;u&=u-1)val=add(val,dp[pos(__builtin_ctz(u),rem)]);
        dp[pos(i,s)]=val;
    }
    Count a=0;
    for(unsigned i=0;i<n;++i)a=add(a,dp[pos(i,all)]);
    return a;
}

int main() {
    const std::vector<Word> core{4,1,130,7,143,159,63,75};
    for(unsigned j=0;j<2;++j) {
        Graph g;
        if(j)g.tape("1010111000111"); else g.tape("1");
        g.append(core,35);
        if(!j)g.tape("000100011001");
        g.tape("0");
        const auto G=g.apex();
        const auto H=count(G);
        if(G.size()!=23||H!=11273425u-2*j)throw std::runtime_error("unexpected literal count");
        std::cout<<"j="<<j<<": "<<H<<" Hamilton paths at order "<<G.size()<<std::endl;
    }
}

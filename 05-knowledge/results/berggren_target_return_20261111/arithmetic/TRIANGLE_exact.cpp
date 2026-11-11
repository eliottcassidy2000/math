// Exact two-birth extension polynomial via ordered three-path covers.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <vector>
using U=std::uint64_t;
using I=std::int64_t;
using State=std::array<int,6>;
U encode(const State& a) {U x=0;for(int j=0;j<6;++j)x|=U(a[j]+1)<<(5*j);return x;}
State decode(U x) {State a;for(int j=0;j<6;++j)a[j]=int((x>>(5*j))&31)-1;return a;}
std::vector<State> singleton(int at) {
    std::vector<State> out;
    for(int bin=0;bin<3;++bin){State a; a.fill(-1);a[bin]=a[bin+3]=at;out.push_back(a);}
    return out;
}
std::vector<State> triangle(int at) {
    std::vector<State> out;
    // Enumerate all ordered partitions into three paths, allowing empty paths.
    for(int assignment=0;assignment<27;++assignment) {
        int code=assignment;
        std::array<std::vector<int>,3> bins;
        for(int v=0;v<3;++v){bins[code%3].push_back(v);code/=3;}
        std::array<std::vector<std::vector<int>>,3> choices;
        for(int b=0;b<3;++b) {
            auto p=bins[b];
            do {
                bool ok=true;
                for(unsigned j=1;j<p.size();++j)if((p[j-1]+1)%3!=p[j])ok=false;
                if(ok)choices[b].push_back(p);
            } while(std::next_permutation(p.begin(),p.end()));
        }
        for(const auto& p:choices[0])for(const auto& q:choices[1])for(const auto& r:choices[2]) {
            std::array<std::vector<int>,3> ps={p,q,r};State s;s.fill(-1);
            for(int b=0;b<3;++b)if(!ps[b].empty()){s[b]=at+ps[b].front();s[b+3]=at+ps[b].back();}
            out.push_back(s);
        }
    }
    if(out.size()!=33)throw std::runtime_error("triangle cover universe");
    return out;
}
int main() {
    int r;std::cin>>r;
    if(r<1 || r>4)throw std::runtime_error("triangle count scope");
    std::vector<int> gaps(r+1);int n=3*r;
    for(auto& g:gaps){std::cin>>g;if(g<0)throw std::runtime_error("negative gap");n+=g;}
    if(n>18)throw std::runtime_error("order scope");
    State empty;empty.fill(-1);
    std::unordered_map<U,I> dp{{encode(empty),1}};
    auto append=[&](const std::vector<State>& options) {
        std::unordered_map<U,I> next;
        for(const auto& [key,value]:dp) {
            State old=decode(key);
            for(const auto& loc:options) {
                State z=old;
                for(int b=0;b<3;++b)if(loc[b]>=0){if(z[b]<0)z[b]=loc[b];z[b+3]=loc[b+3];}
                next[encode(z)]+=value;
            }
        }
        dp.swap(next);
    };
    int at=0;
    for(int block=0;block<=r;++block) {
        for(int j=0;j<gaps[block];++j)append(singleton(at++));
        if(block<r){append(triangle(at));at+=3;}
    }
    std::map<U,I> polynomial;
    I total_covers=0;
    for(const auto& [key,value]:dp) {
        total_covers+=value;
        State s=decode(key);
        for(int first=0;first<2;++first) {
            int second=1-first;
            std::vector<std::pair<int,bool>> factors; // bool means use 1-x.
            if(s[0]>=0)factors.emplace_back(first*n+s[3],true);
            if(s[1]>=0){factors.emplace_back(first*n+s[1],false);factors.emplace_back(second*n+s[4],true);}
            else factors.emplace_back(2*n,first==1); // y means new0 -> new1.
            if(s[2]>=0)factors.emplace_back(second*n+s[2],false);
            std::vector<std::pair<U,I>> terms{{0,value}};
            for(auto [index,negative]:factors) {
                if(negative) {
                    unsigned count=terms.size();
                    for(unsigned j=0;j<count;++j)terms.emplace_back(terms[j].first|(U(1)<<index),-terms[j].second);
                } else for(auto& term:terms)term.first|=U(1)<<index;
            }
            for(auto [monomial,coefficient]:terms)polynomial[monomial]+=coefficient;
        }
    }
    I content=0;std::size_t nonzero=0;
    for(auto [monomial,coefficient]:polynomial)if(coefficient){content=std::gcd(content,coefficient);++nonzero;}
    std::cout<<n<<' '<<dp.size()<<' '<<total_covers<<' '<<nonzero<<' '<<content<<'\n';
    for(auto [monomial,coefficient]:polynomial)if(coefficient)std::cout<<monomial<<' '<<coefficient<<'\n';
}

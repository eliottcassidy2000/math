// Input from fourth.py --board. Each bit is a physical fiber or slot number.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <vector>

using Mask = std::uint64_t;
using Triple = std::array<int,3>;
constexpr int radius = 7;
std::array<std::vector<Triple>,56> domain;
std::vector<int> branches;
Mask required[3]{}, owners[3][56]{};
std::uint64_t tested=0, hits=0;
Mask bit(int i) { return Mask(1)<<i; }

bool fill(Mask rows, const std::array<Mask,3>& available) {
    if (!rows) return true;
    int h=__builtin_ctzll(rows);
    for (const Triple& t : domain[h]) {
        std::array<Mask,3> next=available;
        bool allowed=true;
        for (int k=0; k<3; ++k) {
            if (!(available[k] & bit(t[k]))) { allowed=false; break; }
            next[k]^=bit(t[k]);
        }
        if (allowed && fill(rows^bit(h), next)) return true;
    }
    return false;
}

void choose(Mask changed, int left, int first) {
    if (left) {
        for (int i=first; i<=int(branches.size())-left; ++i)
            choose(changed|bit(branches[i]), left-1, i+1);
        return;
    }
    if (!(changed&(bit(6)|bit(17))) || !(changed&(bit(40)|bit(53)))) return;
    std::array<Mask,3> available={required[0],required[1],required[2]};
    for (int k=0; k<3; ++k)
        for (int val=0; val<56; ++val)
            if (owners[k][val] & ~changed) available[k] &= ~bit(val);
    ++tested;
    if (fill(changed,available)) ++hits;
}

int main() {
    int h,a,b,c,reference;
    while(std::scanf("%d%d%d%d%d",&h,&a,&b,&c,&reference)==5) {
        if(domain[h].empty()) branches.push_back(h);
        Triple t={a,b,c};
        domain[h].push_back(t);
        for(int k=0; k<3; ++k) {
            required[k]|=bit(t[k]);
            if(reference) owners[k][t[k]] |= bit(h);
        }
    }
    if (branches.size()!=53) return 2;
    choose(0,radius,0);
    std::printf("radius %d, supports tested %llu, hits %llu\n",radius,
                (unsigned long long)tested,(unsigned long long)hits);
    return tested==8493464 && hits==0 ? 0 : 1;
}

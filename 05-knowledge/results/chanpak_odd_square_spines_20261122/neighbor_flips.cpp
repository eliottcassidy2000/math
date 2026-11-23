// Exact Hamilton-path derivative for an edge of an ordinal sum with an apex.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>
#include <gmpxx.h>

using Table = std::array<std::array<unsigned,12>,12>;
using State = std::array<mpz_class,12>;
struct Segment { int in, out; }; // entry vertex and exit vertex of a contracted word
int state(int prefix_ok, int middle, int suffix_used) {
    return (3*prefix_ok+middle)*2+suffix_used;
}
struct Host {
    std::vector<std::string> arc;
    std::vector<std::vector<int>> blocks;
    std::map<std::array<int,6>,Table> cache;

    bool beats(int u,int v) const { return arc[u][v]=='1'; }
    // Four booleans extracted from a single distinct physical segment order.
    int contacts(int v,Segment a,Segment b) const {
        return int(beats(v,a.in))+2*int(beats(a.out,v))+4*int(beats(v,b.in))+8*int(beats(b.out,v));
    }

    const Table& local(int index, int removed, Segment a, Segment b) {
        const auto& block=blocks[index];
        std::array<int,6> key{index,removed,0,0,0,0};
        std::vector<int> vertices;
        for(int i=0;i<int(block.size());i++) if(!(removed>>i&1)) {
            int v=block[i], flags=contacts(v,a,b);
            vertices.push_back(v);
            for(int j=0;j<4;j++) if(flags>>j&1) key[j+2]|=1<<i;
        }
        auto [it,inserted]=cache.try_emplace(key,Table{});
        if(!inserted) return it->second;
        auto& M=it->second;
        using Description=std::array<int,5>; // prefix end, middle start/end, suffix start, middle present
        std::map<Description,unsigned> desc;
        int m=int(vertices.size());
        do {
            for(int i=0;i<=m;i++) for(int j=i;j<=m;j++) {
                bool valid=true;
                for(int k=1;k<m;k++) if(k!=i&&k!=j&&!beats(vertices[k-1],vertices[k])) valid=false;
                if(!valid) continue;
                int p=i?int(beats(vertices[i-1],a.in)):-1;
                int ms=i<j?int(beats(a.out,vertices[i])):0;
                int me=i<j?int(beats(vertices[j-1],b.in)):0;
                int s=j<m?int(beats(b.out,vertices[j])):-1;
                ++desc[{p,ms,me,s,int(i<j)}];
            }
        } while(std::next_permutation(vertices.begin(),vertices.end()));
        for(auto [z,weight]:desc) {
            auto [p,ms,me,s,used]=z;
            for(int i=0;i<2;i++) for(int j=0;j<3;j++) for(int k=0;k<2;k++) {
                if((used&&j==0&&!ms)||(s==0&&!k)) continue;
                int to=state(p<0?i:p, used?1+me:j, s<0?k:1);
                M[state(i,j,k)][to]+=weight;
            }
        }
        return M;
    }

    mpz_class ordered(int u,int v,Segment a,Segment b) {
        State live{};
        live[state(1,0,0)]=1;
        for(int ix=0;ix<int(blocks.size());ix++) {
            const auto& block=blocks[ix];
            if(block.size()==1) {
                int w=block[0];
                if(w==u||w==v) continue;
                int mask=contacts(w,a,b);
                State next{};
                for(int i=0;i<2;i++) for(int j=0;j<3;j++) for(int k=0;k<2;k++) {
                    const auto& n=live[state(i,j,k)];
                    next[state(mask&1,j,k)]+=n;
                    if(j||(mask&2)) next[state(i,1+int(bool(mask&4)),k)]+=n;
                    if(k||(mask&8)) next[state(i,j,1)]+=n;
                }
                live.swap(next);
            } else {
                int removed=0;
                for(int j=0;j<int(block.size());j++) if(block[j]==u||block[j]==v) removed|=1<<j;
                const auto& M=local(ix,removed,a,b);
                State next{};
                for(int i=0;i<12;i++) for(int j=0;j<12;j++) if(M[i][j])
                    mpz_addmul_ui(next[j].get_mpz_t(),live[i].get_mpz_t(),M[i][j]);
                live.swap(next);
            }
        }
        mpz_class total=0;
        for(int k=0;k<2;k++) {
            total+=live[state(1,2,k)];
            if(beats(a.out,b.in)) total+=live[state(1,0,k)];
        }
        return total;
    }

    // Apex already inside the contracted ordered pair; introduce a dummy
    // sentinel fixed at the beginning. The empty prefix is forced by using
    // a sentinel with no incoming edges, and it points to every real vertex.
    // Equivalently count one physical separator using the original two slots.
    mpz_class endpoint(int u,int v) {
        int apex=int(arc.size())-2, sentinel=apex+1;
        if(u==apex||v==apex) return ordered(u,v,{sentinel,sentinel},{u,v});
        return ordered(u,v,{apex,apex},{u,v})+ordered(u,v,{u,v},{apex,apex});
    }
};

int main() {
    int n,b;
    while(std::cin>>n>>b) {
        Host host;
        int pos=0;
        for(int j=0;j<b;j++) {
            int s;std::cin>>s;
            if(s<1||s>6) return 2;
            std::vector<int> block(s);
            std::iota(block.begin(),block.end(),pos);
            pos+=s;host.blocks.push_back(block);
        }
        if(pos!=n-1) return 2;
        for(int j=0;j<n;j++) {
            std::string row;std::cin>>row;
            if(int(row.size())!=n) return 2;
            host.arc.push_back(row+"0");
        }
        host.arc.push_back(std::string(n,'1')+"0");
        for(int u=0;u<n;u++) for(int v=u+1;v<n;v++) {
            int a=host.beats(u,v)?u:v,b= u+v-a;
            auto loss=host.endpoint(a,b);
            auto gain=host.endpoint(b,a);
            std::cout<<u<<" "<<v<<" "<<gain-loss<<"\n";
        }
    }
}

// Exact unordered triple-event union census. The Python driver supplies
// exact profile weights and restores the factor 3! for factorial moments.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using U64 = std::uint64_t;

static void need(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}

static U64 signature(U64 points) {
    std::array<unsigned,16> adj{};
    std::array<unsigned,16> degree{};
    unsigned occupied=0, ends=0;
    while (points) {
        unsigned p=__builtin_ctzll(points); points &= points-1;
        unsigned r=p>>3, c=8+(p&7);
        if (++degree[r]>2 || ++degree[c]>2) return 0;
        adj[r] |= 1u<<c; adj[c] |= 1u<<r;
        occupied |= (1u<<r)|(1u<<c);
    }
    for (unsigned j=0;j<16;++j) if(degree[j]==1) ends|=1u<<j;
    std::array<unsigned,9> tokens{};
    unsigned count=0;
    while(occupied) {
        unsigned component=0, front=occupied&-occupied;
        while(front) {
            component |= front;
            unsigned neighbors=0;
            while(front) {
                unsigned j=__builtin_ctz(front); front &= front-1;
                neighbors |= adj[j];
            }
            front=neighbors & ~component;
        }
        occupied &= ~component;
        unsigned r=__builtin_popcount(component&255), c=__builtin_popcount(component>>8);
        bool cycle=(component&ends)==0;
        unsigned m=r+c-(cycle?0:1);
        need(count<9 && m<=9,"signature edge budget");
        tokens[count++]=3*m+(r>c?1:(r<c?2:0));
    }
    std::sort(tokens.begin(),tokens.begin()+count);
    U64 key=0;
    for(unsigned j=0;j<count;++j) key |= U64(tokens[j])<<(5*j);
    return key;
}

static std::vector<U64> triples(int n) {
    std::vector<U64> result;
    for(int x1=0;x1<n;++x1) for(int x2=x1+1;x2<n;++x2) for(int x3=x2+1;x3<n;++x3)
    for(int y1=0;y1<n;++y1) for(int y2=0;y2<n;++y2) {
        if(y1==y2) continue;
        int numerator=(y2-y1)*(x3-x1);
        if(numerator%(x2-x1)) continue;
        int y3=y1+numerator/(x2-x1);
        if(y3<0 || y3>=n) continue;
        result.push_back((U64(1)<<(8*x1+y1))|(U64(1)<<(8*x2+y2))|(U64(1)<<(8*x3+y3)));
    }
    std::sort(result.begin(),result.end());
    need(std::adjacent_find(result.begin(),result.end())==result.end(),"duplicate grid triple");
    return result;
}

int main(int argc,char** argv) {
    need(argc==3 || argc==4,"usage: census n minimum_union_edges [literal]");
    int n=std::atoi(argv[1]), minimum=std::atoi(argv[2]);
    bool literal=argc==4;
    need(!literal || std::string(argv[3])=="literal","literal replay mode");
    need(2<=n && n<=8 && 3<=minimum && minimum<=9,"bounded universe");
    auto ts=triples(n);
    std::unordered_map<U64,U64> counts;
    U64 total=0, rejected_degree=0, rejected_size=0, accepted=0;
    U64 board=0;
    for(int r=0;r<n;++r) board |= U64((1u<<n)-1)<<(8*r);
    for(std::size_t i=0;i<ts.size();++i) for(std::size_t j=i+1;j<ts.size();++j) {
        U64 pair=ts[i]|ts[j];
        unsigned free_cols=(1u<<n)-1, free_rows=(1u<<n)-1;
        std::array<unsigned,8> col_degree{};
        for(int r=0;r<n;++r) {
            unsigned row=(pair>>(8*r))&255;
            if(__builtin_popcount(row)==2) free_rows &= ~(1u<<r);
            while(row) {
                unsigned c=__builtin_ctz(row); row &= row-1;
                if(++col_degree[c]==2) free_cols &= ~(1u<<c);
            }
        }
        // A new edge must avoid each already-saturated row and column.
        // Existing edges remain allowed because the union is a set.
        U64 allowed=pair;
        for(int r=0;r<n;++r) if(free_rows&(1u<<r)) allowed |= U64(free_cols)<<(8*r);
        need((allowed&~board)==0,"allowed grid universe");
        for(std::size_t k=j+1;k<ts.size();++k) {
            ++total;
            if(literal) {
                U64 points=pair|ts[k], key=signature(points);
                if(!key) {++rejected_degree;continue;}
                if(__builtin_popcountll(points)<minimum) {++rejected_size;continue;}
                ++counts[key]; ++accepted;
                continue;
            }
            if(ts[k]&~allowed) {++rejected_degree;continue;}
            U64 points=pair|ts[k];
            if(__builtin_popcountll(points)<minimum) {++rejected_size;continue;}
            U64 key=signature(points);
            need(key!=0,"saturation filter versus literal degrees");
            ++counts[key]; ++accepted;
        }
    }
    need(total==U64(ts.size())*(ts.size()-1)*(ts.size()-2)/6,"complete unordered event universe");
    need(total==rejected_degree+rejected_size+accepted,"census partition");
    std::vector<std::pair<U64,U64>> sorted(counts.begin(),counts.end());
    std::sort(sorted.begin(),sorted.end());
    std::cout<<"META "<<n<<" "<<minimum<<" "<<ts.size()<<" "<<total<<" "<<rejected_degree<<" "<<rejected_size<<" "<<accepted<<" "<<counts.size()<<"\n";
    std::cout<<"GRID";
    for(U64 t:ts) std::cout<<" "<<t;
    std::cout<<"\n";
    for(auto [key,count]:sorted) std::cout<<"TYPE "<<key<<" "<<count<<"\n";
}

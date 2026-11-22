// Exact images of two transitive modules replaced by arbitrary positive widths.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using Wide = unsigned __int128;
using Row = std::array<unsigned,8>;
using Counts = std::array<unsigned,64>;
constexpr unsigned dim=8;
constexpr Wide center=(Wide(1)<<60)*((Wide(1)<<61)-1);
constexpr Wide cap=center+2;
constexpr unsigned max_width=121;
using Readings=std::array<std::array<Wide,dim>,max_width+1>;

bool edge(const Row& rows, unsigned a, unsigned b) { return (rows[a]>>b)&1; }
Wide add(Wide a, Wide b) { return a>=cap-b ? cap : a+b; }
Wide mul(Wide a, Wide b) { return a && b>cap/a ? cap : a*b; }

Readings surjections() {
    Readings ways{};
    ways[0][0]=1;
    for(unsigned w=1;w<=max_width;++w)
        for(unsigned j=1;j<dim;++j)
            ways[w][j]=mul(j,add(ways[w-1][j],ways[w-1][j-1]));
    return ways;
}

bool strong(const Row& rows,unsigned size) {
    unsigned full=(1u<<size)-1;
    for(bool backward:{false,true}) {
        unsigned seen=1;
        for(unsigned prev=0;prev!=seen;) {
            prev=seen;
            for(unsigned i=0;i<size;++i) if((seen>>i)&1)
                seen |= backward ? full & ~(rows[i]|(1u<<i)) : rows[i];
        }
        if(seen!=full) return false;
    }
    return true;
}

// Each gap between ordinary vertices contains one of empty, P, Q, PQ, where P beats Q.
unsigned options(const Row& rows,int before,int after,unsigned p,unsigned q) {
    auto arrive=[&](unsigned u){ return before<0 || edge(rows,before,u); };
    auto depart=[&](unsigned u){ return after<0 || edge(rows,u,after); };
    unsigned nil=(before<0 || after<0 || edge(rows,before,after));
    return nil | ((arrive(p)&&depart(p))<<1) | ((arrive(q)&&depart(q))<<2)
               | ((arrive(p)&&depart(q))<<3);
}

// The commutative product depends only on the multiset of gap types.
const Counts& product(uint32_t key,std::unordered_map<uint32_t,Counts>& cache) {
    auto it=cache.find(key);
    if(it!=cache.end()) return it->second;
    Counts counts{};
    counts[0]=1;
    for(uint32_t remaining=key;remaining;remaining>>=4) {
        unsigned bits=remaining&15;
        Counts next{};
        for(unsigned i=0;i<dim;++i) for(unsigned j=0;j<dim;++j) {
            unsigned v=counts[dim*i+j];
            if(!v) continue;
            if(bits&1) next[dim*i+j]+=v;
            if(bits&2) next[dim*(i+1)+j]+=v;
            if(bits&4) next[dim*i+j+1]+=v;
            if(bits&8) next[dim*(i+1)+j+1]+=v;
        }
        counts=next;
    }
    return cache.emplace(key,counts).first->second;
}

Counts skeletons(const Row& rows,unsigned size,unsigned p,unsigned q,
                 std::unordered_map<uint32_t,Counts>& cache) {
    if(!edge(rows,p,q)) std::swap(p,q);
    std::vector<unsigned> rest;
    for(unsigned i=0;i<size;++i) if(i!=p && i!=q) rest.push_back(i);
    unsigned gaps[9][9]{};
    for(unsigned i=0;i<=size;++i) for(unsigned j=0;j<=size;++j)
        gaps[i][j]=options(rows,i==size ? -1 : int(i),j==size ? -1 : int(j),p,q);
    Counts total{};
    do {
        unsigned previous=size;
        std::array<unsigned,7> kinds{};
        for(unsigned j=0;j<=rest.size();++j) {
            unsigned after=j==rest.size() ? size : rest[j];
            kinds[j]=gaps[previous][after];
            if(!kinds[j]) break;
            previous=after;
            if(j==rest.size()) {
                std::sort(kinds.begin(),kinds.begin()+j+1);
                uint32_t key=0;
                for(unsigned k=0;k<=j;++k) key=(key<<4)|kinds[k];
                auto& contribution=product(key,cache);
                for(unsigned l=0;l<dim*dim;++l) total[l]+=contribution[l];
            }
        }
    } while(std::next_permutation(rest.begin(),rest.end()));
    return total;
}

struct Term { unsigned r,s,c; };
using Terms=std::vector<Term>;
Wide height(const Terms& terms,const Readings& J,unsigned a,unsigned b) {
    Wide sum=0;
    for(const auto& term:terms) {
        Wide first=J[a][term.r],second=J[b][term.s];
        if(!first || !second) continue;
        // These small factors are exact also in the saturated semiring.
        sum=add(sum,mul(mul(first,second),term.c));
        if(sum==cap) break;
    }
    return sum;
}

unsigned literal(const Row& rows,unsigned n) {
    unsigned paths[256][8]{};
    for(unsigned i=0;i<n;++i) paths[1u<<i][i]=1;
    for(unsigned mask=1;mask<(1u<<n);++mask) for(unsigned i=0;i<n;++i) {
        if(!paths[mask][i]) continue;
        unsigned choices=rows[i]&~mask;
        while(choices) {
            unsigned j=__builtin_ctz(choices);
            paths[mask|(1u<<j)][j]+=paths[mask][i];
            choices&=choices-1;
        }
    }
    return std::accumulate(paths[(1u<<n)-1],paths[(1u<<n)-1]+n,0u);
}

unsigned literal_expansion(const Row& rows,unsigned n,unsigned p,unsigned q,
                           unsigned a,unsigned b) {
    std::vector<unsigned> labels;
    for(unsigned i=0;i<n;++i)
        labels.insert(labels.end(),i==p ? a : (i==q ? b : 1),i);
    if(labels.size()>8) throw std::logic_error("literal control too large");
    Row expanded{};
    for(unsigned i=0;i<labels.size();++i) for(unsigned j=0;j<labels.size();++j) {
        bool arc=labels[i]==labels[j] ? i<j : edge(rows,labels[i],labels[j]);
        if(arc) expanded[i] |= 1u<<j;
    }
    return literal(expanded,labels.size());
}

int main(int argc,char** argv) {
    if(argc!=3) { std::cerr<<"usage: two-variable-ports catalogue.rows target-exponent\n";return 2; }
    int exponent=std::stoi(argv[2]);
    if(exponent<2 || exponent>61) return 2;
    Wide value=(Wide(1)<<(exponent-1))*((Wide(1)<<exponent)-1);
    std::ifstream input(argv[1]);
    if(!input) throw std::runtime_error("cannot read catalogue");
    auto J=surjections();
    std::unordered_map<uint32_t,Counts> cache;
    unsigned sources=0, charts=0, low=0, high=0;
    std::string line;
    while(std::getline(input,line)) {
        std::istringstream entries(line);
        Row rows{};
        unsigned n=0,v;
        while(entries>>v) {
            if(n>=8) throw std::runtime_error("maximum supported quotient order is eight");
            rows[n++]=v;
        }
        if(n<3 || !strong(rows,n)) continue;
        ++sources;
        unsigned quotient_height=literal(rows,n);
        for(unsigned i=0;i<n;++i) for(unsigned j=i+1;j<n;++j) {
            unsigned p=i,q=j;
            if(!edge(rows,p,q)) std::swap(p,q);
            auto counts=skeletons(rows,n,p,q,cache);
            if(counts[dim+1]!=quotient_height) throw std::logic_error("singleton control");
            ++charts;
            Terms terms;
            for(unsigned r=1;r<dim;++r) for(unsigned s=1;s<dim;++s)
                if(unsigned count=counts[r*dim+s]) terms.push_back({r,s,count});
            std::reverse(terms.begin(),terms.end());
            if(n<=6) for(auto widths : {std::pair<unsigned,unsigned>{2,2},{1,3},{3,1}}) {
                auto [a,b]=widths;
                if(height(terms,J,a,b)!=literal_expansion(rows,n,p,q,a,b))
                    throw std::logic_error("independent inflated path count");
            }
            unsigned b=2*exponent-1;
            for(unsigned a=1;a<=unsigned(2*exponent-1) && b;++a) {
                while(b && height(terms,J,a,b)>value+1) --b;
                for(unsigned w=b;w;--w) {
                    Wide h=height(terms,J,a,w);
                    if(h<value-1) break;
                    if(h!=value-1 && h!=value+1) throw std::logic_error("even tournament height");
                    if(h<value) ++low; else ++high;
                    if(exponent>10) {
                        std::cout<<"hit "<<(h<value ? "lower":"upper")<<" rows";
                        for(unsigned u=0;u<n;++u) std::cout<<" "<<rows[u];
                        std::cout<<" ports "<<p<<" "<<q<<" widths "<<a<<" "<<w<<"\n";
                    }
                }
            }
        }
    }
    std::cout<<"strong sources "<<sources<<" two-port charts "<<charts
             <<" lower presentations "<<low<<" upper presentations "<<high<<"\n";
}

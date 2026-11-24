// Input: one upper triangle per line, as emitted by nauty gentourng -cq.
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using Count = long long;
using Wide = __int128;
constexpr int max_order = 10;
using Rows = std::array<int,max_order>;
struct Profile { Count h, first, second, both; };
struct Polynomial { Count a,b,c,d; };

Rows decode(int n, const std::string &line) {
    if (int(line.size()) != n*(n-1)/2) throw std::runtime_error("wrong triangle length");
    Rows rows{};
    int at=0;
    for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) {
        if (line[at]=='1') rows[i] |= 1<<j;
        else if (line[at]=='0') rows[j] |= 1<<i;
        else throw std::runtime_error("nonbinary triangle");
        ++at;
    }
    return rows;
}

// For each port intersect the pairs of exterior vertices in all its triangles.
Rows common_vertices(const Rows &rows,int n) {
    Rows common{};
    for (int i=0;i<n;i++) {
        common[i]=((1<<n)-1)^(1<<i);
        for (int out=rows[i];out;out&=out-1) {
            int j=__builtin_ctz(out);
            for (int in=rows[j]&~rows[i]&~(1<<i);in;in&=in-1)
                common[i] &= (1<<j)|(in&-in);
        }
    }
    return common;
}

Profile profile(const Rows &rows,int n,int u,int v) {
    // State is the ordinary subset, the two visit counters and the last label.
    int masks=1<<(n-2), ordinary=((1<<n)-1)^(1<<u)^(1<<v);
    Rows bits{};
    for (int j=0,bit=1;j<n;j++) if ((ordinary>>j)&1) {
        bits[j]=bit;
        bit*=2;
    }
    std::vector<std::array<Count,max_order>> ways(masks*9);
    for (int e=0;e<n;e++) ways[bits[e]*9+3*(e==u)+(e==v)][e]=1;
    Count total[3][3]{};
    for (int mask=0;mask<masks;mask++) {
        int unused=0;
        for (int j=0;j<n;j++) if (bits[j] && !(mask&bits[j])) unused |= 1<<j;
        for (int r=0;r<=2;r++) for (int s=0;s<=2;s++) {
            int state=9*mask+3*r+s;
            for (int e=0;e<n;e++) {
                Count count=ways[state][e];
                if (!count) continue;
                if (mask==masks-1) total[r][s]+=count;
                if (r<2 && (rows[e]>>u&1)) ways[state+3][u]+=count;
                if (s<2 && (rows[e]>>v&1)) ways[state+1][v]+=count;
                for (int t=rows[e]&unused;t;t&=t-1) {
                    int j=__builtin_ctz(t);
                    ways[state+9*bits[j]][j]+=count;
                }
            }
        }
    }
    return {total[1][1],total[2][1],total[1][2],total[2][2]};
}

Polynomial polynomial(Profile p) {
    return {p.both,p.first-2*p.both,p.second-2*p.both,
            p.h-2*(p.first+p.second)+4*p.both};
}

std::vector<int> ray_offsets(Polynomial f) {
    if (f.a<=0 || (f.a&(f.a-1))) return {};
    int sum=1-__builtin_ctzll(f.a);
    std::vector<int> options;
    if (f.b) options.push_back(-__builtin_ctzll(std::abs(f.b)));
    if (f.c) options.push_back(sum+__builtin_ctzll(std::abs(f.c)));
    if (f.b && f.c) {
        int twice=sum+__builtin_ctzll(std::abs(f.c))-__builtin_ctzll(std::abs(f.b));
        if (twice%2==0) options.push_back(twice/2);
    }
    std::sort(options.begin(),options.end());
    options.erase(std::unique(options.begin(),options.end()),options.end());
    std::vector<int> valid;
    for (int i:options) {
        int j=sum-i, shift=-std::min({0,i,j});
        if (f.b*(Wide(1)<<(i+shift))+f.c*(Wide(1)<<(j+shift))==-(Wide(1)<<shift))
            valid.push_back(i);
    }
    return valid;
}

int valuation(Wide positive) {
    auto low=static_cast<unsigned long long>(positive);
    return low ? __builtin_ctzll(low) : 64+__builtin_ctzll(static_cast<unsigned long long>(positive>>64));
}

std::pair<int,int> scalar_widths(Polynomial f,int p,int sign) {
    Wide target=(Wide(1)<<(p-1))*((Wide(1)<<p)-1)+sign;
    Wide difference=target-f.d;
    int limit=2*p;
    if (difference) limit=std::min(limit,valuation(difference>0 ? difference : -difference));
    for (int order=0;order<2;order++) {
        for (int m=1;m<=limit;m++) {
            Wide x=Wide(1)<<m, numerator=difference-f.b*x, divisor=f.a*x+f.c;
            if (divisor<=0) throw std::runtime_error("missing second return");
            if (numerator<=0) continue;
            int k=valuation(numerator)-valuation(divisor);
            if (k<m || k>2*p || (order && k==m)) continue;
            if ((numerator>>k)==divisor)
                return order ? std::make_pair(k,m) : std::make_pair(m,k);
        }
        std::swap(f.b,f.c);
    }
    return {-1,-1};
}

int cutoff(Count bound) {
    int t=0;
    while ((Count(1)<<t)<=12*bound) ++t;
    Count limit=bound*((Count(1)<<(t-1))+1);
    int p=2;
    while ((Wide(1)<<p)<=4*limit+2 || p<t+1 || p<2*(64-__builtin_clzll(bound))+2) ++p;
    return p;
}

int main(int argc,char**argv) {
    if (argc!=2) throw std::runtime_error("usage: port_rays ORDER < triangles");
    int n=std::stoi(argv[1]);
    if (n<3 || n>max_order) throw std::runtime_error("supported orders 3 through 10");
    Count graphs=0,eligible_graphs=0,pairs=0,hits=0,coefficient_bound=0;
    std::map<Count,Count> offsets;
    std::set<std::array<Count,4>> polynomials;
    std::set<std::pair<int,int>> small_hits;
    std::string line;
    while (std::cin>>line) {
        ++graphs;
        Rows rows=decode(n,line), common=common_vertices(rows,n);
        bool eligible=false;
        for (int u=0;u<n;u++) for (int v=u+1;v<n;v++) {
            if (!(common[u]&~(1<<v)) || !(common[v]&~(1<<u))) continue;
            eligible=true;
            ++pairs;
            Polynomial f=polynomial(profile(rows,n,u,v));
            coefficient_bound=std::max({coefficient_bound,f.a,std::abs(f.b),std::abs(f.c),std::abs(f.d)+1});
            auto key=std::array<Count,4>{f.a,std::min(f.b,f.c),std::max(f.b,f.c),f.d};
            if (polynomials.insert(key).second) {
                Count bound=std::max({f.a,std::abs(f.b),std::abs(f.c),std::abs(f.d)+1});
                int stop=cutoff(bound);
                if (stop>61) throw std::runtime_error("cutoff exceeds exponent 61");
                for (int p=2;p<stop;p++) for (int sign:{-1,1}) {
                    if (small_hits.count({p,sign})) continue;
                    auto [m,k]=scalar_widths(f,p,sign);
                    if (m>=1) {
                        small_hits.insert({p,sign});
                        std::cout<<"scalar E"<<p<<(sign>0?"+1":"-1")<<" triangle "<<line<<" ports "
                                 <<u<<" "<<v<<" widths "<<m<<" "<<k<<"\n";
                    }
                }
                for (int sign:{-1,1}) if (scalar_widths(f,61,sign).first>=1) ++hits;
            }
            for (int i:ray_offsets(f)) {
                if (std::abs(f.d)>100) continue;
                if (++offsets[f.d]==1)
                    std::cout<<"ray offset "<<f.d<<" first triangle "<<line<<" ports "<<u<<" "<<v
                             <<" width shifts from p-1 "<<i<<" "<<1-__builtin_ctzll(f.a)-i<<"\n";
            }
        }
        eligible_graphs+=eligible;
    }
    if (graphs==0) throw std::runtime_error("empty input");
    std::cout<<"graphs "<<graphs<<" eligible graphs "<<eligible_graphs<<" eligible pairs "<<pairs
             <<" distinct E61 polynomial/sign hits "<<hits<<" coefficient bound "<<coefficient_bound
             <<" asymptotic cutoff "<<(coefficient_bound?cutoff(coefficient_bound):0)<<"\nRays with absolute offset <=100:";
    for (auto [d,count]:offsets) std::cout<<" "<<d<<"="<<count;
    std::cout<<"\n";
}

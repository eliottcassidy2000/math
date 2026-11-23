// All rooted homogeneous-module forks, without an apex or spine hypothesis.
// Input is the upper-triangular ASCII tournament format from nauty gentourng.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

using Rows = std::vector<int>;
using Edge = std::pair<int,int>;
using Count = std::int64_t;
using Table = std::vector<Count>;

Table starts(const Rows& rows) {
    const int n = rows.size(), full = (1<<n)-1;
    Table dp((full+1)*n);
    for (int i=0; i<n; ++i) dp[(1<<i)*n+i]=1;
    for (int set=1; set<=full; ++set) {
        for (int v=set; v; v&=v-1) {
            int i=__builtin_ctz((unsigned)v), rest=set^(1<<i);
            Count& ways=dp[set*n+i];
            for (int next=rest&rows[i]; next; next&=next-1)
                ways+=dp[rest*n+__builtin_ctz((unsigned)next)];
        }
    }
    return dp;
}

Count height(const Rows& rows) {
    auto dp=starts(rows);
    Count total=0;
    for (size_t i=dp.size()-rows.size(); i<dp.size(); ++i) total+=dp[i];
    return total;
}

Rows reverse_edge(Rows rows, Edge edge) {
    auto [i,j]=edge;
    rows[i]^=1<<j;
    rows[j]^=1<<i;
    return rows;
}

Rows expand_port(Rows rows, int port, int width) {
    int n=rows.size();
    for (int k=1; k<width; ++k) {
        int last=n+k-1;
        rows.push_back(rows[port] & ((1<<n)-1));
        for (int i=0; i<last; ++i)
            if (i==port || i>=n || (rows[i]&(1<<port))) rows[i]|=1<<last;
    }
    return rows;
}

struct Survey {
    int n;
    long classes=0, base_centers=0, marked=0;
    std::map<std::vector<Count>, long> functions;

    std::map<Count,long> parents;
    bool print_centers=false;

    void examine(const Rows& rows) {
        ++classes;
        const int full=(1<<n)-1;
        Table forward=starts(rows);
        Count mid=0;
        for (int i=0; i<n; ++i) mid+=forward[full*n+i];
        if (mid%3) return;

        Rows converse(n);
        for (int i=0; i<n; ++i) converse[i]=full^rows[i]^(1<<i);
        Table backward=starts(converse);
        // An ordered partition ending at i and then starting at j counts
        // the paths gained or lost when the orientation of {i,j} changes.
        Count consecutive[10][10]{};
        for (int left=1; left<full; ++left) {
            int right=full^left;
            for (int v=left; v; v&=v-1) {
                int i=__builtin_ctz((unsigned)v);
                Count ways=backward[left*n+i];
                if (!ways) continue;
                for (int u=right; u; u&=u-1) {
                    int j=__builtin_ctz((unsigned)u);
                    consecutive[i][j]+=ways*forward[right*n+j];
                }
            }
        }
        std::array<std::vector<Edge>,3> choices;
        for (int i=0; i<n; ++i) for (int j=i+1; j<n; ++j) {
            Count diff=consecutive[j][i]-consecutive[i][j];
            if (!(rows[i]&(1<<j))) diff=-diff;
            if (diff==-2*mid/3) choices[0].emplace_back(i,j);
            if (diff==-2) choices[1].emplace_back(i,j);
            if (diff==2) choices[2].emplace_back(i,j);
        }
        for (auto& edges:choices) if (edges.empty()) return;
        ++base_centers;
        ++parents[mid/3];
        if (print_centers) {
            std::cout<<"base "<<mid;
            for (int row:rows) std::cout<<" "<<row;
            std::cout<<std::endl;
        }
        for (int port=0; port<n; ++port) {
            auto edges=choices;
            for (auto& group:edges)
                group.erase(std::remove_if(group.begin(),group.end(),[port](Edge e) {
                    return e.first==port || e.second==port;
                }),group.end());
            int out=__builtin_popcount((unsigned)rows[port]);
            int max_slots=std::min(out,n-1-out)+1;
            if (max_slots<2) continue;
            std::vector<Count> samples{mid/3};
            for (int s=2; s<=max_slots; ++s) {
                if (edges[0].empty() || edges[1].empty() || edges[2].empty()) break;
                Rows blown=expand_port(rows,port,s);
                Count center=height(blown);
                if (center%3) { edges[0].clear(); break; }
                samples.push_back(center/3);
                const Count want[3]={center/3,center-2,center+2};
                for (int k=0;k<3;++k) {
                    auto& group=edges[k];
                    group.erase(std::remove_if(group.begin(),group.end(),[&](Edge e) {
                        return height(reverse_edge(blown,e))!=want[k];
                    }),group.end());
                }
            }
            if (edges[0].empty() || edges[1].empty() || edges[2].empty()) continue;
            // Recover exterior coefficients of ordered nonempty path covers.
            std::vector<Count> coeff(samples.size()), ordered{1};
            for (size_t s=1;s<=samples.size();++s) {
                ordered.push_back(0);
                for (size_t k=s;k>=1;--k) ordered[k]=k*(ordered[k]+ordered[k-1]);
                ordered[0]=0;
                Count residual=samples[s-1];
                for (size_t k=1;k<s;++k) residual-=ordered[k]*coeff[k-1];
                if (residual<0 || residual%ordered[s]) std::abort();
                coeff[s-1]=residual/ordered[s];
            }
            while (!coeff.empty() && !coeff.back()) coeff.pop_back();
            if (coeff.size()<2) continue;
            ++marked;
            if (functions[coeff]++) continue;
            std::cout<<"profile";
            for (Count x:coeff) std::cout<<" "<<x;
            std::cout<<" center";
            for (int row:rows) std::cout<<" "<<row;
            std::cout<<" port "<<port<<" arcs";
            for (auto& group:edges) std::cout<<" "<<group.front().first<<","<<group.front().second;
            std::cout<<std::endl;
        }
    }
};

int main(int argc,char** argv) {
    if (argc!=2 && !(argc==3 && std::string(argv[2])=="--bases")) return 2;
    int n=std::stoi(argv[1]);
    if (n<2 || n>10) return 2;
    Survey survey{n,0,0,0,{},{},argc==3};
    std::string bits;
    while (std::cin>>bits) {
        if ((int)bits.size()!=n*(n-1)/2) return 2;
        Rows rows(n);
        int k=0;
        for (int i=0;i<n;++i) for (int j=i+1;j<n;++j) {
            if (bits[k]=='1') rows[i]|=1<<j;
            else if (bits[k]=='0') rows[j]|=1<<i;
            else return 2;
            ++k;
        }
        survey.examine(rows);
    }
    for (auto [parent,frequency]:survey.parents)
        std::cout<<"parent-height "<<parent<<" ordinary-centers "<<frequency<<"\n";
    std::cout<<"order "<<n<<" classes "<<survey.classes<<" base-centers "<<survey.base_centers
             <<" growing-ports "<<survey.marked<<" profiles "<<survey.functions.size()<<std::endl;
}

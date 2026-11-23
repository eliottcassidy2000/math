// Find marked arc neighborhoods that remain at offsets -2,0,+2 after any dyadic padding.
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using Rows = std::vector<int>;
using Key = std::pair<int,int>; // dyadic coefficient and constant
using Edge = std::pair<int,int>;
struct Witness { Rows rows; int marks; Edge lower,upper; };

using Both = std::array<Key,2>;

std::vector<Both> readings(const Rows& rows, int& base_height) {
    const int n=rows.size(), limit=1<<n, full=limit-1;
    std::vector<std::array<int,8>> starts(limit),ends(limit);
    std::vector<int> paths(limit);
    paths[0]=1;
    for(int subset=1;subset<limit;++subset) {
        for(int i=0;i<n;++i) if(subset&(1<<i)) {
            int rest=subset^(1<<i);
            if(!rest) starts[subset][i]=ends[subset][i]=1;
            for(int j=0;j<n;++j) if(rest&(1<<j)) {
                if(rows[i]&(1<<j)) starts[subset][i]+=starts[rest][j];
                else ends[subset][i]+=ends[rest][j];
            }
            paths[subset]+=ends[subset][i];
        }
    }
    base_height=paths[full];
    std::array<int,8> weight{}, bias{};
    std::array<std::array<int,8>,8> contact{};
    for(int subset=1;subset<=full;++subset) {
        int rest=full^subset;
        for(int j=0;j<n;++j) if(subset&(1<<j)) {
            int contribution=starts[subset][j]*paths[rest];
            weight[j]+=contribution;
            bias[j]+=__builtin_popcount((unsigned)subset)%2 ? contribution : -contribution;
            for(int i=0;i<n;++i) if(rest&(1<<i))
                contact[i][j]+=ends[rest][i]*starts[subset][j];
        }
    }
    std::vector<Both> result(limit);
    for(int mask=0;mask<limit;++mask) {
        int c=0, e=0, height=0;
        for(int j=0;j<n;++j) {
            if(mask&(1<<j)) {
                c+=weight[j]; e+=bias[j];
                height+=starts[full][j];
                for(int i=0;i<n;++i) if(!(mask&(1<<i))) height+=contact[i][j];
            } else height+=ends[full][j];
        }
        result[mask]={Key{c+base_height,height}, Key{c,height+e-c}};
    }
    return result;
}

int read(int x) {
    int total=0;
    while(x) {
        int v=x/2;
        if(!(x&1)) total+=v+(v&1);
        x=v;
    }
    return total;
}

struct Change {Edge arc; int dh, dk;};
struct Fork {
    Witness graph;
    int word, mode;
    Key chart; // actual slope and intercept in the run parameter 2^s
};


// Independent subset DP on full labelled tournaments supplied by the cube replay.
int full_counts() {
    int n;
    while(std::cin>>n) {
        if(n<1 || n>16) return 2;
        Rows rows(n);
        for(int& r:rows) if(!(std::cin>>r)) return 2;
        const int limit=1<<n;
        std::vector<unsigned long long> dp(size_t(limit)*n);
        for(int i=0;i<n;++i) dp[size_t(1<<i)*n+i]=1;
        for(int mask=1;mask<limit;++mask) for(int rem=mask;rem;rem&=rem-1) {
            int i=__builtin_ctz((unsigned)rem);
            auto ways=dp[size_t(mask)*n+i];
            if(!ways) continue;
            for(int next=rows[i]&~mask;next;next&=next-1) {
                int j=__builtin_ctz((unsigned)next);
                dp[size_t(mask|(1<<j))*n+j]+=ways;
            }
        }
        unsigned long long total=0;
        for(int i=0;i<n;++i) total+=dp[size_t(limit-1)*n+i];
        std::cout<<total<<"\n";
    }
    return 0;
}
int main(int argc,char** argv) {
    if(argc==2 && std::string(argv[1])=="--heights") return full_counts();
    if(argc!=4) return 2;
    int largest=std::stoi(argv[2]), max_word=std::stoi(argv[3]);
    if(largest<1 || largest>8 || max_word<1 || max_word>127) return 2;
    std::array<std::vector<Key>,2> tails;
    for(int w=1;w<=max_word;w+=2) {
        tails[0].emplace_back(w,w-1-read(w));
        tails[1].emplace_back(w,-read(w-1));
    }
    std::map<Key,Fork> forks;
    for(int n=1;n<=largest;++n) {
        std::ifstream file(std::string(argv[1])+"/classes"+std::to_string(n)+".rows");
        if(!file) return 1;
        std::string line;
        int tested=0, centers=0, basic=0;
        while(std::getline(file,line)) {
            std::istringstream input(line); Rows rows(n);
            for(int& entry:rows) if(!(input>>entry)) return 1;
            int base;
            auto core=readings(rows,base);
            std::vector<std::array<std::vector<Change>,2>> moves(core.size());
            auto compare=[&](int m,Both other,int h,Edge edge) {
                for(int t=0;t<2;++t) if(core[m][t].first>0 && other[t].first==core[m][t].first)
                    moves[m][t].push_back({edge,h-base,other[t].second-core[m][t].second});
            };
            for(int i=0;i<n;++i) {
                for(int m=0;m<(1<<n);++m) compare(m,core[m^(1<<i)],base,{i,n});
                for(int j=i+1;j<n;++j) {
                    Rows changed=rows;
                    changed[i]^=1<<j; changed[j]^=1<<i;
                    int h; auto other=readings(changed,h);
                    for(int m=0;m<(1<<n);++m) compare(m,other[m],h,{i,j});
                }
            }
            for(int m=0;m<(1<<n);++m) {
                ++tested;
                for(int t=0;t<2;++t) if(moves[m][t].size()>=2) for(auto [w,d]:tails[t]) {
                    Edge lo{-1,-1},hi{-1,-1};
                    for(auto change:moves[m][t]) {
                        int delta=change.dk-d*change.dh;
                        if(delta==-2) lo=change.arc;
                        if(delta==2) hi=change.arc;
                    }
                    if(lo.first<0 || hi.first<0) continue;
                    ++centers;
                    if(w==1 && t==0) ++basic;
                    Key key={core[m][t].first*w,core[m][t].second-base*d};
                    forks.emplace(key,Fork{{rows,m,lo,hi},w,t,key});
                }
            }
        }
        std::cout<<"order "<<n<<" tagged "<<tested<<" base forks "<<basic
                 <<" all forks "<<centers<<" cumulative charts "<<forks.size()<<std::endl;
    }

    auto show=[](const Witness& w) {
        std::cout<<" rows";for(int v:w.rows) std::cout<<" "<<v;
        std::cout<<" mask "<<w.marks<<" arcs "<<w.lower.first<<","<<w.lower.second
                 <<" "<<w.upper.first<<","<<w.upper.second;
    };
    // Two slopes agree up to a shift in s iff their odd parts agree.
    auto unit=[](int L) {while(!(L&1)) L/=2; return L;};
    std::map<Key,std::vector<Key>> wanted;
    for(auto& [key, fork]:forks) {
        if(key.first%3==0 && key.second%3==0)
            wanted[{unit(key.first/3),key.second/3}].push_back(key);
    }
    std::set<Key> found;
    for(int n=1;n<=largest;++n) {
        std::ifstream file(std::string(argv[1])+"/classes"+std::to_string(n)+".rows");
        std::string line;
        while(std::getline(file,line)) {
            std::istringstream input(line);Rows rows(n);
            for(int& x:rows) input>>x;
            int base;auto data=readings(rows,base);
            for(int m=0;m<(1<<n);++m) {
                for(int t=0;t<2;++t) if(data[m][t].first>0) for(auto [w,d]:tails[t]) {
                    Key key={data[m][t].first*w,data[m][t].second-base*d};
                    auto it=wanted.find({unit(key.first),key.second});
                    if(it==wanted.end()) continue;
                    for(auto f:it->second) if(found.insert(f).second) {
                        std::cout<<"parent "<<key.first<<" "<<key.second<<" word "<<w<<" mode "<<t;
                        show({rows,m,{-1,-1},{-1,-1}});
                        std::cout<<" child "<<f.first<<" "<<f.second<<" word "<<forks.at(f).word<<" mode "<<forks.at(f).mode;
                        show(forks.at(f).graph);std::cout<<std::endl;
                    }
                }
            }
        }
    }
    for(auto& [key, f]:forks) {
        std::cout<<"fork "<<key.first<<" "<<key.second<<" word "<<f.word<<" mode "<<f.mode;show(f.graph);std::cout<<"\n";
    }
    std::cout<<"matched charts "<<found.size()<<" of "<<forks.size()<<"\n";
}

// Read the upper triangles emitted by gentourng -cq; count target fibers.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
using Count = std::int64_t;
using Signature = std::pair<int,int>;
using Profile = std::vector<Count>;

int main(int argc,char**argv) {
    if(argc<4 || argc%2) return 2;
    const int n=std::stoi(argv[1]);
    if(n<3 || n>9) return 2;
    std::set<Signature> desired;
    for(int i=2;i<argc;i+=2) desired.emplace(std::stoi(argv[i]),std::stoi(argv[i+1]));
    std::map<std::pair<Signature,Profile>,int> multiplicity;
    const int size=1<<n,full=size-1;
    std::string word;
    int inputs=0;
    while(std::cin>>word) {
        if(int(word.size())!=n*(n-1)/2) return 1;
        ++inputs;
        int incoming[9]={},pos=0;
        for(int i=0;i<n;++i) for(int j=i+1;j<n;++j) {
            char bit=word[pos++];
            if(bit=='1') incoming[j]|=1<<i;
            else if(bit=='0') incoming[i]|=1<<j;
            else return 1;
        }
        int ends[512][9],h[512];
        h[0]=1;
        for(int mask=1;mask<size;++mask) {
            int total=0;
            for(int vs=mask;vs;vs&=vs-1) {
                int v=__builtin_ctz(vs),rest=mask^(1<<v);
                int count=(rest==0);
                for(int us=rest&incoming[v];us;us&=us-1) count+=ends[rest][__builtin_ctz(us)];
                ends[mask][v]=count;
                total+=count;
            }
            h[mask]=total;
        }
        if(desired.lower_bound({h[full],0})==desired.upper_bound({h[full],INT32_MAX})) continue;
        int q=0;
        for(int s=0;s<size;++s) q+=h[s]*h[full^s];
        const Signature signature={h[full],q};
        if(!desired.count(signature)) continue;
        if(n==8 && h[full]==331 && q==5368) {
            std::vector<int> deletions;
            for(int i=0;i<n;++i) deletions.push_back(h[full^(1<<i)]);
            std::sort(deletions.begin(),deletions.end());
            std::cout<<"sample "<<word;
            for(int a:deletions) std::cout<<" "<<a;
            std::cout<<"\n";
        }

        // Unordered path covers. The path containing the least vertex singles
        // out one block of a set partition. Its internal order contributes H.
        Count covers[512][10]={};
        covers[0][0]=1;
        for(int mask=1;mask<size;++mask) {
            int anchor=mask&-mask,other=mask^anchor;
            for(int sub=other;;sub=(sub-1)&other) {
                int chosen=sub|anchor,rest=mask^chosen;
                for(int j=0;j<=__builtin_popcount(unsigned(rest));++j)
                    covers[mask][j+1]+=h[chosen]*covers[rest][j];
                if(!sub) break;
            }
        }
        Profile profile(covers[full]+1,covers[full]+n+1);
        ++multiplicity[{signature,profile}];
    }
    std::cout<<"inputs "<<inputs<<"\n";
    for(auto &[key,m]:multiplicity) {
        auto [s,p]=key;
        std::cout<<s.first<<" "<<s.second<<" "<<m;
        for(auto c:p) std::cout<<" "<<c;
        std::cout<<"\n";
    }
}

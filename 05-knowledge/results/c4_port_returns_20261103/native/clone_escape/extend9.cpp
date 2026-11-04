#include <algorithm>
#include <bit>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
using T = vector<unsigned>;

long long H(const T& rows) {
    const int n=rows.size(), size=1<<n;
    vector<long long> dp(size*n);
    for(int v=0;v<n;v++) dp[(1<<v)*n+v]=1;
    for(int mask=1;mask<size;mask++) for(int v=0;v<n;v++) {
        long long c=dp[mask*n+v];
        if(!c) continue;
        unsigned available=rows[v]&~mask;
        while(available) {
            unsigned bit=available&-available;
            available-=bit;
            dp[(mask|bit)*n+countr_zero(bit)]+=c;
        }
    }
    long long out=0;
    for(int v=0;v<n;v++) out+=dp[(size-1)*n+v];
    return out;
}

bool strong(const T& rows) {
    T reach=rows; int n=rows.size();
    for(int i=0;i<n;i++) reach[i]|=1u<<i;
    for(int k=0;k<n;k++) for(int i=0;i<n;i++)
        if(reach[i]>>k&1) reach[i]|=reach[k];
    return all_of(reach.begin(),reach.end(),[&](unsigned r){return r==(1u<<n)-1;});
}

uint64_t canonical(const T& rows) {
    int n=rows.size(); vector<int> order(n);
    for(int i=0;i<n;i++) order[i]=i;
    sort(order.begin(),order.end(),[&](int a,int b){
        return popcount(rows[a])<popcount(rows[b]) ||
               (popcount(rows[a])==popcount(rows[b]) && a<b);});
    vector<pair<int,int>> groups;
    for(int i=0;i<n;) {
        int j=i+1;
        while(j<n && popcount(rows[order[i]])==popcount(rows[order[j]]))j++;
        groups.emplace_back(i,j); i=j;
    }
    uint64_t best=UINT64_MAX;
    function<void(int)> search=[&](int g){
        if(g==(int)groups.size()) {
            uint64_t code=0; int k=0;
            for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++)
                if(rows[order[i]]>>order[j]&1)code|=uint64_t(1)<<k;
            best=min(best,code); return;
        }
        auto [a,b]=groups[g]; sort(order.begin()+a,order.begin()+b);
        do {search(g+1);} while(next_permutation(order.begin()+a,order.begin()+b));
    };
    search(0); return best;
}

int main() {
    int count; cin>>count;
    map<uint64_t,T> targets;
    long long births=0,strong_births=0,target_births=0;
    for(int p=0;p<count;p++) {
        T parent(8); for(auto& r:parent)cin>>r;
        if(!strong(parent) || H(parent)>97) return 2;
        for(unsigned cut=0;cut<256;cut++) {
            T child=parent;
            for(int i=0;i<8;i++)if(cut>>i&1)child[i]|=1u<<8;
            child.push_back(255u^cut); births++;
            if(!strong(child))continue;
            strong_births++;
            long long h=H(child);
            if(h!=91 && h!=97)continue;
            target_births++;
            targets.emplace(canonical(child),child);
        }
    }
    cout<<"{\"parents\":"<<count<<",\"births\":"<<births
        <<",\"strong_births\":"<<strong_births<<",\"target_births\":"<<target_births
        <<",\"representatives\":[";
    bool comma=false;
    for(auto&[code,rows]:targets) {
        if(comma)cout<<','; comma=true;
        cout<<"{\"n\":9,\"code\":"<<code<<",\"h\":"<<H(rows)
            <<",\"strong\":true,\"rows\":[";
        for(int i=0;i<9;i++){if(i)cout<<',';cout<<rows[i];}
        cout<<"]}";
    }
    cout<<"]}\n";
}

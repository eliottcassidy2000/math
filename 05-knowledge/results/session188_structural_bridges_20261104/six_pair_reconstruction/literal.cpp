#include <algorithm>
#include <cstdint>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
using namespace std;
using Count=uint64_t;
void need(bool ok,const char* why){if(!ok)throw runtime_error(why);}
vector<unsigned> quotient(unsigned code){
    vector<unsigned> a(6);unsigned k=0;
    for(unsigned i=0;i<6;i++)for(unsigned j=i+1;j<6;j++,k++){
        if(code>>k&1)a[i]|=1u<<j;else a[j]|=1u<<i;
    }
    return a;
}
unsigned expand(unsigned mask){unsigned out=0;for(unsigned i=0;i<6;i++)if(mask>>i&1)out|=3u<<(2*i);return out;}
vector<unsigned> core(unsigned code,unsigned flips){
    auto q=quotient(code);vector<unsigned> a(12);
    for(unsigned i=0;i<6;i++){
        if(flips>>i&1)a[2*i+1]|=1u<<(2*i);else a[2*i]|=1u<<(2*i+1);
        for(unsigned j=0;j<6;j++)if(q[i]>>j&1){a[2*i]|=3u<<(2*j);a[2*i+1]|=3u<<(2*j);}
    }
    return a;
}
vector<unsigned> attach(vector<unsigned>a,unsigned am,unsigned cm){
    unsigned n=a.size();a.push_back(1u<<(n+1));a.push_back(0);
    for(unsigned i=0;i<n;i++){
        if(am>>i&1)a[i]|=1u<<n;else a[n]|=1u<<i;
        if(cm>>i&1)a[i]|=1u<<(n+1);else a[n+1]|=1u<<i;
    }
    return a;
}
Count hp(const vector<unsigned>&a){
    unsigned n=a.size(),N=1u<<n;need(n<=16,"uint64 exactness boundary");
    vector<Count>d(size_t(N)*n);
    for(unsigned i=0;i<n;i++)d[size_t(1u<<i)*n+i]=1;
    for(unsigned mask=1;mask<N;mask++){
        unsigned lastbits=mask;
        while(lastbits){unsigned bit=lastbits&-lastbits;lastbits-=bit;unsigned last=__builtin_ctz(bit);Count value=d[size_t(mask)*n+last];if(!value)continue;
            unsigned next=a[last]&~mask;
            while(next){unsigned b=next&-next;next-=b;d[size_t(mask|b)*n+__builtin_ctz(b)]+=value;}
        }
    }
    Count out=0;for(unsigned i=0;i<n;i++)out+=d[size_t(N-1)*n+i];return out;
}
vector<pair<unsigned,unsigned>> battery(){
    set<pair<unsigned,unsigned>>low;
    for(unsigned a=0;a<64;a++)if(__builtin_popcount(a)<=2)low.emplace(a,0);
    for(unsigned j=0;j<6;j++)for(unsigned a=0;a<64;a++)if(__builtin_popcount(a)<=1||(__builtin_popcount(a)==2&&(a>>j&1)))low.emplace(a,1u<<j);
    need(low.size()==94,"low battery");auto all=low;
    for(auto [a,c]:low)all.emplace(63^c,63^a);
    need(all.size()==188,"complete battery");return {all.begin(),all.end()};
}
vector<unsigned> ambient(vector<unsigned>a){
    unsigned n=a.size();a.push_back(1u<<(n+1));a.push_back(0);
    for(unsigned i=0;i<n;i++){a[i]|=1u<<n;a[n+1]|=1u<<i;}
    return a;
}
int main(){
    auto b=battery();cout<<"{\"status\":\"PASS\",\"coarse_battery\":[";
    for(unsigned i=0;i<b.size();i++){if(i)cout<<",";cout<<"["<<b[i].first<<","<<b[i].second<<"]";}cout<<"],\"states\":[";
    vector<pair<unsigned,unsigned>>states={{0,0},{81,0},{32767-81,0},{81,21}};
    vector<vector<Count>> values;
    for(unsigned s=0;s<states.size();s++){
        auto [code,flips]=states[s];if(s)cout<<",";cout<<"{\"quotient_code\":"<<code<<",\"pair_flip_mask\":"<<flips<<",\"exact_values\":[";
        vector<Count>v;
        for(unsigned k=0;k<b.size();k++){auto [a,c]=b[k];auto count=hp(attach(core(code,flips),expand(a),expand(c)));v.push_back(count);if(k)cout<<",";cout<<count;}
        values.push_back(v);cout<<"]}";
    }
    need(values[1]==values[3],"coarse pair flip invariance");
    cout<<"],\"pair_decoder_exact_values\":[";
    auto q=quotient(81);vector<unsigned>masks;
    for(unsigned i=0;i<6;i++){
        unsigned mask=1u<<(2*i);
        for(unsigned j=0;j<6;j++)if(q[j]>>i&1)mask|=3u<<(2*j);
        masks.push_back(mask);
    }
    vector<Count>base;set<vector<unsigned>>signatures;
    for(unsigned flips=0;flips<64;flips++){
        if(flips){cout<<",";}
        cout<<"[";vector<unsigned>sig;
        for(unsigned i=0;i<6;i++){
            Count v=hp(attach(core(81,flips),masks[i],masks[i]));if(i)cout<<",";cout<<v;
            if(flips==0)base.push_back(v);
            need(((v%8+8-base[i]%8)%8)==((flips>>i&1)?4u:0u),"pair isolation gap");sig.push_back(v%8);
        }
        need(signatures.insert(sig).second,"pair decoder signature collision");cout<<"]";
    }
    cout<<"],\"pair_decoder_physical_masks\":[";for(unsigned i=0;i<6;i++){if(i)cout<<",";cout<<masks[i];}
    unsigned witness=0;while(witness<b.size()&&values[1][witness]%8==values[2][witness]%8)witness++;need(witness<b.size(),"root separator exists");auto [a,c]=b[witness];
    Count f=hp(ambient(attach(core(81,0),expand(a),expand(c))));Count g=hp(ambient(attach(core(32767-81,0),expand(a),expand(c))));
    unsigned gap=(values[1][witness]%8+8-values[2][witness]%8)%8;
    need((f%8+8-g%8)%8==(5*gap)%8,"strong ambient multiplier");
    cout<<"],\"strong_ambient_control\":{\"coarse_a\":"<<a<<",\"coarse_c\":"<<c<<",\"core_source\":"<<values[1][witness]<<",\"core_target\":"<<values[2][witness]<<",\"ambient_source\":"<<f<<",\"ambient_target\":"<<g<<",\"multiplier_mod8\":5},\"literal_query_count\":"<<(4*188+64*6+2)<<"}\n";
}

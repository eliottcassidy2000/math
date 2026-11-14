#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <vector>
#include <stdexcept>
using namespace std;
using G=vector<int>;using P=pair<long long,long long>;
P profile_forward(const G&a){int n=a.size(),N=1<<n;vector<long long>E(N*n),H(N);H[0]=1;for(int i=0;i<n;i++)E[(1<<i)*n+i]=1;for(int s=1;s<N;s++)for(int u=s;u;u&=u-1){int i=__builtin_ctz((unsigned)u);long long e=E[s*n+i];H[s]+=e;if(e)for(int z=a[i]&~s;z;z&=z-1){int j=__builtin_ctz((unsigned)z);E[(s|(1<<j))*n+j]+=e;}}long long q=0;for(int s=0;s<N;s++)q+=H[s]*H[N-1-s];return{H[N-1],q};}
P profile_incoming(const G&a) {
    int n=a.size(),N=1<<n;
    vector<long long> E(N*n),H(N); vector<int> in(n); H[0]=1;
    for(int i=0;i<n;i++) { E[(1<<i)*n+i]=1; for(int j=0;j<n;j++) if(a[j]>>i&1) in[i]|=1<<j; }
    for(int s=1;s<N;s++) for(int u=s;u;u&=u-1) {
        int i=__builtin_ctz((unsigned)u),rest=s^(1<<i);
        if(rest) for(int z=in[i]&rest;z;z&=z-1) E[s*n+i]+=E[rest*n+__builtin_ctz((unsigned)z)];
        H[s]+=E[s*n+i];
    }
    long long q=0; for(int s=0;s<N;s++) q+=H[s]*H[N-1-s]; return {H[N-1],q};
}
P profile(const G&a) {
#ifdef INCOMING
    return profile_incoming(a);
#else
    return profile_forward(a);
#endif
}
G canon(const G&a){int n=a.size();map<int,G>gs;for(int i=0;i<n;i++)gs[__builtin_popcount((unsigned)a[i])].push_back(i);vector<G>groups;for(auto&kv:gs)groups.push_back(kv.second);G p,best;function<void(int)>go=[&](int j){if(j==(int)groups.size()){G r(n);for(int i=0;i<n;i++)for(int k=0;k<n;k++)if(a[p[i]]>>p[k]&1)r[i]|=1<<k;if(best.empty()||r<best)best=r;return;}auto&z=groups[j];do{int old=p.size();p.insert(p.end(),z.begin(),z.end());go(j+1);p.resize(old);}while(next_permutation(z.begin(),z.end()));};go(0);return best;}
map<pair<int,int>,vector<G>>memo;
const vector<G>&catalog(int n,int cap){auto key=make_pair(n,cap);if(memo.count(key))return memo.at(key);vector<G>ans;if(cap<(1<<n))return memo.emplace(key,ans).first->second;if(n==0){ans.push_back({});return memo.emplace(key,ans).first->second;}if(cap<3*(1<<(n-1))){G T(n);for(int i=0;i<n;i++)T[i]=((1<<n)-1)^((1<<(i+1))-1);ans.push_back(canon(T));return memo.emplace(key,ans).first->second;}set<G>seen;for(const auto&A:catalog(n-1,cap/2)){for(int cut=0;cut<(1<<(n-1));cut++){G T=A;for(int i=0;i<n-1;i++)if(!(cut>>i&1))T[i]|=1<<(n-1);T.push_back(cut);auto[h,q]=profile(T);if(q<=cap)seen.insert(canon(T));}}ans.assign(seen.begin(),seen.end());cerr<<"catalog "<<n<<' '<<cap<<' '<<ans.size()<<'\n';return memo.emplace(key,ans).first->second;}
int main(int argc,char**argv){int cap=argc>1?atoi(argv[1]):4096;map<long long,map<long long,G>>P;vector<int> levels;for(int n=1;(1<<n)<=cap;n++){const auto&T=catalog(n,cap);levels.push_back(T.size());cerr<<"FINAL "<<n<<' '<<T.size()<<'\n';for(const auto&A:T){auto[h,q]=profile(A);P[q].try_emplace(h,A);}}if(cap==4096) {
    if(levels!=vector<int>({1,1,2,4,12,56,456,2160,716,151,18,1})) throw runtime_error("class totals");
    int total=0; vector<array<long long,3>> phase,sign;
    for(auto&[q,hs]:P) { total+=hs.size(); for(auto i=hs.begin();i!=hs.end();i++) for(auto j=next(i);j!=hs.end();j++) {
        auto a=i->first,b=j->first; if(gcd(a,99LL)!=1||gcd(b,99LL)!=1)continue;
        if((b-a)%396==198) phase.push_back({q,a,b});
        if((a+b)%396==0) sign.push_back({q,a,b});
    }}
    if(total!=919 || phase!=vector<array<long long,3>>({{4096,19,217},{4096,25,223}}) || !sign.empty()) throw runtime_error("phase minima");
}
for(auto&[q,hs]:P){for(auto&[h,A]:hs){cout<<q<<' '<<h<<' '<<A.size();for(int x:A)cout<<' '<<x;cout<<'\n';}for(auto i=hs.begin();i!=hs.end();i++)for(auto j=next(i);j!=hs.end();j++){auto a=i->first,b=j->first;if(gcd(a,99LL)!=1||gcd(b,99LL)!=1)continue;if((b-a)%396==198)cerr<<"PHASE "<<q<<' '<<a<<' '<<b<<'\n';if((a+b)%396==0)cerr<<"SIGN "<<q<<' '<<a<<' '<<b<<'\n';}}}

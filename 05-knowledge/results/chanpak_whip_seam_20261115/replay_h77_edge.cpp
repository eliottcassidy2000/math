#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <unordered_set>
#include <vector>
using namespace std;
using Row=array<int,10>;
using U=uint64_t;
int n;
U canonical(const Row& a, vector<int> p) {
    while (true) {
        vector<pair<U,int>> sig;
        for(int i=0;i<n;i++) {
            U x=0;
            for(int b:p) x=(x<<4)|__builtin_popcount(unsigned(a[i]&b));
            sig.emplace_back(x,i);
        }
        sort(sig.begin(),sig.end());
        vector<int> q;
        for(int b:p) {
            U old=~U(0);
            for(auto [key,i]:sig) if((b>>i)&1) {
                if(key!=old) q.push_back(0);
                q.back()|=1<<i;old=key;
            }
        }
        if(q==p) break;
        p=move(q);
    }
    if((int)p.size()==n) {
        U x=0;
        for(int i=0;i<n;i++) for(int j=i+1;j<n;j++)
            x=(x<<1)|!!(a[__builtin_ctz(unsigned(p[i]))]&p[j]);
        return x;
    }
    int i=0; while(__builtin_popcount(unsigned(p[i]))==1) i++;
    U best=~U(0);
    for(int b=p[i];b;b&=b-1) {
        int bit=b&-b;
        vector<int> q=p;
        q[i]^=bit;q.insert(q.begin()+i,bit);
        best=min(best,canonical(a,move(q)));
    }
    return best;
}
Row decode(U a,int n) {
    Row r={};
    for(int i=n-2;i>=0;i--) for(int j=n-1;j>i;j--) {
        if(a&1) r[i]|=1<<j;else r[j]|=1<<i;
        a>>=1;
    }
    return r;
}
int ham(const Row& a) {
    int paths[512][9]={};
    for(int m=1;m<(1<<n);m++) for(int b=m;b;b&=b-1) {
        int v=__builtin_ctz(unsigned(b));int s=m^(1<<v);
        if(!s) paths[m][v]=1;
        else for(int t=s&a[v];t;t&=t-1)
            paths[m][v]+=paths[s][__builtin_ctz(unsigned(t))];
    }
    int h=0;
    for(int v=0;v<n;v++)h+=paths[(1<<n)-1][v];
    return h;
}
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>
#include <functional>
#include <stdexcept>
using namespace std;
using ll=long long;
using P=array<int,7>; // a,Q,H_with_apex,s,e,t,f
using Key=pair<ll,ll>; // N,bias
vector<P> bank;
map<pair<int,ll>,set<Key>> memo;
bool strong(const vector<int>& r){
    int n=r.size(),f=(1<<n)-1;
    for(int v=0;v<n;v++){
        int b=1<<v;
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(b>>j&1)b|=r[j];
        if(b!=f)return false;
    }
    return true;
}
int cnt(const vector<int>& r){
    int n=r.size(),f=(1<<n)-1, a[1024][10]={};
    for(int s=1;s<=f;s++)for(int i=0;i<n;i++)if(s>>i&1){
        int b=s^(1<<i);
        if(!b)a[s][i]=1;
        else for(int j=b&r[i];j;j&=j-1)a[s][i]+=a[b][__builtin_ctz(unsigned(j))];
    }
    int h=0;for(int i=0;i<n;i++)h+=a[f][i];return h;
}
void enumerate(){
    set<P> seen; int checks=0;
    unordered_set<U> old={0};
    for(n=2;n<=9;n++) {
        unordered_set<U> nx;
        for(U code:old){
            Row a=decode(code,n-1);
            for(int b=0;b<(1<<(n-1));b++){
                Row q=a;q[n-1]=b;
                for(int j=0;j<n-1;j++) if(!((b>>j)&1)) q[j]|=1<<(n-1);
                nx.insert(canonical(q,{(1<<n)-1}));
            }
        }
        old=move(nx);
        for(U code:old){
        Row r0=decode(code,n);vector<int>r(n);for(int j=0;j<n;j++)r[j]=r0[j];
        int hh=ham(r0);if(hh!=3&&hh!=5&&hh!=7&&hh!=11&&hh!=77)continue;
        if(!strong(r))continue;
        if(77%hh)continue;
        int f=(1<<n)-1;
        int h[512]={1},DP[512][9][9]={};
        for(int s=1;s<=f;s++)for(int i=0;i<n;i++)if(s>>i&1){
            int b=s^(1<<i);
            if(!b)DP[s][i][i]=1;
            else for(int j=b&r[i];j;j&=j-1)
                for(int k=0;k<n;k++)DP[s][i][k]+=DP[b][__builtin_ctz(unsigned(j))][k];
            for(int k=0;k<n;k++)h[s]+=DP[s][i][k];
        }
        int S[9]={},D[9]={},T[9]={},F[9]={},K[9][9]={},Q=0;
        for(int s=0;s<=f;s++){
            Q+=h[s]*h[f^s];
            int plus=__builtin_popcount(unsigned(s))%2?1:-1;
            for(int i=0;i<n;i++)for(int j=0;j<n;j++){
                int z=DP[s][i][j]*h[f^s];
                S[i]+=z;D[i]+=plus*z;T[j]+=z;F[j]+=plus*z;
                if(plus<0)K[i][j]+=2*z;
            }
        }
        for(int m=0;m<=f;m++){
            P a={hh,Q,hh,0,0,0,0};
            for(int i=0;i<n;i++){
                if(m>>i&1){
                    a[3]+=S[i];a[4]+=D[i];
                    for(int j=0;j<n;j++)if(!(m>>j&1))a[2]+=K[i][j];
                }else {a[5]+=T[i];a[6]+=F[i];}
            }
            if(seen.insert(a).second){
                vector<int> full=r;
                for(int i=0;i<n;i++)if(!(m>>i&1))full[i]|=1<<n;
                full.push_back(m);
                if(cnt(full)!=a[2])throw runtime_error("literal mismatch");
                checks++;
            }
        }
    }
    }
    seen.insert({1,2,1,0,0,1,1});seen.insert({1,2,1,1,1,0,0});
    bank.assign(seen.begin(),seen.end());
    for(int a:{1,7,11,77}){
        int j=0;for(auto p:bank)j+=p[0]==a;cout<<"local signatures h="<<a<<": "<<j<<endl;
    }
    cout<<"direct literal checks for "<<checks<<" nonsingleton signatures"<<endl;
}
const set<Key>& search(int h,ll S){
    auto z=make_pair(h,S);
    if(memo.count(z))return memo[z];
    auto &out=memo[z];
    if(!S){out.emplace(h,0);return out;}
    for(auto a:bank){
        if(h%a[0])continue;
        int prior_h=h/a[0];
        ll j=S-prior_h*a[3];
        if(j<0||j%a[1])continue;
        for(auto [N,d]:search(prior_h,j/a[1])){
            if(a[5]*j/a[1]+a[6]*d<0)throw runtime_error("negative seam");
            out.emplace(a[0]*N+(a[2]-a[0])*prior_h+a[5]*(j/a[1])+a[6]*d,prior_h*a[4]);
        }
    }
    return out;
}
ll p(ll t){if(!t)return 0;ll a=t/2;return p(a)+(t%2?0:a+(a%2));}
ll Safter(ll S,int h,ll t){
    int L=0;ll w=0;
    for(ll j=t;j;j>>=1){w=2*w+(1-(j&1));L++;}
    return (S<<L)+h*w;
}

int main(){
    enumerate();
    ll D=0,C=0;
    for(auto p:bank)D=max(D,ll(77/p[0])*ll(abs(p[4])));
    set<int> zero;
    map<ll,set<ll>> hits;
    function<void(int,ll,ll,ll)> scan=[&](int h,ll N,ll S,ll d){
        if(N*(77/h)>177)return;
        C++;
        if(h==77&&(N==171||N==173||N==175||N==177)){
            ll x=S+h;hits[N].insert(x/(x&-x));
        }
        if(S>177+D)return;
        if(S==0&&!zero.insert(h).second)return;
        for(auto p:bank)if(77%(h*p[0])==0)
            scan(h*p[0],p[0]*N+(p[2]-p[0])*h+p[5]*S+p[6]*d,p[1]*S+h*p[3],h*p[4]);
    };
    scan(1,1,0,0);
    cout<<"bounded signature visits "<<C<<" Dbound "<<D<<endl;
    for(ll N:{171,173,175,177}){
        cout<<"N "<<N<<" normalized signatures "<<hits[N].size()<<" joins 181 "<<hits[N].count(181)<<endl;
        if(hits[N].count(181)!=(N==173||N==175))throw runtime_error("unexpected lane");
    }
}

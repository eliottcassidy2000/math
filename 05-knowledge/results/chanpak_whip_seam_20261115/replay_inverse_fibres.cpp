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
    int n=r.size(),f=(1<<n)-1, a[128][7]={};
    for(int s=1;s<=f;s++)for(int i=0;i<n;i++)if(s>>i&1){
        int b=s^(1<<i);
        if(!b)a[s][i]=1;
        else for(int j=b&r[i];j;j&=j-1)a[s][i]+=a[b][__builtin_ctz(unsigned(j))];
    }
    int h=0;for(int i=0;i<n;i++)h+=a[f][i];return h;
}
void enumerate(){
    set<P> seen;
    int checks=0;
    for(int n=3;n<=6;n++)for(int spec=0;spec<(1<<(n*(n-1)/2));spec++){
        vector<int> r(n);int v=spec,f=(1<<n)-1;
        for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
            if(v&1)r[i]|=1<<j;else r[j]|=1<<i;v>>=1;
        }
        if(!strong(r))continue;
        int hh=cnt(r);if(15%hh)continue;
        int h[64]={1},DP[64][6][6]={};
        for(int s=1;s<=f;s++)for(int i=0;i<n;i++)if(s>>i&1){
            int b=s^(1<<i);
            if(!b)DP[s][i][i]=1;
            else for(int j=b&r[i];j;j&=j-1)
                for(int k=0;k<n;k++)DP[s][i][k]+=DP[b][__builtin_ctz(unsigned(j))][k];
            for(int k=0;k<n;k++)h[s]+=DP[s][i][k];
        }
        int S[6]={},D[6]={},T[6]={},F[6]={},K[6][6]={},Q=0;
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
    seen.insert({1,2,1,0,0,1,1});seen.insert({1,2,1,1,1,0,0});
    bank.assign(seen.begin(),seen.end());
    for(int a:{1,3,5,15}){
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
void global(){
    const ll cap=150413;
    vector<ll> tg={150401,150403,150407,150409,150413};
    vector<set<ll>> K(tg.size());
    set<int> empty_seen;
    ll D=0,C=0;
    for(auto p:bank)D=max(D,ll(15/p[0])*ll(abs(p[4])));
    function<void(int,ll,ll,ll)> scan=[&](int h,ll N,ll S,ll d){
        if(N*(15/h)>cap)return;
        C++;
        if(h==15){
            auto q=find(tg.begin(),tg.end(),N);
            if(q!=tg.end()) {ll x=S+h;K[q-tg.begin()].insert(x/(x&-x));}
        }
        if(S>cap+D)return;
        if(S==0&&!empty_seen.insert(h).second)return;
        for(auto p:bank)if(15%(h*p[0])==0)
            scan(h*p[0],p[0]*N+(p[2]-p[0])*h+p[5]*S+p[6]*d,p[1]*S+h*p[3],h*p[4]);
    };
    scan(1,1,0,0);
    cout<<"global visits "<<C<<" Dbound "<<D<<endl;
    set<ll> hits=K[0],next;
    for(int i=0;i<(int)tg.size();i++){
        next.clear();
        set_intersection(K[i].begin(),K[i].end(),hits.begin(),hits.end(),inserter(next,next.end()));
        hits=next;
        cout<<"global target "<<tg[i]<<" normalized options "<<K[i].size()<<" common "<<hits.size()<<endl;
    }
    for(ll x:hits)cout<<"common "<<x<<endl;
}
int main(){
    enumerate();
    global();
    for(int j:{3,4,5,7,8,9,11,12,13}){
        int sh=(j%4==3)?2:(j%4==0?4:8);
        ll x=17LL<<j,S=2*(x-15),u=(S-sh)/15;
        if(S-15*u!=sh||u%2)throw runtime_error("fold");
        ll t=2*u, N=27+x+S*t+15*p(t);
        int width=4+sh/2;
        if(j<=9){
            ll lo=Safter(S,15,t+1),hi=Safter(S,15,t);
            int best=0,maxrun=0,ways=0;ll best_S=0;
            set<ll> vals;
            for(ll ss=lo;ss<=hi;ss++){
                auto &o=search(15,ss);
                set<ll> hit;
                for(auto [v,d]:o)if(!d&&v>=N&&v<=N+2*(width-1))hit.insert(v);
                if(hit.size()>best){best=hit.size();best_S=ss;vals=hit;}
                ways+=hit.size()==width;
                for(auto v:hit){int a=1;while(hit.count(v+2*a))a++;maxrun=max(maxrun,a);}
            }
            cout<<"INTERPOLATE "<<j<<" "<<lo<<".."<<hi<<" best "<<best<<" at "<<best_S
                <<" full "<<ways<<" run "<<maxrun<<":";
            for(auto a:vals)cout<<" "<<a;
            cout<<endl;
        }
        for(ll z:{t,t+1}){
            ll ss=Safter(S,15,z);
            auto &o=search(15,ss);
            set<ll> values;
            for(auto [v,d]:o)if(!d)values.insert(v);
            cout<<"height "<<j<<" address "<<z<<" S "<<ss<<" size "<<values.size()<<" interval";
            for(int g=0;g<width;g++)if(values.count(N+2*g))cout<<" "<<N+2*g;
            cout<<endl;
        }
    }
}

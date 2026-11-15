#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using R=std::vector<int>;
using Big=long long;
using S=std::array<Big,6>;

void check(bool b){if(!b)throw std::runtime_error("literal replay mismatch");}
R join(R a,const R& b){
    int n=a.size(),m=b.size(),mask=((1<<m)-1)<<n;
    for(int& x:a)x|=mask;
    for(int x:b)a.push_back(x<<n);
    return a;
}
R apex(R r,int m){
    int n=r.size();
    for(int i=0;i<n;i++)if(!(m>>i&1))r[i]|=1<<n;
    r.push_back(m);return r;
}
R opp(const R& r){
    R b(r.size());
    for(unsigned i=0;i<r.size();i++)for(unsigned j=0;j<r.size();j++)
        if(r[i]>>j&1)b[j]|=1<<i;
    return b;
}
struct Data{
    int n, N;
    std::vector<Big> d,h;
    explicit Data(const R& r):n(r.size()),N(1<<n),d(size_t(n)*N),h(N) {
        h[0]=1;
        for(int s=1;s<N;s++) for(int z=s;z;z&=z-1){
            int i=__builtin_ctz(unsigned(z)),t=s^(1<<i);
            Big zc=t==0;
            for(int u=t&r[i];u;u&=u-1)
                zc+=d[size_t(t)*n+__builtin_ctz(unsigned(u))];
            d[size_t(s)*n+i]=zc;
            h[s]+=zc;
        }
    }
    auto marked(int m) const {
        std::array<Big,2> se{};
        for(int s=1;s<N;s++)for(int b=s&m;b;b&=b-1){
            int i=__builtin_ctz(unsigned(b));Big w=d[size_t(s)*n+i]*h[(N-1)^s];
            se[0]+=w;se[1]+=(__builtin_popcount(unsigned(s))&1)?w:-w;
        }
        return se;
    }
};
S profile(const R& r,int m){
    Data a(r),b(opp(r));Big q=0;
    for(int i=0;i<a.N;i++)q+=a.h[i]*a.h[(a.N-1)^i];
    auto[s,e]=a.marked(m);auto[t,f]=b.marked((a.N-1)^m);
    return {a.h.back(),q,s,e,t,f};
}
int main(){
    std::array<R,2> A={R{126,124,120,48,96,64,8},R{118,124,120,113,96,64,0}};
    std::array<R,2> B={R{254,220,248,240,224,194,128,0},R{254,252,184,240,224,192,132,0}};
    int am[]={38,50},bm[]={54,74};
    for(int a=0;a<2;a++){
        A[a]=join(join(R{0},A[a]),R{0});am[a]=am[a]<<1|1<<8;
        check(profile(A[a],am[a])==S{5,1024,209,5,303,5});
        check(profile(B[a],bm[a])==S{9,672,236,0,545,9});
        check(Data(apex(A[a],am[a])).h.back()==67+2*a);
        check(Data(apex(B[a],bm[a])).h.back()==207+4*a);
    }
    std::cout<<"Module profiles h5: 67,69; h9: 207,211; direct DP passed."<<std::endl;
    for(int a=0;a<2;a++)for(int b=0;b<2;b++){
        auto r=join(A[a],B[b]);int m=am[a]|bm[b]<<9;
        check(profile(r,m)==S{45,688128,141628,0,560807,45});
        Big N=Data(apex(r,m)).h.back();
        check(N==115543+18*a+20*b);
        // Left new slot is marked, right new slot unmarked.
        r=join(join(R{0},r),R{0});m=m<<1|1;
        Big G=Data(apex(r,m)).h.back();
        check(G==N+1390608);
        std::cout<<a<<' '<<b<<": "<<N<<" -> "<<G<<std::endl;
    }
}

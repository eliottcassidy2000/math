// Exact cycle-packing DP for the ordered-spine/two-apex presentation.
// Build with c++ -std=c++17 -O3 two_apex_band.cpp -o two_apex_band.
#include <array>
#include <cstdint>
#include <iostream>
#include <unordered_set>
#include <algorithm>
using I=uint64_t;
struct Pair{I h,q; bool operator==(const Pair&o)const{return h==o.h&&q==o.q;}};
struct Hash{I operator()(Pair p)const{return p.h*0x9e3779b9u+p.q*0xc2b2ae35u;}};
I power(I a,int n){I p=1;while(n-->0)p*=a;return p;}
struct Dp{
    // Per run: 0 = empty, otherwise 1+2*parity+endGood.
    int next[2][2][5];
    Dp(){for(int b=0;b<2;b++)for(int d=0;d<2;d++)for(int s=0;s<5;s++)
         next[b][d][s]=s?1+2*(1-(s-1)/2)+d:(b?3+d:-1);}
    I cycles(int m,I a,I b,int z,bool split){
        std::array<I,25> x{},y{};x[0]=1;
        for(int i=0;i<m;i++){
            int v1=a>>i&1,v2=b>>i&1;
            auto u=next[v1][split?1-v1:1-v2];
            auto v=next[v2][split?1-v2:1-v1];
            for(int s=0;s<25;s++)y[s]=z*x[s];
            for(int s=0;s<5;s++)for(int t=0;t<5;t++)if(I old=x[5*s+t]){
                if(u[s]>=0) y[5*u[s]+t]+=old;
                if(v[t]>=0) y[5*s+v[t]]+=old;
            }
            x=y;
        }
        // A -> B, so the first run may be empty; second cannot be.
        return split? x[5*2+2] : x[5*0+4]+x[5*2+4]+x[5*4+2];
    }
    I value(int m,I a,I b,int z){
        I f=power(z,m+2), full=(I(1)<<m)-1;
        for(int d=1;d<m;d++){
            I cnt=__builtin_popcountll((a<<d)&(full^a))+
                  __builtin_popcountll((b<<d)&(full^b));
            f+=cnt*power(z,m+1-d)*(power(z+1,d-1)+power(z-1,d-1));
        }
        f+=2*z*cycles(m,a,b,z,false)+4*z*z*cycles(m,a,b,z,true);
        return f;
    }
};
int main(int argc,char**argv){
    Dp d;
    if(argc==4){
        int m=std::stoi(argv[1]);I a=std::stoull(argv[2]),b=std::stoull(argv[3]);
        std::cout<<d.value(m,a,b,1)<<","<<d.value(m,a,b,2)<<"\n";
        return 0;
    }
    int maxm=argc==2?std::stoi(argv[1]):10;
    std::cout<<"n,words,distinct_pairs,distinct_F,min_F,max_F,contiguous_odds_from_min,longest_odd_run\n";
    for(int m=2;m<=maxm;m++){
        std::unordered_set<Pair,Hash> pairs;
        for(I x=0;x<(I(1)<<(m-2));x++) for(I b=1;b<(I(1)<<m);b++){
            // A beats first and loses to last; A -> B and B beats at least
            // one spine vertex, hence whole source is strongly connected.
            I a=1|(x<<1);
            pairs.insert({d.value(m,a,b,1),d.value(m,a,b,2)});
        }
        std::unordered_set<I> vals;
        I low=~I(0),high=0;
        for(Pair p:pairs){I f=p.h+p.q;vals.insert(f);low=std::min(low,f);high=std::max(high,f);}
        I initial=0,maxlen=0;
        while(vals.count(low+2*initial))++initial;
        for(I f:vals) if(!vals.count(f-2)){
            I l=1;while(vals.count(f+2*l))++l;maxlen=std::max(maxlen,l);
        }
        std::cout<<m+2<<","<<((I(1)<<(m-2))*((I(1)<<m)-1))<<","<<pairs.size()<<","<<vals.size()<<","<<low<<","<<high<<","<<initial<<","<<maxlen<<std::endl;
    }
}

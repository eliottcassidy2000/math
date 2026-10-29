// Exact least original source of k growing cheaper-ternary returns.
// k<=60 ensures every intermediate below fits unsigned __int128.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using U=unsigned __int128;
struct Block {int j,d,c,r;};
const Block blocks[5]={{1,1,1,3},{4,6,179,125},{7,11,16177,1357},
                       {7,11,12289,2269},{7,11,9697,2877}};
std::string decimal(U n) {
    if (!n) return "0";
    std::string s;
    while(n){s.push_back('0'+n%10);n/=10;}
    std::reverse(s.begin(),s.end());return s;
}
struct Child {int i,J,D;uint64_t r;U h;};
struct Search {
    int k; uint64_t best,visited=0; U endpoint;
    std::array<int,60> path{},answer{};
    std::array<int,421> inv{};
    std::array<U,68> p3{};
    Search(int depth):k(depth),best((uint64_t(1)<<(k+1))-1) {
        if(k<1 || k>60)throw std::runtime_error("horizon outside proved width");
        inv[0]=1;for(int j=1;j<421;j++)inv[j]=inv[j-1]*683%2048;
        p3[0]=1;for(int j=1;j<68;j++)p3[j]=3*p3[j-1];
        endpoint=2*p3[k]-1;
    }
    void visit(int depth,int J,int D,uint64_t r,U h) {
        visited++;
        std::vector<Child> next;
        for(int i=0;i<5;i++) {
            auto b=blocks[i]; int mod=1<<b.d;
            int delta=(b.r-int(h%(2*mod)))/2;
            int t=(delta*inv[J])%mod;if(t<0)t+=mod;
            uint64_t rr=r;U lifted=h;
            if(t) {
                // A nonzero lift is impossible once the source modulus
                // reaches the authenticated incumbent. No huge 2^D needed.
                if(D+1>=64 || uint64_t(t)>((best-1-r)>>(D+1)))continue;
                rr += uint64_t(t)<<(D+1);
                if(J>61)throw std::runtime_error("nonzero lift width guard");
                lifted += 2*p3[J]*t;
            }
            if(rr>=best)continue;
            U numerator=p3[b.j]*lifted+b.c;
            if(numerator%mod)throw std::runtime_error("inexact endpoint");
            U hh=numerator/mod;
            if(!(hh&1) || hh<=lifted)throw std::runtime_error("non-growing guarded word");
            next.push_back({i,J+b.j,D+b.d,rr,hh});
        }
        std::sort(next.begin(),next.end(),[](auto a,auto b){return a.r<b.r;});
        for(auto z:next) {
            if(z.r>=best)continue;
            path[depth]=z.i;
            if(depth+1==k) {best=z.r;answer=path;endpoint=z.h;}
            else visit(depth+1,z.J,z.D,z.r,z.h);
        }
    }
};
int main(int argc,char**argv) {
    int limit=argc>1?std::stoi(argv[1]):50;
    if(limit<1||limit>60)throw std::runtime_error("require 1<=horizon<=60");
    std::cout<<"{\"scope\":\"FINITE-EXACT all positive odd sources at each finite horizon\",\"rows\":[\n";
    for(int k=1;k<=limit;k++) {
        Search s(k);s.visit(0,0,0,1,1);
        if(k>1)std::cout<<",\n";
        std::cout<<"{\"k\":"<<k<<",\"least_source\":"<<s.best
          <<",\"visited\":"<<s.visited<<",\"endpoint\":\""<<decimal(s.endpoint)<<"\",\"word\":[";
        for(int j=0;j<k;j++){if(j)std::cout<<",";std::cout<<s.answer[j];}
        std::cout<<"]}"<<std::flush;
    }
    std::cout<<"\n]}\n";
}

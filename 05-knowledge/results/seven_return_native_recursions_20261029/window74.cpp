// Structured hostile search: source bound <=2^66, horizon<=74.
// Retains every native return branch and original source; unsigned 128-bit arithmetic.
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
struct Child {int i,J,D;U r,h;};
struct Search {
    int k; U best,visited=0,endpoint;
    std::array<int,74> path{},answer{};
    std::array<int,520> inv{};
    std::array<U,75> p3{};
    Search(int depth,U bound):k(depth),best(bound),endpoint(0) {
        if(k<1 || k>74 || bound<=1 || bound>(U(1)<<66))throw std::runtime_error("outside proved width");
        inv[0]=1;for(int j=1;j<520;j++)inv[j]=inv[j-1]*683%2048;
        p3[0]=1;for(int j=1;j<75;j++)p3[j]=3*p3[j-1];

    }
    void visit(int depth,int J,int D,U r,U h) {
        visited++;
        std::vector<Child> next;
        for(int i=0;i<5;i++) {
            auto b=blocks[i]; int mod=1<<b.d;
            int delta=(b.r-int(h%(2*mod)))/2;
            int t=(delta*inv[J])%mod;if(t<0)t+=mod;
            U rr=r;U lifted=h;
            if(t) {
                // A nonzero lift is impossible once the source modulus
                // reaches the authenticated incumbent. No huge 2^D needed.
                if(D+1>=128 || U(t)>((best-1-r)>>(D+1)))continue;
                rr += U(t)<<(D+1);
                if(J>65)throw std::runtime_error("nonzero lift width guard");
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
    if(argc!=3)throw std::runtime_error("supply horizon and exclusive bound");
    int k=std::stoi(argv[1]);U bound=0;for(char c:std::string(argv[2])){if(c<'0'||c>'9')throw std::runtime_error("decimal bound required");if(bound>(U(1)<<66))throw std::runtime_error("bound too wide");bound=10*bound+c-'0';}
    Search s(k,bound);s.visit(0,0,0,1,1);
    std::cout<<"{\"k\":"<<k<<",\"exclusive_bound\":"<<decimal(bound)
      <<",\"visited\":"<<decimal(s.visited)<<",\"counterexample\":"<<(s.best<bound?"true":"false");
    if(s.best<bound){
      std::cout<<",\"source\":"<<decimal(s.best)<<",\"endpoint\":\""<<decimal(s.endpoint)<<"\",\"word\":[";
      for(int j=0;j<k;j++){if(j)std::cout<<",";std::cout<<s.answer[j];}
      std::cout<<"]";
    }
    std::cout<<"}\n";
}

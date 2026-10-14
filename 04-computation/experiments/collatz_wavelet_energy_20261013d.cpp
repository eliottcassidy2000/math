// Source-preserving ternary-wavelet energy census from authenticated parity bits.
// Usage: prog parity26.bits max_M_exponent max_q output.json
// Checks every M in [1,2^max_M_exponent), every r mod 3^q, 1<=q<=max_q.
// t(n)=(-1)^(strict Terras first-hit length); packed format as 20261013c census.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using I = int64_t;
using U = uint64_t;
using Z = __int128_t;
using W = __uint128_t;
void require(bool ok, const char* msg) { if (!ok) throw std::runtime_error(msg); }
std::string decimal(Z z) {
    bool negative=z<0; W u=negative?W(-z):W(z); std::string s;
    do {s+=char('0'+u%10);u/=10;} while(u);
    if(negative)s+='-'; std::reverse(s.begin(),s.end());return s;
}
struct Cursor {
    U n=0; std::vector<std::vector<I>> counts;
    Cursor(const std::vector<int>& mod) {for(auto m:mod)counts.emplace_back(m,0);}
    template<class T> void advance(U target,const std::vector<int>& mod,T t) {
        while(n<target){++n;I value=t(n);for(size_t q=0;q<mod.size();++q)counts[q][n%mod[q]]+=value;}
    }
    I h(int q,int r) const {return 3*counts[q][r]-counts[q-1][r%counts[q-1].size()];}
};
struct Channel {int q,r,e,o;};
struct Record { Z xx=0, yy=0, xy=0, even=0, ez=0, zz=0, energy=0; };
int main(int argc,char**argv) {try {
    require(argc==5,"usage: bits max_M_exponent max_q output");
    int E=std::stoi(argv[2]),Q=std::stoi(argv[3]);
    require(1<=E&&E<=24&&1<=Q&&Q<=3,"declared finite scope");
    const U saved=1ULL<<26,M=1ULL<<E;
    std::vector<U> bits((saved+64)/64);
    std::ifstream input(argv[1],std::ios::binary|std::ios::ate);
    require(bool(input)&&input.tellg()==std::streamoff(bits.size()*8),"expected exact parity26 file size");
    input.seekg(0);input.read(reinterpret_cast<char*>(bits.data()),bits.size()*8);
    require(bool(input),"read bits");
    auto t=[&](U n)->I {require(n>0&&n<=saved,"positive source in bit universe");return ((bits[n/64]>>(n%64))&1)?-1:1;};
    U strict_controls=0;
    auto control=[&](U n){W value=n;U steps=0;while(value!=1){require(value<(W(1)<<120),"strict overflow");value=(value&1)?(3*value+1)/2:value/2;require(++steps<1000000,"strict step cap");}require(t(n)==(steps%2?-1:1),"strict parity mismatch");++strict_controls;};
    for(U n=1;n<=4096;++n)control(n);
    U random=20261013;
    for(int j=0;j<512;++j){random^=random<<13;random^=random>>7;random^=random<<17;control(1+random%(3*M));}
    std::vector<int> mod(Q+2,1);for(int q=1;q<=Q+1;++q)mod[q]=3*mod[q-1];
    std::vector<Channel> channels;
    for(int q=1;q<=Q;++q)for(int r=0;r<mod[q];++r)
        channels.push_back({q,r,((mod[q]+1)/2*r)%mod[q],((mod[q+1]+1)/2*(3*r+1))%mod[q+1]});
    std::vector<std::vector<Record>> rows(E,std::vector<Record>(channels.size()));
    Cursor one(mod),two(mod),three(mod);U identities=0;
    for(U m=1;m<M;++m){
        one.advance(m,mod,t);two.advance(2*m,mod,t);three.advance(3*m,mod,t);
        require(t(2*m)==-t(m)&&t(2*m-1)==-t(3*m-1),"first-step laws");
        int octave=63-__builtin_clzll(m);I next=t(2*m+1);
        for(size_t j=0;j<channels.size();++j){auto c=channels[j];auto& out=rows[octave][j];
            I x=one.h(c.q,c.e),y=three.h(c.q+1,c.o),even=two.h(c.q,c.r);
            require(even==-x-y,"exact source refinement");++identities;
            I z=(3*(int((2*m+1)%mod[c.q])==c.r)-(int((2*m+1)%mod[c.q-1])==c.r%mod[c.q-1]))*next;
            out.xx+=Z(x)*x;out.yy+=Z(y)*y;out.xy+=Z(x)*y;out.even+=Z(even)*even;
            out.ez+=Z(even)*z;out.zz+=Z(z)*z;out.energy+=Z(even)*even+Z(even+z)*(even+z);
        }
        two.advance(2*m+1,mod,t);
        for(auto c:channels){I expected=-one.h(c.q,c.e)-three.h(c.q+1,c.o);
            I z=(3*(int((2*m+1)%mod[c.q])==c.r)-(int((2*m+1)%mod[c.q-1])==c.r%mod[c.q-1]))*next;
            require(two.h(c.q,c.r)==expected+z,"odd prefix carry");++identities;
        }
        if((m+1)%(1ULL<<22)==0)std::cerr<<"checked M through "<<m<<'\n';
    }
    std::ofstream output(argv[4]);
    output<<"{\"status\":\"FINITE-EXACT; no future energy or independence claim\",\"max_M_exponent\":"<<E<<",\"max_q\":"<<Q<<",\"strict_controls\":"<<strict_controls<<",\"refinement_identities\":"<<identities<<",\"rows\":[\n";
    bool comma=false;
    for(int k=0;k<E;++k)for(size_t j=0;j<channels.size();++j){auto c=channels[j];auto r=rows[k][j];
        require(r.even==r.xx+r.yy+2*r.xy,"Gram energy identity");
        require(r.energy==2*r.even+2*r.ez+r.zz,"energy including odd carry");
        require(r.xx<(Z(1)<<62)&&r.yy<(Z(1)<<62)&&r.xy>-(Z(1)<<62)&&r.xy<(Z(1)<<62),"Gram product overflow guard");
        require(r.xx*r.yy>=r.xy*r.xy,"PSD covariance");
        if(comma)output<<",\n";comma=true;
        output<<"{\"k\":"<<k<<",\"X\":"<<(1ULL<<k)<<",\"q\":"<<c.q<<",\"r\":"<<c.r<<",\"even_address\":"<<c.e<<",\"odd_address\":"<<c.o;
        output<<",\"xx\":\""<<decimal(r.xx)<<"\",\"yy\":\""<<decimal(r.yy)<<"\",\"xy\":\""<<decimal(r.xy)<<"\",\"even\":\""<<decimal(r.even)<<"\",\"carry_cross\":\""<<decimal(r.ez)<<"\",\"carry_square\":\""<<decimal(r.zz)<<"\",\"child_energy\":\""<<decimal(r.energy)<<"\"}";
    }
    output<<"]}\n";require(bool(output),"output");
    std::cerr<<"complete: "<<identities<<" exact refinement identities\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}

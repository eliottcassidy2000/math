// Search every high word on the smaller physical flank and every opposite word.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using U = unsigned __int128;
using Word = std::uint64_t;
constexpr int width = 40;
constexpr Word base = Word(1)<<width;
constexpr U target = (U(1)<<60)*((U(1)<<61)-1)-1;
struct Data { Word h,q,c,t,e,r,k; };
constexpr Data data[] = {{2465,2293760,1344270,712944,1740,986,423385},
                          {37265,144834560,84839910,57602544,24360,14906,33838015}};
constexpr Word anchor[2][2] = {{735435556549430,1575929833923500},
                                {4814468855604222,240731609083994}};
std::string decimal(U x) {
    if (x < 10) return std::to_string(unsigned(x));
    return decimal(x/10) + char('0'+x%10);
}
U reader(U z) {
    U f = 0;
    while(z) {
        bool zero = !(z&1);
        z >>= 1;
        if(zero) f += z+(z&1);
    }
    return f;
}
Word reverse_low(Word z) {
    Word result = 0;
    for (int j=0;j<width;++j) { result=2*result+(z&1); z>>=1; }
    return result;
}
Word inverse(Word x, Word h) {
    std::int64_t a=1,b=0;
    Word y=h;
    while(y) {
        Word q=x/y, r=x%y;
        auto c=a-std::int64_t(q)*b;
        a=b;b=c;x=y;y=r;
    }
    assert(x==1);
    return (a%std::int64_t(h)+h)%h;
}
struct Progression { Word residue=0, period=0; };
struct Search {
    Data d;
    U N;
    Word limit;
    std::vector<std::array<Progression,2>> classes;
    std::vector<Word> high_f;
    Word suffixes=0,rows=0,admissible=0,reads=0,hits=0;
    Search(int spine,int side,U goal):d(data[spine]),N(goal) {
        if (side) { std::swap(d.c,d.t); std::swap(d.e,d.r); }
        limit=0;
        for(int i=60;i>=0;--i) {
            Word k=limit+(Word(1)<<i);
            if(U(k)*k <= (N-d.k)/d.q) limit=k;
        }
        classes.resize(2*d.h);
        for(Word a=0;a<2*d.h;++a) for(Word parity=0;parity<2;++parity) {
            Word coef=(d.q*a+d.c)%d.h;
            Word fixed=(d.k+d.t*a+d.r*(a&1)+d.e*parity+coef*parity)%d.h;
            Word rhs=(Word(N%d.h)+d.h-fixed)%d.h;
            coef=2*coef%d.h;
            Word g=std::gcd(coef,d.h);
            if(rhs%g) continue;
            Word reduced=d.h/g;
            Word q=reduced==1?0:(rhs/g)*inverse(coef/g,reduced)%reduced;
            classes[a][parity]={2*q+parity,2*reduced};
        }
        for (Word x=0;x<=limit/base;++x) high_f.push_back(Word(reader(x)));
    }
    // A floating seed only selects a nearby integer; multiplication verifies both
    // quotient endpoints. Outlying unbalanced words use native integer division.
    static U quotient(U value,U divisor) {
        if (divisor < (U(1)<<67)) return value/divisor;
        Word seed=Word(double(value)/double(divisor));
        U multiple=divisor*seed;
        while(multiple>value) { --seed; multiple-=divisor; }
        while(value-multiple>=divisor) { ++seed; multiple+=divisor; }
        return seed;
    }
    void consider(Word low) {
        ++suffixes;
        if(low>limit) return;
        Word f=Word(reader(low)), z=base-1-reverse_low(low);
        Word a=low%(2*d.h),step=base%(2*d.h);
        for (Word high=0,x=low;x<=limit;++high,x+=base,a=(a+step>=2*d.h ? a+step-2*d.h : a+step)) {
            ++rows;
            auto& arr=classes[a];
            if(!(arr[0].period||arr[1].period)) continue;
            ++admissible;
            Word fx=f+z*high+high_f[high]+(low<base/2 && (high&1));
            U R=N-d.k-U(d.t)*x-d.r*(x&1)-U(d.h)*fx, A=U(d.q)*x+d.c;
            U hi=quotient(R,A);
            if(hi<x) continue;
            for(const auto &cl:arr) {
                if(!cl.period) continue;
                Word mod=hi<=UINT64_MAX ? Word(hi)%cl.period : Word(hi%cl.period);
                Word delta=mod>=cl.residue ? mod-cl.residue : mod+cl.period-cl.residue;
                if(delta>hi) continue;
                U y=hi-delta;
                while(y>=x && (A+d.h)*y+d.e>=R) {
                    ++reads;
                    if(A*y+d.h*reader(y)+d.e*(y&1)==R) {
                        ++hits;
                        std::cout<<"hit smaller "<<x<<" opposite "<<decimal(y)<<'\n';
                    }
                    if(y<cl.period) break;
                    y-=cl.period;
                }
            }
        }
    }
    void walk(Word low,int remaining,int first) {
        consider(low);
        if(!remaining) return;
        for(int bit=first;bit<width;++bit) walk(low^(Word(1)<<bit),remaining-1,bit+1);
    }
};
int main(int argc,char**argv) {
    if(argc<5||argc>6) return 2;
    int spine=std::stoi(argv[1]), w=std::stoi(argv[2]), side=std::stoi(argv[3]),r=std::stoi(argv[4]);
    if((spine!=4&&spine!=8)||w<0||w>1||side<0||side>1||r<0||r>8) return 2;
    U n=target;
    if(argc==6) {
        if(std::string(argv[5])!="planted" || spine!=4) return 2;
        n+=w?184814588:-844710;
    }
    Search s(spine==8,side,n);
    s.walk(anchor[w][side]%(base),r,0);
    std::cout<<spine<<' '<<w<<' '<<side<<' '<<r<<" suffixes "<<s.suffixes<<" rows "<<s.rows
             <<" congruent "<<s.admissible<<" reads "<<s.reads<<" hits "<<s.hits<<'\n';
}

// Invert s*y + h*f(y) + K[parity] on a fixed face of a two-flank chart.
// The slope buckets use exact rational inequalities; each output is checked integrally.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
#include <gmpxx.h>

using Word = unsigned long;
static_assert(sizeof(Word) == 8);
using Wide = unsigned __int128;

Wide parse(const std::string& s) {
    Wide n = 0;
    if (s.empty()) throw std::invalid_argument("empty decimal");
    for (char c : s) {
        if (c < '0' || c > '9' || n > (~Wide(0) - (c-'0')) / 10)
            throw std::invalid_argument("invalid decimal");
        n = 10*n + c-'0';
    }
    return n;
}
std::string decimal(Wide n) {
    if (!n) return "0";
    std::string s;
    while (n) { s.push_back('0'+n%10); n/=10; }
    std::reverse(s.begin(),s.end());
    return s;
}
mpz_class big(Wide n) {
    std::array<Word,2> limbs{Word(n),Word(n>>64)};
    mpz_class result;
    mpz_import(result.get_mpz_t(),2,-1,sizeof(Word),0,0,limbs.data());
    return result;
}
Word reverse(Word n,int width) {
    if (!width) return 0;
#if defined(__has_builtin)
#if __has_builtin(__builtin_bitreverse64)
    return __builtin_bitreverse64(n)>>(64-width);
#endif
#endif
    Word out=0;
    for (int i=0;i<width;++i) { out=2*out+(n&1); n>>=1; }
    return out;
}
std::array<Word,65536> small_f;
Wide reading(Wide x) {
    Wide f=0;
    while(x>=65536) {
        Word r=Word(x)&65535;
        x>>=16;
        f+=small_f[r]+x*(65535-reverse(r,16))+((r<32768)&&(x&1));
    }
    return f+small_f[Word(x)];
}
Word inverse(Word a, Word m) {
    if (m==1) return 0;
    mpz_class inv;
    mpz_class aa(a), mm(m);
    if(!mpz_invert(inv.get_mpz_t(),aa.get_mpz_t(),mm.get_mpz_t()))
        throw std::logic_error("noninvertible residue");
    return inv.get_ui();
}
struct Face {
    Wide slope,rem;
    Word h,mod,residue;
    int parity,k,b=0;
    Word W,B,step,inverseB,inverseW;
    std::vector<Word> m_reverse,counts,starts,indices;
    Word blocks=0, candidates=0, hits=0;
    mpz_class numerator;
    Word slack_up,slack_down;

    Face(Wide s, Word h_, Wide M, int j, int requested)
        : slope(s),rem(M),h(h_),parity(j),k(requested) {
        Word g=std::gcd(Word(s%h),h);
        if (M%g) { mod=0; return; }
        mod=h/g;
        residue=Word(((M/g)%mod)*inverse(Word((s/g)%mod),mod)%mod);
        // Select B so the unknown f(z) contributes at most about one bucket.
        while(b<k || big(slope)*(mpz_class(1)<<b)*big(slope)*(mpz_class(1)<<b)*mod
                         < big(rem)*h*(mpz_class(1)<<k)) ++b;
        if(b>62) throw std::range_error("low word exceeds 62 bits");
        W=Word(1)<<k; B=Word(1)<<b; step=B/W;
        if(b==k) { ++b; B*=2; step*=2; } // group never crosses the parity boundary
        if(big(slope+h)*B*W > (mpz_class(1)<<126))
            throw std::range_error("face exceeds supported intermediate width");
        inverseB=inverse(B%mod,mod);
        inverseW=inverse(mod%W,W);
        // Bound (W/mod)*h*(f(z)+boundary+product error)/D in both directions.
        mpz_class Dmin=big(slope)*B;
        auto ceil_ratio=[](mpz_class a,const mpz_class& den) -> Word {
            a=(a+den-1)/den; return a.get_ui();
        };
        slack_up=ceil_ratio(mpz_class(h)*b*W,mpz_class(mod)*Dmin);
        slack_down=ceil_ratio(mpz_class(h)*(big(rem)/Dmin+b+1)*W,mpz_class(mod)*Dmin);
        if(slack_up+slack_down+10>=W) throw std::range_error("choose larger bucket");
        numerator=big(rem)*h*W*W;
        m_reverse.resize(W); counts.resize(W); starts.resize(W+1); indices.resize(W);
        for(Word i=0;i<W;++i) m_reverse[i]=(inverseW*reverse(i,k))&(W-1);
    }
    Wide den(Word v) const { return slope*B+Wide(h)*(B-1-v); }
    Word slope_index(Word v) const {
        mpz_class D=big(den(v));
        mpz_class index=numerator/(mpz_class(mod)*D*D);
        return index.get_ui();
    }
    void bucket(Word slope_num) {
        std::fill(counts.begin(), counts.end(),0);
        auto key=[&](Word i){return (m_reverse[i]-(slope_num*i>>k))&(W-1);};
        for(Word i=0;i<W;++i) ++counts[key(i)];
        starts[0]=0;
        std::partial_sum(counts.begin(),counts.end(),starts.begin()+1);
        std::copy(starts.begin(),starts.end()-1,counts.begin());
        for(Word i=0;i<W;++i) indices[counts[key(i)]++]=i;
    }
    void test(Word R, Word i, Wide D0, Wide z0, Word phase,
              Word slope_num, long unwrapped_bucket) {
        Word j=reverse(i,k), r=R+step*j, E=m_reverse[i];
        long raw=long(E)-long((slope_num*i)>>k);
        long translate=(unwrapped_bucket-raw)/long(W);
        __int128 zsigned=__int128(z0)-phase
            +(__int128(mod)*E-j)/W+__int128(translate)*mod;
        if(zsigned<0) return;
        Wide z=Wide(zsigned), D=D0-Wide(h)*i;
        Wide main=D*z+slope*r;
        if(main>rem || rem-main>Wide(h)*(z+r+1)) return;
        Wide known=main+Wide(h)*reading(r);
        if(known>rem || rem-known>Wide(h)*z) return;
        if(rem-known==Wide(h)*(reading(z)+((r<B/2)&&(z&1)))) {
            Wide value=z*B+r;
            if(slope*value+Wide(h)*reading(value)!=rem) throw std::logic_error("bad split");
            ++hits; std::cout<<"address "<<decimal(value)<<"\n";
        }
    }
    void scan(Word first,Word last) { // inclusive block indices
        Word num=slope_index(first*W);
        if(slope_index((last+1)*W-1)>num+1 && first<last) {
            Word mid=(first+last)/2; scan(first,mid); scan(mid+1,last); return;
        }
        Word rise=slope_index((last+1)*W-1)-num+1;
        bucket(num);
        for(Word block=first;block<=last;++block) {
            ++blocks;
            Word v=block*W, R=reverse(block,b-k);
            Wide D=den(v), z=rem/D, rest=rem-z*D;
            Word off=Word(Wide((residue+mod-R%mod)%mod)*inverseB%mod);
            Word phase=(Word(z%mod)+mod-off)%mod;
            // floor(W*rest/D) is at most W, and rest*W must fit.
            if(rest>(~Wide(0)>>k)) throw std::range_error("fractional phase overflows");
            Word base=Word((Wide(W)*phase+((rest<<k)/D))/mod);
            long lo=long(base)-long(slack_down)-1;
            Word count=slack_down+slack_up+rise+4;
            for(Word j=0;j<std::min(W,count);++j) {
                Word key=(lo+j)&(W-1);
                for(Word a=starts[key];a<starts[key+1];++a) {
                    ++candidates;
                    test(R,indices[a],D,z,phase,num,lo+long(j));
                }
            }
        }
    }
    void run(Word piece,Word pieces) {
        if(!mod) return;
        Word half=step/2;
        Word first=parity*half+Word(Wide(half)*piece/pieces);
        Word end=parity*half+Word(Wide(half)*(piece+1)/pieces);
        if(first<end) scan(first,end-1);
        std::cout<<"parity "<<parity<<" low bits "<<b<<" table bits "<<k
                 <<" blocks "<<blocks<<" candidates "<<candidates<<" hits "<<hits<<std::endl;
    }
};
int main(int argc,char**argv) {
    if(argc!=7 && argc!=9) {
        std::cerr<<"usage: farey_axis slope h even-intercept odd-intercept target table-bits [shard total]\n";
        return 2;
    }
    for(Word i=1;i<65536;++i)
        small_f[i]=small_f[i/2]+((i&1)?0:i/2+((i/2)&1));
    Wide s=parse(argv[1]), hw=parse(argv[2]), N=parse(argv[5]);
    int k=std::stoi(argv[6]);
    Word piece=argc==9?std::stoull(argv[7]):0, pieces=argc==9?std::stoull(argv[8]):1;
    if(s>(Wide(1)<<124) || N>(Wide(1)<<123) || !s||!hw||!(hw&1)||hw>(Word(1)<<48)||k<4||k>22||piece>=pieces) return 2;
    for(int j=0;j<2;++j) {
        Wide K=parse(argv[3+j]);
        if(N>=K) { Face f(s,Word(hw),N-K,j,k); f.run(piece,pieces); }
    }
}

// Direct contact-word certificate, plus optional exhaustive optimality interval.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
using Word=std::uint64_t;
constexpr Word target=Word(1)<<60, u0=102580081, w0=5619616816, best=2758518367;
constexpr Word block=1<<12, mask=block-1;
std::array<Word,block> low{}, coefficient{};
void tables() {
    for (Word x=0; x<block; ++x) {
        if (x) low[x]=low[x/2]+((x&1)?0:x/2+(x/2&1));
        for (int j=0;j<12;++j)
            coefficient[x]=2*coefficient[x]+!(x&(Word(1)<<j));
    }
}
Word append(Word prefix, Word reading, Word suffix) {
    return reading+prefix*coefficient[suffix]+low[suffix]+((prefix&1)&&!(suffix&(block/2)));
}
Word reader(Word x) {
    // All addresses in the proof interval are below 2^36.
    Word first=x>>24, two=x>>12;
    return append(two,append(first,low[first],two&mask),x&mask);
}
std::string binary(Word w) {
    std::string bits;
    while(w) { bits.push_back(char('0'+(w&1))); w>>=1; }
    std::reverse(bits.begin(),bits.end());
    return bits;
}
Word direct(const std::string& bits) {
    Word answer=1;
    for(unsigned i=0;i<bits.size();++i) if(bits[i]=='1')
        for(unsigned j=i+1;j<bits.size();++j) if(bits[j]=='0')
            answer+=(Word(1)<<(j-i-1))+(j==i+1);
    return answer;
}
void certificate() {
    std::string suffix=binary(w0);
    std::reverse(suffix.begin(),suffix.end());
    for(char& b:suffix) b=b=='0'?'1':'0';
    auto left=direct(binary(u0)+"0"+suffix),right=direct(binary(u0)+"1"+suffix);
    auto mersenne=direct("1"+std::string(29,'0')+"101"+std::string(29,'0'));
    assert(left==target-best && right==target+best && mersenne==(Word(1)<<61)-1);
    assert(reader(u0)==45284624 && reader(w0)==5351979710);
    std::cout<<"dyadic counts "<<left<<" "<<right<<", Mersenne factor "<<mersenne<<"\n";
}
int main(int argc,char** argv) {
    tables();
    certificate();
    if (argc==1) return 0;
    if (argc!=3) return 1;
    Word lower=std::stoull(argv[1]), upper=std::stoull(argv[2]);
    if (lower<50000000 || upper>759250125 || lower>upper) return 2;
    for(Word u=lower;u<=upper;++u) {
        Word cu=(u+1)/2, remainder=target-reader(u)-1-cu;
        Word first=(2*remainder-1+4*u+2)/(4*u+3);
        Word last=2*remainder/(4*u+1);
        for(Word w=std::max(first,u);w<=last;++w) {
            Word cw=(w+1)/2;
            if (((cw^cu)&1) && 2*u*w+reader(w)+cw==remainder)
                std::cout<<"hit "<<u<<" "<<w<<" halfgap "<<cw-cu<<"\n";
        }
    }
}

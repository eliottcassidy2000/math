// Complete meet-in-the-middle inverse for the two top-order receiver charts.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using Word = std::uint64_t;
using Wide = unsigned __int128;
constexpr Word BLOCK = 1u<<16;
std::array<unsigned,BLOCK> rev16{},f16{};

std::string decimal(Wide n) {
    if(!n) return "0";
    std::string out;
    while(n) { out+=char('0'+n%10); n/=10; }
    std::reverse(out.begin(),out.end());
    return out;
}

Word reverse_width(Word x,int b) {
    Word rev=(Word(rev16[x&65535])<<48) | (Word(rev16[(x>>16)&65535])<<32)
            | (Word(rev16[(x>>32)&65535])<<16) | rev16[(x>>48)&65535];
    return rev>>(64-b);
}

Word read_word(Word x) {
    Word result=0;
    while(x>=BLOCK) {
        unsigned low=x&(BLOCK-1);
        x>>=16;
        result+=f16[low]+x*(BLOCK-1-rev16[low])+((low<BLOCK/2) && (x&1));
    }
    return result+f16[x];
}

Wide read_wide(Wide x) {
    Wide result=0;
    while(x) {
        bool zero=!(x&1);
        x>>=1;
        if(zero) result+=x+(x&1);
    }
    return result;
}

struct Chart { int n,c,h,k0,k1,ones; };
const Chart SINGLE{1,81,25,43,68,0}, FOUR{4,1021,125,915,940,5};

template<typename T>
std::vector<Wide> inverse(const Chart& chart,T target,int width,bool report,Word shard=0,Word shards=1,Word radius=0) {
    // Strip the terminal ones in the four-vertex chart before choosing the split.
    Word scale=Word(1)<<chart.ones,c=chart.c*scale,h=chart.h;
    std::array<Word,2> K={Word(chart.k0),Word(chart.k1)};
    if(chart.ones) K[0]=K[1]=chart.k1+chart.c*(scale-1);
    int w=width-chart.ones;
    int b=1;
    Wide desired_square=Wide(target)*(radius?h:1);
    while(b<w-1) {
        Wide cB=Wide(c)<<b;
        if(cB>=desired_square/cB+(desired_square%cB!=0)) break;
        ++b;
    }
    Word B=Word(1)<<b,zfirst=Word(1)<<(w-b-1),zlast=2*zfirst-1;
    std::vector<Word> residues(h);
    for(Word z=0;z<h;++z) residues[(B%h)*(c%h)*z%h]=z;
    std::array<Word,2> residue_constant={Word((Wide(target)-K[0])%h),Word((Wide(target)-K[1])%h)};
    std::vector<Wide> hits;
    Word trials=0;
    for(Word r=shard;r<B;r+=shards) {
        T base=K[r&1]+T(c)*r+h*T(read_word(r));
        if(base>target+radius) continue;
        T rem=target+radius-base;
        T min_rem=rem>=2*radius ? rem-2*radius : 0;
        T den=T(c)*B+T(h)*(B-1-reverse_width(r,b));
        Word lo=std::max(zfirst,Word((min_rem+den+h-1)/(den+h)));
        Word hi=std::min(zlast,Word(rem/den));
        if(lo>hi) continue;
        Word rhs=(residue_constant[r&1]+h-(c%h)*(r%h)%h)%h;
        if(!radius) lo+=(residues[rhs]+h-lo%h)%h;
        for(Word z=lo;z<=hi;z+=radius?1:h) {
            ++trials;
            T cost=T(z)*den+h*T(read_word(z)+((r<B/2)&&(z&1)));
            if(cost<min_rem || cost>rem) continue;
            Wide y=(Wide(z)*B+r+1)*scale-1;
            Wide expected=Wide(chart.c)*y+chart.h*read_wide(y)+((y&1)?chart.k1:chart.k0);
            if(expected!=base+cost) throw std::runtime_error("literal binary replay disagrees");
            hits.push_back(y);
        }
    }
    if(report) {
        std::cout<<"width "<<width<<", receiver "<<chart.n<<", low bits "<<b<<", residues tested "
                 <<(shard<B ? 1+(B-1-shard)/shards : 0)<<", quotient trials "<<trials<<", hits "<<hits.size();
        if(!radius) for(auto y:hits) std::cout<<" "<<decimal(y);
        std::cout<<std::endl;
        if(radius) {
            Word down=radius+1,up=radius+1;
            Wide below=0,above=0;
            for(auto y:hits) {
                Wide v=Wide(chart.c)*y+chart.h*read_wide(y)+((y&1)?chart.k1:chart.k0);
                if(v<=target && target-v<down) {down=Word(target-v);below=y;}
                if(v>=target && v-target<up) {up=Word(v-target);above=y;}
            }
            std::cout<<"  below "<<(down<=radius ? std::to_string(down)+" at "+decimal(below) : "none")
                     <<", above "<<(up<=radius ? std::to_string(up)+" at "+decimal(above) : "none")<<std::endl;
        }
    }
    return hits;
}

void controls() {
    for(auto chart:{SINGLE,FOUR}) for(int w:{10,13,16}) {
        int first=1<<(w-1);
        for(int j=0;j<8;++j) {
            int scale=1<<chart.ones;
            int planted=(first+(133*j)%first)|(scale-1);
            Word target=Word(chart.c)*planted+chart.h*read_word(planted)
                       +((planted&1)?chart.k1:chart.k0)+2*(j&1);
            Word radius=j%3 ? 0 : 32;
            auto got=inverse(chart,target,w,false,0,1,radius);
            std::sort(got.begin(),got.end());
            std::vector<Wide> wanted;
            for(int y=first;y<2*first;++y) {
                Word v=Word(chart.c)*y+chart.h*read_wide(y)+((y&1)?chart.k1:chart.k0);
                if((y+1)%scale==0 && v+radius>=target && v<=target+radius) wanted.push_back(y);
            }
            if(wanted!=got) throw std::runtime_error("exhaustive small control failed");
        }
    }
}

int main(int argc,char** argv) {
    if(argc<2) {
        std::cerr<<"usage: top_order_inverse P [--receiver 1|4] [--shard I N] [--radius 0..512]\n";
        return 2;
    }
    int p=std::stoi(argv[1]),which=0;
    Word shard=0,shards=1,radius=0;
    for(int i=2;i<argc;++i) {
        std::string opt=argv[i];
        if(opt=="--receiver"&&i+1<argc) { which=std::stoi(argv[++i]); if(which!=1&&which!=4) return 2; }
        else if(opt=="--shard"&&i+2<argc) { shard=std::stoull(argv[++i]);shards=std::stoull(argv[++i]); }
        else if(opt=="--radius"&&i+1<argc) radius=std::stoull(argv[++i]);
        else return 2;
    }
    if(p<11 || p>61 || !(p&1) || !shards || shard>=shards || shards>(Word(1)<<56) || radius>512) return 2;
    for(Word i=1;i<BLOCK;++i) {
        rev16[i]=(rev16[i>>1]>>1)|((i&1)<<15);
        f16[i]=f16[i>>1]+((i&1)?0:(i>>1)+((i>>1)&1));
    }
    controls();
    Wide target=(Wide(1)<<(p-1))*((Wide(1)<<p)-1)-1;
    std::cout<<"exponent "<<p<<" shard "<<shard<<" / "<<shards<<std::endl;
    for(auto chart:{SINGLE,FOUR}) {
        if(which && chart.n!=which) continue;
        int w=2*p-6-chart.n;
        if(p<=31) inverse(chart,Word(target),w,true,shard,shards,radius);
        else inverse(chart,target,w,true,shard,shards,radius);
    }
}

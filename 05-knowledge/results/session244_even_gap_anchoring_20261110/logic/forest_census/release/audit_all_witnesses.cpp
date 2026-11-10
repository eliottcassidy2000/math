#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

constexpr unsigned MAXN=12;
static std::array<std::array<uint32_t,MAXN>,1u<<MAXN> dp{};

// Literal endpoint DP: no graft identity or census parent-state reuse.
uint32_t count(unsigned n,uint32_t code,unsigned &depth){
    std::array<unsigned,MAXN> choice{},incoming{},dep{};
    for(unsigned v=n;v>1;--v){choice[v-1]=code%v;code/=v;}
    if(code)throw std::runtime_error("invalid code");
    depth=0;
    for(unsigned v=1;v<n;++v){
        dep[v]=choice[v]?dep[choice[v]-1]+1:0;depth=std::max(depth,dep[v]);
        for(unsigned u=0;u<v;++u){
            if(choice[v]==u+1)incoming[u]|=1u<<v;
            else incoming[v]|=1u<<u;
        }
    }
    const unsigned full=(1u<<n)-1u;
    for(unsigned s=1;s<=full;++s){
        unsigned bits=s;
        while(bits){
            const unsigned v=static_cast<unsigned>(__builtin_ctz(bits));bits&=bits-1;
            const unsigned rest=s^(1u<<v);uint32_t val=0;
            if(rest==0)val=1;
            else {
                unsigned pred=rest&incoming[v];
                while(pred){const unsigned u=static_cast<unsigned>(__builtin_ctz(pred));pred&=pred-1;val+=dp[rest][u];}
            }
            dp[s][v]=val;
        }
    }
    uint32_t ans=0;for(unsigned v=0;v<n;++v)ans+=dp[full][v];return ans;
}

unsigned parse_order(const std::string &s){
    if(s.empty()||s.find_first_not_of("0123456789")!=std::string::npos)
        throw std::runtime_error("order must be an integer in 1..12");
    const auto n=std::stoul(s);
    if(n<1||n>MAXN)throw std::runtime_error("order must be 1..12");
    return static_cast<unsigned>(n);
}

int run(int argc,char **argv){
    if(argc!=3)throw std::runtime_error("usage: audit_all_witnesses max_order directory");
    const unsigned maxn=parse_order(argv[1]);
    uint64_t checked=0,factorial=1;
    for(unsigned n=1;n<=maxn;++n){
        factorial*=n;
        const std::string file=std::string(argv[2])+"/n"+std::to_string(n)+"_witnesses.csv";
        std::ifstream f(file);
        if(!f)throw std::runtime_error("missing csv: "+file);
        std::string line;std::getline(f,line);
        if(line!="H,chronological_code,min_forest_depth,min_depth_code")throw std::runtime_error("csv header");
        uint64_t rows=0,last=0;
        while(std::getline(f,line)){
            if(std::count(line.begin(),line.end(),',')!=3)throw std::runtime_error("csv column count");
            if(line.find_first_not_of("0123456789,")!=std::string::npos)throw std::runtime_error("csv character");
            std::replace(line.begin(),line.end(),',',' ');std::istringstream s(line);
            uint64_t h=0,code=0,d=0,mincode=0;
            if(!(s>>h>>code>>d>>mincode))throw std::runtime_error("csv format");
            std::string tail;if(s>>tail)throw std::runtime_error("trailing csv data");
            if(h<=last||h>factorial||!(h&1u)||h==7||h==21)throw std::runtime_error("csv count range/order/exclusion");
            if(code>=factorial||mincode>=factorial||d>=n)throw std::runtime_error("csv witness range");
            unsigned depth=0;
            if(count(n,static_cast<uint32_t>(code),depth)!=h)throw std::runtime_error("full witness mismatch");
            if(d>depth)throw std::runtime_error("minimum depth exceeds first witness depth");
            if(count(n,static_cast<uint32_t>(mincode),depth)!=h||depth!=d)throw std::runtime_error("minimum-depth witness mismatch");
            ++rows;checked+=2;last=h;
        }
        if(!f.eof()||rows==0)throw std::runtime_error("csv read failure or empty spectrum");
        std::cout<<"{\"n\":"<<n<<",\"rows_verified\":"<<rows<<",\"both_codes_recounted\":true}\n";
    }
    std::cerr<<"All "<<checked<<" literal graph witnesses PASS\n";
    return 0;
}
int main(int argc,char **argv){
    try{return run(argc,argv);}
    catch(const std::exception &e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}
}

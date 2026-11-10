#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

constexpr unsigned LIM=12, MASKS=1u<<LIM;
using Row=std::array<uint32_t,LIM>;
static std::array<Row,MASKS> starts{},ends{};
static std::array<uint32_t,MASKS> H{};
static std::array<unsigned,LIM> out{},in{};
static std::array<std::vector<uint32_t>,LIM+1> witness;
static std::array<uint64_t,LIM+1> processed{};
static std::array<unsigned,LIM> depth{};
static std::array<std::vector<unsigned char>,LIM+1> min_depth;
static std::array<std::vector<uint32_t>,LIM+1> min_depth_code;
static std::array<uint32_t,LIM+1> maximum{};
static unsigned target;
static uint64_t leafparents=0;
static auto began=std::chrono::steady_clock::now();

void record(unsigned n,uint32_t h,uint32_t code,unsigned dep){
    if (!(h&1u)) throw std::runtime_error("even full count");
    auto &v=witness[n];
    if(h>=v.size()){
        v.resize(static_cast<size_t>(h)+1024u,0);
        min_depth[n].resize(v.size(),255);
        min_depth_code[n].resize(v.size(),0);
    }
    if(dep<min_depth[n][h]){min_depth[n][h]=static_cast<unsigned char>(dep);min_depth_code[n][h]=code;}
    if(v[h]==0)v[h]=code+1u;
    maximum[n]=std::max(maximum[n],h);++processed[n];
}

void fill_new(unsigned n){
    const unsigned top=1u<<(n-1),full=(1u<<n)-1u;
    for(unsigned s=top;s<=full;++s){
        uint32_t total=0;
        unsigned vs=s;
        while(vs){
            const unsigned r=static_cast<unsigned>(__builtin_ctz(vs));vs&=vs-1;
            const unsigned rem=s^(1u<<r);
            uint32_t a=0,b=0;
            if(rem==0){a=1;b=1;}
            else {
                unsigned choices=rem&out[r];
                while(choices){const unsigned u=static_cast<unsigned>(__builtin_ctz(choices));choices&=choices-1;a+=starts[rem][u];}
                choices=rem&in[r];
                while(choices){const unsigned u=static_cast<unsigned>(__builtin_ctz(choices));choices&=choices-1;b+=ends[rem][u];}
            }
            starts[s][r]=a;ends[s][r]=b;total+=a;
        }
        H[s]=total;
    }
}

void visit(unsigned n,uint32_t code,unsigned dep){
    const unsigned full=(1u<<n)-1u;
    record(n,H[full],code,dep);
    if(n+1u==target){
        std::array<uint32_t,LIM> U{};
        for(unsigned s=1;s<=full;++s){
            const uint32_t complement=H[full^s];
            unsigned vs=s;
            while(vs){const unsigned r=static_cast<unsigned>(__builtin_ctz(vs));vs&=vs-1;U[r]+=starts[s][r]*complement;}
        }
        record(target,H[full],code*target,dep);
        for(unsigned r=0;r<n;++r)record(target,H[full]+U[r]-ends[full][r],code*target+r+1u,std::max(dep,depth[r]+1u));
        ++leafparents;
        if(leafparents%100000u==0u){
            const double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
            std::cerr<<"parents "<<leafparents<<" elapsed "<<sec<<" s\n";
        }
        return;
    }
    const auto savedout=out,savedin=in;
    const unsigned bit=1u<<n;
    for(unsigned choice=0;choice<=n;++choice){
        out=savedout;in=savedin;out[n]=0;in[n]=0;
        for(unsigned r=0;r<n;++r){
            if(choice==r+1u){out[n]|=1u<<r;in[r]|=bit;}
            else {out[r]|=bit;in[n]|=1u<<r;}
        }
        depth[n]=choice==0?0:depth[choice-1u]+1u;
        fill_new(n+1u);visit(n+1u,code*(n+1u)+choice,std::max(dep,depth[n]));
    }
    out=savedout;in=savedin;
}

int main(int argc,char **argv){
    if(argc!=3)throw std::runtime_error("usage: census max_order output_directory");
    const std::string argument=argv[1];
    if(argument.empty()||argument.find_first_not_of("0123456789")!=std::string::npos)
        throw std::runtime_error("order must be an integer in 2..12");
    const auto parsed=std::stoul(argument);
    if(parsed<2||parsed>LIM)throw std::runtime_error("order must be 2..12");
    target=static_cast<unsigned>(parsed);
    H[0]=1;H[1]=1;starts[1][0]=ends[1][0]=1;
    visit(1,0,0);
    uint64_t factorial=1;
    for(unsigned n=1;n<=target;++n){
        factorial*=n;
        if(processed[n]!=factorial)throw std::runtime_error("factorial history count mismatch");
        std::ofstream f(std::string(argv[2])+"/n"+std::to_string(n)+"_witnesses.csv");
        if(!f)throw std::runtime_error("cannot write CSV");
        f<<"H,chronological_code,min_forest_depth,min_depth_code\n";
        uint64_t values=0,depth2=0,depth3=0;uint32_t first=0,first_deep=0;
        for(uint32_t h=1;h<=maximum[n];h+=2){
            if(witness[n][h]){
                ++values;
                depth2+=min_depth[n][h]<=2u;depth3+=min_depth[n][h]<=3u;
                if(first_deep==0&&min_depth[n][h]>2u)first_deep=h;
                f<<h<<','<<witness[n][h]-1u<<','<<static_cast<unsigned>(min_depth[n][h])<<','<<min_depth_code[n][h]<<'\n';
            }
            else if(first==0&&h!=7&&h!=21)first=h;
        }
        if(first==0){first=maximum[n]+2u;while(first==7u||first==21u)first+=2u;}
        std::cout<<"{\"n\":"<<n<<",\"histories\":"<<processed[n]<<",\"distinct_values\":"<<values<<",\"max_H\":"<<maximum[n]<<",\"first_allowed_miss\":"<<first<<",\"depth2_values\":"<<depth2<<",\"depth3_values\":"<<depth3<<",\"first_needing_depth_gt2\":"<<first_deep<<"}\n";
    }
    std::cerr<<"complete elapsed "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<" s\n";
}

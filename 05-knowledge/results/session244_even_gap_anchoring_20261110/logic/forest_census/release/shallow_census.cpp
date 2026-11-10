#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Independent shallow-only traversal. It never visits a depth-three birth,
// stores no minimum-depth array, and checks the complete image against CSV.
constexpr unsigned LIMIT=12;
using Row=std::array<uint32_t,LIMIT>;
static std::array<Row,1u<<LIMIT> ending{},beginning{};
static std::array<uint32_t,1u<<LIMIT> total{};
static std::array<unsigned,LIMIT> outgoing{},incoming{},level{};
static std::vector<unsigned char> seen;
static unsigned order=0;
static bool normalize=false;
static uint64_t histories=0,parents=0;
static uint32_t largest=0;
static auto began=std::chrono::steady_clock::now();

void store(uint32_t h){
    if(!(h&1u))throw std::runtime_error("even count");
    if(h>=seen.size())seen.resize(static_cast<size_t>(h)+1024u,0);
    seen[h]=1;largest=std::max(largest,h);++histories;
}
void update(unsigned k){
    const unsigned first=1u<<(k-1),last=(1u<<k)-1u;
    for(unsigned mask=first;mask<=last;++mask){
        uint32_t sum=0;
        for(unsigned v=0;v<k;++v)if(mask&(1u<<v)){
            const unsigned rest=mask^(1u<<v);uint32_t e=0,b=0;
            if(!rest)e=b=1;
            else {
                unsigned neighbours=incoming[v]&rest;
                while(neighbours){const unsigned u=static_cast<unsigned>(__builtin_ctz(neighbours));neighbours&=neighbours-1;e+=ending[rest][u];}
                neighbours=outgoing[v]&rest;
                while(neighbours){const unsigned u=static_cast<unsigned>(__builtin_ctz(neighbours));neighbours&=neighbours-1;b+=beginning[rest][u];}
            }
            ending[mask][v]=e;beginning[mask][v]=b;sum+=e;
        }
        total[mask]=sum;
    }
}
void search(unsigned k){
    const unsigned full=(1u<<k)-1u;
    if(k+1==order){
        store(total[full]); // sink
        for(unsigned root=0;root<k;++root){
            if(level[root]>=2||(normalize&&root+1==k))continue;
            uint32_t cover=0;
            // This order of summation differs from the all-history census.
            const unsigned rootbit=1u<<root;
            const unsigned other=full^rootbit;
            unsigned subset=other;
            for(;;){
                cover+=beginning[subset|rootbit][root]*total[other^subset];
                if(!subset)break;
                subset=(subset-1u)&other;
            }
            store(total[full]+cover-ending[full][root]);
        }
        ++parents;
        if(parents%100000u==0u)std::cerr<<"parents "<<parents<<" elapsed "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<" s\n";
        return;
    }
    const auto oldout=outgoing,oldin=incoming;
    for(unsigned choice=0;choice<=k;++choice){
        if(choice&&(level[choice-1]>=2||(normalize&&choice==k)))continue;
        outgoing=oldout;incoming=oldin;outgoing[k]=incoming[k]=0;
        level[k]=choice?level[choice-1]+1:0;
        for(unsigned v=0;v<k;++v){
            if(choice==v+1){outgoing[k]|=1u<<v;incoming[v]|=1u<<k;}
            else{outgoing[v]|=1u<<k;incoming[k]|=1u<<v;}
        }
        update(k+1);search(k+1);
    }
    outgoing=oldout;incoming=oldin;
}
uint64_t expected(){
    using Table=std::array<std::array<std::array<uint64_t,3>,LIMIT+1>,LIMIT+1>;
    Table table{};table[1][0][0]=1;
    for(unsigned n=1;n<order;++n){
        Table next{};
        for(unsigned r=0;r<=n;++r)for(unsigned s=0;s<=n-r;++s)for(unsigned last=0;last<3;++last){
            const uint64_t count=table[r][s][last];if(!count)continue;
            next[r+1][s][0]+=count;
            const unsigned roots=r-(normalize&&last==0?1u:0u);
            const unsigned ones=s-(normalize&&last==1?1u:0u);
            next[r][s+1][1]+=count*roots;
            next[r][s][2]+=count*ones;
        }
        table=next;
    }
    uint64_t sum=0;for(const auto&r:table)for(const auto&s:r)for(auto count:s)sum+=count;return sum;
}
int run(int argc,char**argv){
    if(argc!=4)throw std::runtime_error("usage: shallow_census order CSV_directory normalize_0_or_1");
    const std::string arg=argv[1];
    if(arg.empty()||arg.find_first_not_of("0123456789")!=std::string::npos)throw std::runtime_error("bad order");
    const auto parsed=std::stoul(arg);if(parsed<2||parsed>LIMIT)throw std::runtime_error("order must be 2..12");
    order=static_cast<unsigned>(parsed);
    if(std::string(argv[3])!="0"&&std::string(argv[3])!="1")throw std::runtime_error("bad normalization mode");
    normalize=std::string(argv[3])=="1";
    total[0]=total[1]=ending[1][0]=beginning[1][0]=1;
    search(1);
    const auto expected_histories=expected();
    if(histories!=expected_histories)throw std::runtime_error("restricted history count mismatch");
    uint64_t values=0;for(auto b:seen)values+=b;
    std::ifstream input(std::string(argv[2])+"/n"+std::to_string(order)+"_witnesses.csv");
    if(!input)throw std::runtime_error("missing CSV");
    std::string line;std::getline(input,line);
    if(line!="H,chronological_code,min_forest_depth,min_depth_code")throw std::runtime_error("CSV header");
    uint64_t shallow_rows=0;std::vector<uint32_t> excluded;
    while(std::getline(input,line)){
        std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);
        uint32_t h=0,code=0,d=0,mincode=0;if(!(row>>h>>code>>d>>mincode))throw std::runtime_error("CSV row");
        const bool found=h<seen.size()&&seen[h]!=0;
        if(found!=(d<=2))throw std::runtime_error("independent shallow image disagrees with full census");
        if(found)++shallow_rows;else excluded.push_back(h);
    }
    if(shallow_rows!=values)throw std::runtime_error("extra shallow values not in full census");
    std::cout<<"{\"status\":\"PASS\",\"n\":"<<order<<",\"no_immediate_parent_normalization\":"<<(normalize?"true":"false")<<",\"histories\":"<<histories<<",\"expected_histories\":"<<expected_histories<<",\"leaf_parents\":"<<parents<<",\"distinct_shallow_values\":"<<values<<",\"max_H\":"<<largest<<",\"complete_image_equals_full_census_depth_at_most_two\":true,\"full_census_values_absent_shallow\":[";
    bool first=true;for(auto h:excluded){if(!first)std::cout<<',';first=false;std::cout<<h;}
    std::cout<<"]}\n";
    std::cerr<<"complete elapsed "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()<<" s\n";
    return 0;
}
int main(int argc,char**argv){try{return run(argc,argv);}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}}

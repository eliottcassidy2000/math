#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

constexpr unsigned N=9, LIM=10;
using Graph=std::array<unsigned,LIM>;
using Signature=std::vector<uint32_t>;
static Graph graph{},depth{};
static std::array<std::array<uint32_t,LIM>,1u<<LIM> dp{};
static uint64_t shallow_histories=0,fiber_histories=0;
static std::map<Signature,uint32_t> signatures;
static std::set<uint32_t> child_union;

// Recount a literal graph from scratch. No U formula or parent DP reuse.
uint32_t hp(unsigned n,const Graph &a){
    Graph incoming{};
    for(unsigned v=0;v<n;++v)for(unsigned u=0;u<n;++u)
        if(a[u]&(1u<<v))incoming[v]|=1u<<u;
    const unsigned full=(1u<<n)-1u;
    for(unsigned s=1;s<=full;++s){
        unsigned bits=s;
        while(bits){
            const unsigned v=static_cast<unsigned>(__builtin_ctz(bits));bits&=bits-1;
            const unsigned rest=s^(1u<<v);uint32_t value=0;
            if(!rest)value=1;
            else {
                unsigned pred=rest&incoming[v];
                while(pred){const unsigned u=static_cast<unsigned>(__builtin_ctz(pred));pred&=pred-1;value+=dp[rest][u];}
            }
            dp[s][v]=value;
        }
    }
    uint32_t answer=0;for(unsigned v=0;v<n;++v)answer+=dp[full][v];return answer;
}
Graph add(const Graph&a,unsigned n,unsigned choice){
    Graph b=a;b[n]=0;
    for(unsigned u=0;u<n;++u){if(choice==u+1)b[n]|=1u<<u;else b[u]|=1u<<n;}
    return b;
}
Signature children(const Graph&a){
    Signature result;
    for(unsigned choice=0;choice<=N;++choice)result.push_back(hp(N+1,add(a,N,choice)));
    return result;
}
Graph decode(uint32_t code,unsigned &maximum_depth){
    std::array<unsigned,N> choices{};
    for(unsigned n=N;n>1;--n){choices[n-1]=code%n;code/=n;}
    if(code)throw std::runtime_error("code range");
    Graph a{},depths{};maximum_depth=0;
    for(unsigned n=1;n<N;++n){
        a=add(a,n,choices[n]);
        depths[n]=choices[n]?depths[choices[n]-1]+1:0;
        maximum_depth=std::max(maximum_depth,depths[n]);
    }
    return a;
}
void visit(unsigned n,uint32_t code){
    if(n==N){
        ++shallow_histories;
        if(hp(n,graph)!=217)return;
        ++fiber_histories;
        auto sig=children(graph);
        child_union.insert(sig.begin(),sig.end());
        std::sort(sig.begin(),sig.end());sig.erase(std::unique(sig.begin(),sig.end()),sig.end());
        signatures.emplace(sig,code);return;
    }
    const Graph saved=graph;
    for(unsigned choice=0;choice<=n;++choice){
        const unsigned d=choice?depth[choice-1]+1:0;
        if(d>2)continue;
        graph=add(saved,n,choice);depth[n]=d;visit(n+1,code*(n+1)+choice);
    }
    graph=saved;
}
uint64_t independent_shallow_history_count(){
    using Table=std::array<std::array<uint64_t,N+1>,N+1>;
    Table table{};table[1][0]=1;
    for(unsigned n=1;n<N;++n){
        Table next{};
        for(unsigned r=0;r<=n;++r)for(unsigned s=0;s<=n-r;++s){
            next[r+1][s]+=table[r][s];
            next[r][s+1]+=table[r][s]*r;
            next[r][s]+=table[r][s]*s;
        }
        table=next;
    }
    uint64_t result=0;for(const auto &row:table)for(auto x:row)result+=x;return result;
}
void print(const Signature &a){std::cout<<'[';bool first=true;for(auto x:a){if(!first)std::cout<<',';first=false;std::cout<<x;}std::cout<<']';}
int run(){
    unsigned example_depth=0,shallow_depth=0;
    const Graph example=decode(26485,example_depth),shallow=decode(723,shallow_depth);
    if(example_depth!=3||shallow_depth!=1)throw std::runtime_error("example history depths");
    if(hp(N,example)!=217||hp(N,shallow)!=217)throw std::runtime_error("parent count");
    const Signature expected{217,2041,2241,653,769,1197,539,437,501,373};
    const auto example_children=children(example);
    if(example_children!=expected)throw std::runtime_error("literal children");
    visit(1,0);
    if(shallow_histories!=independent_shallow_history_count())throw std::runtime_error("shallow history total");
    Signature absent;
    for(auto h:example_children)if(!child_union.count(h))absent.push_back(h);
    std::sort(absent.begin(),absent.end());
    if(absent!=Signature{373,437,501,539,2041,2241})throw std::runtime_error("absence certificate");
    std::array<unsigned,N> order{};for(unsigned v=0;v<N;++v)order[v]=v;
    std::map<unsigned,unsigned> ordering_depth_histogram;
    do {
        Graph d{};bool valid=true;unsigned maximum=0;
        for(unsigned i=0;i<N&&valid;++i){
            unsigned count=0,parent=0;
            for(unsigned j=0;j<i;++j)if(example[order[i]]&(1u<<order[j])){++count;parent=j;}
            if(count>1){valid=false;break;}
            d[i]=count?d[parent]+1:0;maximum=std::max(maximum,d[i]);
        }
        if(valid)++ordering_depth_histogram[maximum];
    }while(std::next_permutation(order.begin(),order.end()));
    if(ordering_depth_histogram.empty()||ordering_depth_histogram.begin()->first!=3)throw std::runtime_error("graph minimum depth");
    std::cout<<"{\n  \"status\":\"PASS\",\n  \"n\":9,\n  \"parent_code\":26485,\n  \"H\":217,\n  \"depth_one_same_H_code\":723,\n  \"children_in_choice_order\":";
    print(example_children);
    std::cout<<",\n  \"depth_at_most_two_histories\":"<<shallow_histories<<",\n  \"same_H_shallow_histories\":"<<fiber_histories<<",\n  \"same_H_shallow_child_signatures\":"<<signatures.size()<<",\n  \"same_H_shallow_child_union_size\":"<<child_union.size()<<",\n  \"children_absent_from_all_same_H_shallow_parents\":";
    print(absent);
    std::cout<<",\n  \"all_9_factorial_orderings_valid_depth_histogram\":{";
    bool first=true;for(const auto &item:ordering_depth_histogram){if(!first)std::cout<<',';first=false;std::cout<<'"'<<item.first<<"\":"<<item.second;}
    std::cout<<"},\n  \"all_same_H_shallow_signatures\":[\n";
    first=true;for(const auto &item:signatures){if(!first)std::cout<<",\n";first=false;std::cout<<"    {\"code\":"<<item.second<<",\"children\":";print(item.first);std::cout<<'}';}
    std::cout<<"\n  ]\n}\n";return 0;
}
int main(){try{return run();}catch(const std::exception&e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}}

// Enumerate simultaneous one-arc parents and both unit neighbors at every prefix width.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#define main original_fork_census_driver
#include "physical_forks.cpp"
#undef main
#pragma GCC diagnostic pop
#include <cstdint>

enum class Kind { impossible, fixed, arbitrary };
struct Condition {
    Kind kind;
    int value=0;
};
Condition solve(int multiplier, int remainder) {
    if (!multiplier) return {remainder ? Kind::impossible : Kind::arbitrary};
    if (remainder%multiplier) return {Kind::impossible};
    return {Kind::fixed,remainder/multiplier};
}
Condition intersect(Condition a,Condition b) {
    if (a.kind==Kind::impossible || b.kind==Kind::impossible) return {Kind::impossible};
    if (a.kind==Kind::arbitrary) return b;
    if (b.kind==Kind::arbitrary || a.value==b.value) return a;
    return {Kind::impossible};
}
struct Choice {Condition condition; Edge edge;};
struct Options {std::vector<Choice> below, above, parent;};
struct Count {
    long long centers=0, parent_arcs=0, apex_parents=0;
    Rows rows;
    int mask=0, word=0,mode=0;
    Edge lo,hi,parent;
};

const std::vector<int>& prefixes(int mode,int d) {
    static std::map<Key,std::vector<int>> cache;
    auto [loc,fresh]=cache.emplace(Key{mode,d},std::vector<int>{});
    if (fresh && (mode ? d<=0 : d>=0) && d%2==0) {
        for (int v=0;v<=std::abs(d);v++) {
            int w=2*v+1;
            int correction=mode ? -read(w-1) : w-1-read(w);
            if (correction==d) loc->second.push_back(w);
        }
    }
    return loc->second;
}

int main(int argc,char** argv) {
    if (argc!=3) return 2;
    int largest=std::stoi(argv[2]);
    if (largest<1 || largest>8) return 2;
    using LongKey = std::pair<std::int64_t,std::int64_t>;
    std::map<LongKey,Count> sums,small;
    int unrestricted=0;
    for(int n=1;n<=largest;n++) {
        std::ifstream file(std::string(argv[1])+"/classes"+std::to_string(n)+".rows");
        if(!file) return 1;
        std::string line;
        long long tagged=0;
        while(std::getline(file,line)) {
            Rows rows(n); std::istringstream input(line);
            for(auto& e:rows) if(!(input>>e)) return 1;
            int h; auto center=readings(rows,h);
            int limit=1<<n;
            std::vector<std::array<Options,2>> options(limit);
            auto offer=[&](int m,Both other,int h2,Edge edge) {
                for(int t=0;t<2;t++) {
                    auto [a,k]=center[m][t];
                    if(a<=0) continue;
                    auto [q,b]=other[t];
                    auto& opt=options[m][t];
                    auto add=[&](std::vector<Choice>& out,Condition cond) {
                        if(cond.kind!=Kind::impossible) out.push_back({cond,edge});
                    };
                    if(q==a) {
                        add(opt.below,solve(h2-h,b-k+2));
                        add(opt.above,solve(h2-h,b-k-2));
                    }
                    if(3*q==a) add(opt.parent,solve(3*h2-h,3*b-k));
                }
            };
            for(int i=0;i<n;i++) {
                for(int m=0;m<limit;m++) offer(m,center[m^(1<<i)],h,{i,n});
                for(int j=i+1;j<n;j++) {
                    Rows changed=rows;changed[i]^=1<<j;changed[j]^=1<<i;
                    int h2; auto other=readings(changed,h2);
                    for(int m=0;m<limit;m++) offer(m,other[m],h2,{i,j});
                }
            }
            tagged+=limit;
            for(int m=0;m<limit;m++) for(int t=0;t<2;t++) {
                const auto& opts=options[m][t];
                std::set<int> seen;
                std::set<std::pair<int,Edge>> seen_edge;
                for(auto parent:opts.parent) for(auto lo:opts.below) for(auto hi:opts.above) {
                    auto constraint=intersect(parent.condition,intersect(lo.condition,hi.condition));
                    if(constraint.kind==Kind::impossible) continue;
                    if(constraint.kind==Kind::arbitrary) {
                        unrestricted++;
                        std::cout<<"all-prefix n "<<n<<" mask "<<m<<" mode "<<t<<" base "<<h<<std::endl;
                        continue;
                    }
                    auto [a,k]=center[m][t];
                    for(int w:prefixes(t,constraint.value)) {
                        Key key{a*w,k-h*constraint.value};
                        bool first=seen.insert(w).second;
                        bool new_arc=seen_edge.insert({w,parent.edge}).second;
                        if(new_arc) {
                            std::cout<<"presentation center "<<key.first<<" "<<key.second
                                     <<" word "<<w<<" mode "<<t<<" rows";
                            for(int v:rows) std::cout<<" "<<v;
                            std::cout<<" mask "<<m<<" arcs "<<lo.edge.first<<","<<lo.edge.second
                                     <<" "<<hi.edge.first<<","<<hi.edge.second<<" "
                                     <<parent.edge.first<<","<<parent.edge.second<<"\n";
                        }
                        auto record=[&](Count& c) {
                            if(first) c.centers++;
                            if(new_arc) {
                                c.parent_arcs++; c.apex_parents+=parent.edge.second==n;
                                if(c.rows.empty()) {
                                    c.rows=rows;c.mask=m;c.word=w;c.mode=t;
                                    c.lo=lo.edge;c.hi=hi.edge;c.parent=parent.edge;
                                }
                            }
                        };
                        record(sums[key]); if(w<=31) record(small[key]);
                    }
                }
            }
        }
        std::cout<<"order "<<n<<" masks "<<tagged<<" finite curves "<<sums.size()
                 <<" within31 "<<small.size()<<" all-prefix triples "<<unrestricted<<std::endl;
    }
    for(const auto& [key,c]:sums) {
        std::cout<<"center "<<key.first<<" "<<key.second<<" within31 "<<small.count(key)
                 <<" presentations "<<c.centers<<" parent_arcs "<<c.parent_arcs
                 <<" apex_parents "<<c.apex_parents<<" word "<<c.word<<" mode "<<c.mode<<" rows";
        for(int v:c.rows) std::cout<<" "<<v;
        std::cout<<" mask "<<c.mask<<" arcs "<<c.lo.first<<","<<c.lo.second<<" "
                 <<c.hi.first<<","<<c.hi.second<<" "<<c.parent.first<<","<<c.parent.second<<"\n";
    }
    return 0;
}

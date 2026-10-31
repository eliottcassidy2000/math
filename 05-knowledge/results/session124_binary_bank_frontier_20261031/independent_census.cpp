// Independent Session 124 audit. No producer implementation is included.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using I = std::int64_t;
using U = std::uint64_t;
using Graph = std::vector<unsigned>;
using Partition = std::vector<std::vector<int>>;
static void require(bool b, const std::string &s) { if (!b) throw std::runtime_error(s); }
static I factorial(int n) { I a=1; for(int k=2;k<=n;++k)a*=k; return a; }
static Graph decode(U code,int n) {
    require(n>=1 && n<=11,"order range");
    Graph g(n); int k=0;
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++k)
        if((code>>k)&1)g[i]|=1u<<j;else g[j]|=1u<<i;
    require(k==64 || (code>>k)==0,"code range"); return g;
}
static Graph extend(const Graph &g,unsigned mask) {
    int n=g.size(); require(n<11 && mask<(1u<<n),"bad cut");
    Graph q=g; q.push_back(mask);
    for(int v=0;v<n;++v)if(!(mask&(1u<<v)))q[v]|=1u<<n;
    return q;
}
static U encode(const Graph &g,const Partition &p) {
    U code=0;int k=0,n=g.size();
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++k)
        if(g[p[i][0]]&(1u<<p[j][0]))code|=U(1)<<k;
    return code;
}
// Equitable refinement by ordered outgoing-count vectors, then individualization.
// This is not the producer's enumeration within outdegree classes.
static Partition refine(const Graph &g,Partition p) {
    for(;;) {
        Partition q;
        for(const auto &cell:p) {
            std::map<std::vector<int>,std::vector<int>> split;
            for(int v:cell) {
                std::vector<int> key;
                for(const auto &target:p) { int count=0;for(int w:target)count+=(g[v]>>w)&1;key.push_back(count); }
                split[key].push_back(v);
            }
            for(auto &kv:split)q.push_back(std::move(kv.second));
        }
        if(q.size()==p.size())return q;
        p=std::move(q);
    }
}
struct Canon { U code=std::numeric_limits<U>::max(); I aut=0; };
static void canon_rec(const Graph &g,Partition p,Canon &best) {
    p=refine(g,std::move(p)); int cell=-1;
    for(int i=0;i<(int)p.size();++i)if(p[i].size()>1){cell=i;break;}
    if(cell<0) {
        U value=encode(g,p);
        if(value<best.code)best={value,1};else if(value==best.code)++best.aut;
        return;
    }
    for(int v:p[cell]) {
        Partition q;
        for(int i=0;i<(int)p.size();++i) {
            if(i!=cell)q.push_back(p[i]);
            else { q.push_back({v});std::vector<int> rest;for(int w:p[i])if(v!=w)rest.push_back(w);q.push_back(std::move(rest)); }
        }
        canon_rec(g,std::move(q),best);
    }
}
static Canon canonical(const Graph &g) {
    std::vector<int> all(g.size()); std::iota(all.begin(),all.end(),0);
    Canon c;canon_rec(g,Partition{all},c);require(c.aut>0,"no labeling");return c;
}
// Removing the FIRST vertex, rather than propagating endpoint counts.
static I direct_h(const Graph &g) {
    int n=g.size(); unsigned N=1u<<n; std::vector<I>d(N*n);
    for(unsigned S=1;S<N;++S)for(int v=0;v<n;++v)if(S&(1u<<v)) {
        unsigned R=S^(1u<<v); I value=0;
        if(!R)value=1;
        else for(unsigned tail=R&g[v];tail;tail&=tail-1) {
            int w=__builtin_ctz(tail);value+=d[R*n+w];
        }
        d[S*n+v]=value;
    }
    I h=0;for(int v=0;v<n;++v)h+=d[(N-1)*n+v];
    require(h>0 && h<=factorial(n),"Hamilton bound");return h;
}
struct Response { I h;std::vector<I> f; };
// Independent implementation: backward subset recurrences and POSITIVE outgoing
// weighted-cut sums. No producer linear-coefficient/cancellation recurrence.
static Response response(const Graph &g) {
    int n=g.size(); unsigned N=1u<<n; require(n<=10,"response order limit");
    std::vector<I> start(N*n),finish(N*n); std::vector<unsigned> in(n);
    for(int u=0;u<n;++u)for(int v=0;v<n;++v)if((g[u]>>v)&1)in[v]|=1u<<u;
    for(unsigned S=1;S<N;++S)for(int v=0;v<n;++v)if(S&(1u<<v)) {
        unsigned R=S^(1u<<v);I a=0,b=0;
        if(!R)a=b=1;
        else {
            for(unsigned A=R&g[v];A;A&=A-1)a+=start[R*n+__builtin_ctz(A)];
            for(unsigned A=R&in[v];A;A&=A-1)b+=finish[R*n+__builtin_ctz(A)];
        }
        start[S*n+v]=a;finish[S*n+v]=b;
    }
    I h=0;for(int v=0;v<n;++v)h+=start[(N-1)*n+v];
    std::vector<I> W(n*n);
    for(unsigned S=1;S+1<N;++S) {
        unsigned R=(N-1)^S;
        for(unsigned A=S;A;A&=A-1) {
            int u=__builtin_ctz(A);I a=finish[S*n+u]; if(!a)continue;
            for(unsigned B=R;B;B&=B-1) {
                int v=__builtin_ctz(B);W[u*n+v]+=a*start[R*n+v];
            }
        }
    }
    I total=0;
    for(int u=0;u<n;++u)for(int v=0;v<n;++v) {
        require(W[u*n+v]>=0 && W[u*n+v]<=factorial(n-1),"W bound");total+=W[u*n+v];
    }
    require(total<=(n-1)*factorial(n),"two-path cover bound");
    std::vector<I> pair(n*n);
    for(int u=0;u<n;++u)for(int v=u+1;v<n;++v) {
        I weight=W[u*n+v]+W[v*n+u];require(weight>0 && weight%2==0,"pair weight");
        pair[u*n+v]=pair[v*n+u]=weight;
    }
    std::vector<I> f(N,h);
    // Explicit positive arc scan at each cut, avoiding the producer recurrence.
    for(unsigned S=0;S<N;++S)for(unsigned A=S;A;A&=A-1) {
        int u=__builtin_ctz(A);
        for(unsigned B=g[u]&~S;B;B&=B-1)f[S]+=pair[u*n+__builtin_ctz(B)];
    }
    for(I value:f)require(value>=h && value<=factorial(n+1),"response bound");
    require(f[0]==h && f[N-1]==h,"extreme cuts");return {h,std::move(f)};
}
static bool has(const std::vector<I> &v,I x){return std::find(v.begin(),v.end(),x)!=v.end();}
struct Stats {I parents=0,eligible=0,one=0,minus=0,plus=0,neither=0,both=0,childcuts=0,grandcuts=0,directcuts=0;};
static void direct_check(const Graph &g,const Response &r,Stats &s) {
    require(direct_h(g)==r.h,"parent direct mismatch");
    for(unsigned cut=0;cut<r.f.size();++cut) {
        require(direct_h(extend(g,cut))==r.f[cut],"direct cut mismatch");++s.directcuts;
    }
}
static Stats bank_census(const std::set<U> &reps,int n,bool quick) {
    Stats s;
    for(U code:reps) {
        Graph g=decode(code,n);Response r=response(g);++s.parents;
        // All representative parents through order 7; every 257th at orders 8,9.
        if((!quick && (n<=7 || s.parents%257==0)) || s.parents==1)direct_check(g,r,s);
        if(r.h<9)continue;
        ++s.eligible;
        bool first[2]={has(r.f,2*r.h-1),has(r.f,2*r.h+1)};
        if(!first[0] || !first[1])continue;
        ++s.one;
        bool branch[2]={false,false};
        // Deliberately examine EVERY qualifying child, not just a first witness.
        for(unsigned cut=0;cut<r.f.size();++cut)for(int digit=0;digit<2;++digit)
            if(r.f[cut]==2*r.h-1+2*digit) {
                Graph child=extend(g,cut);Response q=response(child);++s.childcuts;s.grandcuts+=q.f.size();
                require(q.h==r.f[cut],"child count mismatch");
                if(!quick)direct_check(child,q,s);
                branch[digit]=branch[digit] || (has(q.f,2*q.h-1)&&has(q.f,2*q.h+1));
            }
        if(branch[0]&&branch[1])++s.both;else if(branch[0])++s.minus;else if(branch[1])++s.plus;else ++s.neither;
    }
    return s;
}
int main(int argc,char **argv)try {
    // A single-graph mode exposes independently computed values, with every cut
    // re-counted on the actual extended tournament by direct subset DP.
    if(argc==4 && std::string(argv[1])=="--graph-code") {
        int n=std::stoi(argv[2]);U code=std::stoull(argv[3]);
        Graph g=decode(code,n);Response r=response(g);Stats s;direct_check(g,r,s);
        std::cout<<"{\"H\":"<<r.h<<",\"direct_cut_checks\":"<<s.directcuts<<",\"cuts\":[";
        for(unsigned c=0;c<r.f.size();++c){if(c)std::cout<<',';std::cout<<r.f[c];}
        std::cout<<"]}\n";return 0;
    }
    int limit=9;std::string input;bool quick=false;
    for(int i=1;i<argc;++i) {
        std::string a=argv[i];if(a=="--max-order" && i+1<argc)limit=std::stoi(argv[++i]);
        else if(a=="--producer-input" && i+1<argc)input=argv[++i];
        else if(a=="--quick")quick=true;else throw std::runtime_error("unknown argument");
    }
    require(limit>=1&&limit<=9,"max order");
    std::set<U> reps={0}; const I expected[]={0,1,1,2,4,12,56,456,6880,191536};
    std::cout<<"{\"independent_method\":\"equitable_refinement_individualization_orbit_mass_and_positive_cut_sums\",\"full_direct_checks\":"<<(quick?"false":"true")<<",\"orders\":["<<std::flush;
    for(int n=1;n<=limit;++n) {
        I mass=0;std::map<I,I> histogram;
        for(U c:reps) {
            Canon a=canonical(decode(c,n));require(a.code==c,"own canonical idempotence");
            require(factorial(n)%a.aut==0,"automorphism divisibility");mass+=factorial(n)/a.aut;++histogram[a.aut];
        }
        require(mass==(I(1)<<(n*(n-1)/2)),"orbit mass does not cover every labeled tournament");
        require((I)reps.size()==expected[n],"independent class count disagreement");
        if(!input.empty()) {
            std::ifstream f(input+"/reps"+std::to_string(n)+".txt");require(bool(f),"missing producer representatives");
            std::set<U> translated;U x;I rows=0;while(f>>x){++rows;translated.insert(canonical(decode(x,n)).code);}
            require(f.eof(),"malformed representative file");require(rows==(I)translated.size(),"producer contains duplicate isomorphism class");
            require(translated==reps,"producer representatives differ from independent universe");
        }
        Stats s=bank_census(reps,n,quick);
        require(!s.both,"COUNTEREXAMPLE to finite claim");
        if(n>1)std::cout<<',';
        std::cout<<"{\"n\":"<<n<<",\"classes\":"<<reps.size()<<",\"labeled_orbit_mass\":"<<mass<<",\"automorphism_histogram\":{";
        bool comma=false;for(auto [a,b]:histogram){if(comma)std::cout<<',';comma=true;std::cout<<'\"'<<a<<"\":"<<b;}
        std::cout<<"},\"H_ge_9\":"<<s.eligible<<",\"one_level\":"<<s.one<<",\"minus_only\":"<<s.minus<<",\"plus_only\":"<<s.plus<<",\"neither\":"<<s.neither<<",\"depth_two\":"<<s.both<<",\"qualifying_children\":"<<s.childcuts<<",\"grandchild_cuts\":"<<s.grandcuts<<",\"direct_child_dp_checks\":"<<s.directcuts<<"}"<<std::flush;
        if(n<limit) {
            std::set<U> next;
            for(U code:reps){Graph g=decode(code,n);for(unsigned cut=0;cut<(1u<<n);++cut)next.insert(canonical(extend(g,cut)).code);}
            reps.swap(next);
        }
    }
    std::cout<<"]}\n";return 0;
} catch(const std::exception &e) { std::cerr<<"AUDIT FAILURE: "<<e.what()<<'\n';return 1; }

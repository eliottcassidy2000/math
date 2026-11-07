// Independent exact endpoint-equality classifier. Standard C++17 only.
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using U = uint64_t;
void need(bool ok, const std::string& why) { if (!ok) throw std::runtime_error(why); }
struct Stats { U tournaments=0, tested_roots=0, equalities=0, degenerate=0; };
std::vector<U> endings(const std::vector<unsigned>& out) {
    int n=out.size(); unsigned full=(1u<<n)-1;
    std::vector<unsigned> in(n);
    for(int i=0;i<n;++i) for(int j=0;j<n;++j) if(out[j]>>i&1) in[i]|=1u<<j;
    std::vector<U> dp((full+1)*n);
    for(int i=0;i<n;++i) dp[(1u<<i)*n+i]=1;
    for(unsigned s=1;s<=full;++s) if(s&(s-1)) {
        unsigned lasts=s;
        while(lasts) { unsigned bit=lasts&-lasts; lasts-=bit; int i=__builtin_ctz(bit);
            unsigned prev=s^bit, pred=prev&in[i]; U sum=0;
            while(pred) { unsigned b=pred&-pred; pred-=b; sum+=dp[prev*n+__builtin_ctz(b)]; }
            dp[s*n+i]=sum;
        }
    }
    return std::vector<U>(dp.begin()+full*n, dp.end());
}
bool transitive(const std::vector<unsigned>& out, unsigned mask) {
    int n=out.size();
    for(int i=0;i<n;++i) if(mask>>i&1) for(int j=i+1;j<n;++j) if(mask>>j&1)
        for(int k=j+1;k<n;++k) if(mask>>k&1)
            if(((out[i]>>j&1)&&(out[j]>>k&1)&&(out[k]>>i&1)) || ((out[j]>>i&1)&&(out[k]>>j&1)&&(out[i]>>k&1))) return false;
    return true;
}
void inspect(const std::vector<unsigned>& out, Stats& stats) {
    int n=out.size(); need(n>=1 && n<=12,"unsupported tournament order"); unsigned full=(1u<<n)-1;
    for(int i=0;i<n;++i) {
        need(!(out[i]&~full) && !(out[i]>>i&1),"invalid row or loop");
        for(int j=i+1;j<n;++j)need(((out[i]>>j&1)+(out[j]>>i&1))==1,"not a tournament");
    }
    auto e=endings(out); U h=0; for(U z:e)h+=z;
    ++stats.tournaments;
    for(int a=0;a<n;++a) {
        if(!out[a])continue;
        if(!e[a]) { if(h==1)++stats.degenerate; continue; }
        ++stats.tested_roots;
        bool actual=(h==2*e[a]+1), predicted=false;
        if(__builtin_popcount(out[a])==1) {
            int r=__builtin_ctz(out[a]); unsigned remain=full^(1u<<a)^(1u<<r);
            predicted=(__builtin_popcount(out[r])==1 && transitive(out,remain));
            if(predicted) {
                int y=__builtin_ctz(out[r]), prefix=0;
                need((remain>>y&1)!=0,"return target outside residual chain");
                for(int x=0;x<n;++x)if((remain>>x&1)&&(out[x]>>y&1))++prefix;
                int core=n-prefix;
                need(core>=3 && h==(U(1)<<(core-2))+1 && e[a]==(U(1)<<(core-3)),"predicted count/ray failure");
            }
        }
        if(actual!=predicted) { std::ostringstream s; s<<"classification fails n="<<n<<" a="<<a<<" H="<<h<<" e="<<e[a]<<" rows";for(auto r:out)s<<' '<<r;throw std::runtime_error(s.str()); }
        if(actual)++stats.equalities;
    }
}
void emit(const std::string& name,int n,const Stats&s) {
    std::cout<<"{\"universe\":\""<<name<<"\",\"n\":"<<n<<",\"tournaments\":"<<s.tournaments<<",\"positive_nonsink_roots\":"<<s.tested_roots<<",\"equalities\":"<<s.equalities<<",\"excluded_H1_e0_equalities\":"<<s.degenerate<<"}\n";
}
int main(int argc,char**argv) {
    for(int n=1;n<=6;++n) {
        Stats s; int edges=n*(n-1)/2;
        for(unsigned c=0;c<(1u<<edges);++c) {
            std::vector<unsigned> out(n); int bit=0;
            for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++bit)if(c>>bit&1)out[i]|=1u<<j;else out[j]|=1u<<i;
            inspect(out,s);
        }
        emit("all_labelled",n,s);
    }
    for(int k=1;k<argc;++k) {
        std::ifstream input(argv[k]);need(bool(input),"cannot read class input");std::string line;Stats s;int n=0;
        while(std::getline(input,line)) {std::istringstream in(line);std::vector<unsigned> out;unsigned row;while(in>>row)out.push_back(row);if(out.empty())continue; n=out.size();inspect(out,s);}
        const U expected[]={0,1,1,2,4,12,56,456,6880};
        need(n<=8 && s.tournaments==expected[n],"unexpected inherited class-file size");
        emit("inherited_complete_unlabelled_rows",n,s);
    }
    U orphans=0;
    for(int m=1;m<=10;++m)for(unsigned mask=0;mask<(1u<<m);++mask)if(__builtin_popcount(mask)>=2) {
        int a=m,r=m+1;std::vector<unsigned> out(m+2);
        for(int i=0;i<m;++i) {for(int j=i+1;j<m;++j)out[i]|=1u<<j;out[i]|=1u<<a;if(mask>>i&1)out[r]|=1u<<i;else out[i]|=1u<<r;}
        out[a]=1u<<r;int i=__builtin_ctz(mask),j=__builtin_ctz(mask^(1u<<i));std::vector<int> p;
        for(int v=0;v<=j;++v) if(v!=i) p.push_back(v);
        p.push_back(a); p.push_back(r); p.push_back(i);
        for(int v=j+1;v<m;++v) p.push_back(v);
        for(size_t v=1;v<p.size();++v)need(out[p[v-1]]>>p[v]&1,"orphan not a path");
        need(p.back()!=a && p.back()!=r && !(out[j]>>r&1),"orphan not outside injection image");
        auto e=endings(out);U h=0;for(U x:e)h+=x;need(h>=2*e[a]+3,"orphan slack failed");++orphans;
    }
    std::cout<<"{\"orphan_witnesses\":"<<orphans<<",\"status\":\"PASS\"}\n";
}

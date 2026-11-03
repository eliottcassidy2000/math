// Compile every strong seven-vertex source chart, retaining its actual path.
// Row-vector coordinates (H,(C-E)/2,(C+E)/2) use boundary vertex zero.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>
using I=std::int64_t;
using V=std::array<I,3>;
void need(bool b,const char* s){if(!b)throw std::runtime_error(s);}

std::vector<unsigned> source(int n,unsigned code){
    std::vector<unsigned> r(n);int bit=0;
    for(int u=0;u<n;++u)for(int v=u+1;v<n;++v){
        bool yes=v==u+1 || ((code>>bit)&1);bit+=v!=u+1;
        r[yes?u:v]|=1u<<(yes?v:u);
    }
    return r;
}
bool strong(const std::vector<unsigned>& r){
    unsigned full=(1u<<r.size())-1;
    for(unsigned s=0;s<r.size();++s){
        unsigned seen=1u<<s,old=0;
        while(seen!=old){old=seen;for(unsigned u=0;u<r.size();++u)if(seen>>u&1)seen|=r[u];}
        if(seen!=full)return false;
    }
    return true;
}
V state(std::vector<unsigned> r,int tail){
    int a=r.size(),n=a+tail;unsigned full=(1u<<n)-1;
    for(auto& row:r)row|=full^((1u<<a)-1);
    for(int j=0;j<tail;++j)r.push_back(full^((1u<<(a+j+1))-1));
    for(int u=0;u<n;++u)r[u]=(full^r[u])&~(1u<<u);
    for(int u=0;u+1<n;++u){r[u]|=1u<<(u+1);r[u+1]&=~(1u<<u);}
    std::vector<unsigned> in(n);
    for(int u=0;u<n;++u)for(int v=0;v<n;++v)if(r[u]>>v&1)in[v]|=1u<<u;
    std::vector<I> dp((full+1)*n),rev(dp.size()),H(full+1);H[0]=1;
    for(int u=0;u<n;++u)dp[(1u<<u)*n+u]=rev[(1u<<u)*n+u]=1;
    for(unsigned mask=1;mask<=full;++mask){
        unsigned todo=mask;
        while(todo){int u=__builtin_ctz(todo);todo&=todo-1;unsigned prev=mask^(1u<<u);
            unsigned x=prev&in[u];while(x){int v=__builtin_ctz(x);x&=x-1;dp[mask*n+u]+=dp[prev*n+v];}
            x=prev&r[u];while(x){int v=__builtin_ctz(x);x&=x-1;rev[mask*n+u]+=rev[prev*n+v];}
            H[mask]+=dp[mask*n+u];
        }
    }
    I C=0,E=dp[full*n];for(unsigned m=1;m<=full;m+=2)C+=rev[m*n]*H[full^m];
    need((C-E)%2==0,"integral chart coordinates");
    return {H[full],(C-E)/2,(C+E)/2};
}
int main(int argc,char**argv){
    need(argc==2||argc==4,"supply output tsv, or --initial8 target-tsv output-tsv");
    int n=argc==2?7:8;std::set<V> targets;
    if(argc==4){std::ifstream input(argv[2]);V v;while(input>>v[0]>>v[1]>>v[2])targets.insert(v);}
    std::ofstream out(argv[argc-1]);need(bool(out),"output open");
    int count=0,hits=0;unsigned end=1u<<((n-1)*(n-2)/2);
    for(unsigned code=0;code<end;++code){
        auto r=source(n,code);if(!strong(r))continue;++count;
        V p0=state(r,0);if(argc==4&&!targets.count(p0))continue;
        V v1=state(r,1),v2=state(r,2);++hits;
        out<<code;for(I x:p0)out<<'\t'<<x;
        for(int j=0;j<3;++j)out<<'\t'<<v2[j]-v1[j];
        for(int j=0;j<3;++j)out<<'\t'<<v1[j]-p0[j];out<<'\n';
    }
    std::cout<<"order "<<n<<" strong_charts "<<count<<" emitted "<<hits<<'\n';
}

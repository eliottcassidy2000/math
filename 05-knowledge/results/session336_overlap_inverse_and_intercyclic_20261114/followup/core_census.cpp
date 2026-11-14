#include <array>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
using V=std::array<long long,36>;
using Pair=std::array<int,2>;
std::vector<int> move[4][36];
std::vector<Pair> lefts,rights;
std::vector<V> suffix;
long long lo,ties,tested;
std::string word,best;
V push(const V& x,int d){
    V y{};
    for(int i=0;i<36;++i) if(x[i]) for(int j:move[d][i]) y[j]+=x[i];
    return y;
}
V pre(const V& x,int d){
    V y{};
    for(int i=0;i<36;++i) for(int j:move[d][i]) y[i]+=x[j];
    return y;
}
bool allowed(int a,int b){return (a&1 && b&2)||(a&2 && b&1);}
void visit(const V& p,int rest) {
    if(rest) {
        for(int d=0;d<4;++d) {
            word.push_back('0'+d); visit(push(p,d),rest-1); word.pop_back();
        }
        return;
    }
    for(unsigned u=0;u<suffix.size();++u) {
        long long count=0;
        for(int i=0;i<36;++i)count+=p[i]*suffix[u][i];
        ++tested;
        if(count<lo) {
            lo=count; ties=0;
            best=word+char('0'+rights[u][0])+char('0'+rights[u][1]);
        }
        if(count==lo)++ties;
    }
}
int main(int argc,char**argv){
    int max=argc>1?atoi(argv[1]):8;
    assert(0<=max&&max<=11);
    V v{},e{};
    v[0]=v[18]=1;
    for(int d=0;d<4;++d)for(int o=0;o<2;++o)
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)for(int c=0;c<2;++c) {
            int u=18*o+6*a+2*b+c, l=d>>(1-o)&1, r=d>>o&1;
            move[d][u].push_back(18*o+6*(l?1:2)+2*b+c);
            if(b||l)move[d][u].push_back(18*o+6*a+2*(r?1:2)+c);
            if(c||r)move[d][u].push_back(18*o+6*a+2*b+1);
            e[u]=(a!=1&&(b==2||(b==0&&o==0)));
        }
    for(int i=0;i<4;++i)for(int j=0;j<4;++j){
        if(allowed(i,j))lefts.push_back({i,j});
        if(allowed(3-i,3-j)){
            rights.push_back({i,j});suffix.push_back(pre(pre(e,j),i));
        }
    }
    assert(lefts.size()==7&&rights.size()==7);
    for(int b=0;b<=max;++b){
        lo=std::numeric_limits<long long>::max();tested=ties=0;
        for(auto a:lefts) {
            word=std::string(1,'0'+a[0])+char('0'+a[1]); visit(push(push(v,a[0]),a[1]),b);
        }
        assert(tested==(49LL<<(2*b)));
        std::cout<<b<<" "<<lo<<" "<<best<<" "<<ties<<" "<<tested<<std::endl;
    }
}

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>
using Mask=std::uint32_t;
using Mat=std::vector<std::vector<int>>;
using Hist=std::map<int,int>;

Mat from_bits(int n,const std::string& bits) {
    Mat a(n,std::vector<int>(n));
    assert(bits.size()==std::size_t(n*(n-1)/2));
    int k=0;
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j) {
        assert(bits[k]=='1'||bits[k]=='0');
        a[i][j]=bits[k++]=='1'?1:-1; a[j][i]=-a[i][j];
    }
    return a;
}
Mat carousel(int n) {
    Mat a(n,std::vector<int>(n));
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)if(i!=j)
        a[i][j]=(j-i+n)%n<=n/2 ? 1 : -1;
    return a;
}
Hist detect(const Mat& a) {
    const int n=a.size();
    Mask neg[32]={};
    int g[32][32]={};
    Hist h;
    for(int i=0;i<n;++i){
        for(int j=0;j<n;++j) if(a[i][j]<0)neg[i]|=1u<<j;
        assert(std::popcount(neg[i])==n/2);
    }
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)
        g[i][j]=g[j][i]=2*std::popcount((neg[i]^neg[j])&~((1u<<i)|(1u<<j)))-(n-2);
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j)for(int k=j+1;k<n;++k){
        const int H=n-3-2*std::popcount((neg[i]^neg[j]^neg[k])&~((1u<<i)|(1u<<j)|(1u<<k)));
        const int J=g[i][j]*(a[k][i]+a[k][j])+g[i][k]*(a[j][i]+a[j][k])+g[j][k]*(a[i][j]+a[i][k]);
        ++h[3*H-J];
    }
    return h;
}
int paths6(const Mat& a,Mask f) {
    const int n=a.size();
    Mask out[32]={};
    for(int i=0;i<n;++i)for(int j=0;j<n;++j)
        if(i!=j && a[i][j]==((((f>>i)^(f>>j))&1u)?-1:1))out[i]|=1u<<j;
    auto rec=[&](auto&& self,int i,Mask used,int left)->int{
        Mask x=out[i]&~used; if(left==1)return std::popcount(x);
        int total=0;
        while(x){Mask y=x&(-x);x^=y;total+=self(self,std::countr_zero(y),used|y,left-1);}
        return total;
    };
    int total=0;
    for(int i=0;i<n;++i)total+=rec(rec,i,1u<<i,5);
    return total;
}
Mat bank(int s) {
    Mat a(15,std::vector<int>(15));
    int base[4]={0,3,8,15};
    for(int b=0;b<3;++b){
        auto c=carousel(base[b+1]-base[b]);
        for(int i=base[b];i<base[b+1];++i){
            for(int j=base[b];j<base[b+1];++j)a[i][j]=(s&(1<<b)?-1:1)*c[i-base[b]][j-base[b]];
            for(int j=base[b+1];j<15;++j){a[i][j]=s&8?-1:1;a[j][i]=-a[i][j];}
        }
    }
    return a;
}
void finite() {
    for(int s=0;s<16;++s){
        const auto a=bank(s);
        const int quartet[4]={0,8,9,11};
        int sum=0;
        for(unsigned f=0;f<16;++f){
            Mask flip=0;
            for(int j=0;j<4;++j)if(f&(1u<<j))flip|=1u<<quartet[j];
            sum+=(std::popcount(f)&1?-1:1)*paths6(a,flip);
        }
        assert(sum==(((s>>2)^(s>>3))&1?-312:-248));
    }
    std::cout<<"(3,5,7) bank: sixteen direct quartic tests, coefficients -15.5 and -19.5\n";
    int reg=0;
    for(unsigned f=0;f<1u<<21;++f){
        int v[7]={},k=0;
        for(int i=0;i<7;++i)for(int j=i+1;j<7;++j)++v[(f>>k++)&1?i:j];
        if(!std::all_of(v,v+7,[](int x){return x==3;}))continue;
        ++reg;
        std::string bits(21,'0');
        for(int i=0;i<21;++i)if((f>>i)&1)bits[i]='1';
        auto h=detect(from_bits(7,bits));
        assert(h[0]==21);
    }
    assert(reg==2640);
    std::cout<<"7 vertices: all 2640 labelled regular tournaments have 14 detecting triples\n";
    for(int p:{7,11,19,23,31}){
        Mat a(p,std::vector<int>(p));bool q[32]={};
        for(int j=1;j<p;++j)q[j*j%p]=true;
        for(int i=0;i<p;++i)for(int j=0;j<p;++j)if(i!=j)a[i][j]=q[(j-i+p)%p]?1:-1;
        auto h=detect(a);int norm=0;
        for(auto [d,c]:h) {assert(d%3==0); norm+=(d/3)*(d/3)*c;}
        assert(6*norm==p*(p-1)*(p-3)*(p+1));
    }
    std::cout<<"Doubly regular sum-of-squares verified in orders 7, 11, 19, 23, 31\n";
}
bool prime_module_decomposition(const Mat& a){
    const unsigned all=(1u<<a.size())-1;
    for(unsigned s=1;s<all;++s){
        if(std::popcount(s)<2)continue;
        bool module=true;
        for(unsigned i=0;i<a.size()&&module;++i)if(!(s&(1u<<i))){
            const int sign=a[i][std::countr_zero(s)];
            for(unsigned j=0;j<a.size();++j)if((s&(1u<<j))&&a[i][j]!=sign){
                module=false;break;
            }
        }
        if(module)return false;
    }
    return true;
}
int main(int argc,char** argv){
    if(argc==1){finite();return 0;}
    int n=std::stoi(argv[1]);assert(n==9||n==11||n==13);
    std::string s,minimizer;int seen=0,mn=n*n*n,mx=0,witness=0,minima=0,first=0;
    while(std::cin>>s){
        ++seen;
        auto h=detect(from_bits(n,s));
        const int k=n*(n-1)*(n-2)/6-h[0];
        if(k>0)++witness;
        if(k<mn) {mn=k;minima=0;first=seen;minimizer=s;}
        if(k==mn)++minima;
        mx=std::max(mx,k);
    }
    assert(seen==(n==9?15:n==11?1223:1495297));
    assert(seen==witness && minima==1 && prime_module_decomposition(from_bits(n,minimizer)));
    std::cout<<n<<" vertices: "<<witness<<"/"<<seen<<" detected, detecting-triple range "<<mn<<".."<<mx<<", minimizers "<<minima<<" (first line "<<first<<"), no nontrivial proper module\n";
}

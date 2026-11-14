#include <array>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "mat.hpp"
using namespace std;
using Z=long long;
using V=array<Z,3>;
Z mx;
vector<bool> seen;
Z g(Z q){Z a=1;for(;q;q/=2)if(!(q&1))a+=q/2+((q/2)&1);return a;}
void reach(V s){
    if(s[0]>mx)return;
    if(s[0]<1||!(s[0]&1)||s[2]<0||s[2]>s[0]||s[1]<s[0])abort();
    seen[s[0]]=true;
    for(const auto& M:MX){
        V t{};
        for(int i=0;i<3;i++)for(int j=0;j<3;j++)t[i]+=M[3*i+j]*s[j];
        reach(t);
    }
}
int main(int argc,char**argv){
    mx=argc>1?atoll(argv[1]):418429;
    seen.resize(mx+1);
    for(int u=0;u<16;u++)reach({g(u),u+1,1-(u&1)});
    for(Z n=1;n<=mx;n+=2)if(!seen[n])cout<<n<<'\n';
}

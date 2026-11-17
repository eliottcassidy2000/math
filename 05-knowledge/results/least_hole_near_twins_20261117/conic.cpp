#include <cmath>
#include <iostream>
#include <initializer_list>
using I=long long;
int main(){
 for(I N: {29LL,89LL,1253LL,10007LL,30011LL,1000000007LL}){
   std::cout<<N<<':';
   for(I K=1;K<=N/29;K+=2){
     I M=13*K*K+(4*N-24)*K+16;
     I W=sqrtl(M);
     while (W*W>M) --W;
     while ((W+1)*(W+1)<=M) ++W;
     if(W*W==M && (W+K-4)%4==0)
       std::cout<<' '<<K<<','<<W;
   }
   std::cout<<'\n';
 }
}

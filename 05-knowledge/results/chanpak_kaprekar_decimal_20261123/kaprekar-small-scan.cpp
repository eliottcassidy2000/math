#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
using namespace std;
template<int n> void run(int x, int y, int count){
 constexpr int F=(1<<n)-1; using T=array<unsigned,n>;
 array<array<unsigned,n>,1<<n>D{};
 auto height=[&](const T&A){
  for(unsigned S=1;S<=F;S++){
   for(unsigned s=S;s;s&=s-1){int i=__builtin_ctz(s);unsigned v=S^(1<<i),z=0;
    if(!v)z=1;
    else for(unsigned m=A[i]&v;m;m&=m-1)z+=D[v][__builtin_ctz(m)];
    D[S][i]=z;
   }
  }
  int h=0;for(auto z:D[F])h+=z;return h;
 };
 constexpr int deg=(n-1)*(n-2)/2;int mx=0,cx=0,cy=0,sz=0;
 for(unsigned k=0;k<1u<<deg;k++){
  T A{};unsigned b=k;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
   if(j==i+1 || (b&1))A[i]|=1<<j;else A[j]|=1<<i;
   if(j>i+1)b>>=1;
  }
  int h=height(A);mx=max(mx,h);if(h==y)++cy;if(h!=x)continue;
  ++cx; bool plus=false,minus=false;
  for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
   auto B=A;B[i]^=1<<j;B[j]^=1<<i;
   int h=height(B);if(h==y){
    if(sz++<1){cout<<n<<" "<<k<<" "<<x<<" "<<y<<" flip "<<i<<" "<<j<<" rows";for(auto z:A)cout<<" "<<z;cout<<endl;}
   }
  }
 }
 cout<<"summary "<<n<<" max "<<mx<<" h counts "<<cx<<" "<<cy<<" seam presentations "<<sz<<endl;
}
int main(int argc,char**argv){run<8>(atoi(argv[1]),atoi(argv[2]),0);}

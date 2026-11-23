#include <cstdint>
#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <cstdlib>
using namespace std;
int main(int argc,char**argv){
 int p=atoi(argv[1]),d=(p-1)/2,m=p*p;
 uint64_t M=1;for(int i=0;i<3*d;i++)M*=p;
 if(M>100000000)return 2;
 vector<bool> seen(M);vector<uint32_t> q;auto test=[&](uint32_t x){if(!seen[x]){seen[x]=true;q.push_back(x);}};
 // t² is the actual even seed T2.
 uint32_t b0=0,z=1;
 for(int i=1;i<=d;i++){b0+=i*i*z;z*=m;}
 for(int i=1;i<=d;i++){b0+=(2*i%p)*z;z*=p;}
 test(b0);
 for(size_t h=0;h<q.size();h++){
  auto v=q[h];array<int,10> a{},b{};
  for(int i=1;i<=d;i++){a[i]=v%m;v/=m;}
  for(int i=1;i<=d;i++){b[i]=v%p;v/=p;}
  a[d+1]=(a[d]-p*b[d]+m)%m;b[d+1]=(p-b[d])%p;
  uint32_t K=0,L=0,z=1;array<int,10> c{},e{};
  for(int i=1;i<=d;i++){
   K+=(i*i*a[i]%m)*z;
   L+=((i*i*a[i]+i*(a[i+1]-a[i-1])+m*m)%m)*z;
   c[i]=(i*i*b[i]+2*i*a[i])%p;
   e[i]=(i*i*b[i]+2*i*a[i]+a[i+1]-a[i-1]+i*(b[i+1]-b[i-1])+m*m)%p;
   z*=m;
  }
  for(int i=1;i<=d;i++){K+=c[i]*z;L+=e[i]*z;z*=p;}
  test(K);test(L);
 }
 uint64_t n=1;for(int i=0;i<d;i++)n*=p;
 cout<<p<<" orbit "<<q.size()<<" maximal "<<(n-1)*n*n<<endl;
}

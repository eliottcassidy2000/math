// Complete action-monoid enumeration over F7; includes source charge.
#include <array>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using A=std::array<int,7>; // Q00,Q01,Q10,Q11,t0,t1,nonzero charge multiplier
constexpr int POW6=117649,N=6*POW6;
int code(A a){int v=0,p=1;for(int i=0;i<6;++i){v+=p*a[i];p*=7;}return v+(a[6]-1)*POW6;}
A decode(int x){A a;a[6]=x/POW6+1;x%=POW6;for(int i=0;i<6;++i){a[i]=x%7;x/=7;}return a;}
A mul(A a,A b){return {(a[0]*b[0]+a[1]*b[2])%7,(a[0]*b[1]+a[1]*b[3])%7,
 (a[2]*b[0]+a[3]*b[2])%7,(a[2]*b[1]+a[3]*b[3])%7,
 (a[0]*b[4]+a[1]*b[5]+a[4]*b[6])%7,(a[2]*b[4]+a[3]*b[5]+a[5]*b[6])%7,a[6]*b[6]%7};}
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
int main(int argc,char**argv){need(argc==2,"generator input");std::ifstream f(argv[1]);
 std::vector<A>G(4);for(auto&a:G)for(int&x:a)f>>x;need(bool(f),"four complete generators");
 for(int k:{2,3,4}){std::vector<int>dist(N,-1),todo;int root=code({1,0,0,1,0,0,1});todo.push_back(root);dist[root]=0;
  int maxd=0;for(size_t i=0;i<todo.size();++i){A a=decode(todo[i]);for(int j=0;j<k;++j){int z=code(mul(a,G[j]));
   if(dist[z]<0){dist[z]=dist[todo[i]]+1;maxd=std::max(maxd,dist[z]);todo.push_back(z);}}}
  std::array<int,3>rank{};
  for(int x=0;x<N;++x){A a=decode(x);int d=(a[0]*a[3]-a[1]*a[2]+49)%7;
   int r=d?2:((a[0]||a[1]||a[2]||a[3])?1:0);
   bool transverse=(a[0]*a[5]-a[2]*a[4])%7 || (a[1]*a[5]-a[3]*a[4])%7;
   bool singular=(r==1&&transverse)||(r==0&&(a[4]||a[5]));
   bool expected=k==4 || ((a[6]==1||a[6]==6) && ((r==2&&!a[4]&&!a[5])||(k==3&&singular)));
   need((dist[x]>=0)==expected,"whole-universe action membership");if(dist[x]>=0)++rank[r];}
  std::cout<<"generators "<<k<<" states "<<todo.size()<<" rank0 "<<rank[0]<<" rank1 "<<rank[1]<<" rank2 "<<rank[2]<<" maximum length "<<maxd<<'\n';
 }
}

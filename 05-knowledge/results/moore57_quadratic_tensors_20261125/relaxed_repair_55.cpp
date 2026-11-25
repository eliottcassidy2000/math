#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;using U=uint64_t;
U bit(int i){return 1ULL<<i;}
struct E{int y;array<int,3> c;U blockers=0;};
vector<E> tab[53];int old[53],sol[53],n;
U init;uint64_t visits=0, nsln=0; int lim; int minextra=53; int best[53];
bool dfs(U open,U active,array<U,3> used){
 visits++;
 if(!open){
   nsln++;U seen=0;for(int i=0;i<n;i++)seen|=bit(sol[i]/56);
   int dup=n-__builtin_popcountll(seen);
   if(dup<minextra){minextra=dup;copy(sol,sol+n,best);}return false;
 }
 int chosen=-1, mn=100000, rem=lim-__builtin_popcountll(init^active);
 vector<const E*> opts;
 for(U f=open;f;f&=f-1){
    int h=__builtin_ctzll(f); vector<const E*> v;
    for(const E&e:tab[h])
      if(!(used[0]&bit(e.c[0]))&&!(used[1]&bit(e.c[1]))&&!(used[2]&bit(e.c[2]))&&
         __builtin_popcountll(e.blockers&active)<=rem)v.push_back(&e);
    if(v.empty())return false;
    if((int)v.size()<mn){chosen=h;mn=v.size();opts=move(v);}
 }
 sort(opts.begin(),opts.end(),[&](auto x,auto y){return __builtin_popcountll(x->blockers&active)<__builtin_popcountll(y->blockers&active);});
 for(auto ep:opts){auto&e=*ep;U drop=e.blockers&active;
   auto V=used;for(int k=0;k<3;k++)V[k]|=bit(e.c[k]);
   sol[chosen]=e.y;
   if(dfs((open^bit(chosen))|drop,active^drop,V))return true;
   sol[chosen]=old[chosen];
 }
 return false;
}
int main(int argc,char**argv){
 cin>>n;int own[3][56];for(auto&rr:own)fill(rr,rr+56,-1);
 for(int r=0;r<n;r++){int num;cin>>old[r]>>num;
  if(old[r]>=0)init|=bit(r);
  while(num--){E e;cin>>e.y>>e.c[0]>>e.c[1]>>e.c[2];
    if(e.y==old[r])for(int k=0;k<3;k++)own[k][e.c[k]]=r;
    tab[r].push_back(e);
  }
 }
 for(int r=0;r<n;r++)for(auto&e:tab[r])for(int k=0;k<3;k++)if(own[k][e.c[k]]>=0)e.blockers|=bit(own[k][e.c[k]]);
 int bound=argc>1?stoi(argv[1]):7;
 for(lim=0;lim<=bound;lim++){copy(old,old+n,sol);visits=0;
   if(dfs(((1ULL<<n)-1)^init,init,{0,0,0})){
      cerr<<"solved limit "<<lim<<" visits "<<visits<<endl;
      for(int i=0;i<n;i++)cout<<sol[i]<<" ";cout<<endl;return 0;
   }
   cerr<<"budget "<<lim<<" visits "<<visits<<" solutions "<<nsln<<" min repeats "<<minextra<<endl;
   if(nsln){for(int i=0;i<n;i++)cout<<best[i]<<" ";cout<<endl;return 0;}
 }
 return 1;
}


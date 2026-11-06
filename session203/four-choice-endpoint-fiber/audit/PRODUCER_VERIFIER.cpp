#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;using I=long long;using G=vector<unsigned>;
void need(bool ok,const char*s){if(!ok)throw runtime_error(s);}
struct Data{vector<I>h;vector<vector<I>>start,end;};
Data profile(const G&a){int n=a.size(),N=1<<n;Data z{vector<I>(N),vector<vector<I>>(N,vector<I>(n)),vector<vector<I>>(N,vector<I>(n))};z.h[0]=1;for(int i=0;i<n;i++)z.start[1<<i][i]=z.end[1<<i][i]=1;for(int S=1;S<N;S++)for(int i=0;i<n;i++)if(S>>i&1){int T=S^(1<<i);if(T)for(int j=0;j<n;j++)if(T>>j&1){if(a[i]>>j&1)z.start[S][i]+=z.start[T][j];else z.end[S][i]+=z.end[T][j];}z.h[S]+=z.start[S][i];}return z;}
array<I,8> profile_data(const G&a){auto z=profile(a);int N=1<<a.size(),F=N-1,T=F^1;I L=0,R=0,beta=0;for(int S=1;S<N;S++){L+=z.end[S][0]*z.h[F^S];R+=z.start[S][0]*z.h[F^S];}for(int S=T;;S=(S-1)&T){beta+=z.h[S]*z.h[T^S];if(!S)break;}I h=z.h[F],s=z.start[F][0],e=z.end[F][0],d=z.h[T];return {h+L-s,e,R+beta,d,s,h,L,R};}
G inflate(G a,int v,G x){int n=a.size(),m=x.size();vector<vector<int>>slot(n);int z=0;for(int i=0;i<n;i++)for(int j=0;j<(i==v?m:1);j++)slot[i].push_back(z++);G b(z);for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(a[i]>>j&1)for(int q:slot[i])for(int r:slot[j])b[q]|=1u<<r;for(int i=0;i<m;i++)for(int j=0;j<m;j++)if(x[i]>>j&1)b[slot[v][i]]|=1u<<slot[v][j];return b;}
I power(I a,int k){I z=1;while(k--)z*=a;return z;}
I expected(int q,int k){I z=0;switch(q){
case 0: z+=(1)*power(0,k);z+=(17)*power(1,k);z+=(-693)*power(2,k);z+=(5376)*power(4,k);z+=(-13200)*power(6,k);z+=(23328)*power(9,k);break;
case 1: break;
case 2: z+=(16)*power(0,k);z+=(-322)*power(1,k);z+=(-396)*power(2,k);z+=(15072)*power(4,k);z+=(-66000)*power(6,k);z+=(165816)*power(9,k);z+=(-182448)*power(12,k);z+=(221184)*power(16,k);break;
case 3: z+=(11)*power(1,k);z+=(-132)*power(2,k);z+=(648)*power(4,k);z+=(-1100)*power(6,k);z+=(1944)*power(9,k);break;
case 4: z+=(8)*power(1,k);z+=(-110)*power(2,k);z+=(456)*power(4,k);z+=(-1100)*power(6,k);z+=(1944)*power(9,k);break;
case 5: z+=(1)*power(0,k);z+=(9)*power(1,k);z+=(-341)*power(2,k);z+=(2480)*power(4,k);z+=(-6600)*power(6,k);z+=(11664)*power(9,k);break;
case 6: z+=(16)*power(1,k);z+=(-462)*power(2,k);z+=(3352)*power(4,k);z+=(-7700)*power(6,k);z+=(13608)*power(9,k);break;
case 7: z+=(8)*power(0,k);z+=(-166)*power(1,k);z+=(132)*power(2,k);z+=(4832)*power(4,k);z+=(-26400)*power(6,k);z+=(71244)*power(9,k);z+=(-91224)*power(12,k);z+=(110592)*power(16,k);break;
}return z;}
int main(){G core={102,0,10,83,135,30,182,47};int total=0;for(int k=0;k<=3;k++)for(int z=0;z<4;z++){int seed[4]={10,14,17,21},offset[4]={0,8,14,22};int m=5+2*k;G x(m);for(int i=0;i<m;i++)for(int j=i+1;j<m;j++)x[i]|=1u<<j;auto a=inflate(core,5,x);unsigned word=(seed[z]<<k)|((1u<<k)-1);for(int j=0;j<m;j++){int v=5+j;if(word>>j&1){a[3]|=1u<<v;a[v]&=~8u;}else{a[3]&=~(1u<<v);a[v]|=8u;}}auto p=profile_data(a);for(int q=0;q<8;q++){I want=expected(q,k)+((q==4||q==5)?offset[z]*(1LL<<k):0);need(p[q]==want,"wrapped profile");}cout<<"k="<<k<<" seed="<<seed[z]<<" vertices="<<a.size()<<" profile=";for(auto v:p)cout<<v<<",";cout<<"\n";total++;}cout<<"PASS graphs="<<total<<"\n";}

#include <iostream>
#include <vector>
#include <array>
#include <cstdint>
using namespace std;
using A=array<int,16>;
int step(const A&a,int x){
 int z0=x&31,z1=(x>>5)&31,z2=(x>>10)&31,z3=x>>15;
 int r=0;for(int i=0;i<4;i++)r|=((a[4*i]*z0+a[4*i+1]*z1+a[4*i+2]*z2+a[4*i+3]*z3)&31)<<(5*i);
 return r;
}
int main(){
 int n,b=0;cin>>n;
 array<int,4> f;
 for(int i=0,x;i<4;i++){cin>>x;b|=(x&31)<<(5*i);}
 for(auto &x:f){cin>>x;x&=31;}
 auto h=[&](int x){int z=0;for(int i=0;i<4;i++)z+=f[i]*((x>>(5*i))&31);z&=31;if(!(z&1)){cerr<<"even\n";exit(1);}return 1<<(z>>1);};
 vector<A>A0(n);
 for(auto&a:A0)for(auto&x:a){cin>>x;x&=31;}
 vector<bool> seen(1<<20,false);vector<int>q={b};seen[b]=true;
 vector<uint16_t>cov(1<<20,0);
 array<int,17> hist1={0},hist2={0};
 for(size_t i=0;i<q.size();i++){
  int x=q[i];
  for(auto&a:A0){
   int y=step(a,x);cov[x]|=h(y);
   if(!seen[y]){seen[y]=true;q.push_back(y);}
  }
  hist1[__builtin_popcount(cov[x])]++;
 }
 for(auto x:q){
  int v=0;
  for(auto&a:A0){v|=cov[step(a,x)];if(v==65535)break;}
  hist2[__builtin_popcount(v)]++;
 }
 if(q.size()!=524288||hist1[15]!=1280||hist1[16]!=523008||hist2[16]!=524288) return 2;
 cout<<"closed "<<q.size()<<" states\n";
 for(int i=0;i<=16;i++)if(hist1[i]||hist2[i])cout<<i<<" "<<hist1[i]<<" "<<hist2[i]<<"\n";
}

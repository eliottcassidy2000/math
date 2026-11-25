#include <array>
#include <bit>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
const int n=8,M=1<<n;
using Graph=array<int,n>;
Graph graph(int base,int ext){
 Graph a{};
 int k=0;
 for(int i=0;i<6;i++)for(int j=i+1;j<6;j++,k++)
  if(base>>k&1)a[j]|=1<<i;else a[i]|=1<<j;
 for(int i=0;i<6;i++)for(int j=6;j<8;j++)
  if(ext>>(2*i+j-6)&1)a[j]|=1<<i;else a[i]|=1<<j;
 a[6]|=1<<7;
 return a;
}
struct Shape{
 int p[M]={},e[M][n]={},s[M][n]={};
 int h,Q,R[n]={},J[n]={},C[n][n]={};
 Shape(const Graph& a){
  Graph back{};
  for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(a[i]>>j&1)back[j]|=1<<i;
  p[0]=1;
  for(int i=0;i<n;i++)e[1<<i][i]=s[1<<i][i]=1;
  for(int S=1;S<M;S++)for(int v=S;v;v&=v-1){
   int j=countr_zero((unsigned)v),T=S^(1<<j);
   for(int x=back[j]&T;x;x&=x-1)e[S][j]+=e[T][countr_zero((unsigned)x)];
   for(int x=a[j]&T;x;x&=x-1)s[S][j]+=s[T][countr_zero((unsigned)x)];
   p[S]+=e[S][j];
  }
  h=p[M-1];Q=0;
  for(int S=0;S<M;S++)Q+=p[S]*p[M-1-S];
 }
 void prep(){
  for(int S=1;S<M;S++)for(int v=S;v;v&=v-1){
   int j=countr_zero((unsigned)v),T=M-1-S;
   R[j]+=p[T]*s[S][j];J[j]+=e[S][j]*p[T];
   for(int w=T;w;w&=w-1){int k=countr_zero((unsigned)w);C[j][k]+=e[S][j]*s[T][k];}
  }
 }
 array<int,4> data(int o)const{
  int r=h,j=0,c=0,e0=0;
  for(int v=0;v<n;v++)if(o>>v&1){r+=R[v];c+=s[M-1][v];}
  else{
    j+=J[v];c+=e[M-1][v];e0+=e[M-1][v];
    for(int w=0;w<n;w++)if(o>>w&1)c+=C[v][w];
  }
  return {r,j,c,e0};
 }
};
int main(){
 int stable=0,seams=0,newh=0;
 map<int,int> types;
 for(int y=0;y<(1<<12);y++){
   auto a=graph(177,y),b=graph(181,y);
   Shape x(a), z(b);
   if(x.h!=z.h||x.Q!=z.Q)continue;
   stable++;
   x.prep();z.prep();
   for(int o=0;o<M;o++){
     auto [r,j,c,e]=x.data(o);auto [R,J,C,E]=z.data(o);
     if(r!=R||j!=J||abs(C-c)!=2)continue;
     seams++;++types[x.h];
     if(x.h!=29&&newh++<40)cout<<"primitive? "<<y<<" "<<o<<" "<<x.h<<" "<<x.Q<<" "<<r<<" "<<j<<" "<<c<<" "<<C<<" e "<<e<<" "<<E<<endl;
   }
 }
 cout<<"stable "<<stable<<" seams "<<seams<<" with different h "<<newh<<endl;
 for(auto [h,q]:types)cout<<h<<":"<<q<<" ";cout<<endl;
}
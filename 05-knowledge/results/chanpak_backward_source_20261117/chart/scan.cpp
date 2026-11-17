#include <array>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <numeric>
using namespace std;
using ll=long long;
using V=array<ll,2>;
struct Hash {size_t operator()(V v)const{return (size_t)v[0]*0x9e3779b97f4a7c15ULL ^ ((size_t)v[1]+0x85ebca6b);}};
int main(){
 for(int m=2;m<=30;m+=2){
  vector<vector<V>> K(m,vector<V>(m));
  for(int q=1;q<=2;q++){
   vector<ll> p(m+1,1),P(m+1,1),M(m+1,1);
   for(int j=1;j<=m;j++)p[j]=q*p[j-1],P[j]=(q+1)*P[j-1],M[j]=(q-1)*M[j-1];
   for(int i=0;i<m;i++)for(int j=i;j<m;j++)K[i][j][q-1]=(i==j?p[m-1]:p[m-1-j+i]*(P[j-i-1]-M[j-i-1])/2);
  }
  int ct=0;
  for(int x=0;x<m;x++)for(int y=0;y<m;y++)if(x!=y){
   V base{}, C[m];
   for(int q=0;q<2;q++){
    for(int z=0;z<m;z++){
      C[z][q]=K[z][x][q]-K[x][z][q]-K[z][y][q]+K[y][z][q];
      base[q]+=K[x][z][q]-K[y][z][q];
    }
    base[q]+=C[y][q];
   }
   const V target={5-base[0],-2-base[1]};
   vector<int> free;for(int z=0;z<m;z++)if(z!=x&&z!=y)free.push_back(z);
   int n=m-2,a=n/2,b=n-a;
   unordered_map<V,int,Hash> T;T.reserve(1<<a);
   vector<V>P(1<<a),Q(1<<b);
   for(int z=0;z<(1<<a);z++){
    if(z){int j=__builtin_ctz(z),u=z&(z-1);P[z]={P[u][0]+C[free[j]][0],P[u][1]+C[free[j]][1]};}
    T.try_emplace(P[z],z);
   }
   for(int z=0;z<(1<<b);z++){
    if(z){int j=__builtin_ctz(z),u=z&(z-1);Q[z]={Q[u][0]+C[free[a+j]][0],Q[u][1]+C[free[a+j]][1]};}
    auto it=T.find({target[0]-Q[z][0],target[1]-Q[z][1]});
    if(it!=T.end()){
       int Y=1<<y;
       int s=it->second|(z<<a);
       for(int j=0;j<n;j++)if(s>>j&1)Y|=1<<free[j];
       int X=Y^(1<<x)^(1<<y);
       if(!ct)cout<<"HIT "<<m<<" "<<x<<" "<<y<<" "<<X<<" "<<Y<<endl;
       ++ct;
    }
   }
  }
  cout<<"size "<<m<<" hits "<<ct<<endl;
 }
}

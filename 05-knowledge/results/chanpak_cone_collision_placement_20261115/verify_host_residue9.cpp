#include <cassert>
#include <cstdint>
#include <array>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>
using namespace std;
const int N=120, MOD=9;
using M=vector<vector<int>>;
M mat(int a,int b){
 M m(N);
 int bits[2]={a,b};
 int src[4][2]={{0,0},{1,0},{0,1},{0,1}};
 int sink[4][2]={{0,0},{1,0},{0,1},{1,0}};
 for(int k=1;k<=2;k++){
  int off=(k-1)*60;
  for(int e=0;e<4;e++){
   int len=e<2?1:2; int sz=e<2?5:25;
   for(int code=0;code<sz;code++){
    int i=off+code;
    for(int z=0;z<k;z++)m[i].push_back(i);
    for(int j=0;j<len;j++){
     int place=(j==0 && len==2)?5:1;int r=(code/place)%5;
     if(r>=3 || (r==0 && !bits[src[e][j]]))continue;
     int q=r==1?2:1;
     m[i].push_back(i+(q-r)*place);
     if(!bits[sink[e][j]])m[i].push_back(i+(q+2-r)*place);
    }
   }
   off+=sz;
  }
 }
 return m;
}
string go(const string&s,const M&m){
 array<uint16_t,N> w={};
 for(int i=0;i<N;i++){int z=s[i];if(z)for(int j:m[i])w[j]+=z;}
 string t(N,0);for(int i=0;i<N;i++)t[i]=w[i]%MOD;
 return t;
}
string back(const string&s,const M&m){
 string t(N,0);
 for(int i=0;i<N;i++){int z=0;for(int j:m[i])z+=s[j];t[i]=z%MOD;}
 return t;
}
string accept(){string v(N,0);
 for(int k=1;k<=2;k++){
   int z=k==1?(-209%MOD+MOD)%MOD:259%MOD;int base=(k-1)*60;
   v[base+4]=v[base+9]=2*k*k*z%MOD;v[base+34]=4*k*k*z%MOD;
   for(int a=0;a<5;a++)for(int b=0;b<5;b++)if((a>=3||a==0)&&b>=3&&((a==3)+(b==3))==1)v[base+35+5*a+b]=2*k*z%MOD;
 }
 return v;
}
int main(){
 auto m0=mat(0,1),m1=mat(1,1),mz=mat(0,0);
 string init(N,0);for(int i:{0,5,10,35})init[i]=init[60+i]=1;
 vector<string> s={init};unordered_map<string,int> ix={{init,0}};vector<array<int,2>> tr;
 for(int i=0;i<(int)s.size();i++){
   auto a=go(s[i],m0),b=go(s[i],m1); 
   int k0=ix.emplace(a,s.size()).first->second;if(k0==(int)s.size())s.push_back(a);
   int k1=ix.emplace(b,s.size()).first->second;if(k1==(int)s.size())s.push_back(b);
   tr.push_back({k0,k1});
 }

 assert(s.size()==209955);
 auto even=[](const string&z){return z[60]==1||z[60]==4||z[60]==7;};
 int even_count=0;
 for(auto &z:s){assert(z[60]==1||z[60]==2||z[60]==4||z[60]==8||z[60]==7||z[60]==5);even_count+=even(z);}
 assert(even_count==104977);

 vector<string> vs={accept()};
 while(true){auto a=back(vs.back(),mz);auto it=find(vs.begin(),vs.end(),a);if(it!=vs.end()){assert((it-vs.begin())==2 && vs.size()==8);break;}vs.push_back(a);}

 // Positive even suffix lengths give the acceptance rows at 2, 4, and 6.
 for(int j: {2,4,6}){
   vector<int> r(s.size());
   for(int i=0;i<(int)s.size();i++){
     int x=0;
     for(int k=0;k<N;k++)x+=s[i][k]*vs[j][k];
     r[i]=1<<(x%MOD);
   }
   for(int d=1;d<=8;d++){
     vector<int> R(s.size());
     for(int i=0;i<(int)s.size();i++)R[i]=r[tr[i][0]]|r[tr[i][1]];
     r.swap(R);
   }
   for(int i=0;i<(int)s.size();i++)if(even(s[i]))assert(r[i]==511);
 }
 cout<<"209955 reachable states; all 104977 even-prefix states admit every residue mod 9 with eight repair bits.\n";
}

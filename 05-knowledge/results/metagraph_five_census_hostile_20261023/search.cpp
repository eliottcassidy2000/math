#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>
using namespace std;
int canon[1024],kind[1024];
vector<pair<int,int>> pairs;
void prepare(){
 for(int i=0;i<5;i++)for(int j=i+1;j<5;j++)pairs.push_back({i,j});
 for(int c=0;c<1024;c++){
  bool a[5][5]={};for(int b=0;b<10;b++){auto [i,j]=pairs[b];a[i][j]=(c>>b)&1;a[j][i]=!a[i][j];}
  array<int,5> p={0,1,2,3,4};int best=1024;
  do{int q=0;for(int b=0;b<10;b++){auto [i,j]=pairs[b];q|=int(a[p[i]][p[j]])<<b;}best=min(best,q);}while(next_permutation(p.begin(),p.end()));
  canon[c]=best;
 }
 map<int,int> ids;for(int c=0;c<1024;c++)ids.emplace(canon[c],0);
 int id=0;for(auto &v:ids)v.second=id++;
 if(id!=12)throw runtime_error("five-class count");
 for(int c=0;c<1024;c++)kind[c]=ids[canon[c]];
}
uint64_t mult(uint64_t s,int a,int p){uint64_t t=0;for(int i=1;i<p;i++)if(s>>i&1)t|=1ULL<<((a*i)%p);return t;}
int main(int argc,char**argv){
 prepare();vector<int> primes={13,17,19,23,31};if(argc>1)primes={atoi(argv[1])};
 for(int p:primes){
  auto start=chrono::steady_clock::now();int d=(p-1)/2;long reps=0;map<array<int,12>,uint64_t> seen;long collisions=0;
  for(uint64_t mask=0;mask<(1ULL<<d);mask++){
   uint64_t s=0;for(int i=1;i<=d;i++)s|=1ULL<<((mask>>(i-1)&1)?i:p-i);
   uint64_t smallest=s;for(int a=1;a<p;a++)smallest=min(smallest,mult(s,a,p));if(smallest!=s)continue;
   reps++;bool a[31][31]={};for(int i=0;i<p;i++)for(int j=0;j<p;j++)if(i!=j)a[i][j]=(s>>((j-i+p)%p))&1;
   array<int,12> counts={};
   for(int v0=0;v0<p;v0++)for(int v1=v0+1;v1<p;v1++)for(int v2=v1+1;v2<p;v2++)for(int v3=v2+1;v3<p;v3++)for(int v4=v3+1;v4<p;v4++){
    int c=a[v0][v1]|(a[v0][v2]<<1)|(a[v0][v3]<<2)|(a[v0][v4]<<3)|(a[v1][v2]<<4)|(a[v1][v3]<<5)|(a[v1][v4]<<6)|(a[v2][v3]<<7)|(a[v2][v4]<<8)|(a[v3][v4]<<9);
    counts[kind[c]]++;
   }
   auto [it,inserted]=seen.emplace(counts,s);
   if(!inserted){collisions++;if(collisions<=5){cout<<"COLLISION p="<<p<<" masks="<<it->second<<","<<s<<" counts=";for(int x:counts)cout<<x<<",";cout<<"\n"<<flush;}}
  }
  double sec=chrono::duration<double>(chrono::steady_clock::now()-start).count();cout<<"SUMMARY p="<<p<<" reps="<<reps<<" keys="<<seen.size()<<" collisions="<<collisions<<" seconds="<<sec<<"\n"<<flush;
 }
}

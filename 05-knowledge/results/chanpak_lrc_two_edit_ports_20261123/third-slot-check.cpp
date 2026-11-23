#include <algorithm>
#include <array>
#include <cassert>
#include <cstdio>
#include <numeric>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>
using namespace std;
using I=long long;using V=vector<int>;
vector<I> allheights(const V& r){
 const int n=r.size(),N=1<<n;
 vector<I> dp(size_t(N)*n),H(N); H[0]=1;
 for(int i=0;i<n;i++) {
  assert(!(r[i]>>i&1));
  for(int j=0;j<i;j++)assert(((r[i]>>j&1)+(r[j]>>i&1))==1);
  dp[(1<<i)*n+i]=1;
 }
 for(int s=1;s<N;s++)for(int z=s;z;z&=z-1){
  const int i=__builtin_ctz(z);I d=dp[size_t(s)*n+i];
  if(d){H[s]+=d; for(int b=r[i]&~s;b;b&=b-1){
    const int j=__builtin_ctz(b);dp[size_t(s|1<<j)*n+j]+=d;
  }}
 }
 return H;
}
array<I,3> profile(const V& r){
 const auto h=allheights(r);const int u=h.size()-1;
 I q=0,t=0;
 for(int s=0;s<=u;s++){
  const int c=u^s;q+=h[s]*h[c];
  for(int j=c;;j=(j-1)&c) {t+=h[s]*h[j]*h[c^j];if(j==0)break;}
 }
 return {h[u],q,t};
}
V literal(string tags,int block,const V& M,int e,int f,int corner) {
 if(corner&1)tags[e/2]=(tags[e/2]-'0' ^ (1<<(e%2)))+'0';
 if(corner&2)tags[f/2]=(tags[f/2]-'0' ^ (1<<(f%2)))+'0';
 vector<pair<int,int>> z;
 for(int i=0;i<int(tags.size());i++)
  for(int j=0;j<(i==block?int(M.size()):1);j++)z.push_back({i,j});
 const int n=z.size(); V rows(n+2);
 for(int i=0;i<n+2;i++)for(int j=i+1;j<n+2;j++){
  const bool forward=j<n?(z[i].first!=z[j].first || (M[z[i].second]>>z[j].second&1)):
   i<n ? !((tags[z[i].first]-'0')>>(j-n)&1) : false;
  if(forward)rows[i]|=1<<j;else rows[j]|=1<<i;
 }
 return rows;
}
bool prime(I n){if(n<2)return false;for(I d=2;d*d<=n;d++)if(n%d==0)return false;return true;}

int main(int argc,char**argv){
 const array<V,4> fillers={V{86,284,392,401,484,79,398,291,33},
   V{6,244,344,51,417,69,153,301,107},V{94,60,184,112,32,193,150,27},
   V{30,92,216,48,32,199,153,27}};
 const array<array<I,3>,4> p={{{2393,38272,317403},{2393,38272,317643},{345,5640,45225},{345,5640,44505}}};
 const array<array<I,4>,2> h={{{7981807,46862827,17465251,55872649},{7986607,46890667,17475811,55906249}}};
 for(int i=0;i<4;i++){
  assert(profile(fillers[i])==p[i]);
  if(i<2){
   assert(prime(h[i][0]));assert(h[i][0]+h[i][3]-h[i][1]-h[i][2]==-2*236811);
   for(int c=0;c<4;c++)assert(allheights(literal("000210230",6,fillers[i],1,5,c)).back()==h[i][c]);
  }
 }
 assert(317403%113!=0 && 317643%113==0 && 45225%113!=0 && 44505%113!=0);
 assert(45225%67==0 && 44505%67!=0);
 // Optional proof of the two small absolute exclusions. Feed the full scan
 // output containing k=7,8,9. Six-site and earlier forms are contained at 7.
 if(argc==2){
  ifstream src(argv[1]);assert(src);
  string line;int size=0;set<pair<int,int>> rows8,rows9;
  while(getline(src,line)){
   if(line.rfind("k=",0)==0)size=stoi(line.substr(2));
   else if(line.rfind("  ",0)==0 && size>=7&&size<=9){
    int a,b,c;assert(sscanf(line.c_str(),"%d,%d,%d",&a,&b,&c)==3);
    if(c==0){rows9.insert({a,b}); if(size<=8)rows8.insert({a,b});}
   }
  }
  assert(!rows8.empty()&&!rows9.empty());
  set<pair<I,I>> profiles;
  for(int n=1;n<=8;n++)profiles.insert({1,1LL<<n});
  // At cap 1723 and slope at least six, nontransitive modules have order <=7.
  for(int n=1;n<=7;n++){
   vector<pair<int,int>> non;
   for(int i=0;i<n;i++)for(int j=i+2;j<n;j++)non.push_back({i,j});
   for(unsigned w=0;w<(1u<<non.size());w++){
    V r(n); for(int i=0;i+1<n;i++)r[i]|=1<<(i+1);
    for(size_t k=0;k<non.size();k++){
     auto [i,j]=non[k];if(w>>k&1)r[i]|=1<<j;else r[j]|=1<<i;
    }
    auto H=allheights(r);I q=0;
    for(size_t s=0;s<H.size();s++)q+=H[s]*H[H.size()-1-s];
    if(q<=287)profiles.insert({H.back(),q});
   }
  }
  for(auto [bound,p0]: {pair<int,I>{8,359},{9,887},{9,1723}}){
   assert(prime(p0));auto& forms=bound==8?rows8:rows9;
   for(auto [a,b]:forms){
    if(b==0){assert(p0%a!=0);continue;}
    const int g=gcd(abs(a),b);
    if(g>1 && p0%g!=0)continue;
    assert(2*b+min(a,0)>=12);
    for(auto [d,q]:profiles)assert(a*d+b*q!=p0);
   }
  }
 }
 cout<<"Fixture and requested exclusion checks passed\n";
}

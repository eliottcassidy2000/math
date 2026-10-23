#include <array>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <random>
#include <string>
#include <cstring>
using namespace std;
static constexpr int n=8,N=1<<n;
using Poly=array<int,n+1>;
struct Key{array<int,2*(n+1)> a; bool operator==(const Key&o)const{return a==o.a;}};
struct KH{size_t operator()(const Key&k)const{size_t z=0x534ae23; for(int x:k.a)z=z*1000003+unsigned(x);return z;}};
struct Val{uint32_t mask;int h;};
int main(int argc,char**argv){int count=argc>1?stoi(argv[1]):1000000;mt19937 rng(370021);unordered_map<Key,Val,KH> seen;seen.reserve(count/10);for(int it=0;it<count;++it){uint32_t mask=0; int eb=0, mb=0; for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++mb) {bool v=(j==i+1)?true:((it>>eb++)&1);if(v)mask|=1u<<mb;}bool A[n][n]={};int b=0;for(int i=0;i<n;++i)for(int j=i+1;j<n;++j,++b)A[i][j]=(mask>>b)&1,A[j][i]=!A[i][j];Poly p[N]={},d[N]={};p[0][0]=d[0][0]=1;int hp[N][n]={};for(int i=0;i<n;++i)hp[1<<i][i]=1;
for(int s=0;s<N;++s){int row=__builtin_popcount((unsigned)s);if(row<n)for(int j=0;j<n;++j)if(!(s>>j&1)&&(row==j||A[row][j])){int shift=(row==j),sign=(__builtin_popcount((unsigned)(s>>(j+1)))&1)?-1:1;for(int k=0;k<=row;++k){p[s|1<<j][k+shift]+=p[s][k];d[s|1<<j][k+shift]+=sign*d[s][k];}}for(int i=0;i<n;++i)if(hp[s][i])for(int j=0;j<n;++j)if(!(s>>j&1)&&A[i][j])hp[s|1<<j][j]+=hp[s][i];}
Key key;for(int k=0;k<=n;++k){key.a[k]=p[N-1][k];key.a[n+1+k]=d[N-1][k];}int h=0;for(int i=0;i<n;++i)h+=hp[N-1][i];auto [pos,added]=seen.emplace(key,Val{mask,h});if(!added&&pos->second.h!=h){cout<<"iteration "<<it<<" mask1 "<<pos->second.mask<<" H1 "<<pos->second.h<<" mask2 "<<mask<<" H2 "<<h<<"\nper(xI+A) low to high:";for(auto x:p[N-1])cout<<' '<<x;cout<<"\ndet(xI+A) low to high:";for(auto x:d[N-1])cout<<' '<<x;cout<<'\n';return 0;}if(it%10000==0)cerr<<it<<" distinct "<<seen.size()<<'\n';}cout<<"NO COLLISION count "<<count<<" classes "<<seen.size()<<'\n';}

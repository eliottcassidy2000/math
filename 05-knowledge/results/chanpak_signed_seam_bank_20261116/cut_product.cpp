#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <tuple>
#include <vector>
using I=long long;
struct R {I p,c;int a,b; bool operator<(const R&r)const{return std::tie(p,c)<std::tie(r.p,r.c);}};
int main(int argc,char**argv){
 int nmax=argc>1?atoi(argv[1]):18,best=0;
 std::vector<int> J(1<<nmax);
 for(int x=1;x<(1<<nmax);x++){
  int i=__builtin_ctz(x), y=x&(x-1);J[x]=J[y]+(y>>(i+1))+((y>>(i+1))&1);
 }
 std::vector<std::vector<int>> rev(nmax);
 for(int k=1;k<nmax;k++){
  rev[k].resize(1<<k);
  for(int x=1;x<(1<<k);x++)rev[k][x]=(rev[k][x>>1]>>1)+((x&1)<<(k-1));
 }
 for(int n=2;n<=nmax;n++)for(int k=1;k<n;k++){
  const int r=n-k;std::vector<R> rows;
  rows.reserve((size_t)1<<(n-2));
  for(int b=0;b<(1<<(k-1));b++)for(int a=0;a<(1<<(r-1));a++)
   rows.push_back({(I)((1<<k)-b)*((1<<r)+a),1+rev[k][b]-J[b]-rev[r][a]-J[a],a,b});
  std::sort(rows.begin(),rows.end());
  size_t w=0;
  for(size_t i=0;i<rows.size();i++)
   if(!w||rows[i].p!=rows[w-1].p||rows[i].c!=rows[w-1].c) rows[w++]=rows[i];
  rows.resize(w);
  int l=0;
  for(size_t i=0;i<rows.size();i++){
   if(i&&rows[i].p==rows[i-1].p&&rows[i].c==rows[i-1].c+2) l++;else l=1;
   if(l>best){
    best=l;auto t=rows[i];printf("N=%d k=%d product=%lld width=%d\n",n,k,t.p,l);
    for(size_t j=i+1-l;j<=i;j++)printf("a=%d b=%d C=%lld\n",rows[j].a,rows[j].b,rows[j].c);
    fflush(stdout);
   }
  }
 }
}

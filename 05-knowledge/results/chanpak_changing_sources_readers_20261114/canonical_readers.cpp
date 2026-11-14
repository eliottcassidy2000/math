#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <set>
#include <vector>
using namespace std;
using U=uint64_t;
int n,a[8],d[8],p[8];U least;
void rec(int j,unsigned left,U mask){
 if(!left){
  U code=0;
  for(int v=0;v<n;v++){
   unsigned row=0;
   for(int w=0;w<n;w++)row|=((a[p[v]]>>p[w])&1)<<w;
   code=(code<<n)|row;
  }
  least=min(least,code);return;
 }
 int m=d[__builtin_ctz(left)];
 for(unsigned r=left;r;r&=r-1)m=min(m,d[__builtin_ctz(r)]);
 for(unsigned r=left;r;r&=r-1){
  int i=__builtin_ctz(r);
  if(d[i]==m){p[j]=i;rec(j+1,left^(1<<i),mask);}
 }
}
void census(const set<U>& v){
 int all=(1<<n)-1;
 set<array<int,3>> seen;
 for(U x:v){
  int out[8]={};for(int i=n-1;i>=0;i--)out[i]=x&all,x>>=n;
  int p[256][8]={},h[256]={1}, odd[8]={},even[8]={};
  for(int m=1;m<=all;m++)for(int rem=m;rem;rem&=rem-1){
   int i=__builtin_ctz(rem), x=m^(1<<i);
   if(!x)p[m][i]=1;else for(int r=(all^out[i])&x;r;r&=r-1) p[m][i]+=p[x][__builtin_ctz(r)];
   h[m]+=p[m][i];
  }
  for(int m=1;m<=all;m++)for(int rem=m;rem;rem&=rem-1){
   int i=__builtin_ctz(rem);
   (__builtin_parity(m)?odd:even)[i]+=p[m][i]*h[all^m];
  }
  for(int i=0;i<n;i++){
   array<int,3> row={h[all],2*odd[i],2*even[i]};
   if(seen.insert(row).second){cout<<n<<" "<<i;for(auto z:row)cout<<" "<<z;for(int j=0;j<n;j++)cout<<" "<<out[j];cout<<"\n";}
  }
 }
}
int main(int argc,char**argv){
 int limit=argc>1?stoi(argv[1]):8;
 set<U> v={0};
 for(n=1;n<=limit;n++){
  if(n>1){
   set<U> q;int small=(1<<(n-1))-1;
   for(U x:v){
    int old[8];for(int i=n-2;i>=0;i--)old[i]=x&small,x>>=(n-1);
    for(int j=0;j<=small;j++){
     a[n-1]=j;
     for(int i=0;i<n-1;i++)a[i]=old[i]|((1-(j>>i&1))<<(n-1));
     for(int i=0;i<n;i++)d[i]=__builtin_popcount(a[i]);
     least=~U(0);rec(0,(1<<n)-1,0);q.insert(least);
    }
   }
   v.swap(q);
  }
  cerr<<n<<" "<<v.size()<<endl;
  census(v);
 }
}

#include <array>
#include <algorithm>
#include <bit>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
constexpr int N=9,F=(1<<N)-1;
array<int,N>A,deg;
unsigned long long cnt=0,d1=0,full=0,chiral=0;
unsigned dp[2][1<<N][N];
unsigned long long cover[N][1<<N];
bool anti(int v,int used,array<int,N>&p) {
 if(v==N)return true;
 for(int j=0;j<N;j++)if(!(used&(1<<j))){
  bool ok=true;
  for(int i=0;i<v;i++)if(((A[i]>>v)&1)==((A[p[i]]>>j)&1)){ok=false;break;}
  if(ok){p[v]=j;if(anti(v+1,used|(1<<j),p))return true;}
 }
 return false;
}
void test(){
 ++cnt;
 for(int s=1;s<=F;s++)for(int v=0;v<N;v++)if(s&(1<<v)){
  int t=s^(1<<v);
  if(!t){dp[0][s][v]=dp[1][s][v]=1;continue;}
  unsigned a=0,b=0;
  for(int w=t;w;w&=w-1){
   int u=countr_zero((unsigned)w);
   if(A[v]&(1<<u)){a+=dp[0][t][u];}
   else {b+=dp[1][t][u];}
  }
  dp[0][s][v]=a;dp[1][s][v]=b;
 }
 for(int i=0;i<N;i++)if(dp[0][F][i]!=dp[1][F][i])return;
 ++d1;
 array<int,N>P{};
 if(anti(0,0,P))return;
 ++chiral;
 cover[0][0]=1;
 for(int s=1;s<=F;s++){
  cover[1][s]=0;
  for(int t=s;t;t&=t-1)cover[1][s]+=dp[0][s][countr_zero((unsigned)t)];
 }
 cover[1][0]=1;
 for(int j=2;j<N;j++)for(int s=0;s<=F;s++){
  auto val=cover[j-1][s];
  for(int t=s;t;t=(t-1)&s)val+=cover[1][t]*cover[j-1][s^t];
  cover[j][s]=val;
 }
 for(int j=1;j<N-1;j++)for(int v=0;v<N;v++){
  long long d=0;
  for(int s=1<<v;s<=F;s++)if(s&(1<<v))
   d+=((long long)dp[0][s][v]-(long long)dp[1][s][v])* (long long)cover[j][F^s];
  if(d)return;
 }
 ++full;
 cout<<"HIT";for(int a:A)cout<<" "<<a;cout<<endl;
}
void rec(int i,int j){
 if(i==N-1){test();return;}
 if(j==N){ if(deg[i]==N/2)rec(i+1,i+2);return;}
 if(deg[i]<N/2&&deg[i]+(N-j)>=N/2&&deg[j]+(N-i-2)>=N/2){
  A[i]^=1<<j;deg[i]++;rec(i,j+1);deg[i]--;A[i]^=1<<j;
 }
 if(deg[j]<N/2&&deg[i]+(N-j-1)>=N/2&&deg[j]+1+(N-i-2)>=N/2){
  A[j]^=1<<i;deg[j]++;rec(i,j+1);deg[j]--;A[j]^=1<<i;
 }
}
int main(){
 for(int i=1;i<N;i++)if(i<=N/2){A[0]|=1<<i;deg[0]++;}else{A[i]|=1;deg[i]++;}
 rec(1,2);
 cout<<cnt<<" normalized regular "<<d1<<" D1=0 "<<chiral<<" D1=0 chiral "<<full<<" all D=0 chiral\n";
}

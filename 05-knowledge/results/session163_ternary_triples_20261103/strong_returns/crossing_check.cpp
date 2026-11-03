// Exact labelled one-cross-edge control; no assert or NDEBUG dependence.
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
using namespace std;
static long dp[8192][13];
void need(bool p,const char *m){if(!p){fprintf(stderr,"FAIL: %s\n",m);exit(1);}}
long paths(const int *a,int n){
 memset(dp,0,sizeof dp); int all=1<<n;
 for(int i=0;i<n;i++)dp[1<<i][i]=1;
 for(int s=1;s<all;s++)for(int t=s;t;t&=t-1){
  int i=__builtin_ctz((unsigned)t);long x=dp[s][i];
  for(int z=a[i]&~s;z;z&=z-1){int j=__builtin_ctz((unsigned)z);dp[s|(1<<j)][j]+=x;}
 }
 long h=0;for(int i=0;i<n;i++)h+=dp[all-1][i];return h;
}
bool strong(const int *a,int n){
 for(int root=0;root<n;root++){
  int seen=1<<root,old=0;
  while(old!=seen){old=seen;for(int i=0;i<n;i++)if((seen>>i)&1)seen|=a[i];}
  if(seen!=(1<<n)-1)return false;
 }
 return true;
}
int main(){
 const int src[2][10]={{1008,1017,1011,997,776,912,816,592,128,256},{996,993,994,999,975,912,800,576,128,256}};
 const long basewant[2]={855,873};
 const long deletion[2][10]={{163,175,95,15,27,123,87,27,81,165},{97,97,97,159,27,15,81,51,81,159}};
 const long delta[2][10]={{2162,1000,760,1080,2790,4426,3670,4470,10572,12570},{1668,1668,1668,822,846,1524,4770,6288,10050,11508}};
 for(int k=0;k<2;k++){
  for(int drop=0;drop<10;drop++){
   int x[9]={},ii=0;
   for(int i=0;i<10;i++)if(i!=drop){int jj=0;for(int j=0;j<10;j++)if(j!=drop){if((src[k][i]>>j)&1)x[ii]|=1<<jj;jj++;}ii++;}
   need(paths(x,9)==deletion[k][drop],"deleted vertex path count");
  }
  long mindel2=1000000;
  for(int p=0;p<10;p++)for(int q=p+1;q<10;q++){
   int x[8]={},ii=0;
   for(int i=0;i<10;i++)if(i!=p && i!=q){int jj=0;for(int j=0;j<10;j++)if(j!=p && j!=q){if((src[k][i]>>j)&1)x[ii]|=1<<jj;jj++;}ii++;}
   mindel2=min(mindel2,paths(x,8));
  }
  need(mindel2==(k?3:1),"two-vertex deletion minimum");
  int a[13]={};for(int i=0;i<3;i++)a[i]=(1<<((i+1)%3))|(((1<<10)-1)<<3);
  for(int j=0;j<10;j++)a[3+j]=src[k][j]<<3;
  long h=paths(a,13);need(h==basewant[k],"base h");
  for(int i=0;i<3;i++)for(int j=0;j<10;j++){
   a[i]^=1<<(j+3);a[j+3]^=1<<i;
   need(strong(a,13),"repair is strong");need(paths(a,13)-h==delta[k][j],"cross reversal delta");
   a[i]^=1<<(j+3);a[j+3]^=1<<i;
  }
  printf("R%d base %ld, minimum jump %ld, all 30 cross flips checked\n",k?4:2,h,*min_element(delta[k],delta[k]+10));
 }
 puts("PASS");
}

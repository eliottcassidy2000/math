// Literal hypotheses and one-step controls for the moving-guard theorem.
#include <cstdio>
#include <cstdlib>
#include <cstring>
using namespace std;
static long dp[262144][18];
void need(bool v,const char*m){if(!v){fprintf(stderr,"FAIL: %s\n",m);exit(1);}}
long paths(const int*a,int n,int start=-1){int N=1<<n;memset(dp,0,N*sizeof(dp[0]));for(int i=0;i<n;i++)if(start<0||start==i)dp[1<<i][i]=1;for(int m=1;m<N;m++)for(int t=m;t;t&=t-1){int i=__builtin_ctz((unsigned)t);long x=dp[m][i];if(x)for(int z=a[i]&~m;z;z&=z-1){int j=__builtin_ctz((unsigned)z);dp[m|1<<j][j]+=x;}}long r=0;for(int i=0;i<n;i++)r+=dp[N-1][i];return r;}
void query(int*a,int n,int p){a[n]=1<<p;for(int i=0;i<n;i++)if(i!=p)a[i]|=1<<n;}
bool tournament(const int*a,int n){for(int i=0;i<n;i++){if(a[i]>>i&1)return false;for(int j=i+1;j<n;j++)if(((a[i]>>j)&1)+((a[j]>>i)&1)!=1)return false;}return true;}
bool strong(const int*a,int n){for(int v=0;v<n;v++){int seen=1<<v,old=0;while(old!=seen){old=seen;for(int i=0;i<n;i++)if(seen>>i&1)seen|=a[i];}if(seen!=(1<<n)-1)return false;}return true;}
int del(const int*a,int n,int q,int*r,int mark){int ii=0,newmark=-1;for(int i=0;i<n;i++)if(i!=q){if(i==mark)newmark=ii;int jj=0;r[ii]=0;for(int j=0;j<n;j++)if(j!=q){if(a[i]>>j&1)r[ii]|=1<<jj;jj++;}ii++;}return newmark;}
void guard(const int*a,int n,int q,int r){int pred=0;for(int i=0;i<n;i++)if(a[i]>>q&1){need(i==r,"unique guard predecessor");pred++;}need(pred==1,"one predecessor");int b[18];int rr=del(a,n,q,b,r);need(rr>=0 && paths(b,n-1,rr)==0,"no Hamilton path starts at predecessor after guard deletion");}
void extend(const int*x,int n,int q,int*y){for(int i=0;i<n;i++)y[i]=x[i]|((i==q)?(1<<n):0);y[n]=((1<<n)-1)^(1<<q);y[n]|=1<<(n+1);y[n+1]=(1<<n)-1;}
int main(){const int seed[2][14]={{1008,4089,8179,997,776,912,816,592,128,256,3065,1017,12283,4095},{996,993,994,999,1999,912,800,576,128,256,15343,1023,11263,3071}};const int q[2]={12,10},r[2]={2,4};
for(int id=0;id<2;id++){int x[18],y[18],n=14,p=9;for(int i=0;i<n;i++)x[i]=seed[id][i];need(tournament(x,n) && strong(x,n),"strong tournament seed");
 for(int k=0;k<3;k++){
  guard(x,n,q[id],r[id]);long h=paths(x,n);extend(x,n,q[id],y);need(tournament(y,n+2) && strong(y,n+2),"strong tournament extension");need(paths(y,n+2)==3*h,"factor 3 in query graph");guard(y,n+2,n,q[id]);
  printf("seed%d query%d H=%ld, sole predecessor %d, factor 3 and moving guard verified\n",id,k,h,r[id]);
  if(k<2){query(x,n,p);p=n;n++;}
 }
}puts("PASS");}

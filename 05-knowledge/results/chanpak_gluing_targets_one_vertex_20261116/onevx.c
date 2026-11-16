/* minimal excess H(T) - p - q over one-vertex overlaps of induced prime witnesses (cover form, T = A u B) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXN 10
static int n, N; static uint32_t out_[MAXN], H[1<<MAXN], dp[1<<MAXN][MAXN]; static char isprime[1<<14];
static int best[1<<12][64];  /* best[p][qi] keyed loosely: use p<q<512 */
static int mn[512][512];
static void allH(void){H[0]=1;for(int S=1;S<N;S++){uint32_t t=0;for(int v=0;v<n;v++){if(!(S>>v&1))continue;int R=S&~(1<<v);uint32_t c;if(!R)c=1;else{c=0;for(int u=0;u<n;u++)if((R>>u&1)&&(out_[u]>>v&1))c+=dp[R][u];}dp[S][v]=c;t+=c;}H[S]=t;}}
int main(int argc,char**argv){n=atoi(argv[1]);N=1<<n;memset(isprime,1,sizeof isprime);isprime[0]=isprime[1]=0;for(int i=2;i*i<(1<<14);i++)if(isprime[i])for(int j=i*i;j<(1<<14);j+=i)isprime[j]=0;
 for(int i=0;i<512;i++)for(int j=0;j<512;j++)mn[i][j]=1<<30;
 char line[256];int V=N-1;
 while(fgets(line,sizeof line,stdin)){int len=strlen(line);while(len&&(line[len-1]=='\n'))line[--len]=0;if(len!=n*(n-1)/2)continue;memset(out_,0,sizeof out_);int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++){if(line[k]=='1')out_[i]|=1u<<j;else out_[j]|=1u<<i;}
  allH();uint32_t HV=H[V];
  for(int v=0;v<n;v++){uint32_t rest=V&~(1u<<v);
   /* split rest into A' and B' (nonempty), A = A' + v, B = B' + v */
   for(uint32_t Ap=rest;Ap;Ap=(Ap-1)&rest){uint32_t Bp=rest&~Ap;if(!Bp)continue;if(Ap>Bp)continue;
    uint32_t A=Ap|(1u<<v),B=Bp|(1u<<v);uint32_t p=H[A],q=H[B];if(p<3||q<3||!isprime[p]||!isprime[q])continue;
    uint32_t lo=p<q?p:q,hi=p<q?q:p;if(hi>=512)continue;int d=(int)HV-(int)(p+q);if(d<mn[lo][hi])mn[lo][hi]=d;}}}
 for(int p=3;p<512;p++)for(int q=p;q<512;q++)if(mn[p][q]<(1<<30))printf("%d %d %d\n",p,q,mn[p][q]);
 return 0;}

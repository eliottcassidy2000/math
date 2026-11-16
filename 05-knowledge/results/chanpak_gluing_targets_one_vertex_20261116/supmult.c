/* test H(A u B) >= H(A) H(B) for induced A,B with |A n B| <= 1 (all tournaments of order n, all such pairs with A u B = V) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXN 10
static int n,N;static uint32_t out_[MAXN],H[1<<MAXN],dp[1<<MAXN][MAXN];
static void allH(void){H[0]=1;for(int S=1;S<N;S++){uint32_t t=0;for(int v=0;v<n;v++){if(!(S>>v&1))continue;int R=S&~(1<<v);uint32_t c;if(!R)c=1;else{c=0;for(int u=0;u<n;u++)if((R>>u&1)&&(out_[u]>>v&1))c+=dp[R][u];}dp[S][v]=c;t+=c;}H[S]=t;}}
int main(int argc,char**argv){n=atoi(argv[1]);N=1<<n;char line[256];int V=N-1;long long tested=0,viol=0,eq=0,tested2=0,viol2=0;
 while(fgets(line,sizeof line,stdin)){int len=strlen(line);while(len&&line[len-1]=='\n')line[--len]=0;if(len!=n*(n-1)/2)continue;memset(out_,0,sizeof out_);int k=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,k++){if(line[k]=='1')out_[i]|=1u<<j;else out_[j]|=1u<<i;}
  allH();uint32_t HV=H[V];
  for(int v=0;v<n;v++){uint32_t rest=V&~(1u<<v);
   for(uint32_t Ap=rest;Ap;Ap=(Ap-1)&rest){uint32_t Bp=rest&~Ap;if(!Bp||Ap>Bp)continue;uint32_t A=Ap|(1u<<v),B=Bp|(1u<<v);
    unsigned long long pr=(unsigned long long)H[A]*H[B];tested++;if(HV<pr){viol++;if(viol<=3)printf("VIOL HV=%u HA=%u HB=%u T=%s A=%u B=%u\n",HV,H[A],H[B],line,A,B);}else if(HV==pr)eq++;}}
  /* two-vertex overlaps, for comparison: how often does H(A u B) < H(A)H(B)? */
  for(int v=0;v<n;v++)for(int w=v+1;w<n;w++){uint32_t rest=V&~(1u<<v)&~(1u<<w);
   for(uint32_t Ap=rest;Ap;Ap=(Ap-1)&rest){uint32_t Bp=rest&~Ap;if(!Bp||Ap>Bp)continue;uint32_t A=Ap|(1u<<v)|(1u<<w),B=Bp|(1u<<v)|(1u<<w);
    tested2++;if((unsigned long long)HV<(unsigned long long)H[A]*H[B])viol2++;}}}
 printf("n=%d one-vertex overlaps tested=%lld violations=%lld equality=%lld | two-vertex overlaps tested=%lld with H(AuB)<H(A)H(B): %lld\n",n,tested,viol,eq,tested2,viol2);return 0;}

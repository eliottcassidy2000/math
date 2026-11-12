/* Auditor: H(K[TT_a at x, TT_b at y]) by literal subset DP. stdin: n x y a b mask_0..mask_{n-1}; prints H. */
#include <stdio.h>
typedef long long ll;
static ll dp[1<<16][16]; static int A[16][16];
int main(void){ int n,x,y,a,b; while(scanf("%d %d %d %d %d",&n,&x,&y,&a,&b)==5){ int K[16]; for(int i=0;i<n;i++) scanf("%d",&K[i]);
  int map[16],cp[16],m=0; for(int v=0;v<n;v++){ int c=(v==x)?a:(v==y)?b:1; for(int t=0;t<c;t++){map[m]=v;cp[m]=t;m++;} }
  for(int i=0;i<m;i++) for(int j=0;j<m;j++){ A[i][j]=0; if(i==j) continue; if(map[i]==map[j]) A[i][j]=cp[i]<cp[j]; else A[i][j]=(K[map[i]]>>map[j])&1; }
  for(int S=0;S<(1<<m);S++) for(int v=0;v<m;v++) dp[S][v]=0; for(int v=0;v<m;v++) dp[1<<v][v]=1;
  for(int S=1;S<(1<<m);S++) for(int v=0;v<m;v++){ ll c=dp[S][v]; if(!c) continue; for(int w=0;w<m;w++) if(!((S>>w)&1)&&A[v][w]) dp[S|1<<w][w]+=c; }
  ll h=0; for(int v=0;v<m;v++) h+=dp[(1<<m)-1][v]; printf("%lld\n",h); } return 0; }

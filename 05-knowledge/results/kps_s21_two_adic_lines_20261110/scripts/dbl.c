/* dbl.c - doubling vertices. For each class K (order m) and vertex x: h_a = H(K[TT_a at x]), a = 1..4.
 * x is DOUBLING if h_a = 1 + (h_1 - 1) 2^(a-1) for a = 2,3,4.
 * output: code H x c2 doubling   (c2 = (h_2 - h_1)/2 = number of Hamiltonian walks through x twice) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MAXV 13
static int adj[MAXV][MAXV], K[MAXV][MAXV];
static u64 dp[1<<MAXV][MAXV];
static u64 hp(int n){
  int full=(1<<n)-1; unsigned out[MAXV];
  for (int v=0;v<n;v++){ out[v]=0; for (int u=0;u<n;u++) if (adj[v][u]) out[v]|=1u<<u; }
  for (int m=1;m<=full;m++) for (int v=0;v<n;v++) dp[m][v]=0;
  for (int v=0;v<n;v++) dp[1<<v][v]=1;
  for (int m=1;m<=full;m++) for (int v=0;v<n;v++){ u64 c=dp[m][v]; if(!c) continue;
    unsigned nx=out[v]&~m; while(nx){ int u=__builtin_ctz(nx); nx&=nx-1; dp[m|1u<<u][u]+=c; } }
  u64 s=0; for (int v=0;v<n;v++) s+=dp[full][v]; return s; }
static void decode(u64 code,int m){ int b=0; memset(K,0,sizeof K);
  for(int i=0;i<m;i++) for(int j=i+1;j<m;j++){ if(code>>b&1) K[i][j]=1; else K[j][i]=1; b++; } }
static int build(int m,int x,int a){ int orig[MAXV],idx[MAXV],n=0;
  for(int v=0;v<m;v++){ int rep=(v==x)?a:1; for(int i=0;i<rep;i++){orig[n]=v;idx[n]=i;n++;} }
  for(int p=0;p<n;p++) for(int q=0;q<n;q++){ if(p==q){adj[p][q]=0;continue;}
    adj[p][q]=(orig[p]==orig[q])?(idx[p]<idx[q]):K[orig[p]][orig[q]]; } return n; }
int main(int argc,char**argv){ int m=atoi(argv[1]); char line[4096];
  while(fgets(line,sizeof line,stdin)){ u64 code=strtoull(line,NULL,10); decode(code,m);
    for(int x=0;x<m;x++){ u64 h[5]; for(int a=1;a<=4;a++){ int n=build(m,x,a); h[a]=hp(n); }
      int d = (h[2]==2*h[1]-1) && (h[3]==4*h[1]-3) && (h[4]==8*h[1]-7);
      printf("%llu %llu %d %llu %d\n", code, h[1], x, (h[2]-h[1])/2, d); } }
  return 0; }

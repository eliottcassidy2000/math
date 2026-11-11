/* twinscan.c - slope rays of seeds. Input: unit-offset roots "m code v H d E U M" (seeds.c output).
 * For each root (L, v) and each vertex u != v replace u by a transitive module TT_a (a = 2, 3, 4) and recompute the
 * rooted state at v. If the offset stays 1 for a = 2, 3, 4, print  m code v u M1 M2 M3 M4  (M_a = slope with module size a).
 * Then M_a = ga + de 2^a when the family is exactly affine (checked: M3 - M2 = 2 (M2 - M1), M4 - M3 = 2 (M3 - M2)). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXV 13
static u64 en[1<<MAXV][MAXV], st[1<<MAXV][MAXV], hs[1<<MAXV];
static int A[MAXV][MAXV];
static void state(int n,int v,i64*H,i64*d,i64*E,i64*U){ int full=(1<<n)-1; unsigned out[MAXV],inn[MAXV];
  for(int i=0;i<n;i++){out[i]=0;inn[i]=0;for(int j=0;j<n;j++){if(A[i][j])out[i]|=1u<<j; if(A[j][i])inn[i]|=1u<<j;}}
  for(int s=0;s<=full;s++) for(int x=0;x<n;x++){en[s][x]=0;st[s][x]=0;}
  for(int x=0;x<n;x++){en[1<<x][x]=1;st[1<<x][x]=1;}
  for(int s=1;s<=full;s++) for(int x=0;x<n;x++){ u64 c=en[s][x]; if(c){unsigned nx=out[x]&~s; while(nx){int u=__builtin_ctz(nx);nx&=nx-1;en[s|1u<<u][u]+=c;}}
    c=st[s][x]; if(c){unsigned nx=inn[x]&~s; while(nx){int u=__builtin_ctz(nx);nx&=nx-1;st[s|1u<<u][u]+=c;}} }
  hs[0]=1; for(int s=1;s<=full;s++){u64 t=0; for(int x=0;x<n;x++) t+=en[s][x]; hs[s]=t;}
  *H=hs[full]; *d=hs[full^(1<<v)]; *E=en[full][v]; i64 t=0;
  for(int X=1;X<=full;X++) if((X>>v&1)&&st[X][v]) t+=(i64)st[X][v]*(i64)hs[full^X]; *U=t; }
int main(void){ char line[1<<12]; int K[MAXV][MAXV];
  while(fgets(line,sizeof line,stdin)){ int m,v; u64 code; i64 h0,d0,e0,u0,M0;
    if(sscanf(line,"%d %llu %d %lld %lld %lld %lld %lld",&m,&code,&v,&h0,&d0,&e0,&u0,&M0)!=8) continue;
    if(M0<=0||m>9) continue; int b=0; memset(K,0,sizeof K);
    for(int i=0;i<m;i++) for(int j=i+1;j<m;j++){ if(code>>b&1) K[i][j]=1; else K[j][i]=1; b++; }
    for(int u=0;u<m;u++){ if(u==v) continue; i64 Ms[5]; int ok=1; Ms[1]=M0;
      for(int a=2;a<=4&&ok;a++){ int orig[MAXV],idx[MAXV],n=0,vv=-1;
        for(int x=0;x<m;x++){ int rep=(x==u)?a:1; for(int i=0;i<rep;i++){ if(x==v) vv=n; orig[n]=x; idx[n]=i; n++; } }
        if(n>MAXV){ok=0;break;}
        for(int p=0;p<n;p++) for(int q=0;q<n;q++){ if(p==q){A[p][q]=0;continue;} A[p][q]=(orig[p]==orig[q])?(idx[p]<idx[q]):K[orig[p]][orig[q]]; }
        i64 H,d,E,U; state(n,vv,&H,&d,&E,&U); if(H-E+2*d-U!=1) ok=0; else Ms[a]=U-d; }
      if(ok) printf("%d %llu %d %d %lld %lld %lld %lld\n",m,code,v,u,Ms[1],Ms[2],Ms[3],Ms[4]); } }
  return 0; }

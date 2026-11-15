// For unit cuts: classify Hamiltonian paths by label-run shape (runs of X / X^c); test simple formulas for delta.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long i64;
static int n, out[16], path[16];
static i64 cntshape[64]; // shape code: sequence of run labels, start label + number of runs
static void dfs(int len,int used,int X){
  if(len==n){ int runs=1; for(int i=0;i+1<n;i++) if(((X>>path[i])&1)!=((X>>path[i+1])&1)) runs++;
    int first=(X>>path[0])&1; if(runs>15) runs=15; cntshape[first*16+runs]++; return; }
  int last=path[len-1]; int nb=out[last]&~used; while(nb){int u=__builtin_ctz(nb); nb&=nb-1; path[len]=u; dfs(len+1,used|1<<u,X);} }
static i64 dpE[1<<15][15], dpS[1<<15][15];
int main(int argc,char**argv){ n=atoi(argv[1]); char line[4096]; int N=1<<n;
  long long ok1=0, ok2=0, tot=0;
  while(fgets(line,sizeof line,stdin)){ if(line[0]=='>') continue;
    int k=0; for(int i=0;i<n;i++) out[i]=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(line[k]=='1') out[i]|=1<<j; else out[j]|=1<<i; k++; }
    memset(dpE,0,sizeof(i64)*N*15); memset(dpS,0,sizeof(i64)*N*15);
    int in[16]; for(int v=0;v<n;v++) in[v]=0; for(int v=0;v<n;v++) for(int u=0;u<n;u++) if((out[v]>>u)&1) in[u]|=1<<v;
    for(int v=0;v<n;v++){ dpE[1<<v][v]=1; dpS[1<<v][v]=1; }
    for(int S=1;S<N;S++) for(int v=0;v<n;v++){ i64 x=dpE[S][v]; if(x){int nb=out[v]&~S; while(nb){int u=__builtin_ctz(nb); nb&=nb-1; dpE[S|1<<u][u]+=x;}}
      i64 y=dpS[S][v]; if(y){int nb=in[v]&~S; while(nb){int u=__builtin_ctz(nb); nb&=nb-1; dpS[S|1<<u][u]+=y;}} }
    i64 H=0; for(int v=0;v<n;v++) H+=dpE[N-1][v];
    for(int X=1;X<N-1;X++){
      i64 tot2=0; for(int v=0;v<n;v++){ if((X>>v)&1) tot2+=dpS[N-1][v]; else tot2+=dpE[N-1][v]; }
      for(int S=1;S<N-1;S++){ int C=(N-1)^S; for(int u=0;u<n;u++) if(((S>>u)&1)&&!((X>>u)&1)&&dpE[S][u]) for(int w=0;w<n;w++) if(((C>>w)&1)&&((X>>w)&1)&&dpS[C][w]) tot2+=dpE[S][u]*dpS[C][w]; }
      i64 delta=tot2-2*H; if(delta!=1&&delta!=-1) continue;
      memset(cntshape,0,sizeof cntshape); for(int s=0;s<n;s++){ path[0]=s; dfs(1,1<<s,X); }
      // shapes: first label f (1 = X), runs r. [X][Xc]: f=1,r=2 ; [Xc][X]: f=0,r=2 ; [X]: impossible (X proper) ...
      i64 A=cntshape[16+2], Z=cntshape[0+2], XcXXc=cntshape[0+3], XXcX=cntshape[16+3], four=cntshape[0+4]+cntshape[16+4];
      tot++;
      if(delta==A-Z) ok1++;
      if(delta==A-Z+XXcX-XcXXc) ok2++;
      if(tot<=12) printf("n=%d X=%d delta=%lld H=%lld  [X][Xc]=%lld [Xc][X]=%lld [Xc][X][Xc]=%lld [X][Xc][X]=%lld 4runs=%lld\n",n,X,delta,H,A,Z,XcXXc,XXcX,four);
    }
  }
  printf("order %d: unit cuts %lld; delta == A-Z in %lld; delta == A-Z+[X][Xc][X]-[Xc][X][Xc] in %lld\n",n,tot,ok1,ok2);
  return 0;}

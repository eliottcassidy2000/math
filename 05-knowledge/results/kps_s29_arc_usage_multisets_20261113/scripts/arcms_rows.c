/* arcms.c: for each class code (classgen file, order n) print: code H zeros | sorted nonzero arc path-counts
   W(u->v) = sum over S containing u, not v, of F(S,u) * B(V\S, v)   (F: paths of T[S] ending at u, B: starting at v) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
#define MAXN 13
static u64 F[(1<<MAXN)*MAXN], Bk[(1<<MAXN)*MAXN];
static int cmp(const void*a,const void*b){ u64 x=*(const u64*)a,y=*(const u64*)b; return x<y?-1:x>y; }
int main(int argc,char**argv){
  char line[8192];
  while(fgets(line,sizeof line,stdin)){
    char*pp=line; int n=strtol(pp,&pp,10); if(n<=0) continue; unsigned N=1u<<n, FULL=N-1; unsigned out[MAXN]={0}; for(int i=0;i<n;i++) out[i]=strtoul(pp,&pp,10); u64 code=0;


    memset(F,0,sizeof(u64)*N*n); memset(Bk,0,sizeof(u64)*N*n);
    for(int v=0;v<n;v++){ F[(size_t)(1u<<v)*n+v]=1; Bk[(size_t)(1u<<v)*n+v]=1; }
    for(unsigned S=1;S<N;S++) for(int v=0;v<n;v++){
      u64 a=F[(size_t)S*n+v]; if(a) for(int w=0;w<n;w++) if(!(S>>w&1)&&(out[v]>>w&1)) F[(size_t)(S|1u<<w)*n+w]+=a;
      u64 b=Bk[(size_t)S*n+v]; if(b) for(int w=0;w<n;w++) if(!(S>>w&1)&&(out[w]>>v&1)) Bk[(size_t)(S|1u<<w)*n+w]+=b; }
    u64 H=0; for(int v=0;v<n;v++) H+=F[(size_t)FULL*n+v];
    u64 Wl[64]; int m=0, zeros=0;
    for(int u=0;u<n;u++) for(int v=0;v<n;v++) if(u!=v && (out[u]>>v&1)){
      u64 W=0;
      for(unsigned S=1;S<N;S++) if((S>>u&1) && !(S>>v&1)){ u64 a=F[(size_t)S*n+u]; if(!a) continue; unsigned R=FULL&~S; W+=a*Bk[(size_t)R*n+v]; }
      if(W==0) zeros++; else Wl[m++]=W; }
    qsort(Wl,m,sizeof(u64),cmp);
    printf("%llu %llu %d |",code,H,zeros); for(int i=m-1;i>=0;i--) printf(" %llu",Wl[i]); putchar(10);
  }
  return 0;
}

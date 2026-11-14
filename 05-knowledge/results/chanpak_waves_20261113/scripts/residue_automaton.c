// BFS over two-hub prefix vectors modulo m: size of the residue automaton.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define NS 36
static int idx(int o,int a,int b,int d){return ((o*3+a)*3+b)*2+d;}
static int T[4][NS][3], Tn[4][NS];
static void build(void){
  for(int let=0;let<4;let++){ int xi=let>>1, yi=let&1;
    for(int o=0;o<2;o++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){ int s=idx(o,a,b,d), k=0;
      int h1=o?yi:xi, h2=o?xi:yi;
      T[let][s][k++]=idx(o,h1?2:1,b,d);
      if(b||h1) T[let][s][k++]=idx(o,a,h2?2:1,d);
      if(d||h2) T[let][s][k++]=idx(o,a,b,1);
      Tn[let][s]=k; } } }
typedef struct { unsigned char v[NS]; } V;
static V *tab; static int *hnext, *hhead; static int nst=0, cap; static unsigned HS;
static unsigned hsh(const V*x){ unsigned h=2166136261u; for(int i=0;i<NS;i++){ h^=x->v[i]; h*=16777619u;} return h; }
static int find_or_add(const V*x){ unsigned h=hsh(x)&(HS-1); for(int i=hhead[h];i>=0;i=hnext[i]) if(!memcmp(tab[i].v,x->v,NS)) return i;
  if(nst>=cap) return -2; tab[nst]=*x; hnext[nst]=hhead[h]; hhead[h]=nst; return nst++; }
int main(int argc,char**argv){
  int m=atoi(argv[1]); cap=atoi(argv[2]); build(); HS=1; while(HS<2u*cap) HS<<=1;
  tab=malloc(sizeof(V)*cap); hnext=malloc(sizeof(int)*cap); hhead=malloc(sizeof(int)*HS); memset(hhead,-1,sizeof(int)*HS);
  V s0; memset(&s0,0,sizeof s0); s0.v[idx(0,0,0,0)]=1%m; s0.v[idx(1,0,0,0)]=1%m; find_or_add(&s0);
  int head=0, depth=0, levelend=1;
  while(head<nst){
    V x=tab[head++];
    for(int let=0;let<4;let++){ int w[NS]; memset(w,0,sizeof w);
      for(int s=0;s<NS;s++) if(x.v[s]) for(int t=0;t<Tn[let][s];t++) w[T[let][s][t]]+=x.v[s];
      V y; for(int s=0;s<NS;s++) y.v[s]=w[s]%m; int r=find_or_add(&y); if(r==-2){ printf("m=%d: more than %d states (depth %d)\n",m,cap,depth); return 0; } }
    if(head==levelend){ depth++; levelend=nst; }
  }
  printf("m=%d: %d reachable prefix vectors mod m (BFS depth %d)\n",m,nst,depth);
  return 0;
}

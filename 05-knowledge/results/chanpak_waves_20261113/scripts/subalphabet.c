// Coverage of two-hub words restricted to a sub-alphabet (bitmask of allowed letters 0..3, letter = 2x+y).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define NS 36
typedef unsigned long long u64;
static int idx(int o,int a,int b,int d){return ((o*3+a)*3+b)*2+d;}
static int T[4][NS][3], Tn[4][NS], fin[2][NS];
static void build(void){
  for(int let=0;let<4;let++){ int xi=let>>1, yi=let&1;
    for(int o=0;o<2;o++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){ int s=idx(o,a,b,d), k=0;
      int h1=o?yi:xi, h2=o?xi:yi;
      T[let][s][k++]=idx(o,h1?2:1,b,d);
      if(b||h1) T[let][s][k++]=idx(o,a,h2?2:1,d);
      if(d||h2) T[let][s][k++]=idx(o,a,b,1);
      Tn[let][s]=k; } }
  for(int e=0;e<2;e++) for(int o=0;o<2;o++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){
      int h12=o?!e:e; fin[e][idx(o,a,b,d)] = !(a==2||b==2||(b==0&&!h12)); } }
static int L, mask; static u64 LIM; static unsigned char *cov; static long long *mult;
static void dfs(int k,u64*v){
  for(int e=0;e<2;e++){ u64 H=0; for(int s=0;s<NS;s++) if(fin[e][s]) H+=v[s]; if(H<=LIM){ cov[H>>1]=1; mult[H>>1]++; } }
  if(k==L) return;
  for(int let=0;let<4;let++) if(mask>>let&1){ u64 w[NS]; memset(w,0,sizeof w);
    for(int s=0;s<NS;s++) if(v[s]) for(int t=0;t<Tn[let][s];t++) w[T[let][s][t]]+=v[s];
    dfs(k+1,w); }
}
int main(int argc,char**argv){
  mask=atoi(argv[1]); L=atoi(argv[2]); LIM=strtoull(argv[3],0,10); build();
  cov=calloc(LIM/2+2,1); mult=calloc(LIM/2+2,sizeof(long long));
  u64 v[NS]; memset(v,0,sizeof v); v[idx(0,0,0,0)]=1; v[idx(1,0,0,0)]=1; dfs(0,v);
  u64 X=23; while(X<=LIM && cov[X>>1]) X+=2;
  int shown=0; printf("alphabet mask %d (letters:",mask); for(int l=0;l<4;l++) if(mask>>l&1) printf(" %d%d",l>>1,l&1); printf("), L<=%d: first miss above 21 = %llu (X/3^L=%.3f); misses <= %llu:",L,X,(double)X/__builtin_powi(3.0,L),LIM);
  u64 cnt=0; for(u64 n=1;n<=LIM;n+=2) if(!cov[n>>1]){ cnt++; if(shown<25){printf(" %llu",n);shown++;} }
  printf(" ... total %llu\n",cnt);
  return 0;
}

// Branch-and-bound absolute-target decoder for two-hub words: H(x,y,e) = N exactly.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef unsigned __int128 u128;
#define NS 36
static int idx(int o,int a,int b,int d){return ((o*3+a)*3+b)*2+d;}
static int T[4][NS][3]; static int Tn[4][NS];      // transitions: letter -> list of successor states (multiplicity via repeats)
static int fin[2][NS];                              // finishing indicator per e
static u128 U[80][NS], Lo[80][NS];
static int L; static u128 N; static long long nodes; static int word[80]; static int foundE;
static void build(void){
  for(int let=0;let<4;let++){ int xi=let>>1, yi=let&1;
    for(int o=0;o<2;o++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){ int s=idx(o,a,b,d), k=0;
      int h1=o?yi:xi, h2=o?xi:yi;
      T[let][s][k++]=idx(o,h1?2:1,b,d);
      if(b||h1) T[let][s][k++]=idx(o,a,h2?2:1,d);
      if(d||h2) T[let][s][k++]=idx(o,a,b,1);
      Tn[let][s]=k; } }
  for(int e=0;e<2;e++) for(int o=0;o<2;o++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){
      int h12=o?!e:e; fin[e][idx(o,a,b,d)] = !(a==2||b==2||(b==0&&!h12)); }
}
static void bounds(int maxr){
  for(int s=0;s<NS;s++){ U[0][s]= (fin[0][s]||fin[1][s]); Lo[0][s]=(fin[0][s]&&fin[1][s]); }
  for(int r=1;r<=maxr;r++) for(int s=0;s<NS;s++){ u128 mx=0, mn=~(u128)0;
      for(int let=0;let<4;let++){ u128 su=0, sl=0; for(int k=0;k<Tn[let][s];k++){ su+=U[r-1][T[let][s][k]]; sl+=Lo[r-1][T[let][s][k]]; }
        if(su>mx) mx=su; if(sl<mn) mn=sl; }
      U[r][s]=mx; Lo[r][s]=mn; }
}
static int dfs(int k, u128 *v){
  nodes++;
  int r=L-k; u128 hi=0, lo=0;
  for(int s=0;s<NS;s++) if(v[s]){ hi+=v[s]*U[r][s]; lo+=v[s]*Lo[r][s]; }
  if(N>hi || N<lo) return 0;
  if(r==0){ for(int e=0;e<2;e++){ u128 H=0; for(int s=0;s<NS;s++) if(fin[e][s]) H+=v[s]; if(H==N){foundE=e; return 1;} } return 0; }
  u128 ws[4][NS]; u128 his[4]; int ord[4];
  for(int let=0;let<4;let++){ memset(ws[let],0,sizeof ws[let]);
    for(int s=0;s<NS;s++) if(v[s]) for(int t=0;t<Tn[let][s];t++) ws[let][T[let][s][t]]+=v[s];
    u128 h=0; for(int s=0;s<NS;s++) if(ws[let][s]) h+=ws[let][s]*U[r-1][s]; his[let]=h; ord[let]=let; }
  for(int i=0;i<4;i++) for(int j=i+1;j<4;j++) if(his[ord[j]]<his[ord[i]]){int t=ord[i];ord[i]=ord[j];ord[j]=t;}
  for(int i=0;i<4;i++){ int let=ord[i]; if(his[let]<N) continue; word[k]=let; if(dfs(k+1,ws[let])) return 1; }
  return 0;
}
static unsigned long long rng=0x9E3779B97F4A7C15ULL;
static unsigned long long xr(void){ rng^=rng<<13; rng^=rng>>7; rng^=rng<<17; return rng; }
int main(int argc,char**argv){
  if(argc==2 && !strcmp(argv[1],"stdin")){ build(); bounds(79); unsigned long long n; while(scanf("%llu",&n)==1){ N=n; int found=0; for(L=1;L<=70&&!found;L++){ u128 v[NS]; memset(v,0,sizeof v); v[idx(0,0,0,0)]=1; v[idx(1,0,0,0)]=1; found=dfs(0,v);} printf("%llu %d\n",n,found?L-1:-1);} return 0; }
  int trials=atoi(argv[1]); double scale=atof(argv[2]); build(); bounds(79);
  int ok=0; long long totn=0; double maxlen=0;
  for(int t=0;t<trials;t++){
    unsigned long long n=(unsigned long long)(scale*(0.5+(xr()>>11)*(1.0/9007199254740992.0))); n|=1; N=n;
    int found=0; nodes=0;
    for(L=1; L<=70 && !found; L++){ u128 v[NS]; memset(v,0,sizeof v); v[idx(0,0,0,0)]=1; v[idx(1,0,0,0)]=1; found=dfs(0,v); }
    if(found){ ok++; totn+=nodes; if(L-1>maxlen) maxlen=L-1; }
    else printf("  FAILED %llu\n",n);
    if(t==0 && found){ printf("  example N=%llu L=%d e=%d word=",n,L-1,foundE); for(int i=0;i<L-1;i++) printf("%d",word[i]); printf("  (letters: 2x+y)\n"); }
  }
  printf("scale %.0e: solved %d/%d, mean search nodes %.0f, max word length %.0f\n",scale,ok,trials,(double)totn/(ok?ok:1),maxlen);
  return 0;
}

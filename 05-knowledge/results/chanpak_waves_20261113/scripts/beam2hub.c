// Coarse-to-fine beam solver: find two-hub words with H(x,y,e) = N.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef unsigned __int128 u128;
#define NS 36
static int idx(int o,int a,int b,int d){return ((o*3+a)*3+b)*2+d;}
static int T[4][NS][3], Tn[4][NS], fin[2][NS];
static u128 U[90][NS];
static void build(void){
  for(int let=0;let<4;let++){ int xi=let&1, yi=let>>1;
    for(int o=0;o<2;o++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){ int s=idx(o,a,b,d),k=0; int h1=o?yi:xi,h2=o?xi:yi;
      T[let][s][k++]=idx(o,h1?2:1,b,d); if(b||h1) T[let][s][k++]=idx(o,a,h2?2:1,d); if(d||h2) T[let][s][k++]=idx(o,a,b,1); Tn[let][s]=k; } }
  for(int e=0;e<2;e++)for(int o=0;o<2;o++)for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int d=0;d<2;d++){ int h12=o?!e:e; fin[e][idx(o,a,b,d)]=!(a==2||b==2||(b==0&&!h12)); }
  for(int s=0;s<NS;s++) U[0][s]=(fin[0][s]||fin[1][s]);
  for(int r=1;r<90;r++) for(int s=0;s<NS;s++){ u128 mx=0; for(int l=0;l<4;l++){ u128 su=0; for(int k=0;k<Tn[l][s];k++) su+=U[r-1][T[l][s][k]]; if(su>mx) mx=su; } U[r][s]=mx; }
}
static void stepv(const u128*v,u128*w,int let){ memset(w,0,sizeof(u128)*NS); for(int s=0;s<NS;s++) if(v[s]) for(int k=0;k<Tn[let][s];k++) w[T[let][s][k]]+=v[s]; }
static u128 finish(const u128*v,int e){ u128 h=0; for(int s=0;s<NS;s++) if(fin[e][s]) h+=v[s]; return h; }
static unsigned long long rs=0x9E3779B97F4A7C15ULL; static unsigned long long xr(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }
typedef struct { u128 v[NS]; unsigned char w[90]; double score; } Node;
static int cmp(const void*a,const void*b){ double x=((Node*)a)->score, y=((Node*)b)->score; return x<y?-1:x>y; }
int solve(u128 N,int L,int B,int S,unsigned char*out,int*oute){
  Node *beam=malloc(sizeof(Node)*B), *cand=malloc(sizeof(Node)*4*B); int nb=1;
  memset(&beam[0],0,sizeof(Node)); beam[0].v[idx(0,0,0,0)]=1; beam[0].v[idx(1,0,0,0)]=1;
  double lN=log((double)N);
  for(int k=0;k<L;k++){
    int nc=0; int r=L-k-1;
    for(int i=0;i<nb;i++) for(int let=0;let<4;let++){
      Node *c=&cand[nc]; stepv(beam[i].v,c->v,let); memcpy(c->w,beam[i].w,k); c->w[k]=let;
      u128 hi=0; for(int s=0;s<NS;s++) if(c->v[s]) hi+=c->v[s]*U[r][s]; if(hi<N) continue;
      // sample random completions
      int above=0, below=0; double dsum=0;
      for(int t=0;t<S;t++){ u128 v[NS],w2[NS]; memcpy(v,c->v,sizeof v); for(int q=0;q<r;q++){ stepv(v,w2,xr()&3); memcpy(v,w2,sizeof v);} u128 h=finish(v,xr()&1);
        if(h>N) above++; else if(h<N) below++; else { memcpy(out,c->w,k+1); /* random suffix hit: reconstruct impossible cheaply */ }
        dsum+= fabs(log((double)h+1.0)-lN); }
      c->score = fabs((double)above-(double)below)/(S) + 0.1*dsum/S + 1e-9*(xr()%1000);
      nc++; }
    if(!nc){ free(beam); free(cand); return 0; }
    qsort(cand,nc,sizeof(Node),cmp); nb= nc<B?nc:B; memcpy(beam,cand,sizeof(Node)*nb);
  }
  for(int i=0;i<nb;i++) for(int e=0;e<2;e++) if(finish(beam[i].v,e)==N){ memcpy(out,beam[i].w,L); *oute=e; free(beam); free(cand); return 1; }
  free(beam); free(cand); return 0;
}
int main(int argc,char**argv){
  int trials=atoi(argv[1]); double scale=atof(argv[2]); int B=atoi(argv[3]), S=atoi(argv[4]); build();
  int ok=0; unsigned char w[90]; int e;
  for(int t=0;t<trials;t++){
    unsigned long long n=(unsigned long long)(scale*(0.5+(xr()>>11)*(1.0/9007199254740992.0))); n|=1;
    int L0=(int)ceil(log(2.0*(double)n)/log(3.0))+1, found=0;
    for(int L=L0; L<=L0+3 && !found; L++) found=solve(n,L,B,S,w,&e);
    if(found){ ok++; if(ok==1){ printf("  example N=%llu e=%d word=",n,e); int L=0; while(L<90 && 0) L++; } }
  }
  printf("scale %.0e beam %d samples %d: solved %d/%d\n",scale,B,S,ok,trials); return 0; }

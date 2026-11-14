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
static int L; static u128 N; static long long nodes; static int firstfail=-1, maxk=0; static int descending=1, fpath[80], fplen=0; static int word[80]; static int foundE;
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
  nodes++; if(k>maxk) maxk=k;
  int r=L-k; u128 hi=0, lo=0;
  for(int s=0;s<NS;s++) if(v[s]){ hi+=v[s]*U[r][s]; lo+=v[s]*Lo[r][s]; }
  if(N>hi || N<lo){ descending=0; return 0; }
  if(r==0){ for(int e=0;e<2;e++){ u128 H=0; for(int s=0;s<NS;s++) if(fin[e][s]) H+=v[s]; if(H==N){foundE=e; return 1;} } return 0; }
  u128 ws[4][NS]; u128 his[4]; int ord[4]; double sc[4];
  for(int let=0;let<4;let++){ memset(ws[let],0,sizeof ws[let]);
    for(int s=0;s<NS;s++) if(v[s]) for(int t=0;t<Tn[let][s];t++) ws[let][T[let][s][t]]+=v[s];
    u128 h=0,l=0; for(int s=0;s<NS;s++) if(ws[let][s]){ h+=ws[let][s]*U[r-1][s]; l+=ws[let][s]*Lo[r-1][s]; } his[let]=h; ord[let]=let; double lc=0.5*(log((double)(l>0?l:1))+log((double)(h>0?h:1))); sc[let]=fabs(log((double)N)-lc); if(h<N||l>N) sc[let]=1e300; }
  for(int i=0;i<4;i++) for(int j=i+1;j<4;j++) if(sc[ord[j]]<sc[ord[i]]){int t=ord[i];ord[i]=ord[j];ord[j]=t;}
  for(int i=0;i<4;i++){ int let=ord[i]; if(his[let]<N) continue; word[k]=let; if(descending){fpath[k]=let; fplen=k+1;} if(dfs(k+1,ws[let])) return 1; }
  descending=0; return 0;
}
int main(int argc,char**argv){
  build(); bounds(79);
  unsigned long long N0=strtoull(argv[1],0,10); int cnt=atoi(argv[2]);
  long long hist_len[80]={0}; int descent=0, tot=0; double nodesum=0; long long worst=0; unsigned long long worstN=0;
  int extraL=0; long long tail[80]={0};
  for(int i=0;i<cnt;i++){ unsigned long long n=N0+2ULL*i; if(n==7||n==21) continue; N=n; int found=0; long long nf=0; int Lmin=0;
    for(L=1;L<=70&&!found;L++){ u128 v[NS]; memset(v,0,sizeof v); v[idx(0,0,0,0)]=1; v[idx(1,0,0,0)]=1; nodes=0; maxk=0; descending=1; fplen=0; found=dfs(0,v); nf=nodes; if(maxk>=L && !Lmin) Lmin=L; }
    tot++; int Lf=L-1; { int cp=0; while(cp<Lf && cp<fplen && fpath[cp]==word[cp]) cp++; tail[Lf-cp]++; } hist_len[Lf]++; if(nf==Lf+1) descent++; nodesum+=nf; if(nf>worst){worst=nf;worstN=n;}
  }
  printf("N in [%llu,%llu): %d targets, first-descent success %d (%.1f%%), mean final-length nodes %.1f, worst %lld at N=%llu; lengths:",N0,N0+2ULL*cnt,tot,descent,100.0*descent/tot,nodesum/tot,worst,worstN);
  for(int l=0;l<80;l++) if(hist_len[l]) printf(" L%d:%lld",l,hist_len[l]); printf("\n  tail (letters after greedy head):"); for(int t=0;t<80;t++) if(tail[t]) printf(" %d:%lld",t,tail[t]); printf("\n");
  return 0;
}

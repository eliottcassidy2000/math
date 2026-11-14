// Moment-steered (and optionally residue-filtered) absolute-target decoder for two-hub words.
// Children are ordered by the normal density of N under the EXACT mean and variance of H over the child's subtree
// (mean = v.mu[r], second moment = v^T Q[r] v, with Q[r] = (1/4) sum_l M_l Q[r-1] M_l^T).
// Optional bottom filter: exact attainable residues of H mod MOD for the last R0 letters (suffix response table).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef unsigned __int128 u128; typedef unsigned long long u64;
#define NS 36
#define RMAX 80
static int idx(int o,int a,int b,int d){return ((o*3+a)*3+b)*2+d;}
static int T[4][NS][3], Tn[4][NS], fin[2][NS];
static u128 U[RMAX][NS], Lo[RMAX][NS];
static double MU[RMAX][NS], Q[RMAX][NS][NS];
static int L; static u128 N; static long long nodes; static int word[RMAX], foundE; static int descending, fplen, fpath[RMAX];
static int R0=0, MOD=1; static unsigned char *RES[16]; static long NRES[16];   // distinct suffix response vectors mod MOD
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
static void tables(void){
  for(int s=0;s<NS;s++){ U[0][s]=(fin[0][s]||fin[1][s]); Lo[0][s]=(fin[0][s]&&fin[1][s]); MU[0][s]=0.5*(fin[0][s]+fin[1][s]);
    for(int t=0;t<NS;t++) Q[0][s][t]=0.5*(fin[0][s]*fin[0][t]+fin[1][s]*fin[1][t]); }
  for(int r=1;r<RMAX;r++){
    for(int s=0;s<NS;s++){ u128 mx=0, mn=~(u128)0; double m=0;
      for(int let=0;let<4;let++){ u128 su=0, sl=0; double sm=0; for(int k=0;k<Tn[let][s];k++){ su+=U[r-1][T[let][s][k]]; sl+=Lo[r-1][T[let][s][k]]; sm+=MU[r-1][T[let][s][k]]; }
        if(su>mx) mx=su; if(sl<mn) mn=sl; m+=0.25*sm; }
      U[r][s]=mx; Lo[r][s]=mn; MU[r][s]=m; }
    for(int s=0;s<NS;s++) for(int t=0;t<NS;t++){ double q=0;
      for(int let=0;let<4;let++) for(int i=0;i<Tn[let][s];i++) for(int j=0;j<Tn[let][t];j++) q+=Q[r-1][T[let][s][i]][T[let][t][j]];
      Q[r][s][t]=0.25*q; } }
}
static int cmprow(const void*a,const void*b){ return memcmp(a,b,NS); }
static void restable(void){   // RES[r]: distinct response vectors mod MOD for suffix length r (with both e)
  for(int r=0;r<=R0;r++){
    long cap=2L<<(2*r); unsigned char *buf=malloc((size_t)cap*NS); long n=0;
    if(r==0){ for(int e=0;e<2;e++){ for(int s=0;s<NS;s++) buf[n*NS+s]=fin[e][s]%MOD; n++; } }
    else { for(long i=0;i<NRES[r-1];i++) for(int let=0;let<4;let++){ unsigned char *x=RES[r-1]+i*NS, *y=buf+n*NS;
        for(int s=0;s<NS;s++){ int q=0; for(int k=0;k<Tn[let][s];k++) q+=x[T[let][s][k]]; y[s]=q%MOD; } n++; } }
    qsort(buf,n,NS,cmprow); long u=0; for(long i=0;i<n;i++){ if(u==0 || memcmp(buf+(u-1)*NS,buf+i*NS,NS)) { if(u!=i) memcpy(buf+u*NS,buf+i*NS,NS); u++; } }
    RES[r]=buf; NRES[r]=u; fprintf(stderr,"  residue table r=%d: %ld distinct of %ld\n",r,u,n); }
}
static int resok(const u128*v,int r){ // exists response x with v.x == N mod MOD
  unsigned vm[NS]; int nz[NS],k=0; for(int s=0;s<NS;s++){ vm[s]=(unsigned)(v[s]%MOD); if(vm[s]) nz[k++]=s; }
  unsigned tgt=(unsigned)(N%MOD);
  for(long i=0;i<NRES[r];i++){ unsigned char*x=RES[r]+i*NS; unsigned acc=0; for(int j=0;j<k;j++) acc+=vm[nz[j]]*x[nz[j]]; if(acc%MOD==tgt) return 1; }
  return 0;
}
static int mode=1; // 0 tightest-first, 1 moment-steered
static int dfs(int k,u128*v){
  nodes++;
  int r=L-k; u128 hi=0, lo=0;
  for(int s=0;s<NS;s++) if(v[s]){ hi+=v[s]*U[r][s]; lo+=v[s]*Lo[r][s]; }
  if(N>hi || N<lo){ descending=0; return 0; }
  if(R0 && r<=R0 && r>0 && !resok(v,r)){ descending=0; return 0; }
  if(r==0){ for(int e=0;e<2;e++){ u128 H=0; for(int s=0;s<NS;s++) if(fin[e][s]) H+=v[s]; if(H==N){foundE=e; return 1;} } descending=0; return 0; }
  u128 ws[4][NS]; double sc[4]; int ord[4];
  for(int let=0;let<4;let++){ memset(ws[let],0,sizeof ws[let]); ord[let]=let;
    for(int s=0;s<NS;s++) if(v[s]) for(int t=0;t<Tn[let][s];t++) ws[let][T[let][s][t]]+=v[s];
    double m=0, q=0; double vd[NS]; for(int s=0;s<NS;s++) vd[s]=(double)ws[let][s];
    if(mode==1 && r-1>R0){ for(int s=0;s<NS;s++) if(vd[s]!=0){ m+=vd[s]*MU[r-1][s]; for(int t=0;t<NS;t++) if(vd[t]!=0) q+=vd[s]*vd[t]*Q[r-1][s][t]; }
      double var=q-m*m; if(var<1) var=1; double z=((double)N-m); sc[let]=-(z*z)/(2*var)-0.5*log(var); }
    else { u128 h=0; (void)q; for(int s=0;s<NS;s++) if(ws[let][s]) h+=ws[let][s]*U[r-1][s]; sc[let]= h<N? -1e300 : -(double)h; } }
  for(int i=0;i<4;i++) for(int j=i+1;j<4;j++) if(sc[ord[j]]>sc[ord[i]]){int t=ord[i];ord[i]=ord[j];ord[j]=t;}
  for(int i=0;i<4;i++){ int let=ord[i]; word[k]=let; if(descending){fpath[k]=let; fplen=k+1;} if(dfs(k+1,ws[let])) return 1; }
  return 0;
}
int main(int argc,char**argv){
  if(argc<4){ fprintf(stderr,"usage: %s mode(0 tight,1 moment) N0 count [R0 MOD]\n",argv[0]); return 1; }
  mode=atoi(argv[1]); u64 N0=strtoull(argv[2],0,10); int cnt=atoi(argv[3]); if(argc>5){ R0=atoi(argv[4]); MOD=atoi(argv[5]); }
  build(); tables(); if(R0) restable();
  int descent=0, tot=0; double nodesum=0; long long worst=0; u64 worstN=0; long tail[RMAX]={0}; long lens[RMAX]={0};
  for(int i=0;i<cnt;i++){ u64 n=N0+2ULL*i; if(n==7||n==21) continue; N=n; int found=0; long long nf=0;
    for(L=1;L<=70&&!found;L++){ u128 v[NS]; memset(v,0,sizeof v); v[idx(0,0,0,0)]=1; v[idx(1,0,0,0)]=1; nodes=0; descending=1; fplen=0; found=dfs(0,v); nf=nodes; }
    tot++; int Lf=L-1; lens[Lf]++; if(nf==Lf+1) descent++; nodesum+=nf; if(nf>worst){worst=nf;worstN=n;}
    int cp=0; while(cp<Lf && cp<fplen && fpath[cp]==word[cp]) cp++; tail[Lf-cp]++; }
  printf("mode %d R0=%d MOD=%d  N in [%llu,+%d): %d targets, first-descent %d (%.1f%%), mean nodes %.1f, worst %lld at %llu; tails:",mode,R0,MOD,N0,2*cnt,tot,descent,100.0*descent/tot,nodesum/tot,worst,worstN);
  for(int t=0;t<RMAX;t++) if(tail[t]) printf(" %d:%ld",t,tail[t]); printf("\n");
  return 0;
}

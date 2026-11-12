/* Independent audit program (auditor's own code; does not reuse the packet's programs).
   usage: audit_G n      (3 <= n <= 6: full checks with subtournament data;  n = 7: light mode)
   Encoding: arc bit for pair (i<j) is 1 iff i->j.  T^X flips every arc with exactly one end in X.
   For every labelled tournament T we compute, from scratch, by a start-indexed DP over subsets:
     f_s[W][v] = number of Hamiltonian paths of T[W] starting at s and ending at v.
   From it: s_T, t_T, e_T, H, hc, and for every subset W: e_{T[W]} and hc(T[W]).
   Checks (full mode):
     G1  Lambda_T(u,w) := e_{T^{w}}(u) + e_T(u); e_{T^X} = e_T + Lambda 1_X for ALL X; symmetric; zero diag;
         zero row sums; Lambda_{T^X} = Lambda_T for all X.
     L1  Lemma 1: e_T(u) = 1 + sum_{v in N+(u)} e_{T-u}(v).
     S1  recursion (*):  Lambda_T(u,w) = e_{T-u}(w) + sum_{v in N+(u)} Lambda_{T-u}(v,w)   (n>=3)
     S2  recursion (**): Lambda_T(u,w) = 1 + <a+b, e_Q> + a^T Lambda_Q b                   (n>=4)
     C1  hc(T^X) = hc(T) + sum_{a in X}(1+e_T(a)) + sum_{pairs} Lambda(a,b)  mod 2, for ALL X
     C2  for each l: c_l(T^X) mod 2 = quadratic with polar form sum_{|W|=l, W>=ab} Lambda_{T[W]}(a,b): check
         c_l(T^X) = c_l(T) + sum_a D_a + sum_{pairs} polar  mod 2 for ALL X;
         (H(T^X)-H(T))/2 = sum lambda_a + sum M(a,b) mod 2 for ALL X; triple sums of H = 0 mod 4.
     C3  M_T(a,b) = sum_{v in R, exactly one of a->v,b->v} E_R(v) + iota_a^T M_R iota_b
     C4  n=6: M(a,b) = 1 + d_a + d_b
     histograms: per class #(H=1 mod 4); distinct Lambda over classes; 2-adic valuations of c_U(H), c_U(hc), all U.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n, NP, full;
static int pid[8][8];
static uint32_t flipm[256];          /* arc mask flipped by switching X */
static int out_[8];

static void setT(uint32_t T){
  for(int i=0;i<n;i++) out_[i]=0;
  for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if((T>>pid[i][j])&1) out_[i]|=1<<j; else out_[j]|=1<<i; }
}

/* per tournament tables */
static uint8_t *Ee;      /* e_T bitmask */
static uint16_t *Hh;     /* H */
static uint16_t *HCc;    /* hc */
static uint8_t *ESUB;    /* e_{T[W]} : index T*(1<<n)+W  (full mode) */
static uint8_t *HCSUB;   /* hc(T[W]) mod 2 */
static uint8_t *CL; /* bit l = parity of c_l(T) */

static uint32_t f[8][128][8];

static void compute_full(uint32_t T){
  setT(T);
  for(int s=0;s<n;s++){
    memset(f[s],0,sizeof(f[s]));
    f[s][1<<s][s]=1;
    for(int W=1;W<=full;W++){ if(!((W>>s)&1)) continue;
      for(int v=0;v<n;v++){ uint32_t c=f[s][W][v]; if(!c) continue;
        int nx=out_[v]&~W; for(int w=0;w<n;w++) if((nx>>w)&1) f[s][W|(1<<w)][w]+=c; } }
  }
  /* subtournament data */
  for(int W=1;W<=full;W++){
    uint8_t ev=0; int k=__builtin_popcount(W);
    for(int v=0;v<n;v++) if((W>>v)&1){
      uint32_t st=0,en=0;
      for(int x=0;x<n;x++){ st+=f[v][W][x]; en+=f[x][W][v]; }
      if((st+en)&1) ev|=1<<v;
    }
    ESUB[(size_t)T*(full+1)+W]=ev;
    uint32_t hc=0;
    if(k>=3){ int m=__builtin_ctz(W); for(int v=0;v<n;v++) if(v!=m && ((W>>v)&1) && ((out_[v]>>m)&1)) hc+=f[m][W][v]; }
    HCSUB[(size_t)T*(full+1)+W]=hc&1;
    if(W==full){ Ee[T]=ev; HCc[T]=hc; }
  }
  uint32_t H=0; for(int s=0;s<n;s++) for(int v=0;v<n;v++) H+=f[s][full][v];
  Hh[T]=H;
}

/* light mode (n=7): end-DP and converse end-DP and cycle DP */
static uint32_t g[128][8];
static void compute_light(uint32_t T){
  setT(T); uint32_t tend[8], sst[8], H=0; uint8_t ev=0;
  memset(g,0,sizeof(g)); for(int v=0;v<n;v++) g[1<<v][v]=1;
  for(int W=1;W<=full;W++) for(int v=0;v<n;v++){ uint32_t c=g[W][v]; if(!c) continue; int nx=out_[v]&~W; while(nx){int w=__builtin_ctz(nx); nx&=nx-1; g[W|(1<<w)][w]+=c;} }
  for(int v=0;v<n;v++){ tend[v]=g[full][v]; H+=tend[v]; }
  memset(g,0,sizeof(g)); for(int v=0;v<n;v++) g[1<<v][v]=1;
  for(int W=1;W<=full;W++) for(int v=0;v<n;v++){ uint32_t c=g[W][v]; if(!c) continue; int nx=(~out_[v])&full&~(1<<v)&~W; while(nx){int w=__builtin_ctz(nx); nx&=nx-1; g[W|(1<<w)][w]+=c;} }
  for(int v=0;v<n;v++){ sst[v]=g[full][v]; if((sst[v]+tend[v])&1) ev|=1<<v; }
  memset(g,0,sizeof(g)); g[1][0]=1;
  for(int W=1;W<=full;W+=2) for(int v=0;v<n;v++){ uint32_t c=g[W][v]; if(!c) continue; int nx=out_[v]&~W; while(nx){int w=__builtin_ctz(nx); nx&=nx-1; g[W|(1<<w)][w]+=c;} }
  uint32_t hc=0; for(int v=1;v<n;v++) if((out_[v]>>0)&1) hc+=g[full][v];
  Ee[T]=ev; Hh[T]=H; HCc[T]=hc;
}

static inline int arc(uint32_t T,int i,int j){ /* 1 iff i->j */ return i<j ? (T>>pid[i][j])&1 : !((T>>pid[j][i])&1); }

/* Lambda of subtournament T[W] at (u,w), via the affine law inside T[W]: (T^{w})[W] = T[W]^{w} */
static int lam_sub(uint32_t T,int W,int u,int w){
  uint32_t T2=T^flipm[1<<w];
  return ((ESUB[(size_t)T*(full+1)+W]^ESUB[(size_t)T2*(full+1)+W])>>u)&1;
}

static int cmp(const void*a,const void*b){ uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b; return x<y?-1:x>y; }
static int v2(long long x){ if(x==0) return 99; int k=0; while(!(x&1)){x>>=1;k++;} return k; }

int main(int argc,char**argv){
  n=atoi(argv[1]); full=(1<<n)-1; NP=0;
  for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ pid[i][j]=pid[j][i]=NP++; }
  for(int X=0;X<=full;X++){ uint32_t m=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) if(((X>>i)&1)!=((X>>j)&1)) m|=1u<<pid[i][j]; flipm[X]=m; }
  uint32_t NT=1u<<NP; int fullmode = (n<=6);
  Ee=malloc(NT); Hh=malloc(NT*2); HCc=malloc(NT*2);
  if(fullmode){ ESUB=malloc((size_t)NT*(full+1)); HCSUB=malloc((size_t)NT*(full+1)); CL=calloc(NT,1); }
  for(uint32_t T=0;T<NT;T++){ if(fullmode) compute_full(T); else compute_light(T); }
  if(fullmode) for(uint32_t T=0;T<NT;T++) for(int W=1;W<=full;W++){ int k=__builtin_popcount(W); if(k>=3 && HCSUB[(size_t)T*(full+1)+W]) CL[T]^=1<<k; }

  long long bad[16]; memset(bad,0,sizeof(bad)); long long nchk[16]; memset(nchk,0,sizeof(nchk));
  for(uint32_t T=0;T<NT;T++){
    int L[8][8];
    for(int u=0;u<n;u++) for(int w=0;w<n;w++) L[u][w]= (u==w)?0: (((Ee[T]^Ee[T^flipm[1<<w]])>>u)&1);
    /* Lemma 2 */
    for(int u=0;u<n;u++){ nchk[0]++; if(((Ee[T]^Ee[T^flipm[1<<u]])>>u)&1) bad[0]++; }
    /* G1 */
    for(int u=0;u<n;u++){ int rs=0; for(int w=0;w<n;w++){ rs^=L[u][w]; if(L[u][w]!=L[w][u]) bad[1]++; } if(rs) bad[1]++; nchk[1]++; }
    for(int X=0;X<=full;X++){
      uint8_t pred=Ee[T]; for(int u=0;u<n;u++){ int s=0; for(int w=0;w<n;w++) if((X>>w)&1) s^=L[u][w]; if(s) pred^=1<<u; }
      nchk[2]++; if(pred!=Ee[T^flipm[X]]) bad[2]++;
      /* class invariance (light mode: single switchings only, which generate the class) */
      uint32_t T2=T^flipm[X];
      if(fullmode || __builtin_popcount(X)==1) for(int w=0;w<n;w++){ uint8_t col=Ee[T2]^Ee[T2^flipm[1<<w]]; uint8_t c0=0; for(int u=0;u<n;u++) if(L[u][w]) c0|=1<<u; nchk[3]++; if(col!=c0) bad[3]++; }
      /* C1 */
      if(n>=3){ int p=HCc[T]&1; for(int a=0;a<n;a++) if((X>>a)&1){ p^=1^((Ee[T]>>a)&1); for(int b=a+1;b<n;b++) if((X>>b)&1) p^=L[a][b]; }
        nchk[4]++; if(p!=(HCc[T2]&1)) bad[4]++; }
    }
    /* triple sums of H mod 4 */
    for(int a=0;a<n;a++) for(int b=a+1;b<n;b++) for(int c=b+1;c<n;c++){ unsigned s=0; for(int m=0;m<8;m++){ int X=((m&1)<<a)|(((m>>1)&1)<<b)|(((m>>2)&1)<<c); s+=Hh[T^flipm[X]]; } nchk[5]++; if(s&3) bad[5]++; }
    if(!fullmode) continue;
    /* L1: Lemma 1 */
    for(int u=0;u<n;u++){ int W=full&~(1<<u); int s=1; for(int v=0;v<n;v++) if(v!=u && arc(T,u,v)) s^=(ESUB[(size_t)T*(full+1)+W]>>v)&1; nchk[6]++; if(s!=((Ee[T]>>u)&1)) bad[6]++; }
    /* S1: (*) */
    for(int u=0;u<n;u++) for(int w=0;w<n;w++) if(u!=w){ int W=full&~(1<<u); int s=(ESUB[(size_t)T*(full+1)+W]>>w)&1;
        for(int v=0;v<n;v++) if(v!=u && v!=w && arc(T,u,v)) s^=lam_sub(T,W,v,w);
        nchk[7]++; if(s!=L[u][w]) bad[7]++; }
    /* S2: (**) */
    if(n>=4) for(int u=0;u<n;u++) for(int w=0;w<n;w++) if(u!=w){ int Q=full&~(1<<u)&~(1<<w); int s=1;
        uint8_t eQ=ESUB[(size_t)T*(full+1)+Q];
        for(int v=0;v<n;v++) if((Q>>v)&1){ if(arc(T,u,v)) s^=(eQ>>v)&1; if(arc(T,w,v)) s^=(eQ>>v)&1; }
        for(int x=0;x<n;x++) if(((Q>>x)&1)&&arc(T,u,x)) for(int y=0;y<n;y++) if(((Q>>y)&1)&&arc(T,w,y)&&x!=y) s^=lam_sub(T,Q,x,y);
        nchk[8]++; if(s!=L[u][w]) bad[8]++; }
    /* C2: c_l parity quadratic with the stated polar form; H mod 4 formula */
    {
      int M[8][8]; memset(M,0,sizeof(M)); int lam[8]; memset(lam,0,sizeof(lam));
      for(int W=1;W<=full;W++){ int k=__builtin_popcount(W); if(k<3) continue;
        for(int a=0;a<n;a++) if((W>>a)&1){ if(k&1) lam[a]^=1^((ESUB[(size_t)T*(full+1)+W]>>a)&1);
          for(int b=0;b<n;b++) if(b!=a && ((W>>b)&1) && (k&1)) M[a][b]^=lam_sub(T,W,a,b); } }
      for(int X=0;X<=full;X++){ uint32_t T2=T^flipm[X]; int p=0;
        for(int a=0;a<n;a++) if((X>>a)&1){ p^=lam[a]; for(int b=a+1;b<n;b++) if((X>>b)&1) p^=M[a][b]; }
        int q=(((int)Hh[T2]-(int)Hh[T])/2)&1; nchk[9]++; if(p!=q) bad[9]++; }
      for(int l=3;l<=n;l++){
        int P[8][8]; memset(P,0,sizeof(P));
        for(int W=1;W<=full;W++) if(__builtin_popcount(W)==l) for(int a=0;a<n;a++) if((W>>a)&1) for(int b=a+1;b<n;b++) if((W>>b)&1) P[a][b]^=lam_sub(T,W,a,b);
        int D[8]; for(int a=0;a<n;a++) D[a]=((CL[T]^CL[T^flipm[1<<a]])>>l)&1;
        for(int X=0;X<=full;X++){ uint32_t T2=T^flipm[X]; int pred=(CL[T]>>l)&1;
          for(int a=0;a<n;a++) if((X>>a)&1){ pred^=D[a]; for(int b=a+1;b<n;b++) if((X>>b)&1) pred^=P[a][b]; }
          nchk[10]++; if(pred!=((CL[T2]>>l)&1)) bad[10]++; }
      }
      /* C3 */
      if(n>=4) for(int a=0;a<n;a++) for(int b=0;b<n;b++) if(a!=b){ int R=full&~(1<<a)&~(1<<b); int s=0;
          for(int v=0;v<n;v++) if(((R>>v)&1) && (arc(T,a,v)^arc(T,b,v))){ int E=0; for(int W=1;W<=full;W++) if((W&R)==W && ((W>>v)&1) && (__builtin_popcount(W)&1) && __builtin_popcount(W)>=3) E^=(ESUB[(size_t)T*(full+1)+W]>>v)&1; s^=E; }
          for(int x=0;x<n;x++) if(((R>>x)&1)&&arc(T,a,x)) for(int y=0;y<n;y++) if(((R>>y)&1)&&arc(T,b,y)&&x!=y){ int mr=0; for(int W=1;W<=full;W++) if((W&R)==W && ((W>>x)&1)&&((W>>y)&1)&&(__builtin_popcount(W)&1)) mr^=lam_sub(T,W,x,y); s^=mr; }
          nchk[11]++; if(s!=M[a][b]) bad[11]++; }
      /* C4 */
      if(n==6){ int d[8]; for(int a=0;a<n;a++){ d[a]=0; for(int b=0;b<n;b++) if(b!=a && arc(T,a,b)) d[a]++; }
        for(int a=0;a<n;a++) for(int b=0;b<n;b++) if(a!=b){ nchk[12]++; if(M[a][b]!=((1+d[a]+d[b])&1)) bad[12]++; } }
    }
  }
  const char *nm[]={"Lemma2","G1 sym/diag/rowsum","G affine law all X","Lambda class-invariant","C1 hc formula all X","H triple sums mod 4","Lemma1","(*)","(**)","C2 H mod 4 formula all X","C2 c_l quadratic w/ polar form","C3 recursion","C4 order-6 closed form"};
  printf("n=%d tournaments=%u mode=%s\n",n,NT,fullmode?"full":"light");
  for(int i=0;i<13;i++) if(nchk[i]) printf("  %-30s checks=%lld violations=%lld\n",nm[i],nchk[i],bad[i]);

  /* class histograms: classes <-> tournaments in which vertex n-1 is a source */
  { long long hist[200]; memset(hist,0,sizeof(hist)); long long ncl=0; long long imb[400]; memset(imb,0,sizeof(imb));
    uint64_t *lams=malloc(sizeof(uint64_t)*(1u<<NP)); long long nl=0; long long lzero=0;
    for(uint32_t T=0;T<NT;T++){ int src=1; for(int v=0;v<n-1;v++) if(!arc(T,n-1,v)) src=0; if(!src) continue; ncl++;
      int c1=0; for(int X=0;X<(1<<(n-1));X++){ if((Hh[T^flipm[X]]&3)==1) c1++; }
      hist[c1]++; imb[c1-((1<<(n-1))-c1)+200]++;
      uint64_t key=0; int z=1; for(int u=0;u<n;u++) for(int w=u+1;w<n;w++){ int l=((Ee[T]^Ee[T^flipm[1<<w]])>>u)&1; if(l){ key|=1ull<<(u*8+w); z=0;} }
      lams[nl++]=key; if(z) lzero++; }
    printf("classes=%lld; histogram #(H=1 mod 4) per class:",ncl); for(int i=0;i<200;i++) if(hist[i]) printf(" %d:%lld",i,hist[i]); printf("\n");
    printf("imbalance #(H=1)-#(H=3) per class:"); for(int i=0;i<400;i++) if(imb[i]) printf(" %+d:%lld",i-200,imb[i]); printf("\n");
    /* distinct Lambda */
    qsort(lams,nl,sizeof(uint64_t),cmp); long long nd=0; for(long long i=0;i<nl;i++) if(i==0||lams[i]!=lams[i-1]) nd++;
    printf("distinct Lambda over classes: %lld ; classes with Lambda=0: %lld\n",nd,lzero); free(lams); }

  /* 2-adic: all T, all U subset V (base T), c_U(H) and c_U(hc) */
  { int minH[9], minhc[9]; for(int k=0;k<9;k++){minH[k]=99;minhc[k]=99;} long long violH=0,violhc=0;
    long long cH[128], cC[128];
    for(uint32_t T=0;T<NT;T++){
      for(int S=0;S<=full;S++){ cH[S]=Hh[T^flipm[S]]; cC[S]=HCc[T^flipm[S]]; }
      for(int i=0;i<n;i++) for(int S=0;S<=full;S++) if((S>>i)&1){ cH[S]-=cH[S^(1<<i)]; cC[S]-=cC[S^(1<<i)]; }
      for(int U=1;U<=full;U++){ int k=__builtin_popcount(U); int a=v2(cH[U]), b=v2(cC[U]);
        if(a<minH[k]) minH[k]=a; if(b<minhc[k]) minhc[k]=b;
        if(a<(k+1)/2) violH++; if(b<(k+1)/2-1) violhc++; }
    }
    printf("2-adic (all T, all nonempty U subset V): violations H %lld, hc %lld; min v2 by |U| H:",violH,violhc);
    for(int k=1;k<=n;k++) printf(" %d:%d",k,minH[k]); printf("  hc:"); for(int k=1;k<=n;k++) printf(" %d:%d",k,minhc[k]); printf("\n");
  }
  return 0;
}

/* aud_engine.c -- independent audit engine for THM-5126 / HYP-9658 (packet kps_s21_two_adic_lines_20261110).
 * Written by the auditor; it shares no code with the packet's scripts. It reads only DATA files of the packet
 * (line_seeds.txt, unit_offset_roots_orders_2_to_9.txt) in the "seeds" and "roots" modes, to compare.
 *
 * build:  gcc -O2 -o aud_engine aud_engine.c
 * modes:
 *   aud_engine census OUTDIR     isomorphism classes of orders 1..9 (own canonical form), spectrum, unit-offset root
 *                                census (orders 2..9, every class, every root), unit sets over every base of order <= 9
 *                                with EVERY out-set O (all sizes; roots of order <= 10), doubling vertices (orders 3..8),
 *                                two-module families (orders 3..8). Writes OUTDIR/summary.txt and data files.
 *   aud_engine seeds FILE        verify every line of line_seeds.txt (states, offset, literal ray points)
 *   aud_engine roots FILE        verify every line of unit_offset_roots_orders_2_to_9.txt; print canonical rows
 *   aud_engine labelled8         spectrum of ALL labelled tournaments of orders 1..8 (no canonical forms)
 *   aud_engine mean              extension-mean law (C) by literal enumeration
 *   aud_engine oneapex LMAX      one-apex counts by literal DP vs g-recursion and trajectory formula (D)
 *
 * Tournament code: for i<j in the order (0,1),(0,2),...,(0,m-1),(1,2),..., bit b of code is 1 iff i -> j.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64; typedef long long i64;
#define MAXN 10
#define BIGN 14
static void die(const char*m){ fprintf(stderr,"FATAL: %s\n",m); fflush(stderr); exit(2); }
#define CHECK(c,msg) do{ if(!(c)) die(msg); }while(0)

static void decode(u64 code,int n,unsigned*out){ int b=0; for(int i=0;i<n;i++) out[i]=0;
  for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(code>>b&1) out[i]|=1u<<j; else out[j]|=1u<<i; b++; } }
static u64 encode(const unsigned*out,int n){ u64 c=0; int b=0;
  for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(out[i]>>j&1) c|=1ULL<<b; b++; } return c; }

/* ------------------------------------------------------------------ canonical form (own design)
 * initial colour of v: (score, sorted scores inside the out-neighbourhood); then colour refinement by
 * (colour, number of out-neighbours of each colour); canonical code = minimum code over all labellings that list
 * the colour cells in order (every permutation inside every cell is tried, no pruning). The number of labellings
 * attaining the minimum is |Aut|. */
static int cf_n; static const unsigned*cf_out; static int cf_col[MAXN], cf_poscell[MAXN], cf_inv[MAXN], cf_bestinv[MAXN], cf_off[MAXN];
static u64 cf_best; static long cf_aut; static unsigned cf_used;
static void cf_rec(int p,u64 partial){
  if(p==cf_n){ if(partial<cf_best){ cf_best=partial; cf_aut=1; memcpy(cf_bestinv,cf_inv,sizeof(int)*cf_n); } else if(partial==cf_best) cf_aut++; return; }
  int c=cf_poscell[p];
  for(int v=0;v<cf_n;v++) if(cf_col[v]==c && !(cf_used>>v&1)){
    u64 add=0; for(int i=0;i<p;i++) if(cf_out[cf_inv[i]]>>v&1) add|=1ULL<<(cf_off[i]+(p-i-1));
    cf_used|=1u<<v; cf_inv[p]=v; cf_rec(p+1,partial|add); cf_used&=~(1u<<v); }
}
static int sigcmp(const int*a,const int*b,int len){ for(int i=0;i<len;i++){ if(a[i]<b[i]) return -1; if(a[i]>b[i]) return 1; } return 0; }
static u64 canon(const unsigned*out,int n,long*aut,int*bestinv){
  cf_n=n; cf_out=out; u64 key[MAXN]; int col[MAXN]; int k;
  for(int v=0;v<n;v++){ int s=__builtin_popcount(out[v]); int deg[MAXN],kk=0;
    for(int u=0;u<n;u++) if(out[v]>>u&1) deg[kk++]=__builtin_popcount(out[u]&out[v]);
    for(int i=1;i<kk;i++){ int x=deg[i],j=i-1; while(j>=0&&deg[j]>x){deg[j+1]=deg[j];j--;} deg[j+1]=x; }
    u64 K=s; for(int i=0;i<kk;i++) K=K*11+deg[i]+1; key[v]=K; }
  { k=0; for(int v=0;v<n;v++){ int r=0; for(int u=0;u<n;u++){ if(key[u]<key[v]){ int dup=0; for(int w=0;w<u;w++) if(key[w]==key[u]){dup=1;break;} if(!dup) r++; } } col[v]=r; if(r+1>k) k=r+1; } }
  for(;;){ int sig[MAXN][MAXN+1]; int ord[MAXN], nc[MAXN];
    for(int v=0;v<n;v++){ sig[v][0]=col[v]; for(int c=0;c<k;c++) sig[v][c+1]=0; for(int u=0;u<n;u++) if(out[v]>>u&1) sig[v][col[u]+1]++; ord[v]=v; }
    for(int i=1;i<n;i++){ int x=ord[i],j=i-1; while(j>=0&&sigcmp(sig[ord[j]],sig[x],k+1)>0){ord[j+1]=ord[j];j--;} ord[j+1]=x; }
    int c=0; nc[ord[0]]=0; for(int i=1;i<n;i++){ if(sigcmp(sig[ord[i]],sig[ord[i-1]],k+1)!=0) c++; nc[ord[i]]=c; }
    int k2=c+1; memcpy(col,nc,sizeof(int)*n); if(k2==k) break; k=k2; }
  memcpy(cf_col,col,sizeof(int)*n);
  { int p=0; for(int c=0;c<k;c++) for(int v=0;v<n;v++) if(col[v]==c) cf_poscell[p++]=c; }
  for(int i=0;i<n;i++) cf_off[i]=i*n-i*(i+1)/2;
  cf_best=~0ULL; cf_aut=0; cf_used=0; cf_rec(0,0);
  if(aut)*aut=cf_aut; if(bestinv) memcpy(bestinv,cf_bestinv,sizeof(int)*n); return cf_best;
}

/* ------------------------------------------------------------------ hash set of u64 */
typedef struct { u64*key; size_t cap,cnt; } hset;
static void hs_init(hset*h,size_t cap){ h->cap=cap; h->cnt=0; h->key=(u64*)malloc(cap*sizeof(u64)); CHECK(h->key,"malloc"); for(size_t i=0;i<cap;i++) h->key[i]=~0ULL; }
static int hs_add(hset*h,u64 k){ size_t i=(size_t)((k*0x9E3779B97F4A7C15ULL)>>17)&(h->cap-1);
  while(h->key[i]!=~0ULL){ if(h->key[i]==k) return 0; i=(i+1)&(h->cap-1); } h->key[i]=k; h->cnt++; CHECK(h->cnt*2<h->cap,"hash full"); return 1; }
static int hs_has(const hset*h,u64 k){ size_t i=(size_t)((k*0x9E3779B97F4A7C15ULL)>>17)&(h->cap-1);
  while(h->key[i]!=~0ULL){ if(h->key[i]==k) return 1; i=(i+1)&(h->cap-1); } return 0; }

/* ------------------------------------------------------------------ path tables, n <= 10 */
static u64 en[1<<MAXN][MAXN], st[1<<MAXN][MAXN], hs_[1<<MAXN];
static void tables(const unsigned*out,int n){ int full=(1<<n)-1; unsigned inn[MAXN];
  for(int v=0;v<n;v++) inn[v]=0; for(int v=0;v<n;v++) for(int u=0;u<n;u++) if(out[v]>>u&1) inn[u]|=1u<<v;
  for(int s=0;s<=full;s++) for(int v=0;v<n;v++){ en[s][v]=0; st[s][v]=0; }
  for(int v=0;v<n;v++){ en[1<<v][v]=1; st[1<<v][v]=1; }
  for(int s=1;s<=full;s++){ unsigned m=s; while(m){ int v=__builtin_ctz(m); m&=m-1;
      u64 c=en[s][v]; if(c){ unsigned nx=out[v]&~s&full; while(nx){ int u=__builtin_ctz(nx); nx&=nx-1; en[s|1u<<u][u]+=c; } }
      c=st[s][v]; if(c){ unsigned nx=inn[v]&~s&full; while(nx){ int u=__builtin_ctz(nx); nx&=nx-1; st[s|1u<<u][u]+=c; } } } }
  hs_[0]=1; for(int s=1;s<=full;s++){ u64 t=0; for(int v=0;v<n;v++) t+=en[s][v]; hs_[s]=t; }
}
/* rooted quantities of (L,v) from the tables of L:
 * H, d = H(L-v), E = # Ham paths ending at v, U = # ordered covers (P1 possibly empty, P2 starting at v),
 * dp = delta' = # crossed covers BY DEFINITION (P1 nonempty, P2 - v nonempty, v beats the last vertex of P1). */
static void rootq(const unsigned*out,int n,int v,i64*H,i64*d,i64*E,i64*U,i64*dp){ int full=(1<<n)-1;
  *H=(i64)hs_[full]; *d=(i64)hs_[full^(1<<v)]; *E=(i64)en[full][v]; i64 u=0,x=0;
  for(int X=1;X<=full;X++) if(X>>v&1){ i64 a=(i64)st[X][v]; if(!a) continue; int Y=full^X; u+=a*(i64)hs_[Y];
    if(X!=(1<<v) && Y){ i64 b=0; unsigned m=out[v]&Y; while(m){ int w=__builtin_ctz(m); m&=m-1; b+=(i64)en[Y][w]; } x+=a*b; } }
  *U=u; *dp=x; }
/* cut data of a base B from the tables of B: UB[y], W[e][f] = Q(e,f)+Q(f,e) */
static void cutdata(int m,i64*UB,i64 W[MAXN][MAXN]){ int full=(1<<m)-1;
  for(int y=0;y<m;y++){ i64 t=0; for(int X=1;X<=full;X++) if((X>>y&1)&&st[X][y]) t+=(i64)st[X][y]*(i64)hs_[full^X]; UB[y]=t; }
  for(int a=0;a<m;a++) for(int c=0;c<m;c++) W[a][c]=0;
  for(int X=1;X<full;X++){ int Y=full^X; unsigned ma=X; while(ma){ int a=__builtin_ctz(ma); ma&=ma-1; i64 ea=(i64)en[X][a]; if(!ea) continue;
      unsigned mc=Y; while(mc){ int c=__builtin_ctz(mc); mc&=mc-1; i64 sc=(i64)st[Y][c]; if(!sc) continue; W[a][c]+=ea*sc; W[c][a]+=ea*sc; } } } }

/* ------------------------------------------------------------------ generic Hamiltonian path count, n <= BIGN */
static u64 bigdp[(size_t)1<<BIGN][BIGN];
static u64 hp(int n,const unsigned*out){ int full=(1<<n)-1;
  for(int s=1;s<=full;s++) for(int v=0;v<n;v++) bigdp[s][v]=0;
  for(int v=0;v<n;v++) bigdp[1<<v][v]=1;
  for(int s=1;s<=full;s++){ unsigned m=s; while(m){ int v=__builtin_ctz(m); m&=m-1; u64 c=bigdp[s][v]; if(!c) continue;
      unsigned nx=out[v]&~s&full; while(nx){ int u=__builtin_ctz(nx); nx&=nx-1; bigdp[s|1u<<u][u]+=c; } } }
  u64 t=0; for(int v=0;v<n;v++) t+=bigdp[full][v]; return t; }
/* L o y : new chain z_1 -> ... -> z_t, every vertex of L - v beats every z_j, v -> z_j iff y_j = 1 */
static int build_tail(int n,const unsigned*out,int v,int t,const int*y,unsigned*B){
  for(int i=0;i<n;i++) B[i]=out[i];
  for(int a=0;a<t;a++){ int z=n+a; B[z]=0; for(int b=a+1;b<t;b++) B[z]|=1u<<(n+b);
    for(int u=0;u<n;u++) if(u!=v) B[u]|=1u<<z;
    if(y[a]) B[v]|=1u<<z; else B[z]|=1u<<v; }
  return n+t; }
/* K[TT_a at x, TT_b at y] (y<0: only x) */
static int build_blow(int m,const unsigned*K,int x,int a,int y,int b,unsigned*B){ int orig[BIGN],idx[BIGN],N=0;
  for(int v=0;v<m;v++){ int rep=(v==x)?a:(v==y)?b:1; for(int i=0;i<rep;i++){ CHECK(N<BIGN,"blow too big"); orig[N]=v; idx[N]=i; N++; } }
  for(int p=0;p<N;p++){ B[p]=0; for(int q=0;q<N;q++){ if(p==q) continue;
      if(orig[p]==orig[q]){ if(idx[p]<idx[q]) B[p]|=1u<<q; } else if(K[orig[p]]>>orig[q]&1) B[p]|=1u<<q; } }
  return N; }

/* ------------------------------------------------------------------ Hamiltonian walks with repeated vertices (n <= 8)
 * c[j] = number of walks of K visiting x exactly j times and every other vertex exactly once (j = 1..5) */
static u64 wf[256][6][8];
static void walk1(int n,const unsigned*out,int x,u64*c){ int full=(1<<n)-1, oth=full&~(1<<x);
  for(int s=0;s<=full;s++) memset(wf[s],0,sizeof wf[s]);
  for(int v=0;v<n;v++){ if(v==x) wf[0][1][x]=1; else wf[1<<v][0][v]=1; }
  for(int s=0;s<=full;s++){ if(s>>x&1) continue; for(int k=0;k<=5;k++) for(int last=0;last<n;last++){ u64 f=wf[s][k][last]; if(!f) continue;
      unsigned nx=out[last]; if((nx>>x&1)&&k<5) wf[s][k+1][x]+=f;
      unsigned m=nx&oth&~s; while(m){ int u=__builtin_ctz(m); m&=m-1; wf[s|1u<<u][k][u]+=f; } } }
  for(int j=0;j<=5;j++){ u64 t=0; for(int last=0;last<n;last++) t+=wf[oth][j][last]; c[j]=t; } }
/* c[i][j] = number of walks visiting x exactly i times, y exactly j times, every other vertex once (i,j = 0..3) */
static u64 wg[256][4][4][8];
static void walk2(int n,const unsigned*out,int x,int y,u64 c[4][4]){ int full=(1<<n)-1, oth=full&~(1<<x)&~(1<<y);
  for(int s=0;s<=full;s++) if(!((s>>x&1)||(s>>y&1))) memset(wg[s],0,sizeof wg[s]);
  for(int v=0;v<n;v++){ if(v==x) wg[0][1][0][x]=1; else if(v==y) wg[0][0][1][y]=1; else wg[1<<v][0][0][v]=1; }
  for(int s=0;s<=full;s++){ if((s>>x&1)||(s>>y&1)) continue;
    for(int t=0;t<=6;t++) for(int kx=0;kx<=3;kx++){ int ky=t-kx; if(ky<0||ky>3) continue;
      for(int last=0;last<n;last++){ u64 f=wg[s][kx][ky][last]; if(!f) continue; unsigned nx=out[last];
        if((nx>>x&1)&&kx<3) wg[s][kx+1][ky][x]+=f;
        if((nx>>y&1)&&ky<3) wg[s][kx][ky+1][y]+=f;
        unsigned m=nx&oth&~s; while(m){ int u=__builtin_ctz(m); m&=m-1; wg[s|1u<<u][kx][ky][u]+=f; } } } }
  for(int i=0;i<=3;i++) for(int j=0;j<=3;j++){ u64 t=0; for(int last=0;last<n;last++) t+=wg[oth][i][j][last]; c[i][j]=t; } }
static const i64 FS[5][5]={{0,0,0,0,0},{0,1,0,0,0},{0,1,2,0,0},{0,1,6,6,0},{0,1,14,36,24}}; /* FS[a][j] = j! S(a,j) = P_j(TT_a) */

/* strong components: comp size of each vertex */
static void compsize(int n,const unsigned*out,int*sz){ unsigned reach[MAXN];
  for(int v=0;v<n;v++) reach[v]=out[v]|(1u<<v);
  for(int k=0;k<n;k++) for(int i=0;i<n;i++) if(reach[i]>>k&1) reach[i]|=reach[k];
  for(int v=0;v<n;v++){ int s=0; for(int u=0;u<n;u++) if((reach[v]>>u&1)&&(reach[u]>>v&1)) s++; sz[v]=s; } }

/* ------------------------------------------------------------------ classes */
static u64*reps[MAXN+1]; static long*auts[MAXN+1]; static int nrep[MAXN+1];
static void gen_classes(int NMAX){
  reps[1]=(u64*)malloc(8); reps[1][0]=0; nrep[1]=1; auts[1]=(long*)malloc(sizeof(long)); auts[1][0]=1;
  for(int n=1;n<NMAX;n++){ int m=n+1; size_t maxc=(m>=9)?200000:8000; hset h; hs_init(&h,(m>=9)?((size_t)1<<20):((size_t)1<<15));
    u64*list=(u64*)malloc(sizeof(u64)*maxc); long*al=(long*)malloc(sizeof(long)*maxc); int cnt=0;
    for(int r=0;r<nrep[n];r++){ unsigned out[MAXN]; decode(reps[n][r],n,out);
      for(unsigned O=0;O<(1u<<n);O++){ unsigned o2[MAXN]; for(int v=0;v<n;v++) o2[v]=out[v]|((O>>v&1)?0u:(1u<<n)); o2[n]=O;
        long aut; u64 c=canon(o2,m,&aut,NULL); if(hs_add(&h,c)){ CHECK((size_t)cnt<maxc,"too many classes"); list[cnt]=c; al[cnt]=aut; cnt++; } } }
    reps[m]=list; auts[m]=al; nrep[m]=cnt; free(h.key); }
}
static int cmp4(const void*a,const void*b){ const i64*x=(const i64*)a,*y=(const i64*)b; for(int i=0;i<4;i++){ if(x[i]<y[i]) return -1; if(x[i]>y[i]) return 1; } return 0; }
static int cmpu64(const void*a,const void*b){ u64 x=*(const u64*)a,y=*(const u64*)b; return x<y?-1:x>y?1:0; }

#define VMAX 4096
static int minord[VMAX]; static u64 witc[VMAX];

static int census(const char*outdir){ char fn[1024];
  snprintf(fn,sizeof fn,"%s/summary.txt",outdir); FILE*fs=fopen(fn,"w"); CHECK(fs,"open summary");
  snprintf(fn,sizeof fn,"%s/roots_unit_offset.txt",outdir); FILE*fr=fopen(fn,"w"); CHECK(fr,"open roots");
  snprintf(fn,sizeof fn,"%s/unitsets_all_sizes.txt",outdir); FILE*fu=fopen(fn,"w"); CHECK(fu,"open unitsets");
  snprintf(fn,sizeof fn,"%s/doubling_vertices.txt",outdir); FILE*fd=fopen(fn,"w"); CHECK(fd,"open dbl");
  snprintf(fn,sizeof fn,"%s/values_witness.txt",outdir); FILE*fv=fopen(fn,"w"); CHECK(fv,"open values");
  snprintf(fn,sizeof fn,"%s/exact_tuples_order_le8.txt",outdir); FILE*ft=fopen(fn,"w"); CHECK(ft,"open tuples");
  gen_classes(9);
  static const i64 A568[10]={1,1,1,2,4,12,56,456,6880,191536};
  for(int m=1;m<=9;m++){ i64 fact=1; for(int i=2;i<=m;i++) fact*=i; i64 mass=0; for(int r=0;r<nrep[m];r++) mass+=fact/auts[m][r];
    fprintf(fs,"classes_%d=%d A000568=%lld mass_sum=%lld two_pow=%lld %s\n",m,nrep[m],A568[m],mass,(i64)1<<(m*(m-1)/2),
      (nrep[m]==A568[m] && mass==((i64)1<<(m*(m-1)/2)))?"OK":"MISMATCH"); }
  fflush(fs);
  for(int i=0;i<VMAX;i++) minord[i]=0;
  i64 nident=0; /* number of roots on which the identities were checked */
  i64 *tup=(i64*)malloc(sizeof(i64)*4*420000); size_t ntup=0; size_t ntup8start=0;
  i64 ntriples=0, nexact=0, nadditive=0, nEviol=0, nlaw1=0, nlaw2=0, ndblflagdiff=0;
  for(int m=1;m<=9;m++){
    i64 nunit=0, nunit_dE=0, nunit_M0=0, outdeg_unit[MAXN]={0}; i64 us_by_size[MAXN+1]={0}; i64 ndbl=0, ndbl_pos=0; i64 maxH=0; i64 nval=0;
    static unsigned char seen[VMAX]; memset(seen,0,sizeof seen);
    if(m==8) ntup8start=ntup;
    for(int r=0;r<nrep[m];r++){ u64 code=reps[m][r]; unsigned out[MAXN]; decode(code,m,out); int full=(1<<m)-1;
      tables(out,m); i64 H=(i64)hs_[full]; CHECK(H<VMAX,"H too large"); if(!seen[H]){ seen[H]=1; nval++; } if(H>maxH) maxH=H;
      if(!minord[H]){ minord[H]=m; witc[H]=code; }
      /* ---- direct root census: every vertex as root */
      if(m>=2) for(int v=0;v<m;v++){ i64 h,d,E,U,dp; rootq(out,m,v,&h,&d,&E,&U,&dp); i64 c0=h-E+2*d-U, M=U-d; nident++;
        CHECK(c0==d-dp,"identity c0 = d - delta' fails"); CHECK(M==h-E+dp,"identity M = H - E + delta' fails"); CHECK((c0&1)!=0,"c0 even");
        if(c0==1){ nunit++; if(d==E) nunit_dE++; if(M==0) nunit_M0++; outdeg_unit[__builtin_popcount(out[v])]++;
          fprintf(fr,"%d %llu %d %lld %lld %lld %lld %lld %d\n",m,code,v,h,d,E,U,M,__builtin_popcount(out[v])); } }
      /* ---- unit sets over the base B = this class: every O (all sizes); roots of order m+1 */
      { i64 UB[MAXN], W[MAXN][MAXN]; cutdata(m,UB,W);
        for(int a=0;a<m;a++) for(int c=a+1;c<m;c++){ CHECK((W[a][c]&1)==0,"W odd"); CHECK(W[a][c]>=2*(i64)hs_[full&~(1<<a)&~(1<<c)],"W < 2H(B-e-f)"); }
        static i64 wO[1<<MAXN], mO[1<<MAXN]; wO[0]=0; mO[0]=0;
        for(int O=1;O<=full;O++){ int e=__builtin_ctz(O); int rest=O&(O-1); i64 w=wO[rest]; unsigned mm=rest; while(mm){ int f=__builtin_ctz(mm); mm&=mm-1; w+=W[e][f]; } wO[O]=w; mO[O]=mO[rest]+UB[e]; }
        int nfound=0; int foundO[1<<MAXN];
        for(int O=0;O<=full;O++) if(wO[O]==H-1){ us_by_size[__builtin_popcount(O)]++; foundO[nfound++]=O;
            fprintf(fu,"%d %llu %d %lld %lld %d\n",m+1,code,O,H,mO[O],__builtin_popcount(O)); }
        /* direct verification of (1') on L = B + v, N+(v) = O: all O for m <= 7, the unit sets found for m = 8, 9 */
        static i64 wsave[1<<MAXN], msave[1<<MAXN]; memcpy(wsave,wO,sizeof(i64)*(full+1)); memcpy(msave,mO,sizeof(i64)*(full+1));
        int nO=(m<=7)?(full+1):nfound;
        for(int q=0;q<nO;q++){ int O=(m<=7)?q:foundO[q]; unsigned o2[MAXN]; for(int v=0;v<m;v++) o2[v]=out[v]|((O>>v&1)?0u:(1u<<m)); o2[m]=(unsigned)O;
          tables(o2,m+1); i64 h,d,E,U,dp; rootq(o2,m+1,m,&h,&d,&E,&U,&dp);
          CHECK(d==H,"d != H(B)"); CHECK(h-E+2*d-U==H-wsave[O],"(1') offset fails"); CHECK(U-d==msave[O],"(1') slope fails"); CHECK(dp==wsave[O],"delta' != W(O)"); nlaw1++; }
        tables(out,m); /* restore tables of the class itself */
      }
      /* ---- doubling vertices, orders 3..8 */
      if(m>=3 && m<=8) for(int x=0;x<m;x++){ u64 c[6]; walk1(m,out,x,c); CHECK((i64)c[1]==H,"c1 != H"); CHECK(c[0]==hs_[full^(1<<x)],"c_0 != H(K - x)");
        i64 hh[5]; for(int a=2;a<=4;a++){ unsigned B[BIGN]; int N=build_blow(m,out,x,a,-1,0,B); hh[a]=(i64)hp(N,B);
          i64 pred=0; for(int j=1;j<=a;j++) pred+=(i64)c[j]*FS[a][j]; CHECK(pred==hh[a],"module law fails"); nlaw2++; }
        int dblA=(2*(i64)c[2]==H-1)&&c[3]==0&&c[4]==0&&c[5]==0; int dblB=(hh[2]==2*H-1)&&(hh[3]==4*H-3)&&(hh[4]==8*H-7);
        if(dblA!=dblB) ndblflagdiff++;
        if(dblA){ ndbl++; if(c[2]>0){ ndbl_pos++; fprintf(fd,"%d %llu %d %lld %llu\n",m,code,x,H,c[2]); } } }
      /* ---- two-module families, orders 3..8 */
      if(m>=3 && m<=8){ int sz[MAXN]; compsize(m,out,sz);
        for(int x=0;x<m;x++) for(int y=0;y<m;y++){ if(x==y) continue; u64 c[4][4]; walk2(m,out,x,y,c); ntriples++;
          i64 h[4][4]; for(int a=1;a<=3;a++) for(int b=1;b<=3;b++){ i64 t=0; for(int i=1;i<=a;i++) for(int j=1;j<=b;j++) t+=(i64)c[i][j]*FS[a][i]*FS[b][j]; h[a][b]=t; }
          CHECK(h[1][1]==H,"c11 != H");
          i64 D=h[2][2]-h[2][1]-h[1][2]+h[1][1]; CHECK(D%4==0,"D mod 4"); i64 de=D/4, be=(h[2][1]-h[1][1])/2-2*de, ga=(h[1][2]-h[1][1])/2-2*de, al=h[1][1]-2*be-2*ga-4*de;
          int exact=1; for(int a=1;a<=3&&exact;a++) for(int b=1;b<=3;b++){ i64 pred=al+be*((i64)1<<a)+ga*((i64)1<<b)+de*((i64)1<<(a+b)); if(pred!=h[a][b]){ exact=0; break; } }
          CHECK(de==(i64)c[2][2],"de != c22");
          int cyc=(sz[x]>=3)&&(sz[y]>=3); if((de>0)!=cyc) nEviol++;
          if(de==0 && be!=0 && ga!=0) nadditive++;
          if(exact){ nexact++; CHECK(ntup<420000,"tuple overflow"); tup[4*ntup]=al; tup[4*ntup+1]=be; tup[4*ntup+2]=ga; tup[4*ntup+3]=de; ntup++; }
          /* literal spot check of the two-module law on the blow-ups (orders <= 6: every triple, a,b <= 3 gives order <= 10; orders 7,8: sparse) */
          if(m<=5 || (ntriples%211)==0){ for(int a=1;a<=3;a++) for(int b=1;b<=3;b++){ if(m+a+b-2>12) continue; unsigned B[BIGN]; int N=build_blow(m,out,x,a,y,b,B); CHECK((i64)hp(N,B)==h[a][b],"two-module law fails"); nlaw2++; } }
        } }
    }
    fprintf(fs,"order=%d distinct_values=%lld max_H=%lld seven=%d twentyone=%d\n",m,nval,maxH,seen[7],seen[21]);
    if(m>=2){ fprintf(fs,"order=%d unit_offset_roots=%lld with_d_eq_E=%lld with_M_zero=%lld outdeg:",m,nunit,nunit_dE,nunit_M0); for(int k=0;k<m;k++) if(outdeg_unit[k]) fprintf(fs," %d:%lld",k,outdeg_unit[k]); fprintf(fs,"\n"); }
    fprintf(fs,"base_order=%d root_order=%d unit_sets_by_size:",m,m+1); for(int k=0;k<=m;k++) if(us_by_size[k]) fprintf(fs," %d:%lld",k,us_by_size[k]); fprintf(fs,"\n");
    if(m>=3&&m<=8) fprintf(fs,"order=%d doubling_vertices=%lld with_c2_positive=%lld\n",m,ndbl,ndbl_pos);
    /* spectrum file of this order */
    snprintf(fn,sizeof fn,"%s/spectrum_order_%d.txt",outdir,m); FILE*fp=fopen(fn,"w"); CHECK(fp,"open spectrum"); for(int i=0;i<VMAX;i++) if(seen[i]) fprintf(fp,"%d\n",i); fclose(fp);
    fflush(fs); fflush(fr); fflush(fu); fflush(fd);
    fprintf(stderr,"order %d done\n",m);
  }
  for(int i=0;i<VMAX;i++) if(minord[i]) fprintf(fv,"%d %d %llu\n",i,minord[i],witc[i]);
  fprintf(fs,"identity_checks_roots=%lld (c0 = d - delta', M = H - E + delta', c0 odd, every root of every class of order 2..9)\n",nident);
  fprintf(fs,"cut_form_checks=%lld (c0 = H(B) - W(O), M = U_B(O), delta' = W(O): every (B,O) with |B| <= 7, and every unit set found over |B| = 8, 9, on the literal L = B + v)\n",nlaw1);
  fprintf(fs,"module_law_literal_checks=%lld\n",nlaw2);
  fprintf(fs,"doubling_flag_definition_vs_h2h3h4_test_differences=%lld\n",ndblflagdiff);
  fprintf(fs,"two_module_triples=%lld exact=%lld additive_de0_be_ga_nonzero=%lld E_violations(de>0 iff both on cycles)=%lld\n",ntriples,nexact,nadditive,nEviol);
  { size_t n8=ntup-ntup8start; i64*t8=(i64*)malloc(sizeof(i64)*4*(n8?n8:1)); memcpy(t8,tup+4*ntup8start,sizeof(i64)*4*n8); qsort(t8,n8,4*sizeof(i64),cmp4);
    size_t d8=0; for(size_t i=0;i<n8;i++) if(i==0||cmp4(t8+4*i,t8+4*(i-1))!=0) d8++;
    qsort(tup,ntup,4*sizeof(i64),cmp4); size_t dall=0; for(size_t i=0;i<ntup;i++) if(i==0||cmp4(tup+4*i,tup+4*(i-1))!=0){ dall++; fprintf(ft,"%lld %lld %lld %lld\n",tup[4*i],tup[4*i+1],tup[4*i+2],tup[4*i+3]); }
    fprintf(fs,"distinct_exact_tuples_order8=%llu distinct_exact_tuples_orders_3_to_8=%llu\n",(u64)d8,(u64)dall); }
  fclose(fs); fclose(fr); fclose(fu); fclose(fd); fclose(fv); fclose(ft); return 0; }

/* ------------------------------------------------------------------ seeds: verify line_seeds.txt */
static int seeds(const char*file){ FILE*f=fopen(file,"r"); CHECK(f,"open seeds file"); char line[1<<12]; long n=0,bad=0,ntail=0,nmod=0,nlit=0; long byord[MAXN+1]={0};
  while(fgets(line,sizeof line,f)){ if(line[0]=='#') continue; long c,s0; char mech[32]; int m,v; u64 code; i64 H,d,E,U;
    if(sscanf(line,"%ld %ld %31s %d %llu %d %lld %lld %lld %lld",&c,&s0,mech,&m,&code,&v,&H,&d,&E,&U)!=10) continue; n++; int ok=1; unsigned out[MAXN]; CHECK(m<=MAXN,"order"); decode(code,m,out);
    if(!strcmp(mech,"tail")){ ntail++; byord[m]++; tables(out,m); i64 h,dd,e,u,dp; rootq(out,m,v,&h,&dd,&e,&u,&dp);
      if(h!=H||dd!=d||e!=E||u!=U) ok=0; if(h-e+2*dd-u!=1) ok=0; if(dp!=dd-1) ok=0; if(u-dd!=(i64)c<<(s0-1)) ok=0; if((c&1)==0) ok=0;
      /* cut form on B = L - v */
      { unsigned ob[MAXN]; int idx[MAXN],k=0; for(int i=0;i<m;i++) if(i!=v) idx[k++]=i; for(int i=0;i<m-1;i++){ ob[i]=0; for(int j=0;j<m-1;j++) if(out[idx[i]]>>idx[j]&1) ob[i]|=1u<<j; }
        tables(ob,m-1); i64 UB[MAXN],W[MAXN][MAXN]; cutdata(m-1,UB,W); i64 w=0,mm=0; for(int i=0;i<m-1;i++) if(out[v]>>idx[i]&1){ mm+=UB[i]; for(int j=i+1;j<m-1;j++) if(out[v]>>idx[j]&1) w+=W[i][j]; }
        if(w!=(i64)hs_[(1<<(m-1))-1]-1) ok=0; if(mm!=u-dd) ok=0; }
      int y[8]={0,0,0,0,0,0,0,0}; for(int t=1;t<=3;t++){ if(m+t>13) break; unsigned B[BIGN]; int N=build_tail(m,out,v,t,y,B); i64 lit=(i64)hp(N,B); nlit++;
        if(lit!=1+((i64)c<<(s0-1+t))) ok=0; }
    } else { nmod++; CHECK(m<=8,"module order"); tables(out,m); i64 h=(i64)hs_[(1<<m)-1]; u64 cj[6]; walk1(m,out,v,cj);
      if(h!=H) ok=0; if(2*(i64)cj[2]!=h-1||cj[3]||cj[4]||cj[5]) ok=0; if(h-1!=(i64)c<<s0) ok=0;
      for(int a=2;a<=5;a++){ if(m+a-1>13) break; unsigned B[BIGN]; int N=build_blow(m,out,v,a,-1,0,B); i64 lit=(i64)hp(N,B); nlit++; if(lit!=1+((h-1)<<(a-1))) ok=0; } }
    if((c==3&&s0<2)||(c==5&&s0<3)) ok=0;
    if(!ok){ bad++; printf("BAD line: %s",line); } }
  fclose(f);
  printf("seeds checked=%ld tail=%ld module=%ld literal_counts=%ld bad=%ld\n",n,ntail,nmod,nlit,bad);
  printf("tail seeds by order:"); for(int m=2;m<=MAXN;m++) printf(" %d:%ld",m,byord[m]); printf("\n");
  return bad?1:0; }

/* ------------------------------------------------------------------ roots: verify the packet's unit-offset root list */
static int roots(const char*file){ FILE*f=fopen(file,"r"); CHECK(f,"open roots file"); char line[1<<12]; long n=0,bad=0;
  while(fgets(line,sizeof line,f)){ int m,v; u64 code; i64 H,d,E,U,M; if(sscanf(line,"%d %llu %d %lld %lld %lld %lld %lld",&m,&code,&v,&H,&d,&E,&U,&M)!=8) continue; n++;
    unsigned out[MAXN]; decode(code,m,out); tables(out,m); i64 h,dd,e,u,dp; rootq(out,m,v,&h,&dd,&e,&u,&dp); int ok=(h==H&&dd==d&&e==E&&u==U&&u-dd==M&&h-e+2*dd-u==1&&dp==dd-1);
    if(!ok){ bad++; fprintf(stderr,"BAD root line: %s",line); }
    u64 cc=canon(out,m,NULL,NULL); printf("%d %llu %lld %lld %lld %lld %lld\n",m,cc,H,d,E,U,M); }
  fclose(f); fprintf(stderr,"roots checked=%ld bad=%ld\n",n,bad); return bad?1:0; }

/* ------------------------------------------------------------------ labelled8: all labelled tournaments, orders 1..8 */
static int labelled8(void){ static unsigned char seen[9][VMAX]; memset(seen,0,sizeof seen); i64 total8=0; i64 ncheck=0;
  for(int n=1;n<=7;n++){ int nb=n*(n-1)/2; int full=(1<<n)-1;
    for(u64 code=0;code<(1ULL<<nb);code++){ unsigned out[MAXN]; decode(code,n,out); tables(out,n); i64 H=(i64)hs_[full]; seen[n][H]=1;
      if(n==7){ /* all 2^7 one-vertex extensions through the cut formula (own derivation: position of the new vertex) */
        i64 Q[7][7]; for(int a=0;a<7;a++) for(int b=0;b<7;b++) Q[a][b]=0;
        for(int X=1;X<full;X++){ int Y=full^X; unsigned ma=X; while(ma){ int a=__builtin_ctz(ma); ma&=ma-1; i64 ea=(i64)en[X][a]; if(!ea) continue;
            unsigned mb=Y; while(mb){ int b=__builtin_ctz(mb); mb&=mb-1; Q[a][b]+=ea*(i64)st[Y][b]; } } }
        for(int O=0;O<=full;O++){ i64 h=0; for(int b=0;b<7;b++) if(O>>b&1) h+=(i64)st[full][b]; else h+=(i64)en[full][b];
          for(int a=0;a<7;a++) if(!(O>>a&1)) for(int b=0;b<7;b++) if(O>>b&1) h+=Q[a][b];
          CHECK(h>0&&h<VMAX,"H range"); seen[8][h]=1; total8+=h;
          if((code*131+O)%1000003==0){ unsigned o2[MAXN]; for(int v=0;v<7;v++) o2[v]=out[v]|((O>>v&1)?0u:(1u<<7)); o2[7]=(unsigned)O; CHECK((i64)hp(8,o2)==h,"cut formula vs DP"); ncheck++; } } } } }
  /* Szele: the sum of H over all labelled tournaments of order 8 is 2^28 * 8!/2^7 */
  printf("sum of H over all 2^28 labelled 8-tournaments = %lld, expected 2^21 * 8! = %lld %s; cut-formula spot checks against DP: %lld\n",total8,(i64)(1<<21)*40320,(total8==(i64)(1<<21)*40320)?"OK":"MISMATCH",ncheck);
  for(int n=1;n<=8;n++){ int cnt=0,mx=0; for(int i=0;i<VMAX;i++) if(seen[n][i]){cnt++;mx=i;} printf("order %d: distinct values %d, max %d, 7 attained: %d, 21 attained: %d\n",n,cnt,mx,seen[n][7],seen[n][21]); }
  printf("values_order8:"); for(int i=0;i<VMAX;i++) if(seen[8][i]) printf(" %d",i); printf("\n");
  { int miss=0; printf("odd numbers below 200 not attained at order <= 8:"); for(int i=1;i<200;i+=2){ int a=0; for(int n=1;n<=8;n++) a|=seen[n][i]; if(!a){ printf(" %d",i); miss++; } } printf("  (count %d)\n",miss); }
  return 0; }

/* ------------------------------------------------------------------ mean: extension-mean law by literal enumeration */
static i64 covP[8][64]; /* covP[j][mask] */
static void Pj(const unsigned*out,int n,i64*P){ int full=(1<<n)-1; tables(out,n);
  for(int s=0;s<=full;s++) covP[0][s]=(s==0); for(int j=1;j<=n&&j<8;j++){ for(int s=0;s<=full;s++){ i64 t=0; for(int sub=s;sub;sub=(sub-1)&s) t+=(i64)hs_[sub]*covP[j-1][s^sub]; covP[j][s]=t; } P[j]=covP[j][full]; }
  for(int j=n+1;j<8;j++) P[j]=0; }
static i64 binom(int n,int k){ if(k<0||k>n) return 0; i64 r=1; for(int i=1;i<=k;i++) r=r*(n-k+i)/i; return r; }
static i64 S2(int n,int k){ if(n==k) return 1; if(k==0||k>n) return 0; return k*S2(n-1,k)+S2(n-1,k-1); }
static int mean_one(int n,u64 code,int k,i64*cnt){ unsigned out[MAXN]; decode(code,n,out); i64 P[8]; for(int j=0;j<8;j++) P[j]=0; Pj(out,n,P);
  /* bound P_j >= j! S(n,j), equality for all j iff transitive */
  { int alleq=1; i64 f=1; for(int j=1;j<=n;j++){ f*=j; if(P[j]<f*S2(n,j)) return 0; if(P[j]!=f*S2(n,j)) alleq=0; } if(alleq!=(P[1]==1)) return 0; }
  int nb=k*n+k*(k-1)/2; i64 tot=0; int N=n+k;
  for(u64 bits=0;bits<(1ULL<<nb);bits++){ unsigned B[BIGN]; for(int i=0;i<n;i++) B[i]=out[i]; for(int i=n;i<N;i++) B[i]=0; int b=0;
    for(int a=0;a<k;a++){ int z=n+a; for(int u=0;u<n;u++){ if(bits>>b&1) B[z]|=1u<<u; else B[u]|=1u<<z; b++; } for(int a2=a+1;a2<k;a2++){ int z2=n+a2; if(bits>>b&1) B[z]|=1u<<z2; else B[z2]|=1u<<z; b++; } }
    tot+=(i64)hp(N,B); (*cnt)++; }
  i64 fk=1; for(int i=2;i<=k;i++) fk*=i; i64 rhs=0; for(int j=1;j<=k+1&&j<8;j++) rhs+=binom(k+1,j)*((i64)1<<(k+1-j))*P[j];
  /* tot * 2^(2k) == 2^nb * k! * sum_j C(k+1,j) 2^(k+1-j) P_j */
  return ( (__int128)tot<<(2*k) ) == ( (__int128)fk*rhs<<nb ); }
static int mean(void){ int allok=1; struct { int k,n; } U[]={{1,1},{1,2},{1,3},{1,4},{1,5},{1,6},{2,1},{2,2},{2,3},{2,4},{2,5},{3,1},{3,2},{3,3},{3,4},{4,1},{4,2},{4,3},{5,1},{5,2},{6,1}};
  for(size_t q=0;q<sizeof U/sizeof U[0];q++){ int k=U[q].k,n=U[q].n; int nb=n*(n-1)/2; i64 cnt=0; int ok=1; long nt=0;
    for(u64 code=0;code<(1ULL<<nb);code++){ if(!mean_one(n,code,k,&cnt)) ok=0; nt++; }
    printf("%s (C) k=%d, all %ld labelled tournaments of order n=%d, %lld extensions counted literally\n",ok?"PASS":"FAIL",k,nt,n,cnt); if(!ok) allok=0; }
  /* k = 3, n = 5 and k = 2, n = 6: sparse deterministic sample */
  { u64 x=88172645463325252ULL; for(int rep=0;rep<6;rep++){ x^=x<<13; x^=x>>7; x^=x<<17; u64 code=x&((1ULL<<10)-1); i64 cnt=0; int ok=mean_one(5,code,3,&cnt); printf("%s (C) k=3, n=5, code %llu, %lld extensions\n",ok?"PASS":"FAIL",code,cnt); if(!ok) allok=0; }
    for(int rep=0;rep<12;rep++){ x^=x<<13; x^=x>>7; x^=x<<17; u64 code=x&((1ULL<<15)-1); i64 cnt=0; int ok=mean_one(6,code,2,&cnt); printf("%s (C) k=2, n=6, code %llu, %lld extensions\n",ok?"PASS":"FAIL",code,cnt); if(!ok) allok=0; } }
  printf(allok?"MEAN ALL PASS\n":"MEAN FAILURES\n"); return allok?0:1; }

/* ------------------------------------------------------------------ oneapex: (D) */
static i64 grec(u64 v){ if(v==0) return 1; if(v&1) return grec(v>>1); return grec(v>>1)+(i64)(v>>1)+(i64)((v>>1)&1); }
static int oneapex(int LMAX){ CHECK(LMAX+1<=BIGN,"LMAX"); int ok_g=1, ok_traj=1, ok_even=1; long nw=0, nm=0;
  i64 *Glit=(i64*)calloc((size_t)1<<(LMAX+1),sizeof(i64)); /* Glit[m] = literal count of the word of length LMAX with value m-1 */
  for(int L=1;L<=LMAX;L++) for(u64 w=0;w<(1ULL<<L);w++){ /* word letters w_1..w_L, w_1 = most significant bit; chain 0 -> 1 -> ... -> L-1; apex L beats chain vertex i iff letter i is 1 */
      unsigned B[BIGN]; for(int i=0;i<=L;i++) B[i]=0; for(int i=0;i<L;i++){ for(int j=i+1;j<L;j++) B[i]|=1u<<j; int bit=(int)(w>>(L-1-i)&1); if(bit) B[L]|=1u<<i; else B[i]|=1u<<L; }
      i64 lit=(i64)hp(L+1,B); nw++; if(lit!=grec(w)) ok_g=0; if(L==LMAX) Glit[w+1]=lit; }
  for(u64 m=1;m<=(1ULL<<LMAX);m++){ nm++; if(m&1){ u64 o=m; i64 tot=1; while(o!=1){ tot+=2*(i64)((o+1)/4); o=o+1; while(!(o&1)) o>>=1; } if(tot!=Glit[m]) ok_traj=0; }
    else if(Glit[m]!=Glit[m/2]) ok_even=0; }
  printf("%s (D) literal one-apex count = g(value) for all %ld words of length <= %d\n",ok_g?"PASS":"FAIL",nw,LMAX);
  printf("%s (D) G(m) = 1 + sum over the trajectory of 2 floor((o+1)/4) for odd m, literal G, all m <= 2^%d (%ld values)\n",ok_traj?"PASS":"FAIL",LMAX,nm);
  printf("%s (D) G(2x) = G(x), literal G\n",ok_even?"PASS":"FAIL");
  /* first odd numbers missed by one-apex counts */
  { static unsigned char hit[1<<16]; memset(hit,0,sizeof hit); for(u64 v=0;v<(1ULL<<20);v++){ i64 g=grec(v); if(g<(1<<16)) hit[g]=1; } printf("odd numbers below 100 that are not one-apex counts (words of length <= 20):"); for(int i=1;i<100;i+=2) if(!hit[i]) printf(" %d",i); printf("\n"); }
  free(Glit); return (ok_g&&ok_traj&&ok_even)?0:1; }

/* ------------------------------------------------------------------ dbl9: doubling vertices at order 9 (audit extension)
 * x is doubling iff c_2 = (H-1)/2 and c_3 = 0 (c_3 = 0 forces c_j = 0 for j >= 4: two x-loops of a walk can be merged
 * into one by Camion's theorem, so a walk with j >= 3 visits yields one with j - 1 visits). */
static u64 wf9[512][4][9];
static void walk1_9(int n,const unsigned*out,int x,u64*c){ int full=(1<<n)-1, oth=full&~(1<<x);
  for(int s=0;s<=full;s++) if(!(s>>x&1)) memset(wf9[s],0,sizeof wf9[s]);
  for(int v=0;v<n;v++){ if(v==x) wf9[0][1][x]=1; else wf9[1<<v][0][v]=1; }
  for(int s=0;s<=full;s++){ if(s>>x&1) continue; for(int k=0;k<=3;k++) for(int last=0;last<n;last++){ u64 f=wf9[s][k][last]; if(!f) continue;
      unsigned nx=out[last]; if((nx>>x&1)&&k<3) wf9[s][k+1][x]+=f;
      unsigned m=nx&oth&~s; while(m){ int u=__builtin_ctz(m); m&=m-1; wf9[s|1u<<u][k][u]+=f; } } }
  for(int j=0;j<=3;j++){ u64 t=0; for(int last=0;last<n;last++) t+=wf9[oth][j][last]; c[j]=t; } }
static int dbl9(void){ gen_classes(9); long ndbl=0,npos=0; static unsigned char val[VMAX]; memset(val,0,sizeof val);
  for(int r=0;r<nrep[9];r++){ unsigned out[MAXN]; decode(reps[9][r],9,out); tables(out,9); i64 H=(i64)hs_[511];
    for(int x=0;x<9;x++){ u64 c[4]; walk1_9(9,out,x,c); CHECK((i64)c[1]==H,"c1 != H (order 9)"); CHECK(c[0]==hs_[511^(1<<x)],"c0 (order 9)");
      if(2*(i64)c[2]==H-1 && c[3]==0){ ndbl++; if(c[2]>0){ npos++; val[H]=1; printf("9 %llu %d %lld %llu\n",reps[9][r],x,H,c[2]); } } } }
  int nv=0; for(int i=0;i<VMAX;i++) nv+=val[i];
  fprintf(stderr,"order 9: doubling vertices %ld, with c_2 > 0: %ld, distinct values H with a doubling realiser of order 9: %d\n",ndbl,npos,nv); return 0; }

/* ------------------------------------------------------------------ exact4: is exactness on the grid a,b <= 3 (the packet's
 * test in modfam.c) the same as true exactness for all a,b? At order <= 8 a walk visits a vertex at most 4 times, so the
 * family H(K[TT_a at x, TT_b at y]) is exactly bilinear in (2^a, 2^b) iff c_ij = 0 whenever i >= 3 or j >= 3 (i,j <= 4). */
static u64 wg4[256][5][5][8];
static int exact4(void){ gen_classes(8); long ntri=0, grid=0, truex=0, mismatch=0; long byord[9]={0};
  for(int m=3;m<=8;m++) for(int r=0;r<nrep[m];r++){ unsigned out[MAXN]; decode(reps[m][r],m,out); int full=(1<<m)-1;
    for(int x=0;x<m;x++) for(int y=0;y<m;y++){ if(x==y) continue; int oth=full&~(1<<x)&~(1<<y); ntri++;
      for(int s=0;s<=full;s++) if(!((s>>x&1)||(s>>y&1))) memset(wg4[s],0,sizeof wg4[s]);
      for(int v=0;v<m;v++){ if(v==x) wg4[0][1][0][x]=1; else if(v==y) wg4[0][0][1][y]=1; else wg4[1<<v][0][0][v]=1; }
      for(int s=0;s<=full;s++){ if((s>>x&1)||(s>>y&1)) continue;
        for(int t=0;t<=8;t++) for(int kx=0;kx<=4;kx++){ int ky=t-kx; if(ky<0||ky>4) continue;
          for(int last=0;last<m;last++){ u64 f=wg4[s][kx][ky][last]; if(!f) continue; unsigned nx=out[last];
            if((nx>>x&1)&&kx<4) wg4[s][kx+1][ky][x]+=f;
            if((nx>>y&1)&&ky<4) wg4[s][kx][ky+1][y]+=f;
            unsigned mm=nx&oth&~s; while(mm){ int u=__builtin_ctz(mm); mm&=mm-1; wg4[s|1u<<u][kx][ky][u]+=f; } } } }
      u64 c[5][5]; for(int i=0;i<=4;i++) for(int j=0;j<=4;j++){ u64 t=0; for(int last=0;last<m;last++) t+=wg4[oth][i][j][last]; c[i][j]=t; }
      int g=1, tr=1; for(int i=1;i<=4;i++) for(int j=1;j<=4;j++) if((i>=3||j>=3)&&c[i][j]){ tr=0; if(i<=3&&j<=3) g=0; }
      grid+=g; truex+=tr; if(g!=tr){ mismatch++; byord[m]++; } } }
  printf("triples (orders 3..8): %ld; exact on the grid a,b <= 3: %ld; truly exact (no walk with 3 or 4 visits of x or y): %ld; grid-exact but not truly exact: %ld (by order:",ntri,grid,truex,mismatch);
  for(int m=3;m<=8;m++) printf(" %d:%ld",m,byord[m]); printf(")\n"); return mismatch?1:0; }

int main(int argc,char**argv){ if(argc<2){ fprintf(stderr,"usage: see header\n"); return 2; }
  if(!strcmp(argv[1],"dbl9")) return dbl9();
  if(!strcmp(argv[1],"exact4")) return exact4();
  if(!strcmp(argv[1],"census")&&argc>=3) return census(argv[2]);
  if(!strcmp(argv[1],"seeds")&&argc>=3) return seeds(argv[2]);
  if(!strcmp(argv[1],"roots")&&argc>=3) return roots(argv[2]);
  if(!strcmp(argv[1],"labelled8")) return labelled8();
  if(!strcmp(argv[1],"mean")) return mean();
  if(!strcmp(argv[1],"oneapex")&&argc>=3) return oneapex(atoi(argv[2]));
  fprintf(stderr,"unknown mode\n"); return 2; }

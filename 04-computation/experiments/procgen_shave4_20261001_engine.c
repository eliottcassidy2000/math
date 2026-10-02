// procgen_shave4_20261001_engine.c  --  shave4 lane (collatz-procgen-20260922), 2026-10-01.
// Parity of embedding counts of spanning oriented graphs S in tournaments T ("Redei graphs").
// Independent of opus S15's code. Single-threaded; memory < 100 MB.
//
// Modes (tournament classes are gentourng -q upper-triangle strings, 1 = i->j for i<j):
//   hp                         < classes   : H, hc, #copies of H_n (= HPs with first->last); THM-4526(A) check
//   til n                      < classes   : tiling table + zeta transform; prints "R <mask> i-j ..." for every
//                                            chord set C with P_n + C Redei, and summary counts
//   chords n [i-j ...]         < classes   : per class, mod-2 counts of HPs containing single chords, chord pairs,
//                                            and every subset of the listed candidate chords
//   dag n classfile            < directg -a -G : parity status of every DAG class via Gray-code completions
//   redei n classfile          < directg -a -G : Redei test by backtracking + proven filters (B4)
//   rigid n classfile          < directg -o -G : parity-rigidity test (constant phi) for every input graph
//   witness n seed tries [i-j ...]            : random tournament with an EVEN number of HPs satisfying the chords
//   unav n classfile e                        : number of forward e-arc sets that embed in every n-tournament (n<=7);
//                                               the masks are printed as "U <mask>"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAXN 16

static int parse_tour(const char *s, int nn, uint32_t *out) {
    int len = (int)strlen(s);
    while (len && (s[len-1]=='\n' || s[len-1]=='\r' || s[len-1]==' ')) len--;
    if (len != nn*(nn-1)/2) return 0;
    for (int i=0;i<nn;i++) out[i]=0;
    int b=0;
    for (int i=0;i<nn;i++) for (int j=i+1;j<nn;j++) { if (s[b]=='1') out[i]|=1u<<j; else out[j]|=1u<<i; b++; }
    return 1;
}
static int read_classes_file(const char *fn, int nn, uint32_t (**REP)[MAXN]) {
    FILE *f = fn ? fopen(fn,"r") : stdin; if (!f) { fprintf(stderr,"cannot open %s\n",fn); exit(1); }
    int cap=1<<16, nc=0; uint32_t (*R)[MAXN] = malloc(sizeof(*R)*cap); char line[256];
    while (fgets(line,sizeof line,f)) { if (nc>=cap) { cap*=2; R=realloc(R,sizeof(*R)*cap); }
        if (!parse_tour(line,nn,R[nc])) { fprintf(stderr,"bad class line\n"); exit(1); } nc++; }
    if (fn) fclose(f);
    *REP=R; return nc;
}

// ======================= hp =======================
static uint64_t dp[1<<MAXN][MAXN];
static uint64_t count_H(const uint32_t *out, int nn) {
    uint32_t full=(1u<<nn)-1;
    for (uint32_t m=0;m<=full;m++) for (int v=0;v<nn;v++) dp[m][v]=0;
    for (int v=0;v<nn;v++) dp[1u<<v][v]=1;
    for (uint32_t m=1;m<=full;m++) for (int v=0;v<nn;v++) { uint64_t x=dp[m][v]; if(!x) continue;
        uint32_t nx = out[v] & ~m; while(nx){ int w=__builtin_ctz(nx); nx&=nx-1; dp[m|(1u<<w)][w]+=x; } }
    uint64_t H=0; for (int v=0;v<nn;v++) H+=dp[full][v]; return H;
}
static void count_from(const uint32_t *out, int nn, int s, uint64_t *P) {
    uint32_t full=(1u<<nn)-1;
    for (uint32_t m=0;m<=full;m++) if (m>>s&1) for (int v=0;v<nn;v++) dp[m][v]=0;
    dp[1u<<s][s]=1;
    for (uint32_t m=1;m<=full;m++) { if(!(m>>s&1)) continue; for (int v=0;v<nn;v++) { uint64_t x=dp[m][v]; if(!x) continue;
        uint32_t nx = out[v] & ~m; while(nx){ int w=__builtin_ctz(nx); nx&=nx-1; dp[m|(1u<<w)][w]+=x; } } }
    for (int t=0;t<nn;t++) P[t]=dp[full][t];
}
static int mode_hp(void) {
    char line[256]; long ncls=0, nzero=0, nodd=0, neven=0; uint64_t mincopies=~0ULL; int nn=-1;
    while (fgets(line,sizeof line,stdin)) {
        int len=(int)strlen(line); while(len&&(line[len-1]=='\n'||line[len-1]=='\r')) line[--len]=0;
        if (nn<0) { nn=2; while(nn*(nn-1)/2<len) nn++; if(nn*(nn-1)/2!=len){fprintf(stderr,"bad len\n");return 1;} }
        uint32_t out[MAXN]; if(!parse_tour(line,nn,out)){fprintf(stderr,"bad line\n");return 1;}
        uint64_t H=count_H(out,nn), P[MAXN];
        count_from(out,nn,0,P); uint64_t hc=0; for(int t=1;t<nn;t++) if(out[t]&1u) hc+=P[t];
        uint64_t direct=0, closing=0, Htot=0;
        for (int s=0;s<nn;s++){ count_from(out,nn,s,P); for(int t=0;t<nn;t++){ if(t==s) continue; Htot+=P[t];
            if(out[s]>>t&1) direct+=P[t]; else closing+=P[t]; } }
        if (Htot!=H) { fprintf(stderr,"H mismatch\n"); return 1; }
        if (closing != (uint64_t)nn*hc) { fprintf(stderr,"closing != n*hc\n"); return 1; }
        ncls++;
        if (direct==0) { nzero++; printf("ZERO %s H=%llu hc=%llu\n", line,(unsigned long long)H,(unsigned long long)hc); }
        else if (direct<mincopies) mincopies=direct;
        if (direct&1) nodd++; else neven++;
    }
    printf("HPSUMMARY n=%d classes=%ld zero=%ld odd=%ld even=%ld min_nonzero=%llu\n", nn, ncls, nzero, nodd, neven,(unsigned long long)mincopies);
    return 0;
}

// ======================= tilings / chords =======================
static int NC, ci_[64], cj_[64], chord_idx[MAXN][MAXN];
static void setup_chords(int nn) {
    NC=0; for(int i=0;i<nn;i++) for(int j=0;j<nn;j++) chord_idx[i][j]=-1;
    for(int i=0;i<nn;i++) for(int j=i+2;j<nn;j++){ chord_idx[i][j]=NC; ci_[NC]=i; cj_[NC]=j; NC++; }
}
static int32_t *tiltab; static uint32_t *hits; static int cur_class; static int hv[MAXN]; static long tbad=0;
static void til_dfs(const uint32_t *out, int nn, int depth, uint32_t used) {
    if (depth==nn) { uint32_t t=0; for (int k=0;k<NC;k++) if (out[hv[ci_[k]]]>>hv[cj_[k]]&1) t|=1u<<k;
        if (tiltab[t]==-1) tiltab[t]=cur_class; else if (tiltab[t]!=cur_class) tbad++; hits[t]++; return; }
    uint32_t cand = depth==0 ? ((1u<<nn)-1) : (out[hv[depth-1]] & ~used);
    while (cand) { int w=__builtin_ctz(cand); cand&=cand-1; hv[depth]=w; til_dfs(out,nn,depth+1,used|(1u<<w)); }
}
static int mode_til(int nn) {
    setup_chords(nn); uint32_t (*REP)[MAXN]; int ncls=read_classes_file(NULL,nn,&REP);
    uint32_t NT=1u<<NC; tiltab=malloc(sizeof(int32_t)*NT); hits=calloc(NT,sizeof(uint32_t));
    for(uint32_t t=0;t<NT;t++) tiltab[t]=-1;
    for (int c=0;c<ncls;c++){ cur_class=c; til_dfs(REP[c],nn,0,0); }
    long unassigned=0; for(uint32_t t=0;t<NT;t++) if(tiltab[t]<0) unassigned++;
    uint32_t *mult=calloc(ncls,sizeof(uint32_t)); long nonuni=0; uint64_t *ntil=calloc(ncls,sizeof(uint64_t));
    for(uint32_t t=0;t<NT;t++){ int c=tiltab[t]; if(c<0) continue; ntil[c]++; if(!mult[c]) mult[c]=hits[t]; else if(mult[c]!=hits[t]) nonuni++; }
    long double fact=1; for(int i=2;i<=nn;i++) fact*=i;
    long double orbit=0; long oddtil=0; for(int c=0;c<ncls;c++){ orbit += fact/mult[c]; if(ntil[c]&1) oddtil++; }
    printf("TIL n=%d chords=%d classes=%d tilings=%u unassigned=%ld conflicts=%ld nonuniform=%ld odd_tiling_classes=%ld orbit_sum=%.0Lf\n",
        nn,NC,ncls,NT,unassigned,tbad,nonuni,oddtil,orbit);
    uint64_t *A=malloc(sizeof(uint64_t)*NT); uint8_t *redei=malloc(NT), *evenr=malloc(NT);
    memset(redei,1,NT); memset(evenr,1,NT);
    for (int b0=0;b0<ncls;b0+=64) {
        int bs = ncls-b0<64 ? ncls-b0 : 64; uint64_t full = bs==64? ~0ULL : ((1ULL<<bs)-1);
        for(uint32_t t=0;t<NT;t++){ int c=tiltab[t]; A[t] = (c>=b0 && c<b0+bs) ? (1ULL<<(c-b0)) : 0; }
        for (int k=0;k<NC;k++){ uint32_t bit=1u<<k; for(uint32_t t=0;t<NT;t++) if(!(t&bit)) A[t]^=A[t|bit]; }
        for(uint32_t t=0;t<NT;t++){ if(A[t]!=full) redei[t]=0; if(A[t]!=0) evenr[t]=0; }
    }
    long nr=0, ne=0;
    for(uint32_t t=0;t<NT;t++){ if(redei[t]){ nr++; printf("R %u",t); for(int k=0;k<NC;k++) if(t>>k&1) printf(" %d-%d",ci_[k],cj_[k]); printf("\n"); } if(evenr[t]) ne++; }
    printf("TILSUMMARY n=%d redei=%ld even_rigid=%ld\n",nn,nr,ne);
    return 0;
}
static int ch_n; static uint32_t ch_out[MAXN]; static int ch_v[MAXN];
static uint64_t acc1, acc2[64]; static uint32_t candmask[256]; static int ncand; static uint8_t candpar[256];
static void ch_dfs(int d, uint32_t used){
    if (d==ch_n){ uint64_t t=0; for(int k=0;k<NC;k++) if(ch_out[ch_v[ci_[k]]]>>ch_v[cj_[k]]&1) t|=1ULL<<k;
        acc1^=t; uint64_t x=t; while(x){ int c=__builtin_ctzll(x); x&=x-1; acc2[c]^=t; }
        for(int q=0;q<ncand;q++) if((t&candmask[q])==candmask[q]) candpar[q]^=1;
        return; }
    uint32_t cand = d==0 ? ((1u<<ch_n)-1) : (ch_out[ch_v[d-1]]&~used);
    while(cand){ int w=__builtin_ctz(cand); cand&=cand-1; ch_v[d]=w; ch_dfs(d+1,used|(1u<<w)); }
}
static int mode_chords(int nn, int argc, char **argv) {
    ch_n=nn; setup_chords(nn); if (NC>64) { fprintf(stderr,"too many chords\n"); return 1; }
    int cc[8], ncc=0; for(int a=0;a<argc;a++){ int i,j; sscanf(argv[a],"%d-%d",&i,&j); cc[ncc++]=chord_idx[i][j]; }
    ncand=1<<ncc; for(int s=0;s<ncand;s++){ candmask[s]=0; for(int k=0;k<ncc;k++) if(s>>k&1) candmask[s]|=1u<<cc[k]; }
    uint64_t all1=~0ULL, all2[64]; for(int c=0;c<64;c++) all2[c]=~0ULL; uint8_t allc[256]; memset(allc,1,sizeof allc);
    char line[256]; long ncls=0;
    while(fgets(line,sizeof line,stdin)){
        if(!parse_tour(line,nn,ch_out)){ fprintf(stderr,"bad line\n"); return 1; }
        acc1=0; memset(acc2,0,sizeof acc2); memset(candpar,0,sizeof candpar);
        ch_dfs(0,0);
        all1&=acc1; for(int c=0;c<NC;c++) all2[c]&=acc2[c]; for(int q=0;q<ncand;q++) if(!candpar[q]) allc[q]=0;
        ncls++;
    }
    printf("CHORDS n=%d classes=%ld\n",nn,ncls);
    printf("SINGLE"); for(int c=0;c<NC;c++) if(all1>>c&1) printf(" %d-%d",ci_[c],cj_[c]); printf("\n");
    printf("PAIRS"); for(int c=0;c<NC;c++) for(int d=c+1;d<NC;d++) if(all2[c]>>d&1) printf(" {%d-%d,%d-%d}",ci_[c],cj_[c],ci_[d],cj_[d]); printf("\n");
    printf("CANDIDATES"); for(int q=0;q<ncand;q++) if(allc[q]){ printf(" {"); int f=1; for(int k=0;k<ncc;k++) if(q>>k&1){ printf("%s%d-%d",f?"":",",ci_[cc[k]],cj_[cc[k]]); f=0;} printf("}"); } printf("\n");
    return 0;
}

// ======================= general graphs: completions (n <= 7) =======================
static int pidx[MAXN][MAXN];
static int16_t *build_class_table(int nn, uint32_t (*REP)[MAXN], int ncls) {
    int P=nn*(nn-1)/2, b=0; for(int i=0;i<nn;i++) for(int j=i+1;j<nn;j++){ pidx[i][j]=pidx[j][i]=b++; }
    uint32_t NT=1u<<P; int16_t *tab=malloc(sizeof(int16_t)*NT); for(uint32_t t=0;t<NT;t++) tab[t]=-1;
    int perm[MAXN]; for(int i=0;i<nn;i++) perm[i]=i; int c_[MAXN]={0}; long conflicts=0;
    #define FILLPERM for(int cc=0;cc<ncls;cc++){ uint32_t code=0; for(int i=0;i<nn;i++) for(int j=i+1;j<nn;j++) if(REP[cc][perm[i]]>>perm[j]&1) code|=1u<<pidx[i][j]; \
        if(tab[code]==-1) tab[code]=(int16_t)cc; else if(tab[code]!=cc) conflicts++; }
    FILLPERM
    int i=0; while(i<nn){ if(c_[i]<i){ if(i%2==0){int t=perm[0];perm[0]=perm[i];perm[i]=t;} else {int t=perm[c_[i]];perm[c_[i]]=perm[i];perm[i]=t;}
            FILLPERM c_[i]++; i=0; } else { c_[i]=0; i++; } }
    long unf=0; for(uint32_t t=0;t<NT;t++) if(tab[t]<0) unf++;
    if (unf||conflicts) { fprintf(stderr,"class table broken: unfilled=%ld conflicts=%ld\n",unf,conflicts); exit(1); }
    return tab;
}
static int mode_dag(int nn, const char *clsfile) {
    uint32_t (*REP)[MAXN]; int ncls=read_classes_file(clsfile,nn,&REP);
    int16_t *tab=build_class_table(nn,REP,ncls); int P=nn*(nn-1)/2;
    char line[2048]; uint8_t *par=malloc(ncls); long ndag=0, nr1=0, nr0=0, nr0odd=0;
    while (fgets(line,sizeof line,stdin)) {
        int vals[400], nv=0; char *p=line,*e; while(1){ long v=strtol(p,&e,10); if(e==p) break; vals[nv++]=(int)v; p=e; }
        if (nv<3) continue; if (vals[0]!=nn) { fprintf(stderr,"wrong n\n"); return 1; }
        int ne=vals[1]; long grp=vals[2]; uint32_t fixed=0, under=0, outS[MAXN]={0}, inS[MAXN]={0};
        for (int k=0;k<ne;k++){ int u=vals[3+2*k], v=vals[4+2*k]; int q=pidx[u][v]; under|=1u<<q; if(u<v) fixed|=1u<<q; outS[u]|=1u<<v; inS[v]|=1u<<u; }
        static uint64_t le[1<<MAXN]; uint32_t full=(1u<<nn)-1; for(uint32_t m=0;m<=full;m++) le[m]=0; le[0]=1;
        for(uint32_t m=0;m<full;m++){ if(!le[m]) continue; for(int v=0;v<nn;v++) if(!(m>>v&1) && (inS[v]&~m)==0) le[m|(1u<<v)]+=le[m]; }
        int freeq[64], nf=0; for(int q=0;q<P;q++) if(!(under>>q&1)) freeq[nf++]=q;
        memset(par,0,ncls); uint32_t code=fixed; par[tab[code]]^=1;
        for (uint64_t g=1; g<(1ULL<<nf); g++){ int k=__builtin_ctzll(g); code^=1u<<freeq[k]; par[tab[code]]^=1; }
        int same=1; for(int c=1;c<ncls;c++) if(par[c]!=par[0]){ same=0; break; }
        int iso=0; for(int v=0;v<nn;v++) if(!(outS[v]|inS[v])) iso++;
        ndag++;
        if (same) { if ((le[full]&1)!=par[0]) { fprintf(stderr,"EXT MISMATCH\n"); return 1; }
            int autodd=(grp%2==1)&&iso<=1;
            if (par[0]) nr1++; else { nr0++; if(autodd) nr0odd++; }
            if (par[0]||autodd) { printf("%s n=%d m=%d grp=%ld iso=%d ext=%llu arcs", par[0]?"REDEI":"EVEN0", nn, ne, grp, iso, (unsigned long long)le[full]);
                for (int k=0;k<ne;k++) printf(" %d>%d", vals[3+2*k], vals[4+2*k]); printf("\n"); } }
    }
    printf("DAGSUMMARY n=%d dags=%ld redei=%ld even_rigid=%ld even_rigid_oddaut=%ld\n",nn,ndag,nr1,nr0,nr0odd);
    return 0;
}

// ======================= backtracking parity (redei / rigid) =======================
static int bn, bncls; static uint16_t (*BTO)[MAXN], (*BTI)[MAXN]; static uint16_t predo[MAXN], predi[MAXN]; static int bimg[MAXN];
static uint32_t emb_parity(int c) {
    uint32_t cnt=0; uint16_t cand[MAXN]; int k=0; uint16_t used=0, full=(uint16_t)((1u<<bn)-1); cand[0]=full;
    while (1) {
        if (!cand[k]) { if (k==0) break; k--; used&=~(1u<<bimg[k]); continue; }
        int y=__builtin_ctz(cand[k]); cand[k]&=cand[k]-1; bimg[k]=y;
        if (k==bn-1) { cnt^=1; continue; }
        used|=1u<<y; k++;
        uint16_t cc=full&~used; uint16_t p=predo[k]; while(p){ int j=__builtin_ctz(p); p&=p-1; cc&=BTO[c][bimg[j]]; }
        p=predi[k]; while(p){ int j=__builtin_ctz(p); p&=p-1; cc&=BTI[c][bimg[j]]; }
        cand[k]=cc;
    }
    return cnt;
}
static uint16_t gadj[MAXN]; static int pinv[MAXN];
static int inv_rec(int start, int nontriv) {   // does the undirected graph gadj have an automorphism of order 2?
    while (start<bn && pinv[start]>=0) start++;
    if (start==bn) { if (!nontriv) return 0;
        for (int u=0;u<bn;u++) for (int v=u+1;v<bn;v++) { int a=(gadj[u]>>v)&1, b=(gadj[pinv[u]]>>pinv[v])&1; if(a!=b) return 0; }
        return 1; }
    pinv[start]=start; if (inv_rec(start+1,nontriv)) return 1; pinv[start]=-1;
    for (int w=start+1; w<bn; w++) if (pinv[w]<0 && __builtin_popcount(gadj[w])==__builtin_popcount(gadj[start])) {
        pinv[start]=w; pinv[w]=start; if (inv_rec(start+1,1)) return 1; pinv[start]=-1; pinv[w]=-1; }
    return 0;
}
static int mode_backtrack(int nn, const char *clsfile, int rigid) {
    bn=nn; uint32_t (*REP)[MAXN]; bncls=read_classes_file(clsfile,nn,&REP);
    BTO=malloc(sizeof(*BTO)*bncls); BTI=malloc(sizeof(*BTI)*bncls);
    for (int c=0;c<bncls;c++) for (int v=0;v<nn;v++) { BTO[c][v]=(uint16_t)REP[c][v]; BTI[c][v]=0; }
    for (int c=0;c<bncls;c++) for (int u=0;u<nn;u++) for (int v=0;v<nn;v++) if (REP[c][u]>>v&1) BTI[c][v]|=1u<<u;
    int *ord=malloc(sizeof(int)*bncls); for(int i=0;i<bncls;i++) ord[i]=i;
    long nread=0, f_aut=0, f_ext=0, f_inv=0, tested=0, hits=0; long long evals=0; char line[4096];
    uint16_t Sout[MAXN], Sin[MAXN];
    while (fgets(line,sizeof line,stdin)) {
        int vals[400], nv=0; char *p=line,*e; while(1){ long v=strtol(p,&e,10); if(e==p) break; vals[nv++]=(int)v; p=e; }
        if (nv<3) continue; nread++;
        int ne=vals[1]; long grp=vals[2]; memset(Sout,0,sizeof Sout); memset(Sin,0,sizeof Sin);
        for(int k=0;k<ne;k++){ int u=vals[3+2*k], v=vals[4+2*k]; Sout[u]|=1u<<v; Sin[v]|=1u<<u; }
        int iso=0; for(int v=0;v<nn;v++) if(!(Sout[v]|Sin[v])) iso++;
        int autodd = (grp%2==1 && iso<=1);
        static uint64_t le[1<<MAXN]; uint32_t full=(1u<<nn)-1; for(uint32_t m=0;m<=full;m++) le[m]=0; le[0]=1;
        for(uint32_t m=0;m<full;m++){ if(!le[m]) continue; for(int v=0;v<nn;v++) if(!(m>>v&1) && (Sin[v]&~m)==0) le[m|(1u<<v)]+=le[m]; }
        for(int v=0;v<nn;v++){ gadj[v]=Sout[v]|Sin[v]; pinv[v]=-1; }
        int hasinv = (ne==0) || inv_rec(0,0);
        if (!rigid) {
            if (!autodd) { f_aut++; continue; }
            if (!(le[full]&1)) { f_ext++; continue; }
            if (!hasinv) { f_inv++; continue; }
        }
        int ordv[MAXN], pos[MAXN]; uint16_t pm=0;
        for (int k=0;k<nn;k++){ int best=-1, bs=-1; for(int v=0;v<nn;v++) if(!(pm>>v&1)){ int s=__builtin_popcount(gadj[v]&pm)*16+__builtin_popcount(gadj[v]); if(s>bs){bs=s;best=v;} } ordv[k]=best; pm|=1u<<best; }
        for(int k=0;k<nn;k++) pos[ordv[k]]=k;
        for (int k=0;k<nn;k++){ predo[k]=0; predi[k]=0; int v=ordv[k];
            for(int w=0;w<nn;w++){ if(pos[w]>=k) continue; if(Sout[w]>>v&1) predo[k]|=1u<<pos[w]; if(Sout[v]>>w&1) predi[k]|=1u<<pos[w]; } }
        tested++;
        uint32_t target = rigid ? (uint32_t)(le[full]&1) : 1u;  // phi(TT) = e(S) mod 2
        int ok=1;
        for (int q=0;q<bncls;q++){ int c=ord[q]; uint32_t par=emb_parity(c); evals++;
            if (par!=target) { ok=0; if(q){ int t=ord[q]; memmove(ord+1,ord,q*sizeof(int)); ord[0]=t; } break; } }
        if (rigid) { printf("%s val=%u autodd=%d Uinv=%d n=%d m=%d grp=%ld ext=%llu arcs", ok?"RIGID":"NONRIGID", target, autodd, hasinv, nn, ne, grp, (unsigned long long)le[full]);
            for(int k=0;k<ne;k++) printf(" %d>%d",vals[3+2*k],vals[4+2*k]); printf("\n"); if (ok) hits++; }
        else if (ok) { hits++; printf("REDEI n=%d m=%d grp=%ld ext=%llu arcs",nn,ne,grp,(unsigned long long)le[full]);
            for(int k=0;k<ne;k++) printf(" %d>%d",vals[3+2*k],vals[4+2*k]); printf("\n"); }
    }
    printf("BTSUMMARY mode=%s n=%d read=%ld filtered_aut=%ld filtered_ext=%ld filtered_noinv=%ld tested=%ld hits=%ld evals=%lld\n",
        rigid?"rigid":"redei", nn, nread, f_aut, f_ext, f_inv, tested, hits, evals);
    return 0;
}

// ======================= witness =======================
static uint64_t rs=88172645463325252ULL; static uint64_t xrnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }
static int w_n, need_at[MAXN][MAXN], nneed[MAXN]; static uint32_t w_out[MAXN]; static int w_v[MAXN]; static uint64_t w_cnt;
static void w_dfs(int d, uint32_t used){
    if(d==w_n){ w_cnt++; return; }
    uint32_t cand = d==0 ? ((1u<<w_n)-1) : (w_out[w_v[d-1]]&~used);
    for(int k=0;k<nneed[d];k++){ int i=need_at[d][k]; cand &= w_out[w_v[i]]; }
    while(cand){ int w=__builtin_ctz(cand); cand&=cand-1; w_v[d]=w; w_dfs(d+1,used|(1u<<w)); }
}
static int mode_witness(int nn, long seed, long tries, int argc, char **argv) {
    w_n=nn; memset(nneed,0,sizeof nneed); rs ^= (uint64_t)seed*0x9E3779B97F4A7C15ULL;
    for(int a=0;a<argc;a++){ int i,j; sscanf(argv[a],"%d-%d",&i,&j); need_at[j][nneed[j]++]=i; }
    for (long s=0;s<tries;s++){
        memset(w_out,0,sizeof w_out); char str[128]; int b=0;
        for(int i=0;i<nn;i++) for(int j=i+1;j<nn;j++){ if(xrnd()&1){ w_out[i]|=1u<<j; str[b]='1'; } else { w_out[j]|=1u<<i; str[b]='0'; } b++; }
        str[b]=0; w_cnt=0; w_dfs(0,0);
        if (!(w_cnt&1)) { printf("WITNESS n=%d count=%llu tour=%s\n",nn,(unsigned long long)w_cnt,str); return 0; }
    }
    printf("NOWITNESS n=%d tries=%ld\n",nn,tries); return 0;
}

// ======================= unavoidability, n <= 7 =======================
static int mode_unav(int nn, const char *clsfile, int e) {
    uint32_t (*REP)[MAXN]; int ncls=read_classes_file(clsfile,nn,&REP);
    int16_t *tab=build_class_table(nn,REP,ncls); int P=nn*(nn-1)/2;
    uint8_t *hit=malloc(ncls); long total=0, good=0;
    // enumerate e-subsets of the P forward pairs (Gosper)
    if (e>P) { printf("UNAV n=%d e=%d sets=0 shavings=0\n",nn,e); return 0; }
    uint32_t m = e ? ((1u<<e)-1) : 0;
    while (1) {
        total++;
        int freeq[64], nf=0; for(int q=0;q<P;q++) if(!(m>>q&1)) freeq[nf++]=q;
        memset(hit,0,ncls); int nh=0; uint32_t code=m; if(!hit[tab[code]]){hit[tab[code]]=1;nh++;}
        for (uint64_t g=1; g<(1ULL<<nf) && nh<ncls; g++){ int k=__builtin_ctzll(g); code^=1u<<freeq[k]; if(!hit[tab[code]]){hit[tab[code]]=1;nh++;} }
        if (nh==ncls) { good++; printf("U %u\n",m); }
        if (e==0) break;
        uint32_t c = m & -m, r = m + c; m = (((r ^ m) >> 2) / c) | r;
        if (m >> P) break;
    }
    printf("UNAV n=%d e=%d sets=%ld shavings=%ld\n",nn,e,total,good);
    return 0;
}

// ======================= embedall: does S embed in every class? (plain recursive search, index order) ======
static int ea_n; static uint32_t ea_T[MAXN]; static uint32_t ea_Sout[MAXN], ea_Sin[MAXN]; static int ea_img[MAXN];
static int ea_rec(int k, uint32_t used) {
    if (k==ea_n) return 1;
    for (int y=0;y<ea_n;y++) { if (used>>y&1) continue; int ok=1;
        for (int j=0;j<k && ok;j++) { if ((ea_Sout[j]>>k&1) && !(ea_T[ea_img[j]]>>y&1)) ok=0; if ((ea_Sout[k]>>j&1) && !(ea_T[y]>>ea_img[j]&1)) ok=0; }
        if (ok) { ea_img[k]=y; if (ea_rec(k+1,used|(1u<<y))) return 1; } }
    return 0;
}
static int mode_embedall(int nn, const char *clsfile, int argc, char **argv) {
    ea_n=nn; memset(ea_Sout,0,sizeof ea_Sout); memset(ea_Sin,0,sizeof ea_Sin);
    for (int a=0;a<argc;a++){ int u,v; if (sscanf(argv[a],"%d>%d",&u,&v)==2){ ea_Sout[u]|=1u<<v; ea_Sin[v]|=1u<<u; } }
    FILE *f=fopen(clsfile,"r"); char line[256]; long ncls=0, avoid=0;
    while (fgets(line,sizeof line,f)) { if(!parse_tour(line,nn,ea_T)) { fprintf(stderr,"bad\n"); return 1; } ncls++;
        if (!ea_rec(0,0)) { avoid++; if (avoid<=5) { int l=(int)strlen(line); if(l&&line[l-1]=='\n') line[l-1]=0; printf("AVOIDER %s\n",line); } } }
    fclose(f);
    printf("EMBEDALL n=%d classes=%ld avoiders=%ld\n",nn,ncls,avoid);
    return 0;
}

// ======================= dcheck: exact identities behind Theorem D1 =======================
// D_n = P_n + {(0,n-2),(1,n-1)}.  A = #HP with v_{n-2}->v_0, B = #HP with v_{n-1}->v_1, AB = both.
// Checks: emb(D_n) = H - A - B + AB; A = sum_w hc(T-w) d^-(w); B = sum_w hc(T-w) d^+(w); AB even.
static int dc_n; static uint32_t dc_out[MAXN]; static int dc_v[MAXN]; static uint64_t dcH, dcA, dcB, dcAB, dcE;
static void dc_dfs(int d, uint32_t used){
    if (d==dc_n){ int n=dc_n; int a = (dc_out[dc_v[n-2]]>>dc_v[0])&1, b=(dc_out[dc_v[n-1]]>>dc_v[1])&1;
        dcH++; if(a) dcA++; if(b) dcB++; if(a&&b) dcAB++; if(!a&&!b) dcE++; return; }
    uint32_t cand = d==0 ? ((1u<<dc_n)-1) : (dc_out[dc_v[d-1]]&~used);
    while(cand){ int w=__builtin_ctz(cand); cand&=cand-1; dc_v[d]=w; dc_dfs(d+1,used|(1u<<w)); }
}
static uint64_t hc_sub(const uint32_t *out, int nn, uint32_t mask) {   // Hamiltonian cycles of T[mask]
    int k=__builtin_popcount(mask); if (k<3) return 0;
    int vs[MAXN], m=0; for(int v=0;v<nn;v++) if(mask>>v&1) vs[m++]=v;
    uint32_t o[MAXN]; for(int i=0;i<m;i++){ o[i]=0; for(int j=0;j<m;j++) if(out[vs[i]]>>vs[j]&1) o[i]|=1u<<j; }
    uint64_t P[MAXN]; count_from(o,m,0,P); uint64_t hc=0; for(int t=1;t<m;t++) if(o[t]&1u) hc+=P[t]; return hc;
}
static int mode_dcheck(int nn) {
    dc_n=nn; char line[256]; long ncls=0, bad=0, nodd=0, formula_ok=0;
    while (fgets(line,sizeof line,stdin)) {
        if(!parse_tour(line,nn,dc_out)){ fprintf(stderr,"bad\n"); return 1; }
        dcH=dcA=dcB=dcAB=dcE=0; dc_dfs(0,0);
        uint64_t sA=0, sB=0, shc=0; uint32_t full=(1u<<nn)-1;
        for (int w=0;w<nn;w++){ uint64_t h=hc_sub(dc_out,nn,full&~(1u<<w)); int dplus=__builtin_popcount(dc_out[w]); int dminus=nn-1-dplus;
            sA+=h*(uint64_t)dminus; sB+=h*(uint64_t)dplus; shc+=h; }
        int ok = (dcE == dcH - dcA - dcB + dcAB) && (dcA==sA) && (dcB==sB) && !(dcAB&1);
        if (!ok) bad++;
        if (dcE&1) nodd++;
        if ((dcE&1) == ((1+shc)&1)) formula_ok++;
        ncls++;
    }
    printf("DCHECK n=%d classes=%ld identity_failures=%ld emb_odd=%ld matches_1+sum_hc(T-w)=%ld\n",nn,ncls,bad,nodd,formula_ok);
    return 0;
}

int main(int argc, char **argv) {
    if (argc<2) { fprintf(stderr,"usage: see header\n"); return 1; }
    if (!strcmp(argv[1],"dcheck")) return mode_dcheck(atoi(argv[2]));
    if (!strcmp(argv[1],"embedall")) return mode_embedall(atoi(argv[2]), argv[3], argc-4, argv+4);
    if (!strcmp(argv[1],"hp")) return mode_hp();
    if (!strcmp(argv[1],"til")) return mode_til(atoi(argv[2]));
    if (!strcmp(argv[1],"chords")) return mode_chords(atoi(argv[2]), argc-3, argv+3);
    if (!strcmp(argv[1],"dag")) return mode_dag(atoi(argv[2]), argv[3]);
    if (!strcmp(argv[1],"redei")) return mode_backtrack(atoi(argv[2]), argv[3], 0);
    if (!strcmp(argv[1],"rigid")) return mode_backtrack(atoi(argv[2]), argv[3], 1);
    if (!strcmp(argv[1],"witness")) return mode_witness(atoi(argv[2]), atol(argv[3]), atol(argv[4]), argc-5, argv+5);
    if (!strcmp(argv[1],"unav")) return mode_unav(atoi(argv[2]), argv[3], atoi(argv[4]));
    fprintf(stderr,"unknown mode\n"); return 1;
}

// procgen_shave4_20261001_u9.c -- shave4 lane, 2026-10-01.  Exact search for K-arc shavings on 9 vertices.
// Every shaving S (spanning oriented graph contained in every 9-tournament) embeds in the host H = C3[C3]
// (|Aut H| = 81).  So it suffices to enumerate arc subsets of H, up to Aut(H), by orderly generation
// (keep a set iff it is lexicographically minimal in its Aut(H)-orbit; children add arcs of larger index).
// Branches are pruned when the set has a directed cycle or fails to embed in some tournament of a POOL
// (failures are inherited by supersets).  At depth K every surviving set is tested against all classes;
// a class that kills it is added to the pool.
// usage: u9 classesN.txt K [host_string|-] [count_depth]  < initial pool (class strings, may be empty)
//   compile with -DN=7 / -DN=8 for other orders (default N = 9); host defaults to C3[C3] (N = 9 only).
//   count_depth: also full-check every surviving set of that size and print the shavings found there.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifndef N
#define N 9
#endif
#define NA (N*(N-1)/2)
#define MAXPOOL 4000
#define MAXAUT 200
static int ncls; static uint16_t (*CO)[N], (*CI)[N];
static uint16_t HO[N];                 // host adjacency (out masks)
static int au[NA], av[NA];             // host arcs, index -> (u,v)
static int naut; static uint8_t gmap[MAXAUT][NA];
static int npool; static int poolidx[MAXPOOL];   // pool = class indices
static int *ord;
static long long nodes[NA+2], fullchecks=0, shavings_at_count=0, searches=0;
static int K, CD;

static int embed_search(int c, const uint16_t *Sout, uint8_t *pi) {
    // find bijection pi: S -> class c with all S arcs mapped to arcs; vertex order by degree
    int ordv[N]; uint16_t deg[N]; uint16_t Sin[N]={0};
    for(int u=0;u<N;u++) for(int v=0;v<N;v++) if(Sout[u]>>v&1) Sin[v]|=1<<u;
    uint16_t placed=0;
    for(int k=0;k<N;k++){ int best=-1,bs=-1; for(int v=0;v<N;v++) if(!(placed>>v&1)){ int s=__builtin_popcount((Sout[v]|Sin[v])&placed)*32+__builtin_popcount(Sout[v]|Sin[v]); if(s>bs){bs=s;best=v;} } ordv[k]=best; placed|=1<<best; (void)deg; }
    int pos[N]; for(int k=0;k<N;k++) pos[ordv[k]]=k;
    uint16_t po[N], pn[N];
    for(int k=0;k<N;k++){ po[k]=pn[k]=0; int v=ordv[k]; for(int w=0;w<N;w++){ if(pos[w]>=k) continue; if(Sout[w]>>v&1) po[k]|=1<<pos[w]; if(Sout[v]>>w&1) pn[k]|=1<<pos[w]; } }
    int img[N]; uint16_t cand[N]; int k=0; uint16_t used=0, full=(1<<N)-1; cand[0]=full; searches++;
    while(1){
        if(!cand[k]){ if(k==0) return 0; k--; used&=~(1<<img[k]); continue; }
        int y=__builtin_ctz(cand[k]); cand[k]&=cand[k]-1; img[k]=y;
        if(k==N-1){ for(int q=0;q<N;q++) pi[ordv[q]]=(uint8_t)img[q]; return 1; }
        used|=1<<y; k++;
        uint16_t cc=full&~used; uint16_t p=po[k]; while(p){ int j=__builtin_ctz(p); p&=p-1; cc&=CO[c][img[j]]; }
        p=pn[k]; while(p){ int j=__builtin_ctz(p); p&=p-1; cc&=CI[c][img[j]]; }
        cand[k]=cc;
    }
}
#define KC 4   // cached embeddings per pool tournament and depth
static int embed_search_multi(int c, const uint16_t *Sout, uint8_t (*outl)[N], int maxk) {
    // collect up to maxk embeddings of S into class c (exact backtracking); returns the number found
    int ordv[N]; uint16_t Sin[N]={0};
    for(int u=0;u<N;u++) for(int v=0;v<N;v++) if(Sout[u]>>v&1) Sin[v]|=1<<u;
    uint16_t placed=0;
    for(int k=0;k<N;k++){ int best=-1,bs=-1; for(int v=0;v<N;v++) if(!(placed>>v&1)){ int s2=__builtin_popcount((Sout[v]|Sin[v])&placed)*32+__builtin_popcount(Sout[v]|Sin[v]); if(s2>bs){bs=s2;best=v;} } ordv[k]=best; placed|=1<<best; }
    int pos[N]; for(int k=0;k<N;k++) pos[ordv[k]]=k;
    uint16_t po[N], pn[N];
    for(int k=0;k<N;k++){ po[k]=pn[k]=0; int v=ordv[k]; for(int w=0;w<N;w++){ if(pos[w]>=k) continue; if(Sout[w]>>v&1) po[k]|=1<<pos[w]; if(Sout[v]>>w&1) pn[k]|=1<<pos[w]; } }
    int img[N]; uint16_t cand[N]; int k=0; uint16_t used=0, full=(1<<N)-1; cand[0]=full; int nf=0; searches++;
    while(1){
        if(!cand[k]){ if(k==0) return nf; k--; used&=~(1<<img[k]); continue; }
        int y=__builtin_ctz(cand[k]); cand[k]&=cand[k]-1; img[k]=y;
        if(k==N-1){ for(int q=0;q<N;q++) outl[nf][ordv[q]]=(uint8_t)img[q]; nf++; if(nf>=maxk) return nf; continue; }
        used|=1<<y; k++;
        uint16_t cc=full&~used; uint16_t p=po[k]; while(p){ int j=__builtin_ctz(p); p&=p-1; cc&=CO[c][img[j]]; }
        p=pn[k]; while(p){ int j=__builtin_ctz(p); p&=p-1; cc&=CI[c][img[j]]; }
        cand[k]=cc;
    }
}
static int PARITY=0;   // if set (env SHAVE4_PARITY=1): at depth K test the Redei property (odd count in every class)
static int emb_parity_count(int c, const uint16_t *Sout) {
    int ordv[N]; uint16_t Sin[N]={0};
    for(int u=0;u<N;u++) for(int v=0;v<N;v++) if(Sout[u]>>v&1) Sin[v]|=1<<u;
    uint16_t placed=0;
    for(int k=0;k<N;k++){ int best=-1,bs=-1; for(int v=0;v<N;v++) if(!(placed>>v&1)){ int s2=__builtin_popcount((Sout[v]|Sin[v])&placed)*32+__builtin_popcount(Sout[v]|Sin[v]); if(s2>bs){bs=s2;best=v;} } ordv[k]=best; placed|=1<<best; }
    int pos[N]; for(int k=0;k<N;k++) pos[ordv[k]]=k;
    uint16_t po[N], pn[N];
    for(int k=0;k<N;k++){ po[k]=pn[k]=0; int v=ordv[k]; for(int w=0;w<N;w++){ if(pos[w]>=k) continue; if(Sout[w]>>v&1) po[k]|=1<<pos[w]; if(Sout[v]>>w&1) pn[k]|=1<<pos[w]; } }
    int img[N]; uint16_t cand[N]; int k=0; uint16_t used=0, full=(1<<N)-1; cand[0]=full; int par=0;
    while(1){
        if(!cand[k]){ if(k==0) return par; k--; used&=~(1<<img[k]); continue; }
        int y=__builtin_ctz(cand[k]); cand[k]&=cand[k]-1; img[k]=y;
        if(k==N-1){ par^=1; continue; }
        used|=1<<y; k++;
        uint16_t cc=full&~used; uint16_t p=po[k]; while(p){ int j=__builtin_ctz(p); p&=p-1; cc&=CO[c][img[j]]; }
        p=pn[k]; while(p){ int j=__builtin_ctz(p); p&=p-1; cc&=CI[c][img[j]]; }
        cand[k]=cc;
    }
}
static long long filtered_ext=0, filtered_inv=0;
static uint16_t ga[N]; static int gp[N];
static int inv_rec9(int start, int nontriv) {   // does the undirected graph ga have an automorphism of order 2?
    while (start<N && gp[start]>=0) start++;
    if (start==N) { if (!nontriv) return 0;
        for (int u=0;u<N;u++) for (int v=u+1;v<N;v++) { int a=(ga[u]>>v)&1, b=(ga[gp[u]]>>gp[v])&1; if(a!=b) return 0; }
        return 1; }
    gp[start]=start; if (inv_rec9(start+1,nontriv)) return 1; gp[start]=-1;
    for (int w=start+1; w<N; w++) if (gp[w]<0 && __builtin_popcount(ga[w])==__builtin_popcount(ga[start])) {
        gp[start]=w; gp[w]=start; if (inv_rec9(start+1,1)) return 1; gp[start]=-1; gp[w]=-1; }
    return 0;
}
static int redei_check(const uint16_t *Sout) {
    // proven necessary conditions first (note B4): e(S) odd, and U(S) has an involution
    { uint16_t Sin[N]={0}; for(int u=0;u<N;u++) for(int v=0;v<N;v++) if(Sout[u]>>v&1) Sin[v]|=1<<u;
      static uint8_t le[1<<N]; uint32_t full=(1u<<N)-1; memset(le,0,sizeof le); le[0]=1;
      for(uint32_t m=0;m<full;m++){ if(!le[m]) continue; for(int v=0;v<N;v++) if(!(m>>v&1) && (Sin[v]&~m)==0) le[m|(1u<<v)]^=1; }
      if(!le[full]) { filtered_ext++; return 0; }
      int ne=0; for(int v=0;v<N;v++){ ga[v]=Sout[v]|Sin[v]; gp[v]=-1; ne+=__builtin_popcount(Sout[v]); }
      if (ne>0 && !inv_rec9(0,0)) { filtered_inv++; return 0; } }
    fullchecks++;
    for(int q=0;q<ncls;q++){ int c=ord[q];
        if(!emb_parity_count(c,Sout)){ if(q){ int t=ord[q]; memmove(ord+1,ord,q*sizeof(int)); ord[0]=t; }
            uint8_t pi[N]; if(!embed_search(c,Sout,pi)){ int in=0; for(int p=0;p<npool;p++) if(poolidx[p]==c) in=1; if(!in && npool<MAXPOOL) poolidx[npool++]=c; }
            return 0; } }
    return 1;
}
static int full_check(const uint16_t *Sout) {
    if (PARITY) return redei_check(Sout);
    uint8_t pi[N]; fullchecks++;
    for(int q=0;q<ncls;q++){ int c=ord[q];
        if(!embed_search(c,Sout,pi)){ if(q){ int t=ord[q]; memmove(ord+1,ord,q*sizeof(int)); ord[0]=t; }
            // add killer to pool
            int in=0; for(int p=0;p<npool;p++) if(poolidx[p]==c) in=1;
            if(!in && npool<MAXPOOL) poolidx[npool++]=c;
            return 0; } }
    return 1;
}
static uint8_t PI[NA+2][MAXPOOL][KC][N];  // cached embeddings per depth (up to KC per pool tournament)
static uint8_t PIN[NA+2][MAXPOOL];        // number of valid cached embeddings (0 = none, search needed)
static uint64_t GS[NA+2][MAXAUT];            // images g(S) per depth
static int lexless(uint64_t a, uint64_t b){ uint64_t d=a^b; if(!d) return 0; return (a & (d & -d))!=0; }
static void print_set(const char *tag, uint64_t S){ printf("%s %d arcs:",tag,__builtin_popcountll(S)); for(int a=0;a<NA;a++) if(S>>a&1) printf(" %d>%d",au[a],av[a]); printf("\n"); fflush(stdout); }
static long found=0;
static int PMIN=-1; static long long redei_at[NA+2];   // env SHAVE4_PARITY_MIN=d: Redei test at every depth >= d
static void dfs(int depth, uint64_t S, int maxidx, uint16_t *Sout) {
    nodes[depth]++;
    if (PMIN>=0 && depth>=PMIN && depth<K) { if (redei_check(Sout)) { redei_at[depth]++; print_set("REDEI",S); } }
    if (depth==CD && CD!=K) { if(full_check(Sout)) { shavings_at_count++; print_set("COUNTED",S); } }
    if (depth==K) { if(full_check(Sout)) { found++; print_set(PARITY?"REDEI":"SHAVING",S); } return; }
    for (int a=maxidx+1;a<NA;a++){
        int u=au[a], v=av[a];
        // cycle check: does v reach u in S?
        uint16_t reach=1<<v, fr=1<<v; while(fr){ uint16_t nf=0; for(int w=0;w<N;w++) if(fr>>w&1) nf|=Sout[w]; nf&=~reach; reach|=nf; fr=nf; }
        if (reach>>u&1) continue;
        uint64_t S2=S|(1ULL<<a);
        // canonicity under Aut(H)
        int canon=1;
        for(int g=1; g<naut; g++){ uint64_t im=GS[depth][g]|(1ULL<<gmap[g][a]); GS[depth+1][g]=im; if(lexless(im,S2)){ canon=0; break; } }
        if(!canon) continue;
        uint16_t So2[N]; memcpy(So2,Sout,sizeof So2); So2[u]|=1<<v;
        // pool check with cached embeddings
        int ok=1;
        for(int p=0;p<npool && ok;p++){ int c=poolidx[p]; int cnt=PIN[depth][p], m=0;
            for(int q=0;q<cnt;q++){ uint8_t *e=PI[depth][p][q]; if (CO[c][e[u]]>>e[v]&1) { memcpy(PI[depth+1][p][m],e,N); m++; } }
            if (m==0) {
                m=embed_search_multi(c,So2,PI[depth+1][p],KC);
                if (m==0) { ok=0; break; }
                // the new embeddings also embed the parent: keep one for the siblings
                if (cnt<KC) { memcpy(PI[depth][p][cnt],PI[depth+1][p][0],N); PIN[depth][p]=(uint8_t)(cnt+1); }
            }
            PIN[depth+1][p]=(uint8_t)m; }
        if(!ok) continue;
        dfs(depth+1,S2,a,So2);
    }
}
int main(int argc,char**argv){
    FILE*f=fopen(argv[1],"r"); char line[64]; int cap=200000; CO=malloc(sizeof(*CO)*cap); CI=malloc(sizeof(*CI)*cap);
    while(fgets(line,sizeof line,f)){ memset(CO[ncls],0,sizeof CO[ncls]); memset(CI[ncls],0,sizeof CI[ncls]); int q=0;
        for(int i=0;i<N;i++) for(int j=i+1;j<N;j++){ if(line[q]=='1'){ CO[ncls][i]|=1<<j; CI[ncls][j]|=1<<i; } else { CO[ncls][j]|=1<<i; CI[ncls][i]|=1<<j; } q++; }
        ncls++; }
    fclose(f);
    K=atoi(argv[2]); CD = argc>4 ? atoi(argv[4]) : -1;
    PARITY = getenv("SHAVE4_PARITY") && atoi(getenv("SHAVE4_PARITY"));
    if (getenv("SHAVE4_PARITY_MIN")) PMIN=atoi(getenv("SHAVE4_PARITY_MIN"));
    ord=malloc(sizeof(int)*ncls); for(int i=0;i<ncls;i++) ord[i]=i;
    if (argc>3 && strcmp(argv[3],"-")) {   // host given as upper-triangle string
        const char *h=argv[3]; int q=0; for(int i=0;i<N;i++) HO[i]=0;
        for(int i=0;i<N;i++) for(int j=i+1;j<N;j++){ if(h[q]=='1') HO[i]|=1<<j; else HO[j]|=1<<i; q++; }
    } else {
        if (N!=9) { fprintf(stderr,"default host needs N=9\n"); return 1; }
        // host C3[C3]: vertex i = (block i/3, pos i%3); (b,x)->(b',x') iff b'=b+1 mod 3, or b'=b and x'=x+1 mod 3
        for(int i=0;i<N;i++){ HO[i]=0; for(int j=0;j<N;j++){ if(i==j) continue; int b=i/3,x=i%3,c=j/3,y=j%3; if(c==(b+1)%3 || (b==c && y==(x+1)%3)) HO[i]|=1<<j; } }
    }
    int na=0; for(int u=0;u<N;u++) for(int v=0;v<N;v++) if(HO[u]>>v&1){ au[na]=u; av[na]=v; na++; }
    if(na!=NA){ fprintf(stderr,"host arcs %d\n",na); return 1; }
    // automorphisms (all 9! permutations)
    int p[N]; for(int i=0;i<N;i++) p[i]=i; naut=0;
    int c_[N]={0}; int it=0;
    while(1){
        int ok=1; for(int u=0;u<N && ok;u++) for(int v=0;v<N;v++) if((HO[u]>>v&1) && !(HO[p[u]]>>p[v]&1)){ ok=0; break; }
        if(ok){ if(naut>=MAXAUT){ fprintf(stderr,"too many aut\n"); return 1; }
            for(int a=0;a<NA;a++){ int gu=p[au[a]], gv=p[av[a]]; int idx=-1; for(int b=0;b<NA;b++) if(au[b]==gu&&av[b]==gv) idx=b; gmap[naut][a]=(uint8_t)idx; }
            naut++; }
        // next permutation (Heap)
        int i=0; while(i<N && c_[i]>=i){ c_[i]=0; i++; }
        if(i==N) break;
        if(i%2==0){ int t=p[0]; p[0]=p[i]; p[i]=t; } else { int t=p[c_[i]]; p[c_[i]]=p[i]; p[i]=t; }
        c_[i]++; it++;
    }
    // make sure identity is element 0
    { int idg=-1; for(int g=0;g<naut;g++){ int isid=1; for(int a=0;a<NA;a++) if(gmap[g][a]!=a) isid=0; if(isid) idg=g; }
      if(idg!=0){ uint8_t t[NA]; memcpy(t,gmap[0],NA); memcpy(gmap[0],gmap[idg],NA); memcpy(gmap[idg],t,NA); } }
    printf("N=%d host arcs=%d |Aut(host)|=%d\n",N,NA,naut);
    // initial pool: read from stdin (class strings) -> find indices by direct match against classes
    npool=0; char buf[128];
    while(fgets(buf,sizeof buf,stdin)){ int l=strlen(buf); while(l&&(buf[l-1]=='\n'||buf[l-1]==' ')) buf[--l]=0; if(!l) continue;
        uint16_t o[N]={0}; int q=0; for(int i=0;i<N;i++) for(int j=i+1;j<N;j++){ if(buf[q]=='1') o[i]|=1<<j; else o[j]|=1<<i; q++; }
        // locate class index by exact string match
        int found_i=-1; for(int c=0;c<ncls && found_i<0;c++){ if(!memcmp(o,CO[c],sizeof o)) found_i=c; }
        if(found_i>=0 && npool<MAXPOOL) poolidx[npool++]=found_i; }
    printf("initial pool size %d\n",npool); fflush(stdout);
    // initial embeddings of the empty graph: identity
    for(int g=0;g<naut;g++) GS[0][g]=0;
    // caches: depth 0 (empty graph) -- KC cyclic shifts are embeddings; all deeper caches start empty (count 0)
    memset(PIN,0,sizeof PIN);
    for(int pp=0;pp<MAXPOOL;pp++){ for(int q=0;q<KC;q++) for(int v=0;v<N;v++) PI[0][pp][q][v]=(uint8_t)((v+q)%N); PIN[0][pp]=KC; }
    uint16_t S0[N]={0};
    dfs(0,0,-1,S0);
    printf("K=%d found=%ld fullchecks=%lld final_pool=%d searches=%lld\n",K,found,fullchecks,npool,searches);
    if (PARITY || PMIN>=0) printf("redei filters: e(S) even %lld, no involution of U(S) %lld\n",filtered_ext,filtered_inv);
    if(CD>=0) printf("count_depth=%d full-check survivors=%lld\n",CD,shavings_at_count);
    for(int d=0;d<=K;d++) printf("depth %d nodes %lld%s\n",d,nodes[d], (PMIN>=0 && d>=PMIN && d<K) ? "" : "");
    if (PMIN>=0) { printf("REDEI_BY_DEPTH"); for(int d=PMIN; d<K; d++) printf(" %d:%lld",d,redei_at[d]); printf("\n"); }
    printf("POOL"); for(int pp=0;pp<npool;pp++) printf(" %d",poolidx[pp]); printf("\n");
    return 0;
}

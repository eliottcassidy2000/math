// procgen_shave4_20261001_closure.c -- shave4 lane (collatz-procgen-20260922), 2026-10-01.
// Horn-closure test: which parity-rigid oriented graphs (n <= NMAX) are derivable from
//   free involutions (|Aut| even => phi = 0), directed paths (Redei), disjoint unions (Lucas),
//   and the deletion-reversal relation phi_S + phi_{rev_f S} = phi_{S-f} (two known => third).
// Reads cls<n>.txt (gentourng -q n) and og<n>.txt (geng -q n | directg -q -o -G) from the cwd;
// optional targets<n>.txt (arc lists like '0>1 1>2') get their derivation printed.
// usage: closure NMAX
// Closure test: is every parity-rigid oriented graph on <= NMAX vertices derivable from
//   (B1) |Aut S| even  => phi_S = 0
//   (B2) directed path P_k => phi = 1 (Redei), plus
//   (U)  disjoint unions: phi_{S1 u S2} = C(n,n1) c1 c2 when S1, S2 rigid with values c1, c2
//   (DR) deletion-reversal: phi_S + phi_{rev_f S} = phi_{S - f}  (two known => third known)
// Usage: clo NMAX  (expects files cls<n>.txt (gentourng) and og<n>.txt (directg -o -G) for n=2..NMAX)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MAXN 8
typedef struct { uint64_t code; int32_t idx; } KV;

static int NN;                       // current n
static int npair; static int pidx[MAXN][MAXN];

static uint64_t mix64(uint64_t x){ x^=x>>33; x*=0xff51afd7ed558ccdULL; x^=x>>33; x*=0xc4ceb9fe1a85ec53ULL; x^=x>>33; return x; }
static int cmpu64(const void*a,const void*b){ uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b; return x<y?-1:x>y; }

// canonical code of oriented graph given by out[] on n vertices: 2 bits per pair (01: low->high, 10: high->low)
static uint64_t best; static int ord_[MAXN]; static int cellstart[MAXN+1], ncell; static int cellv[MAXN];
static const uint8_t *G_out; static int G_n;
static uint64_t code_of_order(const int *ord, int n, const uint8_t *out) {
    uint64_t c=0; int b=0;
    for (int p=0;p<n;p++) for (int q=p+1;q<n;q++) { int u=ord[p], v=ord[q];
        if (out[u]>>v&1) c|=1ULL<<(2*b); else if (out[v]>>u&1) c|=2ULL<<(2*b); b++; }
    return c;
}
static int used_[MAXN];
static void perm_rec(int pos) {
    if (pos==G_n) { uint64_t c=code_of_order(ord_,G_n,G_out); if(c<best) best=c; return; }
    // which cell does pos belong to?
    int cl=0; while (!(cellstart[cl]<=pos && pos<cellstart[cl+1])) cl++;
    for (int k=cellstart[cl];k<cellstart[cl+1];k++){ int v=cellv[k]; if(used_[v]) continue; used_[v]=1; ord_[pos]=v; perm_rec(pos+1); used_[v]=0; }
}
static uint64_t canon(const uint8_t *out, int n) {
    uint8_t in[MAXN]={0}; for(int u=0;u<n;u++) for(int v=0;v<n;v++) if(out[u]>>v&1) in[v]|=1<<u;
    uint64_t inv[MAXN], nw[MAXN];
    for (int v=0;v<n;v++) inv[v]=mix64(0x1234+__builtin_popcount(out[v])*16+__builtin_popcount(in[v]));
    for (int r=0;r<n;r++) {
        for (int v=0;v<n;v++) {
            uint64_t a[MAXN]; int k=0; for(int w=0;w<n;w++) if(out[v]>>w&1) a[k++]=inv[w]; qsort(a,k,8,cmpu64);
            uint64_t h=mix64(inv[v]^0x9e37); for(int i=0;i<k;i++) h=mix64(h^a[i]);
            k=0; for(int w=0;w<n;w++) if(in[v]>>w&1) a[k++]=inv[w]; qsort(a,k,8,cmpu64);
            h=mix64(h^0x7777); for(int i=0;i<k;i++) h=mix64(h^(a[i]*3));
            nw[v]=h;
        }
        memcpy(inv,nw,sizeof inv);
    }
    // cells sorted by invariant
    int idx[MAXN]; for(int v=0;v<n;v++) idx[v]=v;
    for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) if (inv[idx[j]]<inv[idx[i]]) { int t=idx[i]; idx[i]=idx[j]; idx[j]=t; }
    ncell=0; cellstart[0]=0;
    for (int i=0;i<n;i++){ cellv[i]=idx[i]; if (i>0 && inv[idx[i]]!=inv[idx[i-1]]) { cellstart[++ncell]=i; } }
    cellstart[++ncell]=n;
    best=~0ULL; G_out=out; G_n=n; memset(used_,0,sizeof used_); perm_rec(0);
    return best;
}

// ----- per-n data -----
typedef struct {
    int n, ncls, ntour;
    uint64_t *code;      // canonical code per class
    int8_t *status;      // true status: -1 non-rigid, 0 rigid even, 1 rigid odd (Redei)
    int8_t *known;       // closure value or -1
    int8_t *reason;      // 1 aut-even, 2 path/forest, 3 union, 4 deletion-reversal
    int32_t *wit1, *wit2;// witnesses for reason 4 (class indices in same n); wit for union: unused
    int8_t *role;        // which of the triple was derived: 0 = S (from X,Y), 1 = X (from S,Y), 2 = Y (from S,X)
    int16_t *arcf;       // arc index (u*8+v) used in triple
    uint8_t (*out)[MAXN];// representative
    KV *sorted;
} Level;
static Level L[MAXN+1];

static int find_class(Level *lv, uint64_t code) {
    int lo=0, hi=lv->ncls-1;
    while (lo<=hi) { int mid=(lo+hi)/2; if (lv->sorted[mid].code==code) return lv->sorted[mid].idx; if (lv->sorted[mid].code<code) lo=mid+1; else hi=mid-1; }
    return -1;
}
static int cmpkv(const void*a,const void*b){ uint64_t x=((const KV*)a)->code,y=((const KV*)b)->code; return x<y?-1:x>y; }

static long binom_mod2(int n,int k){ return ((k & (n-k))==0) ? 1 : 0; } // Lucas: C(n,k) odd iff k & (n-k) == 0


static int8_t *printed[MAXN+1]; static int gid_counter=0; static int32_t *gid[MAXN+1];
static void print_graph(int n, const uint8_t *out){ int first=1; for(int u=0;u<n;u++) for(int v=0;v<n;v++) if(out[u]>>v&1){ printf("%s%d>%d", first?"":" ", u,v); first=0; } if(first) printf("(no arcs)"); }
static int derive(int n, int i){
    Level *lv=&L[n];
    if (!printed[n]) { printed[n]=calloc(lv->ncls,1); gid[n]=malloc(4*lv->ncls); }
    if (printed[n][i]) return gid[n][i];
    int a=-1,b=-1;
    if (lv->reason[i]==4){ a=derive(n,lv->wit1[i]); b=derive(n,lv->wit2[i]); }
    int me=++gid_counter; printed[n][i]=1; gid[n][i]=me;
    printf("  G%d [n=%d, phi=%d] = ", me, n, lv->known[i]); print_graph(n, lv->out[i]); printf("\n      by ");
    switch(lv->reason[i]){
      case 1: printf("free automorphism of even order (|Aut| even)\n"); break;
      case 2: printf("Redei (directed Hamiltonian path)\n"); break;
      case 3: printf("disjoint union of rigid components (Lucas)\n"); break;
      case 4: { int u=lv->arcf[i]/8, v=lv->arcf[i]%8;
        if (lv->role[i]==0) printf("DR on its arc %d>%d: phi = phi(G%d = minus arc) + phi(G%d = arc reversed)\n",u,v,a,b);
        else if (lv->role[i]==1) printf("DR: it is G%d minus arc %d>%d; phi = phi(G%d) + phi(G%d = that arc reversed)\n",a,u,v,a,b);
        else printf("DR: it is G%d with arc %d>%d reversed; phi = phi(G%d) + phi(G%d = minus that arc)\n",a,u,v,a,b);
        break; }
      default: printf("??\n");
    }
    return me;
}

int main(int argc, char **argv) {
    int NMAX=atoi(argv[1]);
    for (int n=1;n<=NMAX;n++) {
        Level *lv=&L[n]; lv->n=n; NN=n;
        npair=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ pidx[i][j]=pidx[j][i]=npair++; }
        // tournament classes and full class table
        char fn[64]; int ntour=0; uint8_t (*TR)[MAXN]=malloc(sizeof(*TR)*10000);
        if (n>=2) { sprintf(fn,"cls%d.txt",n); FILE*f=fopen(fn,"r"); char line[128];
            while(fgets(line,sizeof line,f)){ memset(TR[ntour],0,MAXN); int b=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ if(line[b]=='1') TR[ntour][i]|=1<<j; else TR[ntour][j]|=1<<i; b++; } ntour++; }
            fclose(f); }
        else { ntour=1; memset(TR[0],0,MAXN); }
        lv->ntour=ntour;
        uint32_t NT=1u<<npair; int16_t *tab=malloc(sizeof(int16_t)*NT); for(uint32_t t=0;t<NT;t++) tab[t]=-1;
        { int perm[MAXN]; for(int i=0;i<n;i++) perm[i]=i; int c_[MAXN]={0};
          #define FILL for(int cc=0;cc<ntour;cc++){ uint32_t code=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) if(TR[cc][perm[i]]>>perm[j]&1) code|=1u<<pidx[i][j]; tab[code]=cc; }
          FILL
          int i=0; while(i<n){ if(c_[i]<i){ if(i%2==0){int t=perm[0];perm[0]=perm[i];perm[i]=t;} else {int t=perm[c_[i]];perm[c_[i]]=perm[i];perm[i]=t;} FILL c_[i]++; i=0; } else { c_[i]=0; i++; } } }
        for(uint32_t t=0;t<NT;t++) if(tab[t]<0){ fprintf(stderr,"table hole\n"); return 1; }
        // read oriented graphs
        sprintf(fn,"og%d.txt",n); FILE*f=fopen(fn,"r"); if(!f){fprintf(stderr,"missing %s\n",fn);return 1;}
        int cap=2200000; lv->code=malloc(8*cap); lv->status=malloc(cap); lv->known=malloc(cap); lv->reason=malloc(cap);
        lv->wit1=malloc(4*cap); lv->wit2=malloc(4*cap); lv->role=malloc(cap); lv->arcf=malloc(2*cap); lv->out=malloc(sizeof(*lv->out)*cap);
        int8_t *autodd=malloc(cap);
        char line[2048]; int nc=0; uint8_t *par=malloc(ntour);
        while (fgets(line,sizeof line,f)) {
            int vals[200], nv=0; char *p=line,*e; while(1){ long v=strtol(p,&e,10); if(e==p) break; vals[nv++]=(int)v; p=e; }
            if (nv<3) continue;
            int ne=vals[1]; long grp=vals[2];
            uint8_t out[MAXN]={0}; uint32_t fixed=0, under=0;
            for(int k=0;k<ne;k++){ int u=vals[3+2*k], v=vals[4+2*k]; out[u]|=1<<v; if(n>=2){ int q=pidx[u][v]; under|=1u<<q; if(u<v) fixed|=1u<<q; } }
            int iso=0; for(int v=0;v<n;v++){ int deg=__builtin_popcount(out[v]); for(int w=0;w<n;w++) if(out[w]>>v&1) deg++; if(!deg) iso++; }
            memcpy(lv->out[nc],out,MAXN);
            lv->code[nc]=canon(out,n);
            autodd[nc] = (grp%2==1 && iso<=1);
            // true status
            if (n==1) lv->status[nc]=1;
            else {
                memset(par,0,ntour); int freeq[64], nf=0; for(int q=0;q<npair;q++) if(!(under>>q&1)) freeq[nf++]=q;
                uint32_t code=fixed; par[tab[code]]^=1;
                for (uint64_t g=1; g<(1ULL<<nf); g++){ int k=__builtin_ctzll(g); code^=1u<<freeq[k]; par[tab[code]]^=1; }
                int same=1; for(int c=1;c<ntour;c++) if(par[c]!=par[0]){same=0;break;}
                lv->status[nc] = same ? par[0] : -1;
            }
            nc++;
        }
        fclose(f); lv->ncls=nc;
        lv->sorted=malloc(sizeof(KV)*nc); for(int i=0;i<nc;i++){ lv->sorted[i].code=lv->code[i]; lv->sorted[i].idx=i; }
        qsort(lv->sorted,nc,sizeof(KV),cmpkv);
        for(int i=1;i<nc;i++) if(lv->sorted[i].code==lv->sorted[i-1].code){ fprintf(stderr,"canon collision n=%d\n",n); return 1; }
        // base cases
        long nrig=0, nred=0;
        for (int i=0;i<nc;i++){ lv->known[i]=-1; lv->reason[i]=0; if(lv->status[i]>=0) nrig++; if(lv->status[i]==1) nred++; }
        for (int i=0;i<nc;i++){
            const uint8_t *out=lv->out[i];
            if (!autodd[i]) { lv->known[i]=0; lv->reason[i]=1; continue; }
            // components (undirected)
            uint8_t adj[MAXN]={0}; for(int u=0;u<n;u++) for(int v=0;v<n;v++) if(out[u]>>v&1){ adj[u]|=1<<v; adj[v]|=1<<u; }
            uint8_t seen=0; int ncomp=0; uint8_t comp[MAXN];
            for(int s=0;s<n;s++) if(!(seen>>s&1)){ uint8_t c=1<<s, fr=1<<s; while(fr){ uint8_t nf2=0; for(int v=0;v<n;v++) if(fr>>v&1) nf2|=adj[v]; nf2&=~c; c|=nf2; fr=nf2; } seen|=c; comp[ncomp++]=c; }
            if (ncomp==1) {
                // directed path?
                int ok=1, ne=0; for(int v=0;v<n;v++){ ne+=__builtin_popcount(out[v]); int indeg=0; for(int w=0;w<n;w++) if(out[w]>>v&1) indeg++; if(__builtin_popcount(out[v])>1||indeg>1) ok=0; }
                if (ok && ne==n-1) { lv->known[i]=1; lv->reason[i]=2; }
                continue;
            }
            // union rule: all components known in their own sizes
            int val=1, rem=n, allk=1;
            for (int k=0;k<ncomp;k++){
                int sz=__builtin_popcount(comp[k]); int map[MAXN], m=0; for(int v=0;v<n;v++) if(comp[k]>>v&1) map[v]=m++;
                uint8_t sub[MAXN]={0}; for(int u=0;u<n;u++) if(comp[k]>>u&1) for(int v=0;v<n;v++) if((comp[k]>>v&1)&&(out[u]>>v&1)) sub[map[u]]|=1<<map[v];
                uint64_t cc=canon(sub,sz); int j=find_class(&L[sz],cc); if(j<0){fprintf(stderr,"comp not found\n");return 1;}
                if (L[sz].known[j]<0) { allk=0; break; }
                val &= L[sz].known[j]; val &= binom_mod2(rem,sz); rem-=sz;
            }
            if (allk) { lv->known[i]=(int8_t)val; lv->reason[i]=3; }
        }
        // closure by deletion-reversal; iterate over truly rigid classes only (sound + complete for all-rigid triples)
        int changed=1, passes=0; long contradictions=0;
        while (changed) { changed=0; passes++;
            for (int i=0;i<nc;i++){ if (lv->status[i]<0) continue;
                const uint8_t *out=lv->out[i];
                for (int u=0;u<n;u++) for (int v=0;v<n;v++) if (out[u]>>v&1) {
                    uint8_t X[MAXN], Y[MAXN]; memcpy(X,out,MAXN); memcpy(Y,out,MAXN);
                    X[u]&=~(1<<v); Y[u]&=~(1<<v); Y[v]|=1<<u;
                    int xi=find_class(lv,canon(X,n)), yi=find_class(lv,canon(Y,n));
                    if (xi<0||yi<0){ fprintf(stderr,"lookup fail\n"); return 1; }
                    int ks=lv->known[i], kx=lv->known[xi], ky=lv->known[yi];
                    int nk=(ks>=0)+(kx>=0)+(ky>=0);
                    if (nk==2) {
                        if (ks<0){ lv->known[i]=(int8_t)(kx^ky); lv->reason[i]=4; lv->role[i]=0; lv->wit1[i]=xi; lv->wit2[i]=yi; lv->arcf[i]=u*8+v; }
                        else if (kx<0){ lv->known[xi]=(int8_t)(ks^ky); lv->reason[xi]=4; lv->role[xi]=1; lv->wit1[xi]=i; lv->wit2[xi]=yi; lv->arcf[xi]=u*8+v; }
                        else { lv->known[yi]=(int8_t)(ks^kx); lv->reason[yi]=4; lv->role[yi]=2; lv->wit1[yi]=i; lv->wit2[yi]=xi; lv->arcf[yi]=u*8+v; }
                        changed=1;
                    } else if (nk==3) { if ((ks^kx^ky)!=0) contradictions++; }
                }
            }
        }
        long nexp=0, nexp_red=0, bad=0, unexp=0;
        for (int i=0;i<nc;i++){ if (lv->known[i]>=0){ if(lv->known[i]!=lv->status[i]) bad++; else { nexp++; if(lv->known[i]==1) nexp_red++; } } else if (lv->status[i]>=0) unexp++; }
        printf("n=%d classes=%d rigid=%ld (redei=%ld) explained=%ld (redei %ld) unexplained=%ld wrong=%ld contradictions=%ld passes=%d\n",
               n,nc,nrig,nred,nexp,nexp_red,unexp,bad,contradictions,passes);
        fflush(stdout);
        for (int i=0;i<nc;i++) if (lv->known[i]<0 && lv->status[i]>=0) {
            printf("  UNEXPLAINED n=%d status=%d arcs",n,lv->status[i]); for(int u=0;u<n;u++) for(int v=0;v<n;v++) if(lv->out[i][u]>>v&1) printf(" %d>%d",u,v); printf("\n"); }
        { char tf[64]; sprintf(tf,"targets%d.txt",n); FILE*ft=fopen(tf,"r");
          if (ft) { char tl[1024]; while(fgets(tl,sizeof tl,ft)){ uint8_t o[MAXN]={0}; char*p=tl; int u,v,k;
                if (tl[0]=='#'||tl[0]=='\n') continue;
                while(sscanf(p,"%d>%d%n",&u,&v,&k)==2){ o[u]|=1<<v; p+=k; }
                int ci=find_class(lv,canon(o,n)); printf("TARGET n=%d: %s",n,tl);
                if (ci<0) { printf("  not found\n"); continue; }
                printf("  true status %d, known %d\n", lv->status[ci], lv->known[ci]);
                if (lv->known[ci]>=0) { gid_counter=0; for(int q=1;q<=n;q++) if(printed[q]) memset(printed[q],0,L[q].ncls); derive(n,ci); } }
            fclose(ft); } }
        free(tab); free(TR); free(par); free(autodd);
    }
    return 0;
}

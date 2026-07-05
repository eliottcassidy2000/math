/* opus-2026-07-05-S84 -- level-3 class census v4: COLUMN-SOLVE-THEN-FILTER.
 *
 * Stage 1: enumerate all assignments (mu_1..mu_12) whose column-1 kills cover all 169
 *          column-1 cells (DFS over column-1 cells only; cores expanded over free
 *          positions).  Occupancy predicts ~13^12 * e^-22.8 ~ 3e3 survivors.
 * Stage 2: filter survivors through columns 2..12 (precomputed 169-bit column masks).
 * Output: the number of FULL level-3 covering refinements of the given level-2 class,
 *         exactly (no caps needed).
 *
 * usage: level3_colfilter k1 ... k12
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MOD3 2197
typedef struct { uint64_t w[3]; } cb;    /* 169-bit column bitset */

static inline void cb_or(cb *d, const cb *s){ d->w[0]|=s->w[0]; d->w[1]|=s->w[1]; d->w[2]|=s->w[2]; }
static inline int cb_get(const cb *b,int i){ return (b->w[i>>6]>>(i&63))&1; }
static inline void cb_set(cb *b,int i){ b->w[i>>6] |= 1ull<<(i&63); }
static cb CFULL;
static inline int cb_full(const cb *b){
    return b->w[0]==CFULL.w[0] && b->w[1]==CFULL.w[1] && b->w[2]==CFULL.w[2];
}
static inline int cb_pop(const cb *b){
    return __builtin_popcountll(b->w[0])+__builtin_popcountll(b->w[1])+__builtin_popcountll(b->w[2]);
}

/* KMC[b][r][mu]: column-b kill mask (169 heights h: cell a = b + 13*h) */
static cb KMC[13][13][13];
static int kap[13];

static long long stage1_nodes=0, survivors=0, covers=0;
static int mu_asn[13];
/* column-1 killer lists per height */
static int k1_r[169][32], k1_mu[169][32], k1_n[169];

static void check_full(void){
    /* all positions assigned: verify columns 2..12 */
    for(int b=2;b<13;b++){
        cb acc; memset(&acc,0,sizeof acc);
        for(int r=1;r<13;r++) cb_or(&acc,&KMC[b][r][mu_asn[r]]);
        if(!cb_full(&acc)) return;
    }
    covers++;
    if(covers<=5){
        printf("  cover: mu=[");
        for(int r=1;r<13;r++) printf("%d%s",mu_asn[r],r<12?",":"");
        printf("]\n");
    }
}

static void expand_free(int r){
    while(r<13 && mu_asn[r]>=0) r++;
    if(r>=13){ survivors++; check_full(); return; }
    for(int mu=0;mu<13;mu++){ mu_asn[r]=mu; expand_free(r+1); }
    mu_asn[r]=-1;
}

static void dfs1(cb covered, int nassigned){
    stage1_nodes++;
    if(cb_full(&covered)){ expand_free(1); return; }
    int nrem=12-nassigned;
    if(169-cb_pop(&covered) > 26*nrem) return;
    int best=-1,bestn=99;
    for(int h=0;h<169;h++){
        if(cb_get(&covered,h)) continue;
        int n=0;
        for(int k=0;k<k1_n[h];k++) if(mu_asn[k1_r[h][k]]<0) n++;
        if(n<bestn){ bestn=n; best=h; if(n<=1) break; }
    }
    if(best<0||bestn==0) return;
    for(int k=0;k<k1_n[best];k++){
        int r=k1_r[best][k], mu=k1_mu[best][k];
        if(mu_asn[r]>=0) continue;
        mu_asn[r]=mu;
        cb nc=covered; cb_or(&nc,&KMC[1][r][mu]);
        dfs1(nc,nassigned+1);
        mu_asn[r]=-1;
    }
}

int main(int argc,char**argv){
    if(argc<13){ fprintf(stderr,"need 12 kappa digits\n"); return 1; }
    for(int r=1;r<13;r++) kap[r]=atoi(argv[r]);
    memset(&CFULL,0,sizeof CFULL);
    for(int i=0;i<169;i++) cb_set(&CFULL,i);
    memset(k1_n,0,sizeof k1_n);
    for(int b=1;b<13;b++)
        for(int r=1;r<13;r++)
            for(int mu=0;mu<13;mu++){
                long long v = r + 13*kap[r] + 169*mu;
                memset(&KMC[b][r][mu],0,sizeof(cb));
                for(int h=0;h<169;h++){
                    long long a = b + 13*h;
                    long long x = (a*v)%MOD3;
                    if(x<=169 || x>=2028){
                        cb_set(&KMC[b][r][mu],h);
                        if(b==1){
                            k1_r[h][k1_n[h]]=r; k1_mu[h][k1_n[h]]=mu; k1_n[h]++;
                        }
                    }
                }
            }
    for(int r=1;r<13;r++) mu_asn[r]=-1;
    cb empty; memset(&empty,0,sizeof empty);
    dfs1(empty,0);
    printf("kappa=[");
    for(int r=1;r<13;r++) printf("%d%s",kap[r],r<12?",":"");
    printf("] col1_solutions=%lld stage1_nodes=%lld LEVEL3_COVERS=%lld\n",
           survivors,stage1_nodes,covers);
    return 0;
}

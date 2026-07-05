/* opus-2026-07-05-S84 -- HYP-4136: the level-3 class census (deciding the S83 conjecture).
 *
 * For a level-2 pattern kappa (12 digits base 13, position r=1..12), count level-3
 * refinements mu in (Z/13)^12 whose family class v_r = r + 13*kappa_r + 169*mu_r
 * COVERS all 2028 unit cells a mod 2197 (i.e. admits NO strict witness
 * dist(v a / 2197) >= 170/2197 <=> for every a some v_r has a*v_r mod 2197 in
 * +-[1,169]).  Conjecture (S83, HYP-4126): above a NON-shadow level-2 cover the count
 * is 0; above a shadow it is >= 1 (the true shadow refines: sanity).
 *
 * usage: level3_census k1 k2 ... k12   (kappa digits 0..12)
 * DFS: branch on the least-options uncovered cell; column-deficit prune.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MOD3 2197
#define NCELLS 2028
#define NW 32              /* 32 * 64 = 2048 >= 2028 bits */

typedef struct { uint64_t w[NW]; } bs;

static inline void bs_or(bs *d, const bs *s){ for(int i=0;i<NW;i++) d->w[i]|=s->w[i]; }
static inline int bs_get(const bs *b,int i){ return (b->w[i>>6]>>(i&63))&1; }
static inline void bs_set(bs *b,int i){ b->w[i>>6] |= 1ull<<(i&63); }
static inline int bs_full(const bs *b,const bs *full){
    for(int i=0;i<NW;i++) if((b->w[i]&full->w[i])!=full->w[i]) return 0;
    return 1;
}
static inline int popcount_and(const bs *a,const bs *b){
    int c=0; for(int i=0;i<NW;i++) c+=__builtin_popcountll(a->w[i]&b->w[i]); return c;
}

static int cellidx[MOD3];          /* a -> compact index, -1 if 13|a */
static int cellval[NCELLS];        /* index -> a */
static bs KM[13][13];              /* KM[r][mu], r=1..12 */
static bs COLMASK[13];             /* by a %% 13 */
static bs FULLMASK;
static int killer_r[NCELLS][26], killer_mu[NCELLS][26], killer_n[NCELLS];

static long long nodes=0, sols=0;
static long long node_cap, sol_cap;
static int assigned_mu[13];        /* -1 = unassigned */
static int opt_cnt[NCELLS];        /* # killers at unassigned positions (incremental) */
static int killers_at_r[13][NCELLS]; /* how many of cell i's killers sit at position r */

static void dfs(bs covered, int nassigned){
    if(nodes>node_cap || sols>=sol_cap) return;
    nodes++;
    if(bs_full(&covered,&FULLMASK)){ sols++; return; }
    int nrem = 12-nassigned;
    for(int b=1;b<13;b++){
        int have = popcount_and(&covered,&COLMASK[b]);
        if(169-have > 26*nrem) return;
    }
    /* least-options uncovered cell (incremental counts) */
    int best=-1, bestn=99;
    for(int i=0;i<NCELLS;i++){
        if(bs_get(&covered,i)) continue;
        int n=opt_cnt[i];
        if(n<bestn){ bestn=n; best=i; if(n<=1) break; }
    }
    if(best<0||bestn==0) return;
    for(int k=0;k<killer_n[best];k++){
        int r=killer_r[best][k], mu=killer_mu[best][k];
        if(assigned_mu[r]>=0) continue;
        assigned_mu[r]=mu;
        for(int i=0;i<NCELLS;i++) opt_cnt[i]-=killers_at_r[r][i];
        bs nc = covered; bs_or(&nc,&KM[r][mu]);
        dfs(nc,nassigned+1);
        for(int i=0;i<NCELLS;i++) opt_cnt[i]+=killers_at_r[r][i];
        assigned_mu[r]=-1;
        if(nodes>node_cap || sols>=sol_cap) return;
    }
}

int main(int argc,char**argv){
    if(argc<13){ fprintf(stderr,"need 12 kappa digits\n"); return 1; }
    int kap[13];
    for(int r=1;r<13;r++) kap[r]=atoi(argv[r]);
    node_cap = (argc>13)? atoll(argv[13]) : 4000000000LL;
    sol_cap  = (argc>14)? atoll(argv[14]) : 1000;

    int n=0;
    for(int a=1;a<MOD3;a++){
        if(a%13==0){ cellidx[a]=-1; continue; }
        cellidx[a]=n; cellval[n]=a; n++;
    }
    memset(&FULLMASK,0,sizeof FULLMASK);
    for(int i=0;i<NCELLS;i++) bs_set(&FULLMASK,i);
    memset(COLMASK,0,sizeof COLMASK);
    for(int i=0;i<NCELLS;i++) bs_set(&COLMASK[cellval[i]%13],i);
    memset(killer_n,0,sizeof killer_n);
    for(int r=1;r<13;r++)
        for(int mu=0;mu<13;mu++){
            long long v = r + 13*kap[r] + 169*mu;
            memset(&KM[r][mu],0,sizeof(bs));
            for(int i=0;i<NCELLS;i++){
                long long x = ((long long)cellval[i]*v)%MOD3;
                if(x<=169 || x>=2028){
                    bs_set(&KM[r][mu],i);
                    int idx=cellidx[cellval[i]];
                    killer_r[idx][killer_n[idx]]=r;
                    killer_mu[idx][killer_n[idx]]=mu;
                    killer_n[idx]++;
                }
            }
        }
    for(int r=1;r<13;r++) assigned_mu[r]=-1;
    memset(killers_at_r,0,sizeof killers_at_r);
    for(int i=0;i<NCELLS;i++){
        opt_cnt[i]=killer_n[i];
        for(int k=0;k<killer_n[i];k++) killers_at_r[killer_r[i][k]][i]++;
    }
    bs empty; memset(&empty,0,sizeof empty);
    dfs(empty,0);
    printf("kappa=[");
    for(int r=1;r<13;r++) printf("%d%s",kap[r],r<12?",":"");
    printf("] level3_covering_refinements=%lld nodes=%lld%s\n",
           sols,nodes,(nodes>node_cap)?" (NODE CAP HIT: undecided)":"");
    return 0;
}

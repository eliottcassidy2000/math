/* lrc_q50_parity_analysis_kps_S11.c -- kind-pasteur-2026-07-05-S11, HYP-4137.
 *
 * THE PARITY-SPLIT REFRAME of the Q50 crux.  Re-runs the mac-mini-S54 gap census
 * (exact profile F1-F5 below), and for EVERY full-profile survivor records WHICH
 * moduli q in [26,50] carry a 2/25-margin witness -- with special focus on the
 * EVEN COMPOSITE family Q2P = {26,34,38,46,50} = {2p : p in 13,17,19,23,25}.
 *
 * QUESTION (ii) of HYP-4137: does the 5-modulus family Q2P alone clear every
 * survivor?  If yes, the Q50 target shrinks from "some q <= 50" to "some q in Q2P".
 *
 * Profile (identical to lrc_gap_census40_S54.c, MISTAKE-108 corrected F4):
 *   F1 covering: every q in 2..12 divides some element
 *   F2 spread:   2*w_max > 23*w_min
 *   F3 24-compression: w_max <= 24*w_2nd
 *   F4 big pair: 2*w_max >= 38
 *   F5 pinning: at every q<=25, some v==0 mod q OR residues hit every unit +-pair
 *   primitive: gcd = 1
 *
 * gcc -O2 -o q50p lrc_q50_parity_analysis_kps_S11.c && ./q50p 40
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int B = 40;
typedef long long ll;
static int gcd_i(int a, int b){ while(b){int t=a%b;a=b;b=t;} return a; }

static ll cnt_f5=0, cnt_hard=0;
/* clearing analysis */
static ll cleared_by_2p=0, not_cleared_by_2p=0;
static ll cleared_2p_only=0;         /* Q2P clears but NO odd q in [26,49] does */
static ll minq_hist[64];             /* min clearing q over [26,50] */
static ll q_clears[64];              /* how many survivors each q in [26,50] clears */
static ll no_clear_2650=0;           /* survivors with NO witness in [26,50] at all */
static int Q2P[5] = {26,34,38,46,50};
/* THE KEY SPLIT: pinned-only moduli (all prime-power factors <= 25, so the
 * witness depends ONLY on residues mod lcm(2..25) = HEIGHT-INDEPENDENT) vs free.
 * pinned-only in [26,50]: */
static int QPIN[16] = {26,28,30,33,34,35,36,38,39,40,42,44,45,46,48,50};
static ll cleared_by_pin=0, not_cleared_by_pin=0;   /* height-independent witness? */
static ll pin_q_clears[64];          /* per pinned-only q clear counts */
static ll pin_minq_hist[64];         /* min pinned-only clearing q */
/* Q0 evidence: smallest pinned-only witness modulus scanning [26,PMAX];
 * pinned-only = all prime-power factors <= 25 (<=> divides lcm(2..25)) */
#define PMAX 200
static int is_pin_only[PMAX+1];
static ll q0_max=0; static int q0_argmax[12];
static void init_pinonly(void){
    for (int q=2;q<=PMAX;q++){
        int n=q, ok=1;
        for (int p=2;(long)p*p<=n;p++){ if(n%p==0){ int pe=1; while(n%p==0){n/=p; pe*=p;} if(pe>25) ok=0; } }
        if (n>1 && n>25) ok=0;   /* leftover prime factor > 25 */
        is_pin_only[q]=ok;
    }
}

static int W[12];

static int pairid23[23], npair23;
static int pairid25[25], npair25;
static void init_pairs(void){
    for (int r=0;r<23;r++) pairid23[r]=-1;
    pairid23[0]=-2; npair23=0;
    for (int b=1;b<=11;b++){ pairid23[b]=npair23; pairid23[23-b]=npair23; npair23++; }
    for (int r=0;r<25;r++) pairid25[r]=-1;
    pairid25[0]=-2; npair25=0;
    for (int b=1;b<13;b++){
        if (gcd_i(b,25)!=1) continue;
        if (pairid25[b]>=0) continue;
        pairid25[b]=npair25; pairid25[25-b]=npair25; npair25++;
    }
}

static int pinning_ok(void){
    for (int q=2;q<=25;q++){
        int has0=0;
        for (int i=0;i<12;i++) if (W[i]%q==0){ has0=1; break; }
        if (has0) continue;
        for (int b=1;b<=q/2;b++){
            if (gcd_i(b,q)!=1) continue;
            int hit=0;
            for (int i=0;i<12;i++){ int r=W[i]%q; if (r==b||r==q-b){ hit=1; break; } }
            if (!hit) return 0;
        }
    }
    return 1;
}

static int dist_q(int x,int q){ int r=x%q; if(r<0)r+=q; return r<q-r?r:q-r; }

/* 2/25-margin witness at modulus q?  returns 1 if some unit dilation makes all safe */
static int witness_at(int q){
    for (int i=0;i<12;i++) if (W[i]%q==0) return 0;
    for (int a=1;a<=q/2;a++){
        if (gcd_i(a,q)!=1) continue;
        int ok=1;
        for (int i=0;i<12;i++){ if (25*dist_q(W[i]*a,q) < 2*q){ ok=0; break; } }
        if (ok) return 1;
    }
    return 0;
}

static void dfs(int pos, int maxnext, int cov, int mask23, int mask25){
    if (pos == 12){
        if (!(mask23 & (1<<11)) && __builtin_popcount(mask23 & 0x7ff) < 11) return;
        if (!(mask25 & (1<<10)) && __builtin_popcount(mask25 & 0x3ff) < npair25) return;
        if (cov != 0x3f) return;
        for (int m=2;m<=6;m++){
            int hit=0;
            for (int i=0;i<12;i++) if (W[i]%m==0){ hit=1; break; }
            if (!hit) return;
        }
        if (2*W[0] <= 23*W[11]) return;      /* F2 */
        if (W[0] > 24*W[1]) return;          /* F3 */
        if (2*W[0] < 38) return;             /* F4 */
        { int g=0; for (int i=0;i<12;i++) g=gcd_i(g,W[i]); if (g!=1) return; }
        if (!pinning_ok()) return;           /* F5 */
        cnt_f5++;
        /* --- clearing analysis over [26,50] --- */
        int c2p = 0;
        for (int k=0;k<5;k++) if (witness_at(Q2P[k])) { c2p=1; }
        int minq=0, odd_clear=0, any_clear=0;
        for (int q=26;q<=50;q++){
            if (witness_at(q)){
                any_clear=1;
                q_clears[q]++;
                if (!minq) minq=q;
                if (q%2==1) odd_clear=1;
            }
        }
        if (minq) minq_hist[minq]++;
        if (!any_clear) no_clear_2650++;
        if (c2p) cleared_by_2p++; else not_cleared_by_2p++;
        if (c2p && !odd_clear) cleared_2p_only++;
        /* --- THE KEY: pinned-only (height-independent) clearing --- */
        {
            int cpin=0, pminq=0;
            for (int k=0;k<16;k++){
                int q=QPIN[k];
                if (witness_at(q)){ cpin=1; pin_q_clears[q]++; if (!pminq) pminq=q; }
            }
            if (cpin){ cleared_by_pin++; if (pminq) pin_minq_hist[pminq]++; }
            else {
                not_cleared_by_pin++;   /* no pinned-only witness in [26,50] */
            }
            /* Q0 evidence: smallest pinned-only witness modulus in [26,PMAX] --
             * this is HEIGHT-INDEPENDENT (q | L). Track the max over survivors. */
            int q0=0;
            for (int q=26;q<=PMAX;q++){ if (is_pin_only[q] && witness_at(q)){ q0=q; break; } }
            if (q0 > q0_max){ q0_max=q0; for(int i=0;i<12;i++) q0_argmax[i]=W[i]; }
        }
        /* HARD = no witness anywhere in [8,60] (the census's original gate) */
        {
            int fq=0; for (int q=8;q<=60 && !fq;q++) if (witness_at(q)) fq=q;
            if (!fq){
                cnt_hard++;
                printf("HARD");
                for (int i=11;i>=0;i--) printf(" %d", W[i]);
                printf("\n"); fflush(stdout);
            }
        }
        return;
    }
    int slots = 12 - pos;
    if (!(mask23 & (1<<11))){
        int missing = 11 - __builtin_popcount(mask23 & 0x7ff);
        if (missing > slots) return;
    }
    if (!(mask25 & (1<<10))){
        int missing = npair25 - __builtin_popcount(mask25 & 0x3ff);
        if (missing > slots) return;
    }
    for (int m=7;m<=12;m++){
        if (!(cov & (1<<(m-7))) && maxnext < m) return;
    }
    for (int v = maxnext; v >= 1; v--){
        if (v < slots) break;
        W[pos] = v;
        int cov2 = cov;
        for (int m=7;m<=12;m++) if (v%m==0) cov2 |= 1<<(m-7);
        int mk23 = mask23, mk25 = mask25;
        { int r=v%23; if (r==0) mk23 |= 1<<11; else if (pairid23[r]>=0) mk23 |= 1<<pairid23[r]; }
        { int r=v%25; if (r==0) mk25 |= 1<<10; else if (pairid25[r]>=0) mk25 |= 1<<pairid25[r]; }
        dfs(pos+1, v-1, cov2, mk23, mk25);
    }
}

int main(int argc, char **argv){
    if (argc > 1) B = atoi(argv[1]);
    init_pairs();
    init_pinonly();
    for (int w0 = B; w0 >= 25; w0--){
        W[0] = w0;
        int cov = 0;
        for (int m=7;m<=12;m++) if (w0%m==0) cov |= 1<<(m-7);
        int mk23=0, mk25=0;
        { int r=w0%23; if (r==0) mk23 |= 1<<11; else if (pairid23[r]>=0) mk23 |= 1<<pairid23[r]; }
        { int r=w0%25; if (r==0) mk25 |= 1<<10; else if (pairid25[r]>=0) mk25 |= 1<<pairid25[r]; }
        dfs(1, w0-1, cov, mk23, mk25);
    }
    printf("=== Q50 PARITY ANALYSIS, census [25..%d] ===\n", B);
    printf("full-profile survivors (F1-F5+prim): %lld\n", cnt_f5);
    printf("HARD (no witness in [8,60]): %lld\n", cnt_hard);
    printf("cleared by Q2P={26,34,38,46,50}: %lld / NOT cleared by Q2P: %lld\n",
           cleared_by_2p, not_cleared_by_2p);
    printf("cleared by Q2P but NO odd q in [26,49] (2p NECESSARY): %lld\n", cleared_2p_only);
    printf("no witness anywhere in [26,50]: %lld\n", no_clear_2650);
    printf("*** HEIGHT-INDEPENDENT (pinned-only) clearing ***\n");
    printf("cleared by SOME pinned-only q in [26,50]: %lld / NEED a free modulus: %lld\n",
           cleared_by_pin, not_cleared_by_pin);
    printf("min pinned-only clearing q histogram:");
    for (int q=26;q<=50;q++) if (pin_minq_hist[q]) printf(" %d:%lld", q, pin_minq_hist[q]);
    printf("\nper pinned-only q clear counts:");
    for (int q=26;q<=50;q++) if (pin_q_clears[q]) printf(" %d:%lld", q, pin_q_clears[q]);
    printf("\n");
    printf("min-clearing-q histogram [26,50]:");
    for (int q=26;q<=50;q++) if (minq_hist[q]) printf(" %d:%lld", q, minq_hist[q]);
    printf("\nper-q clear counts [26,50]:");
    for (int q=26;q<=50;q++) if (q_clears[q]) printf(" %d:%lld", q, q_clears[q]);
    printf("\n*** Q0 EVIDENCE (height-independent bound) ***\n");
    printf("MAX over survivors of [smallest pinned-only witness modulus in 26..%d] = %lld\n", PMAX, q0_max);
    printf("  witnessed by base:");
    for (int i=11;i>=0;i--) printf(" %d", q0_argmax[i]);
    printf("\n  => the height-independent (pinned-only) witness bound is >= %lld over the census\n", q0_max);
    return 0;
}

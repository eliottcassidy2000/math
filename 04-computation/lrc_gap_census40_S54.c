/* lrc_gap_census40_S54.c -- mac-mini-2026-07-05-S54, HYP-4117.
 *
 * THE GAP CENSUS EXTENSION [1,B], B = 40 (kps-S3 did [1,24]: 488,894 filtered,
 * ZERO in the gap).  Enumerate ALL primitive 12-subsets of [1,B] passing the
 * kernel-checked gap-violator profile (kps S2-S6, klein, opus):
 *   F1 covering: every q in 2..12 divides some element        (not_loose_dvd)
 *   F2 spread:   w_max * 2 > 23 * w_min                       (ratio gate)
 *   F3 24-compression: w_max <= 24 * w_2nd                    (kps-S6)
 *   F4 big pair: w_max + w_2nd >= 38                          (merge exclusion)
 *   F5 pinning at every q <= 25, every unit a: some v*a in {0,+-1} mod q
 *      (not_loose_near_unit; the heavy filter -- checked as: for every q,
 *       either 0 in W mod q, or W's residues hit every unit +-pair of q)
 * Survivors: exact M via a rational scan gate (margin >= 2/25 witness search
 * at q <= 60ish) -- those with NO witness printed for python exact-M.
 * ALSO report the residue-shape anatomy mod 13 (full system vs doubled pairs
 * vs 13-multiple) of all filtered survivors: the doubled-residue stratum's
 * first census.
 *
 * Enumeration: DFS over elements in decreasing order with branch pruning:
 *   - remaining-slots feasibility for the covering moduli {7,8,9,10,11,12}
 *   - mod-23 pair-hitting budget (11 pairs must be hit; prune when
 *     slots_left < pairs_missing)  [23 = the most rigid pinning modulus]
 *   - mod-25 same (10 unit pairs unless a 25-multiple appears; B=40 => only 25)
 * gcc -O2 -o gap40 lrc_gap_census40_S54.c && ./gap40 40
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int B = 40;
typedef long long ll;
static int gcd_i(int a, int b){ while(b){int t=a%b;a=b;b=t;} return a; }

static ll cnt_visit=0, cnt_leaf=0, cnt_f1=0, cnt_f2=0, cnt_f3=0, cnt_f4=0, cnt_f5=0, cnt_prim=0, cnt_scan=0, cnt_hard=0;
static ll shape_full=0, shape_mult13=0, shape_doubled=0;

static int W[12];

/* pair index mod q: unit pairs {b, q-b}; map residue -> pair id or -1; 0 -> -2 */
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

/* full pinning check at all q <= 25 (only at leaves; 23/25 pruned in-branch) */
static int pinning_ok(void){
    for (int q=2;q<=25;q++){
        int has0=0;
        for (int i=0;i<12;i++) if (W[i]%q==0){ has0=1; break; }
        if (has0) continue;
        /* need: every unit a: some v*a == +-1 mod q  <=> residues hit every
           unit pair {b, q-b}, b unit (v == +-a^{-1})  */
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

/* rational witness gate: margin >= 2/25 at some t=a/q, q in library */
static int scan_gate(void){
    static int qs[64]; int nq=0;
    for (int q=8;q<=60;q++) qs[nq++]=q;
    for (int qi=0;qi<nq;qi++){
        int q=qs[qi], bad=0;
        for (int i=0;i<12;i++) if (W[i]%q==0){ bad=1; break; }
        if (bad) continue;
        for (int a=1;a<=q/2;a++){
            if (gcd_i(a,q)!=1) continue;
            int ok=1;
            for (int i=0;i<12;i++){ if (25*dist_q(W[i]*a,q) < 2*q){ ok=0; break; } }
            if (ok) return 1;
        }
    }
    return 0;
}

/* DFS: choose elements in DECREASING order. pos = index (0 = largest). */
static int need_mult[6];        /* moduli 7..12 covered? bit per m-7 */
static void dfs(int pos, int maxnext, int cov, int mask23, int mask25){
    cnt_visit++;
    if (pos == 12){
        cnt_leaf++;
        /* cheap mask gates first (23/25 pinning necessary conditions) */
        if (!(mask23 & (1<<11)) && __builtin_popcount(mask23 & 0x7ff) < 11) return;
        if (!(mask25 & (1<<10)) && __builtin_popcount(mask25 & 0x3ff) < npair25) return;
        /* F1 remaining moduli 2..6 (7..12 tracked in cov) */
        if (cov != 0x3f) return;                    /* 7..12 all covered */
        for (int m=2;m<=6;m++){
            int hit=0;
            for (int i=0;i<12;i++) if (W[i]%m==0){ hit=1; break; }
            if (!hit) return;
        }
        cnt_f1++;
        /* F2 spread (W[0]=max, W[11]=min) */
        if (2*W[0] <= 23*W[11]) return;
        cnt_f2++;
        /* F3 24-compression */
        if (W[0] > 24*W[1]) return;
        cnt_f3++;
        /* F4 big pair */
        if (W[0]+W[1] < 38) return;
        cnt_f4++;
        /* primitivity */
        { int g=0; for (int i=0;i<12;i++) g=gcd_i(g,W[i]); if (g!=1) return; }
        cnt_prim++;
        /* F5 full pinning */
        if (!pinning_ok()) return;
        cnt_f5++;
        /* residue shape mod 13 */
        {
            int c13[13]; memset(c13,0,sizeof c13);
            for (int i=0;i<12;i++) c13[W[i]%13]++;
            if (c13[0]) shape_mult13++;
            else {
                int full=1;
                for (int r=1;r<=12;r++) if (c13[r]!=1){ full=0; break; }
                if (full) shape_full++; else shape_doubled++;
            }
        }
        /* witness gate */
        if (scan_gate()){ cnt_scan++; return; }
        cnt_hard++;
        printf("HARD");
        for (int i=11;i>=0;i--) printf(" %d", W[i]);
        printf("\n"); fflush(stdout);
        return;
    }
    int slots = 12 - pos;
    /* prune: mod-23 pairs: 11 - popcount(mask23) still missing (if no 0 yet) */
    if (!(mask23 & (1<<11))){                        /* bit 11 = "has 23-mult" */
        int missing = 11 - __builtin_popcount(mask23 & 0x7ff);
        if (missing > slots) return;
    }
    if (!(mask25 & (1<<10))){
        int missing = npair25 - __builtin_popcount(mask25 & 0x3ff);
        if (missing > slots) return;
    }
    /* prune: covering 7..12: each uncovered m needs a multiple <= maxnext */
    for (int m=7;m<=12;m++){
        if (!(cov & (1<<(m-7))) && maxnext < m) return;
    }
    for (int v = maxnext; v >= 1; v--){
        /* feasibility: enough room below v for remaining slots */
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
    /* w_max >= 25 (census [1,24] done by kps-S3); w_max <= B */
    for (int w0 = B; w0 >= 25; w0--){
        W[0] = w0;
        int cov = 0;
        for (int m=7;m<=12;m++) if (w0%m==0) cov |= 1<<(m-7);
        int mk23=0, mk25=0;
        { int r=w0%23; if (r==0) mk23 |= 1<<11; else if (pairid23[r]>=0) mk23 |= 1<<pairid23[r]; }
        { int r=w0%25; if (r==0) mk25 |= 1<<10; else if (pairid25[r]>=0) mk25 |= 1<<pairid25[r]; }
        dfs(1, w0-1, cov, mk23, mk25);
        fprintf(stderr, "w_max=%d done: visit=%lld leaf=%lld F1=%lld F2=%lld F3=%lld F4=%lld prim=%lld F5=%lld scan-ok=%lld HARD=%lld\n",
                w0, cnt_visit, cnt_leaf, cnt_f1, cnt_f2, cnt_f3, cnt_f4, cnt_prim, cnt_f5, cnt_scan, cnt_hard);
    }
    printf("CENSUS [25..%d]: visited=%lld leaves=%lld | F1(cover)=%lld F2(spread)=%lld F3(24comp)=%lld F4(pair38)=%lld prim=%lld F5(pinning)=%lld | witness-ok=%lld HARD=%lld\n",
           B, cnt_visit, cnt_leaf, cnt_f1, cnt_f2, cnt_f3, cnt_f4, cnt_prim, cnt_f5, cnt_scan, cnt_hard);
    printf("residue shapes mod 13 of fully-filtered: full-system=%lld with-13-mult=%lld DOUBLED-PAIR=%lld\n",
           shape_full, shape_mult13, shape_doubled);
    return 0;
}

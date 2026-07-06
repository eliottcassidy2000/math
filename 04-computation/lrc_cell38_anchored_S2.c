/* lrc_cell38_anchored_S2.c -- mac-mini-2026-07-06-S2, HYP-4242.
 * THE ANCHORED k=1 SWEEP: 12-subsets of [1,B] CONTAINING a both-odd pair
 * (v, 38-v), v in {1,3,...,17}, against the FULL gap profile + the k=1
 * witness-determinism filter (all runners clear >= 3 mod 38 under the FORCED
 * dilation m = 3 v^{-1} mod 38 -- necessary for 3/38-attainment with this
 * binding pair).  Survivors: witness gate >= 2/25 (disproves attainment);
 * HARD rows -> exact-M post-pass.  Verdict expected: ZERO in-cell.
 * Derived from lrc_gap_census40_S54.c (same profile filters + prunes).
 *
 * ORIGINAL HEADER: lrc_gap_census40_S54.c -- mac-mini-2026-07-05-S54, HYP-4117.
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
/* spectroscopy */
static ll firstq_hist[64];              /* first witness q (8..60) */
static ll comp_ok=0, comp_fail=0;       /* composite-CRT witness exists? */
static ll compq_hist[640];              /* first composite witness q */
static ll bot_in12=0, bot_not12=0;      /* bottom-6 inside [1,12]? */
static ll wmin_hist[64];                /* w_min histogram (scale stat) */
static ll botmax_hist[160];             /* w_(6) = bottom-6 max histogram */
static int COMPQ[256]; static int nCOMPQ=0;
static void init_compq(void){
    unsigned char seen[640]; memset(seen,0,sizeof seen);
    for (int a1=13;a1<=25;a1++) for (int b1=a1;b1<=25;b1++){
        int pr=a1*b1;
        for (int q=26;q<=600 && q<640;q++) if (pr%q==0 && !seen[q]){ seen[q]=1; }
    }
    for (int q=26;q<640;q++) if (seen[q]) COMPQ[nCOMPQ++]=q;
}

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

/* rational witness gate: margin >= 2/25 at some t=a/q; returns q or 0 */
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
static int scan_gate(void){
    for (int q=8;q<=60;q++) if (witness_at(q)) return q;
    return 0;
}
static int comp_gate(void){
    for (int i=0;i<nCOMPQ;i++) if (witness_at(COMPQ[i])) return COMPQ[i];
    return 0;
}

static int anchored_leaf_extra(void);   /* fwd decl (defined near main) */
/* DFS: choose elements in DECREASING order. pos = index (0 = largest). */
static int need_mult[6];        /* moduli 7..12 covered? bit per m-7 */
static void dfs(int pos, int maxnext, int cov, int mask23, int mask25){
    cnt_visit++;
    if (pos == 12){
        cnt_leaf++;
        /* cheap mask gates first (23/25 pinning necessary conditions) */
        if (!(mask23 & (1<<11)) && __builtin_popcount(mask23 & 0x7ff) < 11) return;
        if (!(mask25 & (1<<10)) && __builtin_popcount(mask25 & 0x3ff) < npair25) return;
        if (!anchored_leaf_extra()) return;         /* k=1 determinism filter */
        /* order statistics computed locally (anchor preload breaks W ordering) */
        int mx=0, mx2=0, mn=1<<30;
        for (int i=0;i<12;i++){
            int w=W[i];
            if (w>mx){ mx2=mx; mx=w; } else if (w>mx2) mx2=w;
            if (w<mn) mn=w;
        }
        /* F1 remaining moduli 2..6 (7..12 tracked in cov) */
        if (cov != 0x3f) return;                    /* 7..12 all covered */
        for (int m=2;m<=6;m++){
            int hit=0;
            for (int i=0;i<12;i++) if (W[i]%m==0){ hit=1; break; }
            if (!hit) return;
        }
        cnt_f1++;
        /* F2 spread */
        if (2*mx <= 23*mn) return;
        cnt_f2++;
        /* F3 24-compression */
        if (mx > 24*mx2) return;
        cnt_f3++;
        /* F4 big pair (corrected form, i=j allowed) */
        if (2*mx < 38) return;
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
        if (mn < 64) wmin_hist[mn]++;
        /* witness gate */
        { int fq = scan_gate();
          if (fq){ cnt_scan++; if (fq<64) firstq_hist[fq]++; return; } }
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
        if (pos >= 2 && (v == W[0] || v == W[1])) continue;   /* anchored pair dupes */
        W[pos] = v;
        int cov2 = cov;
        for (int m=7;m<=12;m++) if (v%m==0) cov2 |= 1<<(m-7);
        int mk23 = mask23, mk25 = mask25;
        { int r=v%23; if (r==0) mk23 |= 1<<11; else if (pairid23[r]>=0) mk23 |= 1<<pairid23[r]; }
        { int r=v%25; if (r==0) mk25 |= 1<<10; else if (pairid25[r]>=0) mk25 |= 1<<pairid25[r]; }
        dfs(pos+1, v-1, cov2, mk23, mk25);
    }
}

static int ANCV = 1;                 /* the anchor v (odd, pair = (v, 38-v)) */
static int FORB[5];                  /* forbidden residues mod 38 for ALL runners */
static int anchored_leaf_extra(void){
    /* determinism filter: every runner clears >= 3 mod 38 under m = 3*v^{-1}:
       equivalently no runner residue in {0, +-13v, +-26v} mod 38 (precomputed) */
    for (int i=0;i<12;i++){
        int r = W[i] % 38;
        for (int j=0;j<5;j++) if (r == FORB[j]) return 0;
    }
    return 1;
}
int main(int argc, char **argv){
    if (argc > 1) B = atoi(argv[1]);
    int vlo = argc > 2 ? atoi(argv[2]) : 1, vhi = argc > 3 ? atoi(argv[3]) : 17;
    init_pairs();
    init_compq();
    for (int v = vlo; v <= vhi; v += 2){
        if (v == 19) continue;
        ANCV = v;
        /* forbidden residues: 13v * {0, +-1, +-2} mod 38 */
        int m13v = (13 * v) % 38;
        FORB[0] = 0;
        FORB[1] = m13v; FORB[2] = (38 - m13v) % 38;
        FORB[3] = (2 * m13v) % 38; FORB[4] = (38 - 2 * m13v % 38) % 38;
        int p1 = v, p2 = 38 - v;
        /* seed the DFS with the anchor pair as the two largest-of-their-kind:
           enumerate the other TEN over [1,B] \ {p1,p2} in decreasing order.
           Reuse dfs with W[0], W[1] preloaded and pos starting at 2; maxnext = B. */
        W[0] = p2 > p1 ? p2 : p1;    /* keep W sorted-desc-ish: pair first */
        W[1] = p2 > p1 ? p1 : p2;
        int cov = 0;
        for (int m=7;m<=12;m++){ if (W[0]%m==0) cov |= 1<<(m-7); if (W[1]%m==0) cov |= 1<<(m-7); }
        int mk23=0, mk25=0;
        for (int t=0;t<2;t++){
            int r=W[t]%23; if (r==0) mk23 |= 1<<11; else if (pairid23[r]>=0) mk23 |= 1<<pairid23[r];
            r=W[t]%25; if (r==0) mk25 |= 1<<10; else if (pairid25[r]>=0) mk25 |= 1<<pairid25[r];
        }
        dfs(2, B, cov, mk23, mk25);
        fprintf(stderr, "anchor v=%d done: visit=%lld leaf=%lld F5=%lld scan-ok=%lld HARD=%lld\n",
                v, cnt_visit, cnt_leaf, cnt_f5, cnt_scan, cnt_hard);
    }
    printf("CENSUS [25..%d]: visited=%lld leaves=%lld | F1(cover)=%lld F2(spread)=%lld F3(24comp)=%lld F4(pair38)=%lld prim=%lld F5(pinning)=%lld | witness-ok=%lld HARD=%lld\n",
           B, cnt_visit, cnt_leaf, cnt_f1, cnt_f2, cnt_f3, cnt_f4, cnt_prim, cnt_f5, cnt_scan, cnt_hard);
    printf("residue shapes mod 13 of fully-filtered: full-system=%lld with-13-mult=%lld DOUBLED-PAIR=%lld\n",
           shape_full, shape_mult13, shape_doubled);
    printf("SPECTROSCOPY: composite-CRT witness exists: %lld / fails: %lld\n", comp_ok, comp_fail);
    printf("first-q histogram (q: count):");
    for (int q=8;q<64;q++) if (firstq_hist[q]) printf(" %d:%lld", q, firstq_hist[q]);
    printf("\nfirst composite-q histogram:");
    { int shown=0; for (int q=26;q<640 && shown<40;q++) if (compq_hist[q]){ printf(" %d:%lld", q, compq_hist[q]); shown++; } }
    printf("\nBOTTOM-6: inside [1,12]: %lld / not: %lld\n", bot_in12, bot_not12);
    printf("w_min histogram:");
    for (int v=1;v<64;v++) if (wmin_hist[v]) printf(" %d:%lld", v, wmin_hist[v]);
    printf("\nbottom-6-max histogram:");
    for (int v=1;v<160;v++) if (botmax_hist[v]) printf(" %d:%lld", v, botmax_hist[v]);
    printf("\n");
    return 0;
}

/* lrc14_subgap_census_klein_S319.c -- klein-2026-07-19-S319, THM-1290.
 *
 * EXHAUSTIVE census of 13-subsets W of [1,B] with M(W) in the OPEN interval
 * (LO, HI) = (1/14, 3/41) -- the n=14 sub-gap below the attained edge
 * (THM-1230/1235).  Adapts mac-mini S54's n=12 harness (lrc_gap_census40_S54.c)
 * to 13 speeds with the new rigorous filters:
 *
 *   F1 covering: every q in 2..13 divides some element (else M >= 1/q >= 1/13 > HI)
 *   F2 spread (THM-1043): 12*w_min < w_max  (else M >= 1/13)
 *   F4 active-pair sum (THM-1269, D = M*s): some pair v_i+v_j = s with an
 *      integer D, LO < D/s < HI.  Admissible s (for HI=3/41): 55,69,83,96,97,...
 *      In particular w_max >= ceil(s_min/2).
 *   F5 pinning: for every q in [2, QPIN] with depth d(q) = ceil(HI*q)-1 >= 0:
 *      either q | some v, or for every unit a mod q some v has
 *      dist(v*a mod q) <= d(q).  (Witness-failure reformulation; d(q)=0 is F1.)
 *   In-branch mask pruning at q = 23, 25, 27 (all depth-1 for HI = 3/41).
 *   Primitivity: gcd = 1.  (Non-primitive g*W' forces g*28 <= B; for B <= 55
 *      only g=1 possible, so primitive enumeration is complete -- see THM-1290.)
 *
 * Survivors: rational witness scan (margin >= HI at t = a/q, q <= QSCAN,
 * integer test den_hi*dist >= num_hi*q) certifies M >= HI (out of gap, above).
 * Families with NO witness print HARD for the exact-M python referee
 * (lrc14_subgap_referee_klein_S319.py).  The referee decides in-gap / at-floor
 * / above and would expose any M < 1/14 LRC(14) counterexample.
 *
 * Debug/rediscovery gates: -g NUM DEN replaces HI by NUM/DEN.  With HI = 2/27
 * (B = 36) the census MUST find {1..11,13,36} (M = 3/41 strictly inside);
 * with HI = 1/13 (B = 30) it MUST find {1..12,26} (M = 2/27 strictly inside).
 * Boundary families ({1..13} at LO, the HI attainers) must NOT be reported.
 *
 * gcc -O2 -o subgap lrc14_subgap_census_klein_S319.c
 * ./subgap B [-g NUM DEN] [-w W0LO W0HI]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define K 13            /* number of speeds (n = 14 runners) */
#define MAXB 64
#define QSCAN 200
#ifndef QPIN
#define QPIN 48         /* full pinning up to here (S320: 41 -> 48; depth 3 at 42..48 for HI=3/41).
                         * Override with -DQPIN=41 to reproduce v1 pinning (spectroscopy runs). */
#endif
#define QMLO 14         /* in-branch mask moduli range [QMLO, QMHI]: all q with depth d(q)=1 used */
#define QMHI 27
#define NQM  (QMHI-QMLO+1)

static int B = 36;
/* open interval (LO_N/LO_D, HI_N/HI_D); defaults (1/14, 3/41).
 * -lo 0 1 -g 1 14 gives LRC(14)-verification mode: M in (0, 1/14) hunts
 * counterexamples; covering 2..14 and spread factor 13 arise automatically
 * from the generic d(q)/spread machinery. */
static long long LO_N = 1, LO_D = 14, HI_N = 3, HI_D = 41;
static int SPREADF = 12;   /* = floor(HI_D/HI_N) - 1 (THM-1043 rung) */
typedef long long ll;
static int gcd_i(int a,int b){ while(b){int t=a%b;a=b;b=t;} return a; }

static int W[K];
static unsigned char admiss[2*MAXB+2];     /* pair-sum admissible? */
static int cntadm[2*MAXB+2];               /* prefix count of admissible */
static int smin_glob;                       /* least admissible pair sum */

static ll cnt_visit=0, cnt_leaf=0, cnt_f1=0, cnt_f2=0, cnt_f4=0, cnt_prim=0,
          cnt_f5=0, cnt_scan=0, cnt_hard=0;
static ll firstq_hist[QSCAN+1];

/* ---------- admissible pair sums: exists D with LO < D/s < HI ---------- */
static void init_admiss(void){
    memset(admiss,0,sizeof admiss); smin_glob = 0;
    for (int s=1; s<=2*MAXB; s++){
        /* D > s*LO_N/LO_D  and  D < s*HI_N/HI_D, strict both sides.
         * dlo = floor(s*LO)+1 (works for exact division too);
         * dhi = (s*HI_N - 1)/HI_D = largest D < s*HI (works both cases). */
        ll dlo = (ll)s*LO_N / LO_D + 1;
        ll dhi = ((ll)s*HI_N - 1) / HI_D;
        if (dlo <= dhi){ admiss[s]=1; if(!smin_glob) smin_glob=s; }
    }
    cntadm[0]=0;
    for (int s=1; s<=2*MAXB; s++) cntadm[s]=cntadm[s-1]+admiss[s];
}

/* ---------- pinning depth per modulus ---------- */
static int dq[QPIN+1];                     /* d(q) = ceil(HI*q) - 1 */
static int useQ[NQM];                      /* mask prune valid only at depth 1 */
static int covreq;                         /* required bits for moduli 7..14 */
static void init_dq(void){
    for (int q=2;q<=QPIN;q++){
        ll c = (HI_N*q + HI_D - 1)/HI_D;   /* ceil(HI*q) */
        dq[q] = (int)c - 1;
        if (dq[q] < 0) dq[q] = 0;
    }
    for (int qi=0;qi<NQM;qi++) useQ[qi] = (dq[qi+QMLO]==1);
    covreq = (dq[14]==0) ? 0xff : 0x7f;    /* track 7..14; require 14 iff d=0 */
}

/* generalized in-branch mask machinery (S320): one bitmask per modulus
 * q in [QMLO, QMHI], used only when d(q) == 1.  pairidQ[qi][r] = unit-pair
 * index, -1 non-unit, -2 zero; bit npQ[qi] = "has a multiple of q".
 * A family with M < HI and no multiple of q must hit EVERY unit pair mod q
 * (necessary condition; the leaf pinning_ok remains the authority). */
static int pairidQ[NQM][QMHI+1];
static int npQ[NQM];
static unsigned mkQ[K+1][NQM];      /* mask state per DFS depth */
static void init_pairs(void){
    for (int q=QMLO;q<=QMHI;q++){
        int qi=q-QMLO;
        for (int r=0;r<=QMHI;r++) pairidQ[qi][r]=-1;
        pairidQ[qi][0]=-2; npQ[qi]=0;
        for (int b=1;b<=q/2;b++){
            if (gcd_i(b,q)!=1 || pairidQ[qi][b]>=0) continue;
            pairidQ[qi][b]=npQ[qi]; pairidQ[qi][q-b]=npQ[qi]; npQ[qi]++;
        }
    }
}
/* fill mkQ[pos+1] from mkQ[pos] and the pushed element v */
static void apply_masks(int pos, int v){
    for (int qi=0;qi<NQM;qi++){
        int q=qi+QMLO, r=v%q;
        unsigned m=mkQ[pos][qi];
        if (r==0) m |= 1u<<npQ[qi];
        else if (pairidQ[qi][r]>=0) m |= 1u<<pairidQ[qi][r];
        mkQ[pos+1][qi]=m;
    }
}

static int dist_q(int x,int q){ int r=x%q; if(r<0)r+=q; return r<q-r?r:q-r; }

/* full pinning at leaves for all q in [2, QPIN]:
 * q | some v, OR for every unit a: min_v dist(v*a mod q) <= d(q).
 * (d(q) == 0 moduli are the covering filter F1 when q <= 13; for q >= 14
 *  d(q) >= 1 automatically at HI >= 3/41.) */
static int pinning_ok(void){
    for (int q=2;q<=QPIN;q++){
        int has0=0;
        for (int i=0;i<K;i++) if (W[i]%q==0){ has0=1; break; }
        if (has0) continue;
        int d = dq[q];
        for (int a=1;a<=q/2;a++){
            if (gcd_i(a,q)!=1) continue;
            int ok=0;
            for (int i=0;i<K;i++) if (dist_q(W[i]*a,q) <= d){ ok=1; break; }
            if (!ok) return 0;
        }
    }
    return 1;
}

/* rational witness: exists a/q with min_v dist >= HI  <=>  HI_D*dist >= HI_N*q */
static int witness_at(int q){
    for (int i=0;i<K;i++) if (W[i]%q==0) return 0;
    for (int a=1;a<=q/2;a++){
        if (gcd_i(a,q)!=1) continue;
        int ok=1;
        for (int i=0;i<K;i++)
            if ((ll)HI_D*dist_q(W[i]*a,q) < (ll)HI_N*q){ ok=0; break; }
        if (ok) return 1;
    }
    return 0;
}
static int scan_gate(void){
    for (int q=14;q<=QSCAN;q++) if (witness_at(q)) return q;
    return 0;
}

/* ---------- DFS, elements strictly decreasing ---------- */
static void dfs(int pos, int maxnext, int cov, int haspair){
    cnt_visit++;
    if (pos == K){
        cnt_leaf++;
        if ((cov & covreq) != covreq) return;          /* 7..13(..14) covered */
        for (int m=2;m<=6;m++){
            int hit=0;
            for (int i=0;i<K;i++) if (W[i]%m==0){ hit=1; break; }
            if (!hit) return;
        }
        /* covering at 14 (and any q with d(q)=0 beyond 13) is enforced by
         * pinning_ok below; in LRC mode d(14)=0 forces a 14-multiple. */
        cnt_f1++;
        if (SPREADF*W[K-1] >= W[0]) return;            /* F2 spread */
        cnt_f2++;
        if (!haspair) return;                          /* F4 (tracked exact) */
        {   /* re-verify F4 exactly at leaf (belt & braces) */
            int ok=0;
            for (int i=0;i<K && !ok;i++)
                for (int j=i+1;j<K;j++)
                    if (admiss[W[i]+W[j]]){ ok=1; break; }
            if (!ok) return;
        }
        cnt_f4++;
        { int g=0; for (int i=0;i<K;i++) g=gcd_i(g,W[i]); if (g!=1) return; }
        cnt_prim++;
        if (!pinning_ok()) return;
        cnt_f5++;
        { int fq = scan_gate();
          if (fq){ cnt_scan++; firstq_hist[fq]++;
                   printf("SURV q=%d :", fq);
                   for (int i=K-1;i>=0;i--) printf(" %d", W[i]);
                   printf("\n"); fflush(stdout);
                   return; } }
        cnt_hard++;
        printf("HARD");
        for (int i=K-1;i>=0;i--) printf(" %d", W[i]);
        printf("\n"); fflush(stdout);
        return;
    }
    int slots = K - pos;
    /* mask prunes over all depth-1 moduli (only when no multiple present) */
    for (int qi=0;qi<NQM;qi++){
        if (!useQ[qi]) continue;
        unsigned m = mkQ[pos][qi];
        if (m & (1u<<npQ[qi])) continue;
        if (npQ[qi] - __builtin_popcount(m & ((1u<<npQ[qi])-1)) > slots) return;
    }
    /* covering 7..13(..14) feasibility: uncovered required m needs a multiple <= maxnext */
    for (int m=7;m<=14;m++)
        if ((covreq & (1<<(m-7))) && !(cov & (1<<(m-7))) && maxnext < m) return;
    /* F4 feasibility */
    if (!haspair){
        int feasible = 0;
        if (2*maxnext-1 >= smin_glob) feasible = 1;    /* two future elements */
        else {
            for (int i=0;i<pos;i++){
                int hi = W[i]+maxnext; if (hi > 2*MAXB) hi = 2*MAXB;
                if (cntadm[hi] > cntadm[W[i]]) { feasible=1; break; }
            }
        }
        if (!feasible) return;
    }
    for (int v = maxnext; v >= 1; v--){
        if (v < slots) break;                           /* room below */
        if (pos == K-1 && SPREADF*v >= W[0]) continue;  /* F2 on the min */
        W[pos] = v;
        int cov2 = cov;
        for (int m=7;m<=14;m++) if (v%m==0) cov2 |= 1<<(m-7);
        apply_masks(pos, v);
        int hp2 = haspair;
        if (!hp2) for (int i=0;i<pos;i++) if (admiss[W[i]+v]){ hp2=1; break; }
        dfs(pos+1, v-1, cov2, hp2);
    }
}

int main(int argc, char **argv){
    int w0lo = -1, w0hi = -1;
    for (int i=1;i<argc;i++){
        if (!strcmp(argv[i],"-g") && i+2<argc){
            HI_N = atoll(argv[i+1]); HI_D = atoll(argv[i+2]); i+=2;
        } else if (!strcmp(argv[i],"-lo") && i+2<argc){
            LO_N = atoll(argv[i+1]); LO_D = atoll(argv[i+2]); i+=2;
        } else if (!strcmp(argv[i],"-w") && i+2<argc){
            w0lo = atoi(argv[i+1]); w0hi = atoi(argv[i+2]); i+=2;
        } else B = atoi(argv[i]);
    }
    if (B > MAXB){ fprintf(stderr,"B > %d unsupported\n", MAXB); return 1; }
    if (HI_N*13 > HI_D){ fprintf(stderr,"HI must be <= 1/13 (F1/F2 soundness)\n"); return 1; }
    SPREADF = (int)(HI_D/HI_N) - 1;            /* THM-1043 rung: M < HI <= 1/n => sigma > n-1 */
    init_admiss(); init_dq(); init_pairs();
    if (!smin_glob){ printf("no admissible pair sums at all: interval empty a priori\n"); return 0; }
    int w0min = (smin_glob+1)/2;              /* v_max >= ceil(s_min/2) */
    if (w0min < K) w0min = K;
    if (w0min < SPREADF+1) w0min = SPREADF+1; /* v_max > SPREADF*v_min >= SPREADF */
    if (w0lo < 0){ w0lo = w0min; w0hi = B; }
    fprintf(stderr, "interval (%lld/%lld, %lld/%lld), B=%d, s_min=%d, w0 in [%d,%d]\n",
            LO_N,LO_D,HI_N,HI_D,B,smin_glob,w0lo,w0hi);
    fprintf(stderr, "admissible s up to %d:", 2*B);
    for (int s=1;s<=2*B;s++) if (admiss[s]) fprintf(stderr, " %d", s);
    fprintf(stderr, "\n");
    for (int w0 = w0hi; w0 >= w0lo; w0--){
        W[0] = w0;
        int cov=0;
        for (int m=7;m<=14;m++) if (w0%m==0) cov |= 1<<(m-7);
        memset(mkQ[0], 0, sizeof mkQ[0]);
        apply_masks(0, w0);
        dfs(1, w0-1, cov, 0);
        fprintf(stderr, "w_max=%d done: visit=%lld leaf=%lld F1=%lld F2=%lld F4=%lld prim=%lld F5=%lld scan=%lld HARD=%lld\n",
                w0, cnt_visit, cnt_leaf, cnt_f1, cnt_f2, cnt_f4, cnt_prim, cnt_f5, cnt_scan, cnt_hard);
    }
    printf("CENSUS 13-subsets of [1,%d], w_max in [%d,%d], interval (%lld/%lld,%lld/%lld):\n",
           B, w0lo, w0hi, LO_N,LO_D,HI_N,HI_D);
    printf("visited=%lld leaves=%lld F1=%lld F2=%lld F4=%lld prim=%lld F5=%lld witness-cleared=%lld HARD=%lld\n",
           cnt_visit, cnt_leaf, cnt_f1, cnt_f2, cnt_f4, cnt_prim, cnt_f5, cnt_scan, cnt_hard);
    printf("first-witness-q histogram:");
    for (int q=14;q<=QSCAN;q++) if (firstq_hist[q]) printf(" %d:%lld", q, firstq_hist[q]);
    printf("\n");
    return 0;
}

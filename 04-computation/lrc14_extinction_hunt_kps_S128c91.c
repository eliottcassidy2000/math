/* kind-pasteur-2026-07-19-S128c91 -- HYP-7975: the extinction prime hunt.
 *
 * c(p) := minimum number of folded danger-APs A(w) = {fold(j*w mod p) : j<=dk},
 * dk = p/n (integer division; n = 14 for the k=13 problem, threshold 1/14),
 * needed to cover the folded line P = [1..(p-1)/2].  If c(p) > 13 then the
 * level-1 improper set I(13,p,1) is EMPTY (13 danger-APs cannot cover), and
 * prime p verifies in the Rosenfeld/S-T sieve for free (Lemma-6 vacuous at
 * level 1).
 * Modes:  (no args)          gate + full run
 *         exact  LO HI       exact c(p) table on primes in [LO,HI]
 *         decide LO HI [BUD] cap-13 decision sweep (extinct iff no <=13 cover)
 *         control            k=12 twin (n=13, cap 12) over the S-T range [167,733]
 * Search: DFS, MRV pivot (least candidates among up to 24 lowest uncovered),
 * include-exclude banning, count cut |U| > dk*rem.  WLOG 1 in W (scaling).
 * Node budget per (p,cap); exceeding = UNKNOWN (never treated as "no").
 * Fractional note (proved one line): weight 1/dk on every w covers each point
 * exactly once (each x lies in exactly dk sets), so the FRACTIONAL covering
 * number is h/dk -> 7 exactly; all excess of c(p) is integrality gap, and
 * greedy gives c(p) <= (h/dk)(1 + ln dk) = O(log p).
 * Build: gcc -O3 -march=native -o lrc14_extinction_hunt_kps_S128c91.exe \
 *            lrc14_extinction_hunt_kps_S128c91.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXP 3072
#define MAXH (MAXP/2 + 2)
#define W 26

static int h, dk, cap;
static uint64_t maskA[MAXH][W];
static int cand[MAXH][256];
static int ncand;
static unsigned char used[MAXH], banned[MAXH];
static long long nodes, budget;
static int found_flag;
static int nw;

static inline int popU(const uint64_t *U){
    int s = 0;
    for (int i = 0; i < nw; i++) s += __builtin_popcountll(U[i]);
    return s;
}

static int fold(long long x, int p){
    x %= p; if (x < 0) x += p;
    return (x <= p/2) ? (int)x : (int)(p - x);
}

static long long modpow(long long b, long long e, long long m){
    long long r = 1; b %= m;
    while (e) { if (e & 1) r = r*b % m; b = b*b % m; e >>= 1; }
    return r;
}

static void build(int p, int n){
    h = (p-1)/2; dk = p/n; ncand = dk;
    nw = (h + 63)/64 + 1; if (nw > W) { fprintf(stderr, "p too big\n"); exit(1); }
    memset(maskA, 0, sizeof(uint64_t)*MAXH*W);
    for (int w = 1; w <= h; w++)
        for (int j = 1; j <= dk; j++){
            int x = fold((long long)j*w, p);
            maskA[w][(x-1)>>6] |= 1ULL << ((x-1)&63);
        }
    for (int x = 1; x <= h; x++)
        for (int j = 1; j <= dk; j++){
            long long inv = modpow(j, p-2, p);
            cand[x][j-1] = fold(inv % p * x, p);
        }
}

static void rec(uint64_t *U, int depth){
    if (found_flag) return;
    if (++nodes > budget) { found_flag = -1; return; }
    int first = -1;
    for (int i = 0; i < nw; i++) if (U[i]) { first = i; break; }
    if (first < 0) { found_flag = 1; return; }
    int rem = cap - depth;
    if (rem <= 0) return;
    if (popU(U) > dk*rem) return;
    int bestx = -1, bestn = 1<<30, seen = 0;
    for (int i = first; i < nw && seen < 24 && bestn > 1; i++){
        uint64_t v = U[i];
        while (v && seen < 24){
            int b = __builtin_ctzll(v); v &= v-1;
            int x = (i<<6) + b + 1;
            int cnt = 0;
            for (int j = 0; j < ncand; j++){
                int w = cand[x][j];
                if (!used[w] && !banned[w]) cnt++;
            }
            if (cnt == 0) return;
            if (cnt < bestn) { bestn = cnt; bestx = x; }
            seen++;
        }
    }
    int local[256], nl = 0;
    for (int j = 0; j < ncand; j++){
        int w = cand[bestx][j];
        if (used[w] || banned[w]) continue;
        used[w] = 1;
        uint64_t U2[W];
        for (int i = 0; i < nw; i++) U2[i] = U[i] & ~maskA[w][i];
        rec(U2, depth+1);
        used[w] = 0;
        if (found_flag) { for (int t = 0; t < nl; t++) banned[local[t]] = 0; return; }
        banned[w] = 1; local[nl++] = w;
    }
    for (int t = 0; t < nl; t++) banned[local[t]] = 0;
}

static int cover_le(int p, int n, int size_cap, long long bud){
    build(p, n);
    cap = size_cap; budget = bud; nodes = 0; found_flag = 0;
    memset(used, 0, h+1); memset(banned, 0, h+1);
    uint64_t U[W];
    for (int i = 0; i < nw; i++) U[i] = 0;
    for (int x = 1; x <= h; x++) U[(x-1)>>6] |= 1ULL << ((x-1)&63);
    used[1] = 1;
    for (int i = 0; i < nw; i++) U[i] &= ~maskA[1][i];
    rec(U, 1);
    return (found_flag == 1) ? 1 : (found_flag == -1 ? -1 : 0);
}

static int isprime(int x){
    if (x < 2) return 0;
    for (int d = 2; (long long)d*d <= x; d++) if (x % d == 0) return 0;
    return 1;
}

static void run_exact(int lo, int hi, long long BUD){
    printf("== exact c(p), p = %d..%d (n=14) ==\n", lo, hi);
    for (int p = lo|1; p <= hi; p += 2){
        if (!isprime(p) || p % 7 == 0) continue;
        int lb = ((p-1)/2 + p/14 - 1) / (p/14);
        int c = -1; long long tot = 0;
        for (int s = lb; s <= 14; s++){
            int r = cover_le(p, 14, s, BUD);
            tot += nodes;
            if (r == 1) { c = s; break; }
            if (r == -1) { c = -2; break; }
        }
        if (c == -2)      printf("  p=%d dk=%d bound=%d c=UNKNOWN(budget) nodes=%lld\n", p, p/14, lb, tot);
        else if (c == -1) printf("  p=%d dk=%d bound=%d c>14 nodes=%lld  <-- EXTINCT++\n", p, p/14, lb, tot);
        else              printf("  p=%d dk=%d bound=%d c=%d nodes=%lld%s\n", p, p/14, lb, c, tot,
                                  c >= 14 ? "  <-- EXTINCT (I(13,p,1) empty)" : "");
    }
}

static void run_decide(int lo, int hi, long long bud2){
    printf("== cap-13 decision, p = %d..%d (budget %lld) ==\n", lo, hi, bud2);
    int run = 0, first = 0;
    for (int p = lo|1; p <= hi; p += 2){
        if (!isprime(p) || p % 7 == 0) continue;
        int r = cover_le(p, 14, 13, bud2);
        if (r == 0){
            printf("  p=%d EXTINCT (no <=13 cover; I(13,p,1) EMPTY) nodes=%lld\n", p, nodes);
            if (!first) first = p; run++;
        } else if (r == 1){ printf("  p=%d alive nodes=%lld\n", p, nodes); run = 0; }
        else { printf("  p=%d UNKNOWN(budget) nodes=%lld\n", p, nodes); run = 0; }
        if (run >= 8){ printf("  8 consecutive extinct -- stopping\n"); break; }
    }
    if (first) printf("  FIRST EXTINCT PRIME (k=13, this range): %d\n", first);
}

static void run_control(void){
    printf("== k=12 control (n=13, cap 12), p in [167,733] ==\n");
    int nex = 0;
    for (int p = 167; p <= 733; p += 2){
        if (!isprime(p) || p % 13 == 0) continue;
        int r = cover_le(p, 13, 12, 4000000000LL);
        if (r == 0){ printf("  p=%d EXTINCT for k=12 (I(12,p,1) EMPTY -- S-T retrodiction)\n", p); nex++; }
        else if (r == -1) printf("  p=%d UNKNOWN(budget)\n", p);
        else printf("  p=%d alive nodes=%lld\n", p, nodes);
    }
    printf("  k=12 extinct in [167,733]: %d\n", nex);
}

int main(int argc, char **argv){
    long long BUD = 4000000000LL;
    setvbuf(stdout, NULL, _IONBF, 0);
    if (argc >= 4 && !strcmp(argv[1], "exact")){
        run_exact(atoi(argv[2]), atoi(argv[3]), BUD); return 0;
    }
    if (argc >= 4 && !strcmp(argv[1], "decide")){
        long long b = (argc >= 5) ? atoll(argv[4]) : BUD;
        run_decide(atoi(argv[2]), atoi(argv[3]), b); return 0;
    }
    if (argc >= 2 && !strcmp(argv[1], "control")){ run_control(); return 0; }

    int gp[9]  = {29,43,61,71,101,113,127,151,173};
    int gc[9]  = { 7, 9, 9, 9, 10, 10, 10, 11, 11};
    printf("== gate: reproduce HYP-7955 exact c(p) ==\n");
    int ok = 1;
    for (int i = 0; i < 9; i++){
        int p = gp[i], c = -1;
        int lb = ((p-1)/2 + p/14 - 1) / (p/14);
        for (int s = lb; s <= 13; s++){
            int r = cover_le(p, 14, s, BUD);
            if (r == 1) { c = s; break; }
            if (r == -1) { c = -2; break; }
        }
        printf("  p=%d c=%d want=%d %s\n", p, c, gc[i], c==gc[i] ? "OK" : "FAIL");
        if (c != gc[i]) ok = 0;
    }
    if (!ok) { printf("!! GATES FAILED\n"); return 1; }
    run_exact(179, 460, BUD);
    run_decide(461, 1200, BUD);
    run_control();
    return 0;
}

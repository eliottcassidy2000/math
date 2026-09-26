/* procgen_price_20260925 -- the sequential certification construction G_L (explicit infinite pairing members).
 *
 * Pairing family.  Offset 0 (Collatz sheet): pairs {2i-1, 2i}, pair(v) = ceil(v/2).  Offset 1 (3n-1 sheet): pairs
 * {2i, 2i+1}, pair(v) = floor(v/2).  A pair has one bit eps_i; v goes UP iff (v odd) XOR eps_{pair(v)}.
 *   q = 3: up/down by the shared length i = pair index:  v -> v + i  or  v - i   (T: eps = 0, offset 0;  U = 3n-1: eps = 0,
 *          offset 1).
 *   q = 5 (drift control, offset 0 only): up v -> (5v + [v odd])/2, down v -> floor(v/2)  (5x+1: eps = 0).
 * P_L: every n >= n0 falls strictly below itself within L steps.
 *
 * Construction (processing n = n0, n0+1, ... in increasing order; every bit is either FREE or FROZEN):
 *   default path of n = path with frozen bits, free bits read as 0 (the base map).
 *   If it descends within L: freeze the bits it used ("certificate").  Otherwise rescue, trying in order
 *     odd n:  q = 3: [A if n = 3 mod 8 (offset 0) / 5 mod 8 (offset 1)] F, B, A;   q = 5: A, then DFS.
 *       A = flip pair(n) (n goes down at once);  B = flip pair(up(n)) (n -> up(n) -> down: 2 steps);
 *       F = flip the first odd up-image y < 2n on the default path whose partner's up-image is even.
 *     even n (a flipped partner): P = flip pair(up(n)) if up(n) is odd;  then F.
 *   then a depth-first search over free bits on the path (<= MAXF flips).  If nothing works: STUCK (counted).
 *   Every accepted certificate is frozen.  A frozen certificate is never changed, so every n <= N that is not STUCK
 *   descends within L steps in the final member, and the bits of pairs <= N/2 are final (rescues of n > N only touch
 *   pairs of points >= n).
 * Output: density of flipped pairs among pairs <= N/2, STUCK count, rescue statistics, descent-time histogram of the
 * final member (independent re-run), and (q = 3) LEMMA CHECKS of the theorem in
 * 05-knowledge/results/procgen_price_20260925_provability_price.md section 2: for L >= 8 every rescued n is odd, in the
 * bad residue mod 4, Collatz-bad; B is used only for n = 3 mod 8 (5 mod 8 on the 3n-1 sheet) with parity prefix 110110110;
 * the flipped pair index lies in (n/2, n]; partners of A/F flips descend in 2 steps and partners of B flips in <= 8.
 * (Proved in the note; P and DFS are then never used and the construction never gets stuck.)
 * Build: cc -O2 -o greedy procgen_price_20260925_greedy.c
 * Usage: greedy L N offset q [MAXF]      Memory: dense table 16N bytes + hash table (<= 2^24 entries).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef unsigned long long u64;
static int L, OFF, Q, MAXF = 3;
static u64 N, DENSE;
/* rescue records (lemma checks) */
static u64 *recN, *recI; static int *recT; static u64 recK = 0, recCap = 0;
static void rec(u64 n, u64 i, int t) { if (recK == recCap) { recCap = recCap ? 2 * recCap : 1 << 20; recN = realloc(recN, recCap * 8); recI = realloc(recI, recCap * 8); recT = realloc(recT, recCap * sizeof(int)); } recN[recK] = n; recI[recK] = i; recT[recK] = t; recK++; }
static signed char *dense;              /* -1 free, 0/1 frozen */
static u64 *hkey; static signed char *hval; static u64 HCAP, HCNT;

static inline u64 hsh(u64 k) { k ^= k >> 33; k *= 0xff51afd7ed558ccdULL; k ^= k >> 33; k *= 0xc4ceb9fe1a85ec53ULL; k ^= k >> 33; return k; }
static int getb(u64 i) {
    if (i < DENSE) return dense[i];
    u64 h = hsh(i) & (HCAP - 1);
    while (hkey[h]) { if (hkey[h] == i) return hval[h]; h = (h + 1) & (HCAP - 1); }
    return -1;
}
static void setb(u64 i, int b) {
    if (i < DENSE) { dense[i] = (signed char)b; return; }
    u64 h = hsh(i) & (HCAP - 1);
    while (hkey[h]) { if (hkey[h] == i) { hval[h] = (signed char)b; return; } h = (h + 1) & (HCAP - 1); }
    if (HCNT * 10 > HCAP * 7) { fprintf(stderr, "hash full\n"); exit(2); }
    hkey[h] = i; hval[h] = (signed char)b; HCNT++;
}
static inline u64 pairof(u64 v) { return OFF ? (v >> 1) : ((v + 1) >> 1); }
static inline u64 upof(u64 v) {
    if (Q == 3) { u64 i = pairof(v); return v + i; }
    return (v & 1) ? (5 * v + 1) / 2 : 5 * (v / 2);
}
static inline u64 downof(u64 v) {
    if (Q == 3) { u64 i = pairof(v); return v - i; }
    return v >> 1;
}
static inline u64 stepb(u64 v, int b) { return (((v & 1) ^ (u64)b) & 1) ? upof(v) : downof(v); }

/* extension: small list of tentative assignments */
typedef struct { u64 i[16]; int b[16]; int k; } Ext;
static int extget(const Ext *e, u64 i) { for (int j = 0; j < e->k; j++) if (e->i[j] == i) return e->b[j]; return -1; }
static int bitat(const Ext *e, u64 i) { int b = getb(i); if (b >= 0) return b; b = e ? extget(e, i) : -1; return b >= 0 ? b : 0; }
static int isfree(const Ext *e, u64 i) { return getb(i) < 0 && (e == NULL || extget(e, i) < 0); }

static u64 usedI[64]; static int usedB[64]; static int usedK;
static int run(u64 n, const Ext *e) {  /* returns steps to descend (1..L) or 0 */
    u64 v = n; usedK = 0;
    for (int s = 1; s <= L; s++) {
        u64 i = pairof(v); int b = bitat(e, i);
        usedI[usedK] = i; usedB[usedK] = b; usedK++;
        v = stepb(v, b);
        if (v < n) return s;
    }
    return 0;
}
static void freeze_used(void) { for (int j = 0; j < usedK; j++) setb(usedI[j], usedB[j]); }

/* DFS over free bits */
static Ext best; static int bestf;
static void dfs(u64 n, u64 v, int depth, Ext *e, int nf) {
    if (bestf >= 0 && nf >= bestf) return;
    if (depth == L) return;
    u64 i = pairof(v); int fb = getb(i), eb = extget(e, i);
    int opts[2], no = 0;
    if (fb >= 0) opts[no++] = fb; else if (eb >= 0) opts[no++] = eb; else { opts[no++] = 0; opts[no++] = 1; }
    for (int t = 0; t < no; t++) {
        int b = opts[t]; int fresh = (fb < 0 && eb < 0);
        int nf2 = nf + (fresh && b == 1);
        if (nf2 > MAXF) continue;
        if (fresh) { e->i[e->k] = i; e->b[e->k] = b; e->k++; }
        u64 w = stepb(v, b);
        if (w < n) { if (bestf < 0 || nf2 < bestf) { bestf = nf2; best = *e; } }
        else if (e->k < 15) dfs(n, w, depth + 1, e, nf2);
        if (fresh) e->k--;
    }
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: greedy L N offset q [MAXF]\n"); return 1; }
    L = atoi(argv[1]); N = strtoull(argv[2], 0, 10); OFF = atoi(argv[3]); Q = atoi(argv[4]);
    if (argc > 5) MAXF = atoi(argv[5]);
    DENSE = 16 * N + 64;
    dense = malloc(DENSE); memset(dense, -1, DENSE);
    HCAP = 1ULL << 24; hkey = calloc(HCAP, sizeof(u64)); hval = calloc(HCAP, 1);
    u64 n0 = OFF ? 2 : 3;
    long long stuck = 0, cnt[8] = {0}; const char *nm[8] = {"A", "B", "F", "P", "dfs", "", "", ""};
    u64 firststuck[8]; int nfs = 0;
    for (u64 n = n0; n <= N; n++) {
        if (run(n, NULL)) { freeze_used(); continue; }
        Ext cand[4]; int nc = 0; int tag[4];
        u64 un = upof(n);
        /* option A */
        Ext eA = {{0}, {0}, 0}; int okA = 0;
        if (n & 1) { u64 i = pairof(n); if (isfree(NULL, i)) { eA.i[0] = i; eA.b[0] = 1; eA.k = 1; okA = 1; } }
        /* option B / P : flip pair of up(n) when up(n) is odd, n goes up by default */
        Ext eB = {{0}, {0}, 0}; int okB = 0;
        if (Q == 3) {
            int nb = bitat(NULL, pairof(n)); int goesup = ((n & 1) ^ nb) & 1;
            u64 r = goesup ? un : 0;
            if (goesup && (r & 1)) { u64 i = pairof(r); if (isfree(NULL, i)) { eB.i[0] = i; eB.b[0] = 1; eB.k = 1; okB = 1; } }
        }
        /* option F: first odd up-image y (reached from odd point by an up move), y < 2n, partner's up-image even */
        Ext eF = {{0}, {0}, 0}; int okF = 0;
        if (Q == 3) {
            u64 v = n;
            for (int s = 0; s < L; s++) {
                u64 i = pairof(v); int b = bitat(NULL, i); u64 w = stepb(v, b);
                int upmove = ((v & 1) ^ b) & 1;
                if (upmove && (v & 1) && (w & 1) && w < 2 * n) {
                    u64 pw = OFF ? w - 1 : w + 1;       /* partner of odd w */
                    u64 pu = upof(pw);                  /* partner goes up when the pair is flipped */
                    u64 iw = pairof(w);
                    if (!(pu & 1) && isfree(NULL, iw)) { eF.i[0] = iw; eF.b[0] = 1; eF.k = 1; okF = 1; break; }
                }
                v = w;
            }
        }
        int preferA = (n & 1) && ((OFF == 0 && n % 8 == 3) || (OFF == 1 && n % 8 == 5));
        if (Q == 5) { if (okA) { cand[nc] = eA; tag[nc++] = 0; } }
        else if (n & 1) {
            if (preferA && okA) { cand[nc] = eA; tag[nc++] = 0; }
            if (okF) { cand[nc] = eF; tag[nc++] = 2; }
            if (okB) { cand[nc] = eB; tag[nc++] = 1; }
            if (!preferA && okA) { cand[nc] = eA; tag[nc++] = 0; }
        } else {
            if (okB) { cand[nc] = eB; tag[nc++] = 3; }
            if (okF) { cand[nc] = eF; tag[nc++] = 2; }
        }
        int done = 0;
        for (int c = 0; c < nc && !done; c++) {
            if (run(n, &cand[c])) { freeze_used(); for (int j = 0; j < cand[c].k; j++) setb(cand[c].i[j], cand[c].b[j]); cnt[tag[c]]++; done = 1; rec(n, cand[c].i[0], tag[c]); }
        }
        if (done) continue;
        Ext e = {{0}, {0}, 0}; bestf = -1; dfs(n, n, 0, &e, 0);
        if (bestf >= 0) {
            if (!run(n, &best)) { fprintf(stderr, "internal\n"); return 3; }
            freeze_used(); for (int j = 0; j < best.k; j++) setb(best.i[j], best.b[j]);
            cnt[4]++; continue;
        }
        stuck++; if (nfs < 8) firststuck[nfs++] = n;
    }
    /* density on pairs <= N/2 */
    u64 H = N / 2, fl = 0;
    for (u64 i = 1; i <= H; i++) if (getb(i) == 1) fl++;
    /* independent re-run: descent times of n in [n0, N] under the final member (free bits read 0) */
    long long hist[80] = {0}, fail = 0; int maxt = 0;
    for (u64 n = n0; n <= N; n++) {
        u64 v = n; int t = 0, ok = 0;
        for (t = 1; t <= 3 * L + 20; t++) { v = stepb(v, bitat(NULL, pairof(v))); if (v < n) { ok = 1; break; } }
        if (!ok) { fail++; continue; }
        if (t < 80) hist[t]++;
        if (t > maxt) maxt = t;
    }
    /* small n: orbits of 1..n0-1 (cycles) */
    printf("L=%d N=%llu offset=%d q=%d MAXF=%d: flipped pairs <= N/2: %llu  density %.6f  stuck %lld",
           L, N, OFF, Q, MAXF, fl, (double)fl / H, stuck);
    for (int j = 0; j < nfs; j++) printf(" %llu", firststuck[j]);
    printf("\n  rescues:");
    for (int j = 0; j < 5; j++) printf(" %s=%lld", nm[j], cnt[j]);
    printf("  hash entries %llu\n  re-run: n in [%llu,N] not descending within %d steps: %lld; max descent time %d; histogram:",
           HCNT, n0, 3 * L + 20, fail, maxt);
    for (int t = 1; t <= maxt && t < 80; t++) if (hist[t]) printf(" %d:%lld", t, hist[t]);

    if (Q == 3) {
        /* Lemma checks (offset 0; offset 1 with U = 3n-1 and the residues 1 mod 4 / 5 mod 8): every rescued n is odd, = 3 mod 4, has a non-descending T-path of length L;
           B only for n = 3 mod 8 with T-word prefix 110110110; flipped pair index i satisfies n/2 < i <= n;
           partner 2i of each flip descends: A/F partners in exactly 2 steps, B partners in <= 8 steps. */
        long long bad_even = 0, bad_mod4 = 0, bad_tdesc = 0, bad_B8 = 0, bad_Bprefix = 0, bad_idx = 0;
        long long pAF[70] = {0}, pB[70] = {0}; long long npart = 0;
        for (u64 r = 0; r < recK; r++) {
            u64 n = recN[r], i = recI[r]; int t = recT[r];
            if (!(n & 1)) bad_even++;
            if ((n & 3) != (OFF ? 1u : 3u)) bad_mod4++;
            { u64 v = n; int d = 0; for (int s = 1; s <= L; s++) { v = (v & 1) ? (OFF ? (3 * v - 1) / 2 : (3 * v + 1) / 2) : v / 2; if (v < n) { d = 1; break; } } if (d) bad_tdesc++; }
            if (t == 1) {
                if (n % 8 != (OFF ? 5u : 3u)) bad_B8++;
                u64 v = n; const char *w = "110110110"; for (int s = 0; s < 9; s++) { if ((int)(v & 1) != w[s] - '0') { bad_Bprefix++; break; } v = (v & 1) ? (OFF ? (3 * v - 1) / 2 : (3 * v + 1) / 2) : v / 2; }
            }
            if (!(2 * i + OFF > n && i <= n)) bad_idx++;
            u64 p = 2 * i; if (p > N) continue;          /* partner = even member of the flipped pair */
            u64 v = p; int tt; for (tt = 1; tt < 69; tt++) { v = stepb(v, bitat(NULL, pairof(v))); if (v < p) break; }
            npart++; if (t == 1) pB[tt]++; else pAF[tt]++;
        }
        printf("\n  LEMMA CHECKS over %llu rescues: rescued even %lld; rescued n off the bad residue mod 4 (3; 1 for 3n-1) %lld; rescued n whose base-map path descends within L %lld;"
               " B with n off 3 mod 8 (5 for 3n-1) %lld; B without prefix 110110110 %lld; flipped index outside (n/2, n] %lld\n",
               recK, bad_even, bad_mod4, bad_tdesc, bad_B8, bad_Bprefix, bad_idx);
        printf("  partner descent times (partners <= N: %lld):  A/F:", npart);
        for (int t = 1; t < 70; t++) if (pAF[t]) printf(" %d:%lld", t, pAF[t]);
        printf("   B:");
        for (int t = 1; t < 70; t++) if (pB[t]) printf(" %d:%lld", t, pB[t]);
        printf("\n");
    }
    printf("\n  pair(1) bit = %d (offset 0: 0 means {1,2} is the root cycle)\n", getb(1));
    return 0;
}

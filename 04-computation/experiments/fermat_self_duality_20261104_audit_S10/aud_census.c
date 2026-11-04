/* aud_census.c -- independent audit code (fermat_self_duality_20261104_audit_S10).
 * Written from scratch for the audit; shares no code with selfdual3_fewcycle_20261104.c.
 *
 * Part 1 (orders 3..8): enumerate ALL labelled tournaments of order n (2^C(n,2) codes) and split them into
 *   isomorphism classes by orbit marking: scan codes in increasing order; an unvisited code is the minimum of a new
 *   orbit; apply all n! relabellings and mark the images.  |Aut| = n!/|orbit|.  No canonical-form heuristic is used.
 *   Self-checks: class counts 2,4,12,56,456,6880 (A000568) and sum of orbit sizes = 2^C(n,2).
 *   Per class: H (Hamiltonian paths), hc (directed Hamiltonian cycles), scores, |Aut|, self-converse flag,
 *   {-3}-self-dual flag (n >= 6; every T - X, |X| = 3, self-converse, via the orbit tables of order n - 3),
 *   decomposability and the types of all modules, and for min out-degree >= 2 with hc in {2,..,5} the number of
 *   v with hc(T - v) = 1.
 * Part 2 (orders 7..10): a pruned labelled search for {-3}-self-dual tournaments.  Start from a representative of
 *   every self-converse class of order m = N - 3 on vertices 0..m-1, add vertices m..N-1 one at a time with every
 *   arc pattern, and keep a partial tournament only if every m-subset containing the new vertex is self-converse.
 *   Complete at the class level: every {-3}-self-dual T of order N has T - X self-converse, so some relabelling of T
 *   starts with a class representative.  Each survivor is tested against explicit TT_N / AT_N recognisers;
 *   at N = 7, 8 the survivors are also reduced to classes (class table at 7, min over all 8! relabellings at 8).
 * Encoding: pair (i,j), i<j, has bit index j(j-1)/2 + i; bit = 1 iff i -> j.  The converse is the bitwise complement.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

typedef uint32_t u32;
typedef unsigned long long u64;

static inline int pidx(int i, int j) { return j * (j - 1) / 2 + i; }   /* i < j */
static int E_(int n) { return n * (n - 1) / 2; }

static void code_to_adj(int n, u32 code, u32 *out) {
    for (int v = 0; v < n; v++) out[v] = 0;
    for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) {
        if ((code >> pidx(i, j)) & 1u) out[i] |= 1u << j; else out[j] |= 1u << i;
    }
}
/* code of the subtournament induced on the sorted vertex list s[0..k-1] */
static u32 sub_code(const u32 *out, const int *s, int k) {
    u32 c = 0;
    for (int q = 1; q < k; q++) for (int p = 0; p < q; p++) if ((out[s[p]] >> s[q]) & 1u) c |= 1u << pidx(p, q);
    return c;
}

/* ---------------- permutations and their action on codes ---------------- */
typedef struct { u32 flip; u32 tab[4][128]; } PT;
static int nperm; static PT *ptab;
static void build_perms(int n) {
    int fact = 1; for (int i = 2; i <= n; i++) fact *= i;
    nperm = fact; ptab = malloc(sizeof(PT) * (size_t)fact);
    int p[16]; for (int i = 0; i < n; i++) p[i] = i;
    /* lexicographic enumeration of permutations */
    for (int k = 0; k < fact; k++) {
        PT *t = &ptab[k]; memset(t, 0, sizeof(PT));
        for (int j = 1; j < n; j++) for (int i = 0; i < j; i++) {
            int s = pidx(i, j), a = p[i], b = p[j], tgt, fl;
            if (a < b) { tgt = pidx(a, b); fl = 0; } else { tgt = pidx(b, a); fl = 1; }
            if (fl) t->flip |= 1u << tgt;
            int ch = s / 7, bit = s % 7;
            for (int x = 0; x < 128; x++) if ((x >> bit) & 1) t->tab[ch][x] |= 1u << tgt;
        }
        /* next permutation */
        int i = n - 2; while (i >= 0 && p[i] > p[i + 1]) i--;
        if (i < 0) break;
        int j = n - 1; while (p[j] < p[i]) j--;
        int tmp = p[i]; p[i] = p[j]; p[j] = tmp;
        for (int a = i + 1, b = n - 1; a < b; a++, b--) { tmp = p[a]; p[a] = p[b]; p[b] = tmp; }
    }
}
static inline u32 apply_perm(const PT *t, u32 c) {
    return t->flip ^ t->tab[0][c & 127] ^ t->tab[1][(c >> 7) & 127] ^ t->tab[2][(c >> 14) & 127] ^ t->tab[3][(c >> 21) & 127];
}

/* ---------------- Hamiltonian paths and cycles (own subset DP) ---------------- */
static u64 dpbuf[(1 << 11) * 11];
static u64 count_H(int n, const u32 *out) {
    int full = (1 << n) - 1;
    memset(dpbuf, 0, sizeof(u64) * (size_t)(full + 1) * n);
    for (int v = 0; v < n; v++) dpbuf[(1 << v) * n + v] = 1;
    for (int S = 1; S <= full; S++) for (int v = 0; v < n; v++) {
        u64 c = dpbuf[S * n + v]; if (!c) continue;
        for (int w = 0; w < n; w++) if (!((S >> w) & 1) && ((out[v] >> w) & 1u)) dpbuf[(S | (1 << w)) * n + w] += c;
    }
    u64 h = 0; for (int v = 0; v < n; v++) h += dpbuf[full * n + v];
    return h;
}
static u64 count_hc(int n, const u32 *out) {
    if (n < 3) return 0;
    int full = (1 << n) - 1;
    memset(dpbuf, 0, sizeof(u64) * (size_t)(full + 1) * n);
    dpbuf[1 * n + 0] = 1;                       /* paths starting at vertex 0 */
    for (int S = 1; S <= full; S += 2) for (int v = 0; v < n; v++) {
        u64 c = dpbuf[S * n + v]; if (!c) continue;
        for (int w = 0; w < n; w++) if (!((S >> w) & 1) && ((out[v] >> w) & 1u)) dpbuf[(S | (1 << w)) * n + w] += c;
    }
    u64 hc = 0; for (int v = 1; v < n; v++) if ((out[v] >> 0) & 1u) hc += dpbuf[full * n + v];
    return hc;
}
static void delete_v(int n, const u32 *out, int x, u32 *res) {
    int k = 0;
    for (int i = 0; i < n; i++) {
        if (i == x) continue;
        u32 r = 0; int kk = 0;
        for (int j = 0; j < n; j++) { if (j == x) continue; if ((out[i] >> j) & 1u) r |= 1u << kk; kk++; }
        res[k++] = r;
    }
}

/* ---------------- class tables for orders 3..7 ---------------- */
static uint16_t *cls[8];     /* class id of every labelled code, orders 3..7 */
static uint8_t *scT[8];      /* self-converse flag of every labelled code, orders 1..7 */
static int nclass[9];
static uint8_t *scbits8;            /* self-converse flag of every labelled 8-tournament, one bit per code */
static u32 *screps8; static int nscreps8;   /* one representative per self-converse class of order 8 */
static inline int sc_lookup(int m, u32 code) {
    if (m <= 7) return scT[m][code];
    return (scbits8[code >> 3] >> (code & 7)) & 1;
}

static int is_module(int n, const u32 *out, u32 M) {
    for (int x = 0; x < n; x++) if (!((M >> x) & 1u)) {
        u32 o = out[x] & M, in = M & ~out[x];
        if (o && in) return 0;
    }
    return 1;
}
static int is_transitive_sub(const u32 *out, u32 M) {
    /* induced subtournament on M transitive iff its in-M scores are pairwise distinct */
    int seen[16] = {0};
    for (int v = 0; v < 16; v++) if ((M >> v) & 1u) {
        int s = __builtin_popcount(out[v] & M);
        if (seen[s]) return 0;
        seen[s] = 1;
    }
    return 1;
}
static int is_sd3(int n, const u32 *out) {
    int m = n - 3, s[16];
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) for (int c = b + 1; c < n; c++) {
        int k = 0;
        for (int v = 0; v < n; v++) if (v != a && v != b && v != c) s[k++] = v;
        if (!sc_lookup(m, sub_code(out, s, m))) return 0;
    }
    return 1;
}
static int is_TT(int n, const u32 *out) {
    int seen[16] = {0};
    for (int v = 0; v < n; v++) { int s = __builtin_popcount(out[v]); if (seen[s]) return 0; seen[s] = 1; }
    return 1;
}
/* AT_n: order a, b, middle..., c, d with every arc forward except d -> a */
static int is_AT(int n, const u32 *out) {
    if (n < 4) return 0;
    int sc[16], top[2], nt = 0, bot[2], nb = 0, mid[16], bys[16];
    for (int s = 0; s < 16; s++) bys[s] = -1;
    for (int v = 0; v < n; v++) {
        sc[v] = __builtin_popcount(out[v]);
        if (sc[v] == n - 2) { if (nt == 2) return 0; top[nt++] = v; }
        else if (sc[v] == 1) { if (nb == 2) return 0; bot[nb++] = v; }
        else { if (sc[v] < 2 || sc[v] > n - 3 || bys[sc[v]] >= 0) return 0; bys[sc[v]] = v; }
    }
    if (nt != 2 || nb != 2) return 0;
    int ord[16], k = 0;
    if ((out[top[0]] >> top[1]) & 1u) { ord[k++] = top[0]; ord[k++] = top[1]; } else { ord[k++] = top[1]; ord[k++] = top[0]; }
    for (int s = n - 3; s >= 2; s--) { if (bys[s] < 0) return 0; ord[k++] = bys[s]; }
    if ((out[bot[0]] >> bot[1]) & 1u) { ord[k++] = bot[0]; ord[k++] = bot[1]; } else { ord[k++] = bot[1]; ord[k++] = bot[0]; }
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        int fwd = (out[ord[i]] >> ord[j]) & 1u;
        if (i == 0 && j == n - 1) { if (fwd) return 0; } else if (!fwd) return 0;
    }
    (void)mid;
    return 1;
}

/* ---------------- Part 1: the census ---------------- */
static long Hdist8[5000];
static void census(int n) {
    clock_t t0 = clock();
    int E = E_(n); u64 total = 1ull << E;
    build_perms(n);
    uint8_t *vis = calloc(total / 8 + 1, 1);
    if (n <= 7) { cls[n] = malloc(sizeof(uint16_t) * total); scT[n] = malloc(total); }
    if (n == 8) { scbits8 = calloc(total / 8 + 1, 1); screps8 = malloc(sizeof(u32) * 7000); nscreps8 = 0; }
    u32 *reps = malloc(sizeof(u32) * 7000); int *aut = malloc(sizeof(int) * 7000);
    int nc = 0; u64 tot_orbit = 0;
    for (u64 c = 0; c < total; c++) {
        if ((vis[c >> 3] >> (c & 7)) & 1) continue;
        int cnt = 0, sc8 = 0;
        if (n == 8) {   /* c is the minimum of its orbit; the class is self-converse iff the converse orbit has the same minimum */
            u32 cc = (u32)c ^ ((1u << E) - 1), mn = 0xffffffffu;
            for (int k = 0; k < nperm; k++) { u32 img = apply_perm(&ptab[k], cc); if (img < mn) mn = img; }
            sc8 = (mn == (u32)c);
            if (sc8) screps8[nscreps8++] = (u32)c;
        }
        for (int k = 0; k < nperm; k++) {
            u32 img = apply_perm(&ptab[k], (u32)c);
            if (!((vis[img >> 3] >> (img & 7)) & 1)) { vis[img >> 3] |= (uint8_t)(1 << (img & 7)); cnt++; if (n <= 7) cls[n][img] = (uint16_t)nc;
                if (sc8) scbits8[img >> 3] |= (uint8_t)(1 << (img & 7)); }
        }
        reps[nc] = (u32)c; aut[nc] = nperm / cnt; tot_orbit += cnt; nc++;
    }
    nclass[n] = nc;
    u32 mask = (E == 32) ? 0xffffffffu : ((1u << E) - 1);
    /* self-converse: is the complement code in the same orbit? */
    int *sc = malloc(sizeof(int) * nc); int nsc = 0;
    for (int i = 0; i < nc; i++) {
        u32 cc = reps[i] ^ mask;
        if (n <= 7) sc[i] = cls[n][cc] == i;
        else sc[i] = (scbits8[reps[i] >> 3] >> (reps[i] & 7)) & 1;
        nsc += sc[i];
    }
    if (n <= 7) for (u64 c = 0; c < total; c++) scT[n][c] = (uint8_t)sc[cls[n][c]];
    printf("\n=== order %d: %d classes (sum of orbit sizes %llu = 2^%d: %s), %d self-converse classes [%.1f s]\n",
           n, nc, tot_orbit, E, tot_orbit == total ? "ok" : "MISMATCH", nsc, (double)(clock() - t0) / CLOCKS_PER_SEC);
    /* per-class properties */
    u64 maxH = 0; long nsd3 = 0;
    long del_hist[6][16]; memset(del_hist, 0, sizeof del_hist); long del_cls[6]; memset(del_cls, 0, sizeof del_cls);
    long Hcount[700]; memset(Hcount, 0, sizeof Hcount);
    for (int i = 0; i < nc; i++) {
        u32 out[16]; code_to_adj(n, reps[i], out);
        u64 H = count_H(n, out), hc = count_hc(n, out);
        if (H > maxH) maxH = H;
        if (n == 8 && H < 5000) Hdist8[H]++;
        if (H < 700) Hcount[H]++;
        int sco[16], mo = 99; for (int v = 0; v < n; v++) { sco[v] = __builtin_popcount(out[v]); if (sco[v] < mo) mo = sco[v]; }
        /* deletion law */
        if (mo >= 2 && hc >= 2 && hc <= 5) {
            int k = 0;
            for (int v = 0; v < n; v++) { u32 r[16]; delete_v(n, out, v, r); if (count_hc(n - 1, r) == 1) k++; }
            del_hist[hc][k]++; del_cls[hc]++;
        }
        if (n >= 6 && is_sd3(n, out)) {
            nsd3++;
            int dec = 0, nmod = 0, nontr = 0, maxtt = 0;
            for (u32 M = 1; M < (1u << n) - 1; M++) {
                int pc = __builtin_popcount(M);
                if (pc < 2 || pc > n - 1 || !is_module(n, out, M)) continue;
                dec = 1; nmod++;
                if (!is_transitive_sub(out, M)) nontr++; else if (pc > maxtt) maxtt = pc;
            }
            int ss[16]; memcpy(ss, sco, sizeof(int) * n);
            for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if (ss[b] < ss[a]) { int t = ss[a]; ss[a] = ss[b]; ss[b] = t; }
            printf("SD3 n=%d dec=%d sc=%d H=%llu hc=%llu aut=%d TT=%d AT=%d modules=%d nontransitive_modules=%d largest_TT_module=%d scores:",
                   n, dec, sc[i], H, hc, aut[i], is_TT(n, out), is_AT(n, out), nmod, nontr, maxtt);
            for (int v = 0; v < n; v++) printf(" %d", ss[v]);
            printf(" code=%u\n", reps[i]);
        }
    }
    printf("order %d: max H = %llu; {-3}-self-dual classes: %ld\n", n, maxH, nsd3);
    for (int hc = 2; hc <= 5; hc++) if (del_cls[hc]) {
        printf("DEL n=%d hc=%d sole-exit-free classes=%ld  #(v: hc(T-v)=1) ->", n, hc, del_cls[hc]);
        for (int k = 0; k < 16; k++) if (del_hist[hc][k]) printf(" %d:%ld", k, del_hist[hc][k]);
        printf("\n");
    }
    printf("N(%d,h) for h = 3,5,9,15,17,31,33,63,65,127,129,255,257,511,513:", n);
    int hs[] = {3, 5, 9, 15, 17, 31, 33, 63, 65, 127, 129, 255, 257, 511, 513};
    for (int q = 0; q < 15; q++) printf(" %d:%ld", hs[q], hs[q] < 700 ? Hcount[hs[q]] : -1L);
    printf("\n");
    if (n == 6) { /* the two classes with H = 17 */
        for (int i = 0; i < nc; i++) { u32 out[16]; code_to_adj(n, reps[i], out); if (count_H(n, out) == 17) {
            int sco[16]; for (int v = 0; v < n; v++) sco[v] = __builtin_popcount(out[v]);
            printf("H=17 at order 6: hc=%llu AT=%d scores(unsorted):", count_hc(n, out), is_AT(n, out)); for (int v = 0; v < n; v++) printf(" %d", sco[v]); printf("\n"); } }
    }
    fflush(stdout);
    free(vis); free(ptab); free(reps); free(aut); free(sc);
}

/* ---------------- Part 2: pruned labelled search for {-3}-self-dual tournaments of order N ---------------- */
static int N_, M_;
static u32 sout[16];
static long long nsurv, nTT, nAT, nother, nodes[16];
static u32 *surv_codes; static long long nsurv_cap;
static int subsets[200][16]; static int nsubs[16];   /* (m-1)-subsets of {0..k-1}, per k */
static int comb_buf[16];
static void gen_subsets(int k, int r, int start, int depth, int *cnt) {
    if (depth == r) { memcpy(subsets[*cnt], comb_buf, sizeof(int) * r); (*cnt)++; return; }
    for (int v = start; v < k; v++) { comb_buf[depth] = v; gen_subsets(k, r, v + 1, depth + 1, cnt); }
}
static int (*subs_k[16])[16]; static int nsub_k[16];
static void extend(int k) {
    nodes[k]++;
    if (k == N_) {
        nsurv++;
        if (is_sd3(N_, sout) != 1) { printf("internal error: survivor not {-3}-self-dual\n"); exit(1); }
        int tt = is_TT(N_, sout), at = is_AT(N_, sout);
        if (tt) nTT++; else if (at) nAT++; else {
            nother++;
            if (nother <= 5) { printf("  OTHER survivor rows:"); for (int v = 0; v < N_; v++) printf(" %x", sout[v]); printf("\n"); }
        }
        if (surv_codes && nsurv <= nsurv_cap) { u32 c = 0; for (int j = 1; j < N_; j++) for (int i = 0; i < j; i++) if ((sout[i] >> j) & 1u) c |= 1u << pidx(i, j); surv_codes[nsurv - 1] = c; }
        return;
    }
    u32 save[16]; memcpy(save, sout, sizeof(u32) * k);
    for (u32 mask = 0; mask < (1u << k); mask++) {
        /* mask bit i = 1: new vertex k beats i */
        for (int i = 0; i < k; i++) sout[i] = save[i] | (((mask >> i) & 1u) ? 0u : (1u << k));
        sout[k] = mask;
        int ok = 1, s[16];
        for (int q = 0; q < nsub_k[k] && ok; q++) {
            for (int t = 0; t < M_ - 1; t++) s[t] = subs_k[k][q][t];
            s[M_ - 1] = k;
            if (!sc_lookup(M_, sub_code(sout, s, M_))) ok = 0;
        }
        if (ok) extend(k + 1);
    }
    memcpy(sout, save, sizeof(u32) * k);
}
static void search(int N) {
    clock_t t0 = clock();
    N_ = N; M_ = N - 3;
    nsurv = nTT = nAT = nother = 0; memset(nodes, 0, sizeof nodes);
    for (int k = M_; k < N; k++) {
        int cnt = 0; static int store[16][200][16];
        int tmpcnt = 0; gen_subsets(k, M_ - 1, 0, 0, &tmpcnt);
        for (int q = 0; q < tmpcnt; q++) memcpy(store[k][q], subsets[q], sizeof(int) * 16);
        cnt = tmpcnt; subs_k[k] = store[k]; nsub_k[k] = cnt;
    }
    /* one representative per self-converse class of order M_ */
    int nbase = 0;
    if (M_ == 8) {
        for (int i = 0; i < nscreps8; i++) { nbase++; code_to_adj(8, screps8[i], sout); extend(8); }
    } else {
        int E = E_(M_); u64 total = 1ull << E;
        int *seen = calloc(nclass[M_] + 1, sizeof(int));
        for (u64 c = 0; c < total; c++) {
            int id = cls[M_][c];
            if (seen[id] || !scT[M_][c]) continue;
            seen[id] = 1; nbase++;
            code_to_adj(M_, (u32)c, sout);
            extend(M_);
        }
        free(seen);
    }
    printf("SEARCH N=%d (m=%d): %d self-converse base classes; nodes per level:", N, M_, nbase);
    for (int k = M_; k <= N; k++) printf(" %lld", nodes[k]);
    printf("\n  labelled survivors (all {-3}-self-dual, re-tested): %lld = TT_%d: %lld, AT_%d: %lld, other: %lld  [%.1f s]\n",
           nsurv, N, nTT, N, nAT, nother, (double)(clock() - t0) / CLOCKS_PER_SEC);
    fflush(stdout);
}

static int cmpu(const void *a, const void *b) { u32 x = *(const u32 *)a, y = *(const u32 *)b; return x < y ? -1 : x > y; }

int main(int argc, char **argv) {
    int nmax = argc > 1 ? atoi(argv[1]) : 8;
    int smax = argc > 2 ? atoi(argv[2]) : 10;
    printf("aud_census (independent audit code): census to order %d, {-3}-self-dual search to order %d\n", nmax, smax);
    /* trivial orders 1, 2 */
    scT[1] = malloc(1); scT[1][0] = 1; scT[2] = malloc(2); scT[2][0] = scT[2][1] = 1; nclass[1] = nclass[2] = 1;
    for (int n = 3; n <= nmax; n++) census(n);
    /* Part 2.  Orders 7 and 8 validate the search against the census; 9 and 10 are the audit targets. */
    for (int N = 7; N <= smax; N++) {
        if (N == 7 || N == 8) { nsurv_cap = 2000000; surv_codes = malloc(sizeof(u32) * nsurv_cap); } else surv_codes = NULL;
        search(N);
        if (N == 7) {
            int *hit = calloc(nclass[7], sizeof(int)); int d = 0;
            for (long long i = 0; i < nsurv && i < nsurv_cap; i++) { int id = cls[7][surv_codes[i]]; if (!hit[id]) { hit[id] = 1; d++; } }
            printf("  order 7: survivors fall into %d distinct classes (class table)\n", d); free(hit);
        }
        if (N == 8) {
            build_perms(8);
            long long m = nsurv < nsurv_cap ? nsurv : nsurv_cap;
            u32 *canon = malloc(sizeof(u32) * m);
            for (long long i = 0; i < m; i++) { u32 mn = 0xffffffffu; for (int k = 0; k < nperm; k++) { u32 img = apply_perm(&ptab[k], surv_codes[i]); if (img < mn) mn = img; } canon[i] = mn; }
            qsort(canon, m, sizeof(u32), cmpu);
            long long d = 0; for (long long i = 0; i < m; i++) if (i == 0 || canon[i] != canon[i - 1]) d++;
            printf("  order 8: survivors fall into %lld distinct classes (minimum code over all 8! relabellings)\n", d);
            free(canon); free(ptab);
        }
        if (surv_codes) { free(surv_codes); surv_codes = NULL; }
    }
    return 0;
}

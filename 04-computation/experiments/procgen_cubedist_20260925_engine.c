/*
 * procgen_cubedist_20260925_engine.c -- fast exact routines for the strategy cube of the
 * m n +- 1 maps (m = 3: Collatz / 3n-1; m = 5: 5n+-1), cube-distance lane, 2026-09-25.
 *
 * A level-k sign strategy sigma : odd residues mod 2^k -> {+1,-1};
 *   T(n) = n/2 (n even),  (m n + sigma(n mod 2^k))/2 (n odd).
 * Parity graph G_sigma: nodes Z/2^k, s -> the two lifts mod 2^k of T(s) mod 2^(k-1).
 * A cycle with a odd nodes and length p is expanding iff m^a > 2^p.
 *
 * Input (stdin):  line 1  "<mode> k m [args]"
 *                 line 2  flip mask as a hex string, least significant nibble LAST
 *                         (Python format(mask,'x')); bit i <-> residue 2i+1; bit = 1 means sigma = -1.
 * Modes:
 *   karp  k m          exact maximum cycle odd-density (Karp, O(n) memory, two passes);
 *                      prints "D a p" (a/p reduced).  Weights are integer odd counts, exact.
 *   cert  k m q r      threshold q/r: odd weight r-q, even weight -q (a cycle is positive iff a/p > q/r).
 *                      Computes psi(s) = sup over walks from s of partial sums (value iteration on the
 *                      reversed graph with a worklist).  If finite prints "OK maxpsi" and then psi(s) for
 *                      s = 0..2^k-1 (one line, space separated); this is an integer certificate:
 *                      psi(t) <= psi(s) - w(s) on every edge s->t.  Otherwise prints "CYCLE p a" and the cycle.
 *   cycles k m q r M   up to M node-disjoint cycles with a/p > q/r (repeated detection + deletion of the
 *                      nodes of each cycle found); prints "C v1,v2,..." per cycle, then "END".
 *   short k m q r P C F   every simple cycle of length <= P with a/p > q/r (DFS rooted at its minimum node,
 *                      each cycle once), at most C of them; F = 1: only cycles through a node with sigma = -.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int K, N, H, MUL;
static unsigned char *flip;   /* by residue */
static int *succ0, *succ1;    /* the two successors */
static unsigned char *alive;

static void build(void) {
    succ0 = malloc(sizeof(int) * N); succ1 = malloc(sizeof(int) * N);
    for (int s = 0; s < N; s++) {
        long long t;
        if (!(s & 1)) t = s >> 1;
        else {
            long long v = (long long)MUL * s + (flip[s] ? -1 : 1);
            t = v >> 1;      /* v even and positive */
        }
        int t0 = (int)(t % H);
        succ0[s] = t0; succ1[s] = t0 + H;
    }
}

static void read_mask(void) {
    flip = calloc(N, 1);
    int cap = (N / 8) + 64;
    char *buf = malloc(cap + 2);
    if (!fgets(buf, cap + 2, stdin)) { buf[0] = 0; }
    int len = (int)strlen(buf);
    while (len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r' || buf[len-1] == ' ')) len--;
    /* hex, most significant first */
    for (int i = 0; i < len; i++) {
        char c = buf[len - 1 - i];
        int v;
        if (c >= '0' && c <= '9') v = c - '0';
        else if (c >= 'a' && c <= 'f') v = c - 'a' + 10;
        else if (c >= 'A' && c <= 'F') v = c - 'A' + 10;
        else { fprintf(stderr, "bad hex char\n"); exit(2); }
        for (int b = 0; b < 4; b++) if ((v >> b) & 1) {
            long long idx = 4LL * i + b;       /* bit index -> residue 2 idx + 1 */
            if (2 * idx + 1 >= N) { fprintf(stderr, "mask too long\n"); exit(2); }
            flip[2 * idx + 1] = 1;
        }
    }
    free(buf);
}

static long long gcdll(long long a, long long b) { if (a < 0) a = -a; if (b < 0) b = -b; while (b) { long long t = a % b; a = b; b = t; } return a; }

/* ---------------------------------------------------------------- Karp (vertex weight = oddness of source) */
static void karp(void) {
    const int NEG = -1000000000;
    int *Dp = malloc(sizeof(int) * N), *Dc = malloc(sizeof(int) * N), *Dn = malloc(sizeof(int) * N);
    for (int v = 0; v < N; v++) Dp[v] = 0;
    for (int j = 1; j <= N; j++) {
        for (int v = 0; v < N; v++) Dc[v] = NEG;
        for (int u = 0; u < N; u++) {
            if (Dp[u] == NEG) continue;
            int val = Dp[u] + (u & 1);
            if (val > Dc[succ0[u]]) Dc[succ0[u]] = val;
            if (val > Dc[succ1[u]]) Dc[succ1[u]] = val;
        }
        int *t = Dp; Dp = Dc; Dc = t;
    }
    memcpy(Dn, Dp, sizeof(int) * N);
    /* best[v] = min_j (Dn[v]-Dj[v])/(N-j) as fraction bn/bd */
    long long *bn = malloc(sizeof(long long) * N), *bd = malloc(sizeof(long long) * N);
    for (int v = 0; v < N; v++) { bn[v] = 1; bd[v] = 0; }   /* +inf */
    for (int v = 0; v < N; v++) Dp[v] = 0;
    for (int j = 0; j < N; j++) {
        if (j > 0) {
            for (int v = 0; v < N; v++) Dc[v] = NEG;
            for (int u = 0; u < N; u++) {
                if (Dp[u] == NEG) continue;
                int val = Dp[u] + (u & 1);
                if (val > Dc[succ0[u]]) Dc[succ0[u]] = val;
                if (val > Dc[succ1[u]]) Dc[succ1[u]] = val;
            }
            int *t = Dp; Dp = Dc; Dc = t;
        }
        for (int v = 0; v < N; v++) {
            if (Dn[v] == NEG || Dp[v] == NEG) continue;
            long long num = (long long)Dn[v] - Dp[v], den = N - j;
            /* compare num/den < bn/bd */
            if (bd[v] == 0 || num * bd[v] < bn[v] * den) { bn[v] = num; bd[v] = den; }
        }
    }
    long long An = -1, Ad = 1;
    for (int v = 0; v < N; v++) {
        if (Dn[v] == NEG || bd[v] == 0) continue;
        if (An < 0 || bn[v] * Ad > An * bd[v]) { An = bn[v]; Ad = bd[v]; }
    }
    long long g = gcdll(An, Ad); if (g == 0) g = 1;
    printf("D %lld %lld\n", An / g, Ad / g);
    free(Dp); free(Dc); free(Dn); free(bn); free(bd);
}

/* ---------------------------------------------------------------- value iteration psi(s)=max(0, w(s)+max psi(succ)) */
static int *rhead, *rnext, *rsrc;   /* reversed adjacency: for each t list of (s) */
static void build_reverse(void) {
    rhead = malloc(sizeof(int) * N); rnext = malloc(sizeof(int) * 2 * N); rsrc = malloc(sizeof(int) * 2 * N);
    for (int t = 0; t < N; t++) rhead[t] = -1;
    int e = 0;
    for (int s = 0; s < N; s++) {
        int ts[2] = { succ0[s], succ1[s] };
        for (int i = 0; i < 2; i++) { rsrc[e] = s; rnext[e] = rhead[ts[i]]; rhead[ts[i]] = e; e++; }
    }
}

/* scan the pointer graph s -> best[s] for a cycle; returns length and fills cyc, or 0 */
static int pointer_cycle(const int *best, int *cyc) {
    int *mark = malloc(sizeof(int) * N);
    for (int s = 0; s < N; s++) mark[s] = -1;
    int found = 0;
    for (int s0 = 0; s0 < N && !found; s0++) {
        if (mark[s0] != -1) continue;
        int v = s0;
        while (v >= 0 && mark[v] == -1) { mark[v] = s0; v = best[v]; }
        if (v >= 0 && mark[v] == s0) {          /* new cycle through v */
            int u = v, L = 0;
            do { cyc[L++] = u; u = best[u]; } while (u != v);
            found = L;
        }
    }
    free(mark);
    return found;
}

/* returns 1 if finite (psi filled), 0 if a positive cycle exists (cycle written to cyc, length *clen) */
static int value_iter(long long wo, long long we, long long *psi, int *cyc, int *clen) {
    int *best = malloc(sizeof(int) * N);       /* argmax successor */
    int *queue = malloc(sizeof(int) * (N + 1));
    unsigned char *inq = calloc(N, 1);
    for (int s = 0; s < N; s++) { psi[s] = 0; best[s] = -1; }
    int qh = 0, qt = 0, qn = 0;
    for (int t = 0; t < N; t++) if (alive[t]) { queue[qt] = t; qt = (qt + 1) % (N + 1); qn++; inq[t] = 1; }
    long long relax = 0, nextscan = N;
    int ok = 1;
    while (qn > 0) {
        int t = queue[qh]; qh = (qh + 1) % (N + 1); qn--; inq[t] = 0;
        for (int e = rhead[t]; e != -1; e = rnext[e]) {
            int s = rsrc[e];
            if (!alive[s]) continue;
            long long ws = (s & 1) ? wo : we;
            long long cand = ws + psi[t];
            if (cand > psi[s]) {
                psi[s] = cand; best[s] = t; relax++;
                if (!inq[s]) { queue[qt] = s; qt = (qt + 1) % (N + 1); qn++; inq[s] = 1; }
            }
        }
        if (relax >= nextscan) {
            nextscan = relax + N;
            int L = pointer_cycle(best, cyc);
            if (L) { *clen = L; ok = 0; break; }
        }
    }
    if (ok) {  /* final sanity: pointer graph must be acyclic */
        int L = pointer_cycle(best, cyc);
        if (L) { *clen = L; ok = 0; }
    }
    if (!ok) {
        long long W = 0; for (int i = 0; i < *clen; i++) W += (cyc[i] & 1) ? wo : we;
        if (W <= 0) { fprintf(stderr, "extracted cycle not positive\n"); exit(3); }
    }
    free(best); free(queue); free(inq);
    return ok;
}

/* ---------------------------------------------------------------- short positive cycles (DFS, rooted at min node) */
static int SP_P, SP_cap, SP_count, SP_needflip;
static long long SP_wo, SP_we;
static int *SP_path;
static unsigned char *SP_on;
static void sp_dfs(int root, int v, int depth, long long w) {
    if (SP_count >= SP_cap) return;
    long long wv = (v & 1) ? SP_wo : SP_we;
    long long w2 = w + wv;
    int ts[2] = { succ0[v], succ1[v] };
    for (int i = 0; i < 2; i++) {
        int t = ts[i];
        if (i == 1 && ts[1] == ts[0]) break;
        if (t == root) {
            if (w2 > 0) {
                int hasflip = 0;
                for (int j = 0; j <= depth; j++) if ((SP_path[j] & 1) && flip[SP_path[j]]) hasflip = 1;
                if (!SP_needflip || hasflip) {
                    printf("C ");
                    for (int j = 0; j <= depth; j++) printf("%d%c", SP_path[j], j == depth ? '\n' : ',');
                    SP_count++;
                    if (SP_count >= SP_cap) return;
                }
            }
            continue;
        }
        if (t < root || SP_on[t] || depth + 1 >= SP_P) continue;
        /* prune: remaining steps (including t) can add at most (SP_P - depth - 1) * wo */
        if (w2 + (long long)(SP_P - depth - 1) * SP_wo <= 0) continue;
        SP_on[t] = 1; SP_path[depth + 1] = t;
        sp_dfs(root, t, depth + 1, w2);
        SP_on[t] = 0;
    }
}

int main(void) {
    char mode[32];
    if (scanf("%31s %d %d", mode, &K, &MUL) != 3) { fprintf(stderr, "bad header\n"); return 2; }
    long long q = 0, r = 0; int M = 0;
    if (!strcmp(mode, "cert")) { if (scanf("%lld %lld", &q, &r) != 2) return 2; }
    if (!strcmp(mode, "cycles")) { if (scanf("%lld %lld %d", &q, &r, &M) != 3) return 2; }
    int SPP = 0, SPC = 0, SPF = 0;
    if (!strcmp(mode, "short")) { if (scanf("%lld %lld %d %d %d", &q, &r, &SPP, &SPC, &SPF) != 5) return 2; }
    int c; while ((c = getchar()) != '\n' && c != EOF) {}
    N = 1 << K; H = N >> 1;
    read_mask();
    build();
    alive = malloc(N); memset(alive, 1, N);
    if (!strcmp(mode, "karp")) { karp(); return 0; }
    build_reverse();
    long long wo = r - q, we = -q;
    long long *psi = malloc(sizeof(long long) * N);
    int *cyc = malloc(sizeof(int) * (N + 1)); int clen = 0;
    if (!strcmp(mode, "cert")) {
        if (value_iter(wo, we, psi, cyc, &clen)) {
            long long mx = 0; for (int s = 0; s < N; s++) if (psi[s] > mx) mx = psi[s];
            printf("OK %lld\n", mx);
            for (int s = 0; s < N; s++) printf("%lld%c", psi[s], s == N - 1 ? '\n' : ' ');
        } else {
            int a = 0; for (int i = 0; i < clen; i++) a += cyc[i] & 1;
            printf("CYCLE %d %d\n", clen, a);
            for (int i = 0; i < clen; i++) printf("%d%c", cyc[i], i == clen - 1 ? '\n' : ',');
        }
        return 0;
    }
    if (!strcmp(mode, "cycles")) {
        for (int it = 0; it < M; it++) {
            if (value_iter(wo, we, psi, cyc, &clen)) break;
            /* the extracted cycle: best pointers go s -> best[s], i.e. cyc[i] -> cyc[i+1] */
            printf("C ");
            for (int i = 0; i < clen; i++) printf("%d%c", cyc[i], i == clen - 1 ? '\n' : ',');
            for (int i = 0; i < clen; i++) alive[cyc[i]] = 0;
        }
        printf("END\n");
        return 0;
    }
    if (!strcmp(mode, "short")) {
        SP_P = SPP; SP_cap = SPC; SP_needflip = SPF; SP_count = 0; SP_wo = wo; SP_we = we;
        SP_path = malloc(sizeof(int) * (SPP + 2)); SP_on = calloc(N, 1);
        for (int root = 0; root < N && SP_count < SP_cap; root++) {
            SP_path[0] = root; SP_on[root] = 1;
            sp_dfs(root, root, 0, 0);
            SP_on[root] = 0;
        }
        printf("END\n");
        return 0;
    }
    fprintf(stderr, "unknown mode\n");
    return 2;
}

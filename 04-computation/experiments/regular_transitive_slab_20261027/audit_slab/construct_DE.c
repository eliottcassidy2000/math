/* construct_DE.c -- independent audit of THM-4819 Theorem D (odd n = 2m+1, T regular, both native roles)
 * and Theorem E (even n = 2m, T balanced, the dipole-forced role), by EXECUTING the proof's construction
 * with randomised choices wherever the proof says "any" and checking every claim it makes on the way.
 *
 * For a target (v, role) the construction runs in the reversal gauge (R_P reverses the n-1 path arcs) on
 * B = T or B = T^op (B = T^op for the "other" role in the odd case, and for s(v) = m-1 in the even case):
 *   step 1: a RANDOM Hamilton cycle of B - v; among its rotations (c_t,...,c_{t-1},v), c_{t-1} -> v, pick
 *           a black one (Lemma B / B' says one exists for EVERY cycle); P1, a = first(P1);
 *   step 2: a random Hamilton cycle of U1 - v (claim: U1 - v strong); a random rotation at an in-neighbour of v
 *           whose start differs from a (claim: >= 2 in-neighbours, so one exists); claim: it is black;
 *   steps j >= 3: a random Hamilton path of U_{j-1} ending at v (claim: U_{j-1} strong); claim: it is black,
 *           and indeed every start u has d(u) + d(v) > n-1 at that time;
 *   end (K = m odd, K = m-1 even): claim: v is a source of U_K, U_K - v strong.
 * The R route is converted to the native gauge (T_j = U_j or U_j^op, Q_j = P or rev P by parity), and when
 * B = T^op the native word is conjugated by the converse (Q -> rev Q).  The resulting native word is then
 * REPLAYED from T itself with an independent native-move routine (reverse every arc not on Q): every Q_j must
 * be a directed Hamilton path of the current state, black there, and the final state must have v as a
 * source (score n-1) or sink (score 0) as requested, with T - v strong.
 * Input: gentourng lines on stdin.  Usage: ./construct_DE n trials seed
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef struct { uint32_t o[32]; } Tour;
static int n, M, m, odd;
static uint32_t FULL;

static uint64_t rs = 88172645463325252ull;
static uint32_t rnd(uint32_t k) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (uint32_t)(rs % k); }

static inline int arc(const Tour *T, int a, int b) { return (int)((T->o[a] >> b) & 1u); }
static int score(const Tour *T, int v) { return __builtin_popcount(T->o[v]); }
static Tour converse(const Tour *T) { Tour R; for (int u = 0; u < n; u++) R.o[u] = (~T->o[u]) & FULL & ~(1u << u); return R; }
static Tour rmove(const Tour *T, const int *P) {   /* reverse the path arcs only */
    Tour R = *T;
    for (int i = 0; i + 1 < n; i++) { int a = P[i], b = P[i + 1]; R.o[a] &= ~(1u << b); R.o[b] |= 1u << a; }
    return R;
}
static Tour kmove(const Tour *T, const int *Q) {   /* native move: reverse every arc NOT on Q */
    Tour R;
    for (int u = 0; u < n; u++) R.o[u] = (~T->o[u]) & FULL & ~(1u << u);
    for (int i = 0; i + 1 < n; i++) { int a = Q[i], b = Q[i + 1]; R.o[a] |= 1u << b; R.o[b] &= ~(1u << a); }
    return R;
}
static int is_ham_path(const Tour *T, const int *P) {
    uint32_t seen = 0;
    for (int i = 0; i < n; i++) { if (P[i] < 0 || P[i] >= n || (seen >> P[i] & 1)) return 0; seen |= 1u << P[i]; }
    for (int i = 0; i + 1 < n; i++) if (!arc(T, P[i], P[i + 1])) return 0;
    return seen == FULL;
}
static int is_black(const Tour *T, const int *P) {
    int sig[32];
    for (int i = 0; i < n; i++) sig[P[i]] = P[n - 1 - i];
    for (int x = 0; x < n; x++) for (int y = x + 1; y < n; y++) if (arc(T, x, y) != arc(T, sig[y], sig[x])) return 1;
    return 0;
}
static int strong_on(const Tour *T, uint32_t W) {  /* is the sub-tournament on W strong? */
    if (__builtin_popcount(W) <= 1) return 1;
    int s = __builtin_ctz(W);
    uint32_t fw = 1u << s, fr = 1u << s;
    for (;;) { uint32_t nf = fw; for (int u = 0; u < n; u++) if (fw >> u & 1) nf |= T->o[u] & W; if (nf == fw) break; fw = nf; }
    for (;;) { uint32_t nf = fr; for (int u = 0; u < n; u++) if ((W >> u & 1) && (T->o[u] & fr)) nf |= 1u << u; if (nf == fr) break; fr = nf; }
    return fw == W && fr == W;
}
/* random Hamilton cycle of the sub-tournament on W (bitmask), randomised DFS; returns length or 0 */
static int cyc[32], cl;
static int hc_dfs(const Tour *T, uint32_t W, uint32_t used, int depth) {
    if (depth == cl) return arc(T, cyc[cl - 1], cyc[0]);
    int cand[32], k = 0;
    for (int w = 0; w < n; w++) if ((W >> w & 1) && !(used >> w & 1) && arc(T, cyc[depth - 1], w)) cand[k++] = w;
    for (int i = k - 1; i > 0; i--) { int j = rnd(i + 1); int t = cand[i]; cand[i] = cand[j]; cand[j] = t; }
    for (int i = 0; i < k; i++) { cyc[depth] = cand[i]; if (hc_dfs(T, W, used | (1u << cand[i]), depth + 1)) return 1; }
    return 0;
}
static int ham_cycle(const Tour *T, uint32_t W) {
    cl = __builtin_popcount(W);
    int vs[32], k = 0; for (int w = 0; w < n; w++) if (W >> w & 1) vs[k++] = w;
    cyc[0] = vs[rnd(k)];
    return hc_dfs(T, W, 1u << cyc[0], 1) ? cl : 0;
}
/* random Hamilton path ending at v: backward randomised DFS */
static int pth[32];
static int hp_dfs(const Tour *T, uint32_t used, int depth) {
    if (depth == n) return 1;
    int cand[32], k = 0;
    for (int w = 0; w < n; w++) if (!(used >> w & 1) && arc(T, w, pth[depth - 1])) cand[k++] = w;
    for (int i = k - 1; i > 0; i--) { int j = rnd(i + 1); int t = cand[i]; cand[i] = cand[j]; cand[j] = t; }
    for (int i = 0; i < k; i++) { pth[depth] = cand[i]; if (hp_dfs(T, used | (1u << cand[i]), depth + 1)) return 1; }
    return 0;
}
static int ham_path_ending_at(const Tour *T, int v, int *P) {
    pth[0] = v;
    if (!hp_dfs(T, 1u << v, 1)) return 0;
    for (int i = 0; i < n; i++) P[i] = pth[n - 1 - i];
    return 1;
}

static long long fails[16];
static const char *failname[16] = {
    "step1: Hamilton cycle of B-v missing", "step1: chosen cycle has no black rotation (Lemma B/B' strong form)",
    "step2: U1-v not strong", "step2: fewer than 2 in-neighbours of v / no rotation with start != a",
    "step2: chosen rotation is blue", "later: state not strong before step K",
    "later: endpoint-score sum <= n-1 for some start", "later: chosen path blue",
    "end: v not a source of U_K", "end: U_K - v not strong",
    "replay: Q_j not a Hamilton path of the native state", "replay: Q_j blue in the native state",
    "replay: wrong final role / v not singleton", "replay: final core not strong", "", "" };

/* run the construction for (T, v, role); role 0 = source, 1 = sink; returns 1 on success */
static int run(const Tour *T, int v, int role, int *forced_role_out) {
    int K = odd ? m : m - 1;
    int base_is_op;
    if (odd) {
        int rho0 = (m % 2 == 0) ? 0 : 1;       /* R-gauge source becomes native source iff K = m even */
        base_is_op = (role != rho0);
    } else {
        base_is_op = (score(T, v) == m - 1);
    }
    Tour B = base_is_op ? converse(T) : *T;
    if (score(&B, v) != m) { fprintf(stderr, "internal: v must have score m in B\n"); exit(5); }
    Tour U[32]; int P[32][32];
    U[0] = B;
    /* step 1 */
    uint32_t W = FULL & ~(1u << v);
    if (!ham_cycle(&U[0], W)) { fails[0]++; return 0; }
    int C[32]; memcpy(C, cyc, sizeof(int) * M);
    int blacks[32], nb = 0;
    for (int t = 0; t < M; t++) {
        if (!arc(&U[0], C[(t - 1 + M) % M], v)) continue;
        int Pt[32]; for (int i = 0; i < M; i++) Pt[i] = C[(t + i) % M]; Pt[M] = v;
        if (!is_ham_path(&U[0], Pt)) { fprintf(stderr, "internal path\n"); exit(6); }
        if (is_black(&U[0], Pt)) blacks[nb++] = t;
    }
    if (!nb) { fails[1]++; return 0; }
    { int t = blacks[rnd(nb)]; for (int i = 0; i < M; i++) P[1][i] = C[(t + i) % M]; P[1][M] = v; }
    int a = P[1][0];
    U[1] = rmove(&U[0], P[1]);
    /* step 2 */
    if (K >= 2) {
        if (!strong_on(&U[1], W)) { fails[2]++; return 0; }
        if (!ham_cycle(&U[1], W)) { fails[2]++; return 0; }
        memcpy(C, cyc, sizeof(int) * M);
        int opts[32], no = 0, nin = 0;
        for (int t = 0; t < M; t++) {
            if (!arc(&U[1], C[(t - 1 + M) % M], v)) continue;
            nin++;
            if (C[t] != a) opts[no++] = t;
        }
        if (nin < 2 || !no) { fails[3]++; return 0; }
        int t = opts[rnd(no)];
        for (int i = 0; i < M; i++) P[2][i] = C[(t + i) % M]; P[2][M] = v;
        if (!is_ham_path(&U[1], P[2])) { fprintf(stderr, "internal path2\n"); exit(6); }
        if (!is_black(&U[1], P[2])) { fails[4]++; return 0; }
        U[2] = rmove(&U[1], P[2]);
    }
    /* later steps */
    for (int j = 3; j <= K; j++) {
        Tour *S = &U[j - 1];
        if (!strong_on(S, FULL)) { fails[5]++; return 0; }
        int dv = score(S, v);
        for (int u = 0; u < n; u++) if (u != v && score(S, u) + dv <= n - 1) { fails[6]++; return 0; }
        if (!ham_path_ending_at(S, v, P[j])) { fails[5]++; return 0; }
        if (!is_ham_path(S, P[j])) { fprintf(stderr, "internal path3\n"); exit(6); }
        if (!is_black(S, P[j])) { fails[7]++; return 0; }
        U[j] = rmove(S, P[j]);
    }
    if (score(&U[K], v) != n - 1) { fails[8]++; return 0; }
    if (!strong_on(&U[K], W)) { fails[9]++; return 0; }
    /* native conversion: native move j (j = 0..K-1) from T_j uses R-path P[j+1] */
    int Q[32][32];
    for (int j = 0; j < K; j++) {
        for (int i = 0; i < n; i++) Q[j][i] = (j % 2 == 0) ? P[j + 1][i] : P[j + 1][n - 1 - i];
        if (base_is_op) { int tmp[32]; for (int i = 0; i < n; i++) tmp[i] = Q[j][n - 1 - i]; memcpy(Q[j], tmp, sizeof tmp); }
    }
    /* replay from T */
    Tour S = *T;
    for (int j = 0; j < K; j++) {
        if (!is_ham_path(&S, Q[j])) { fails[10]++; return 0; }
        if (!is_black(&S, Q[j])) { fails[11]++; return 0; }
        S = kmove(&S, Q[j]);
    }
    int sv = score(&S, v);
    int got = (sv == n - 1) ? 0 : (sv == 0 ? 1 : -1);
    if (forced_role_out) *forced_role_out = got;
    if (got != role) { fails[12]++; return 0; }
    if (!strong_on(&S, W)) { fails[13]++; return 0; }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s n trials seed\n", argv[0]); return 1; }
    n = atoi(argv[1]); int trials = atoi(argv[2]); rs ^= (uint64_t)atoll(argv[3]) * 0x9E3779B97F4A7C15ull;
    M = n - 1; odd = n & 1; m = odd ? (n - 1) / 2 : n / 2; FULL = (n == 32) ? 0xffffffffu : ((1u << n) - 1);
    char line[2048]; long ntour = 0; long long inst = 0, ok = 0;
    long long rolecheck_bad = 0;
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] != '0' && line[0] != '1') continue;
        Tour T; memset(&T, 0, sizeof T);
        int k = 0;
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') T.o[i] |= 1u << j; else T.o[j] |= 1u << i; k++; }
        for (int v = 0; v < n; v++) {
            int s = score(&T, v);
            if (odd ? s != m : (s != m && s != m - 1)) { fprintf(stderr, "input not regular/balanced\n"); return 2; }
        }
        ntour++;
        for (int v = 0; v < n; v++) {
            int roles[2], nr = 0;
            if (odd) { roles[nr++] = 0; roles[nr++] = 1; }
            else {
                /* dipole-forced role: score m -> source iff m odd; score m-1 -> the opposite */
                int s = score(&T, v);
                int r = (s == m) ? ((m % 2 == 1) ? 0 : 1) : ((m % 2 == 1) ? 1 : 0);
                roles[nr++] = r;
            }
            for (int r = 0; r < nr; r++) {
                for (int tr = 0; tr < trials; tr++) {
                    inst++;
                    int got = -2;
                    if (run(&T, v, roles[r], &got)) ok++;
                    else if (got != -2 && got != roles[r]) rolecheck_bad++;
                }
            }
        }
    }
    printf("n=%d m=%d: tournaments %ld; (vertex, role, random trial) runs %lld; construction succeeded and replayed: %lld; failed: %lld\n",
           n, m, ntour, inst, ok, inst - ok);
    for (int i = 0; i < 14; i++) if (fails[i]) printf("   failure kind '%s': %lld\n", failname[i], fails[i]);
    return 0;
}

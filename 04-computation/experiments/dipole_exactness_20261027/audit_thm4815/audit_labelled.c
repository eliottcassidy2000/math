/* audit_labelled.c -- independent audit of THM-4815 section C (forward-arc parity law, P_7 distances, the 6-move word).
 * Written from scratch for the audit.  Labelled tournaments on {0..n-1} encoded by C(n,2) bits: bit k = 1 iff i -> j
 * for the k-th pair (i<j) in row order (the gentourng order), so F(T) = popcount(code).
 * Moves: tau_Q(T) = T^op xor E(Q) for every Hamiltonian path Q of T.
 * Mode "parity n": over every labelled tournament and every Hamiltonian path check
 *     F' = C(n,2) + n - 1 - F - 2 back(Q)      (back(Q) = number of path arcs q_i -> q_{i+1} with q_i > q_{i+1})
 *   and A1 identities: tau_Q(tau_Q T) = T; (tau_Q T)^op = tau_{Q^rev}(T^op); tau_Q(T) = rho_{Q^rev}(T^op).
 *   Then BFS components of the labelled move graph and bipartiteness.
 * Mode "p7": BFS from the Paley labelling of Z_7 (i -> i+1, i+2, i+4) over all 2^21 states: distance to TT_id,
 *   distance to the transitive class (minimum over all 5,040 transitive labellings) and its distribution; replay of the
 *   six-move word of THM-4815 (C).
 * Build: cc -O3 -march=native audit_labelled.c -o audit_labelled ; Run: ./audit_labelled parity 6 ; ./audit_labelled p7
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int n, NP, pidx[16][16];
static void decode(uint32_t code, uint32_t *o) { for (int i = 0; i < n; i++) o[i] = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (code >> pidx[i][j] & 1) o[i] |= 1u << j; else o[j] |= 1u << i; } }
static uint32_t encode(const uint32_t *o) { uint32_t c = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) if (o[i] >> j & 1) c |= 1u << pidx[i][j]; return c; }
static uint32_t conv(uint32_t code) { return (~code) & ((NP == 32) ? 0xffffffffu : ((1u << NP) - 1)); }
static uint32_t tau_code(uint32_t code, const int *Q) {   /* converse, then restore the arcs of Q */
    uint32_t c = conv(code);
    for (int i = 0; i + 1 < n; i++) { int u = Q[i], v = Q[i + 1]; int a = u < v ? u : v, b = u < v ? v : u; int k = pidx[a][b];
        if (u < v) c |= 1u << k; else c &= ~(1u << k); }
    return c;
}
static uint32_t rho_code(uint32_t code, const int *P) {   /* reverse the arcs of the Hamiltonian path P of code */
    uint32_t c = code;
    for (int i = 0; i + 1 < n; i++) { int u = P[i], v = P[i + 1]; int a = u < v ? u : v, b = u < v ? v : u; c ^= 1u << pidx[a][b]; }
    return c;
}
static int isHP(uint32_t code, const int *Q) {
    uint32_t o[16]; decode(code, o); uint32_t seen = 0;
    for (int i = 0; i < n; i++) { if (Q[i] < 0 || Q[i] >= n || (seen >> Q[i] & 1)) return 0; seen |= 1u << Q[i]; }
    for (int i = 0; i + 1 < n; i++) if (!(o[Q[i]] >> Q[i + 1] & 1)) return 0;
    return 1;
}
/* enumerate Hamiltonian paths; callback */
typedef void (*pathcb)(uint32_t code, const int *Q, void *ctx);
static void paths(uint32_t code, pathcb cb, void *ctx) {
    uint32_t o[16]; decode(code, o);
    int Q[16]; uint32_t cand[16];
    for (int st = 0; st < n; st++) {
        int depth = 0; Q[0] = st; uint32_t used = 1u << st; cand[0] = o[st] & ~used;
        while (depth >= 0) {
            if (depth == n - 1) { cb(code, Q, ctx); used &= ~(1u << Q[depth]); depth--; continue; }
            if (!cand[depth]) { used &= ~(1u << Q[depth]); depth--; continue; }
            int v = __builtin_ctz(cand[depth]); cand[depth] &= cand[depth] - 1;
            depth++; Q[depth] = v; used |= 1u << v; cand[depth] = o[v] & ~used;
        }
    }
}
static long parity_fail, invol_fail, conv_fail, rho_fail, nmoves;
static void cb_parity(uint32_t code, const int *Q, void *ctx) {
    (void)ctx; nmoves++;
    uint32_t t = tau_code(code, Q);
    int F = __builtin_popcount(code), Fp = __builtin_popcount(t), back = 0;
    for (int i = 0; i + 1 < n; i++) if (Q[i] > Q[i + 1]) back++;
    if (Fp != NP + n - 1 - F - 2 * back) parity_fail++;
    if (!isHP(t, Q) || tau_code(t, Q) != code) invol_fail++;
    int Qr[16]; for (int i = 0; i < n; i++) Qr[i] = Q[n - 1 - i];
    uint32_t c_op = conv(code);
    if (!isHP(c_op, Qr) || conv(t) != tau_code(c_op, Qr)) conv_fail++;
    if (rho_code(c_op, Qr) != t) rho_fail++;
}
static int8_t *dist; static uint32_t *queue; static long qt;
static void cb_bfs(uint32_t code, const int *Q, void *ctx) { int d = *(int *)ctx; uint32_t t = tau_code(code, Q); if (dist[t] < 0) { dist[t] = d + 1; queue[qt++] = t; } }
static long bip_fail;
static void cb_bip(uint32_t code, const int *Q, void *ctx) { (void)ctx; uint32_t t = tau_code(code, Q); if (((dist[t] ^ dist[code]) & 1) == 0) bip_fail++; }
static int is_transitive(uint32_t code) { uint32_t o[16]; decode(code, o); uint32_t seen = 0; for (int v = 0; v < n; v++) { int s = __builtin_popcount(o[v]); if (seen >> s & 1) return 0; seen |= 1u << s; } return 1; }


static int black_only = 0;
static int isBlueQ(uint32_t code, const int *Q) { uint32_t o[16]; decode(code, o); int r[16]; for (int i = 0; i < n; i++) r[Q[i]] = Q[n - 1 - i];
    for (int u = 0; u < n; u++) for (uint32_t x = o[u]; x; x &= x - 1) { int v = __builtin_ctz(x); if (!(o[r[v]] >> r[u] & 1)) return 0; } return 1; }
static void cb_bfs2(uint32_t code, const int *Q, void *ctx) { int d = *(int *)ctx; if (black_only && isBlueQ(code, Q)) return; uint32_t t = tau_code(code, Q); if (dist[t] < 0) { dist[t] = d + 1; queue[qt++] = t; } }
static int is_regular(uint32_t code) { uint32_t o[16]; decode(code, o); for (int v = 0; v < n; v++) if (2 * __builtin_popcount(o[v]) != n - 1) return 0; return 1; }
int main(int argc, char **argv) {
    if (argc < 2) return 1;

    if (!strcmp(argv[1], "from")) {
        n = atoi(argv[2]); NP = n * (n - 1) / 2; int kk = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) pidx[i][j] = kk++;
        uint32_t NS2 = 1u << NP; black_only = !strcmp(argv[4], "black"); int reg = !strcmp(argv[5], "reg");
        uint32_t start = 0; for (int b = 0; b < NP; b++) if (argv[3][b] == '1') start |= 1u << b;
        dist = malloc(NS2); queue = malloc(sizeof(uint32_t) * NS2); memset(dist, -1, NS2);
        long qh2 = 0; qt = 0; dist[start] = 0; queue[qt++] = start;
        while (qh2 < qt) { uint32_t c = queue[qh2++]; int d = dist[c]; paths(c, cb_bfs2, &d); }
        int best = 99; long cnt = 0; for (uint32_t c = 0; c < NS2; c++) if (dist[c] >= 0 && (reg ? is_regular(c) : is_transitive(c))) { cnt++; if (dist[c] < best) best = dist[c]; }
        uint32_t o2[16]; decode(start, o2); int s[16]; for (int v = 0; v < n; v++) s[v] = __builtin_popcount(o2[v]);
        printf("from %s (n=%d, %s moves): reached %ld states; target %s labellings reached %ld; distance to the %s class = %d; scores", argv[3], n, black_only ? "black" : "any", qt, reg ? "regular" : "transitive", cnt, reg ? "regular" : "transitive", best == 99 ? -1 : best);
        for (int v = 0; v < n; v++) printf(" %d", s[v]); printf("\n");
        return 0;
    }
    int p7 = !strcmp(argv[1], "p7");
    n = p7 ? 7 : atoi(argv[2]);
    NP = n * (n - 1) / 2; int k = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) pidx[i][j] = k++;
    uint32_t NS = 1u << NP;
    if (!p7) {
        for (uint32_t c = 0; c < NS; c++) paths(c, cb_parity, NULL);
        printf("n=%d labelled tournaments %u, moves %ld: parity-formula failures %ld; involution/HP failures %ld; converse-covariance failures %ld; rho failures %ld\n", n, NS, nmoves, parity_fail, invol_fail, conv_fail, rho_fail);
        /* components and bipartiteness */
        dist = malloc(NS); queue = malloc(sizeof(uint32_t) * NS); memset(dist, -1, NS);
        int comps = 0; long sizes[8] = {0};
        for (uint32_t s0 = 0; s0 < NS; s0++) if (dist[s0] < 0) {
            long qh = 0; qt = 0; dist[s0] = 0; queue[qt++] = s0;
            while (qh < qt) { uint32_t c = queue[qh++]; int d = dist[c]; paths(c, cb_bfs, &d); }
            if (comps < 8) sizes[comps] = qt; comps++;
        }
        bip_fail = 0; for (uint32_t c = 0; c < NS; c++) paths(c, cb_bip, NULL);
        printf("   labelled move graph: %d components (sizes:", comps); for (int i = 0; i < comps && i < 8; i++) printf(" %ld", sizes[i]);
        printf("); edges joining same-parity BFS levels (non-bipartite witnesses): %ld; (n-1)(n+2)/2 mod 2 = %d\n", bip_fail, ((n - 1) * (n + 2) / 2) & 1);
        /* does every component have constant F parity (n = 1,2 mod 4)? */
        return 0;
    }
    /* P_7 */
    uint32_t o[16]; for (int i = 0; i < 7; i++) { o[i] = 0; for (int d = 1; d <= 4; d *= 2) o[i] |= 1u << ((i + d) % 7); }
    uint32_t P7 = encode(o), TTid = (1u << NP) - 1;
    printf("P_7 code %u, F(P_7) = %d, F(TT_id) = %d\n", P7, __builtin_popcount(P7), __builtin_popcount(TTid));
    dist = malloc(NS); queue = malloc(sizeof(uint32_t) * NS); memset(dist, -1, NS);
    long qh = 0; qt = 0; dist[P7] = 0; queue[qt++] = P7;
    while (qh < qt) { uint32_t c = queue[qh++]; int d = dist[c]; paths(c, cb_bfs, &d); }
    long reached = qt; int dmin = 99; long hist[64] = {0}; long ntrans = 0;
    for (uint32_t c = 0; c < NS; c++) if (is_transitive(c)) { ntrans++; if (dist[c] >= 0) { hist[(int)dist[c]]++; if (dist[c] < dmin) dmin = dist[c]; } }
    int maxd = 0; for (uint32_t c = 0; c < NS; c++) if (dist[c] > maxd) maxd = dist[c];
    printf("BFS from P_7: reached %ld of %u labelled tournaments; eccentricity %d\n", reached, NS, maxd);
    printf("d_lab(P_7, TT_id) = %d; transitive labellings %ld; distance to the transitive class = %d; distribution:", dist[TTid], ntrans, dmin);
    for (int d = 0; d < 64; d++) if (hist[d]) printf(" d=%d:%ld", d, hist[d]); printf("\n");
    bip_fail = 0; for (uint32_t c = 0; c < NS; c++) paths(c, cb_bip, NULL);
    printf("edges joining same-parity levels (from P_7 BFS): %ld\n", bip_fail);
    /* parity-class check for all transitive labellings at the min distance */
    long parity_mismatch = 0; for (uint32_t c = 0; c < NS; c++) if (dist[c] >= 0 && ((dist[c] + __builtin_popcount(c) + __builtin_popcount(P7)) & 1)) parity_mismatch++;
    printf("states with d_lab(P_7, U) != F(P_7) + F(U) mod 2: %ld\n", parity_mismatch);
    /* replay the six-move word of THM-4815 (C) */
    const char *word[6] = { "0123456", "1203645", "0352641", "2316405", "4053162", "2163405" };
    uint32_t cur = P7; int okall = 1;
    for (int st = 0; st < 6; st++) {
        int Q[16]; for (int i = 0; i < 7; i++) Q[i] = word[st][i] - '0';
        int hp = isHP(cur, Q);
        uint32_t oo[16]; decode(cur, oo); printf("  step %d: F=%2d scores", st, __builtin_popcount(cur)); for (int v = 0; v < 7; v++) printf(" %d", __builtin_popcount(oo[v]));
        printf("  path %s is a Hamiltonian path: %s\n", word[st], hp ? "yes" : "NO"); if (!hp) { okall = 0; break; }
        cur = tau_code(cur, Q);
    }
    uint32_t oo[16]; decode(cur, oo); printf("  end: F=%2d scores", __builtin_popcount(cur)); for (int v = 0; v < 7; v++) printf(" %d", __builtin_popcount(oo[v]));
    printf("  transitive: %s; d_lab(P_7, end) = %d\n", is_transitive(cur) ? "yes" : "NO", dist[cur]);
    printf("six-move word replay: %s\n", (okall && is_transitive(cur)) ? "VALID" : "INVALID");
    return 0;
}

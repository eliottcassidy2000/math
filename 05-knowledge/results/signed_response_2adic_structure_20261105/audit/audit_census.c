/* audit_census.c -- independent adversarial audit of THM-4883 (own code; nothing imported from the packet).
 *
 * Modes (all output is deterministic text on stdout):
 *
 *   audit_census census LOG2N
 *       stdin: lines "a b M cstar L low W1 W2" (M = 2^(s+1) a power of two, cstar = predicted full class or -1,
 *       L = T2 suffix length or 0, low = required low L bits of n - 1, W1 = window for the list of unhit C* targets,
 *       W2 = window for the hit bitmap / density counts).
 *       For every odd n with 3 <= n < 2^LOG2N it evaluates F = a n - b g(n), g by its defining recurrence
 *       g(1) = 1, g(2m) = g(m), g(2m - 1) = g(m) + m - (m mod 2), and reports per residue class mod M:
 *       count, min (first n attaining), max (first n attaining), distinct hits in [500, W2] and [-W2, -500],
 *       distinct hits in [-W1, W1]; the T2-suffix census; and the unhit targets of class cstar in [-W1, W1].
 *
 *   audit_census graph
 *       stdin: "nv root" then nv lines, line i = out-neighbour bitmask of vertex i.
 *       Prints "D E U d s C": Hamilton paths, paths ending at root, ordered covers (P, Q) with Q nonempty starting at
 *       root and P possibly empty, Hamilton paths of T - root, Hamilton paths starting at root, directed Hamilton
 *       cycles (up to rotation).  Subset dynamic programming, nv <= 22.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

typedef long long ll;
typedef unsigned long long ull;

static int32_t *G;

static void build_g(ll N) {
    G = (int32_t *)malloc(sizeof(int32_t) * (size_t)(N + 1));
    if (!G) { fprintf(stderr, "oom\n"); exit(2); }
    G[0] = 0;
    G[1] = 1;
    for (ll n = 2; n <= N; n++) {
        if ((n & 1) == 0) G[n] = G[n >> 1];
        else { ll m = (n + 1) >> 1; G[n] = (int32_t)(G[m] + (m - (m & 1))); }
    }
}

static int census(int log2n) {
    ll N = 1LL << log2n;
    build_g(N);
    /* sanity: g odd, g <= n, 2g(n) >= n + 3 for odd n >= 3 */
    ll bad = 0;
    for (ll n = 1; n < N; n++) {
        if ((G[n] & 1) == 0 || G[n] > n) bad++;
        if ((n & 1) && n >= 3 && 2LL * G[n] < n + 3) bad++;
    }
    printf("LEMMA1_BAD %lld N %lld\n", bad, N);
    ll a, b, M, cstar, L, low, W1, W2;
    while (scanf("%lld %lld %lld %lld %lld %lld %lld %lld", &a, &b, &M, &cstar, &L, &low, &W1, &W2) == 8) {
        ll *cnt = calloc(M, sizeof(ll)), *mn = malloc(M * sizeof(ll)), *mx = malloc(M * sizeof(ll));
        ll *amn = calloc(M, sizeof(ll)), *amx = calloc(M, sizeof(ll));
        for (ll c = 0; c < M; c++) { mn[c] = LLONG_MAX; mx[c] = LLONG_MIN; }
        unsigned char *bm = calloc((size_t)(2 * W2 + 1), 1);
        ll sufc = 0, sufv = 0, lowmask = (L > 0) ? ((1LL << L) - 1) : 0;
        for (ll n = 3; n < N; n += 2) {
            ll F = a * n - b * (ll)G[n];
            ll c = F & (M - 1);
            cnt[c]++;
            if (F < mn[c]) { mn[c] = F; amn[c] = n; }
            if (F > mx[c]) { mx[c] = F; amx[c] = n; }
            if (F >= -W2 && F <= W2) bm[F + W2] = 1;
            if (L > 0 && (((n - 1) & lowmask) == low) && (n - 1) >= (1LL << L)) {
                sufc++;
                if (c != cstar) sufv++;
            }
        }
        printf("SLOPE %lld %lld M %lld N %lld\n", a, b, M, N);
        for (ll c = 0; c < M; c++) {
            if (!cnt[c]) { printf("C %lld 0\n", c); continue; }
            ll hp = 0, hn = 0, hw = 0;
            for (ll r = 500; r <= W2; r++) if ((r & (M - 1)) == c && bm[r + W2]) hp++;
            for (ll r = -W2; r <= -500; r++) if ((r & (M - 1)) == c && bm[r + W2]) hn++;
            for (ll r = -W1; r <= W1; r++) if ((r & (M - 1)) == c && bm[r + W2]) hw++;
            printf("C %lld %lld %lld %lld %lld %lld %lld %lld %lld\n", c, cnt[c], mn[c], amn[c], mx[c], amx[c], hp, hn, hw);
        }
        printf("SUFFIX %lld %lld\n", sufc, sufv);
        printf("MISS");
        if (cstar >= 0)
            for (ll r = -W1; r <= W1; r++) if ((r & (M - 1)) == cstar && !bm[r + W2]) printf(" %lld", r);
        printf("\nEND\n");
        fflush(stdout);
        free(cnt); free(mn); free(mx); free(amn); free(amx); free(bm);
    }
    return 0;
}

/* number of directed paths with vertex set exactly mask ending at v; optional fixed start vertex */
static ull *path_table(int nv, const ull *out, int start) {
    size_t S = (size_t)1 << nv;
    ull *t = calloc(S * (size_t)nv, sizeof(ull));
    if (!t) { fprintf(stderr, "oom table\n"); exit(2); }
    if (start < 0) { for (int v = 0; v < nv; v++) t[((size_t)1 << v) * nv + v] = 1; }
    else t[((size_t)1 << start) * nv + start] = 1;
    for (size_t mask = 1; mask < S; mask++) {
        if (start >= 0 && !((mask >> start) & 1)) continue;
        for (int v = 0; v < nv; v++) {
            ull x = t[mask * nv + v];
            if (!x) continue;
            ull nb = out[v] & ~(ull)mask;
            while (nb) {
                int w = __builtin_ctzll(nb);
                nb &= nb - 1;
                t[(mask | ((size_t)1 << w)) * nv + w] += x;
            }
        }
    }
    return t;
}

static int graph_mode(void) {
    int nv, root;
    while (scanf("%d %d", &nv, &root) == 2) {
        ull out[64];
        for (int i = 0; i < nv; i++) scanf("%llu", &out[i]);
        size_t S = (size_t)1 << nv, full = S - 1;
        ull *e = path_table(nv, out, -1);
        ull *hp = malloc(S * sizeof(ull));
        hp[0] = 1;
        for (size_t m = 1; m < S; m++) { ull s = 0; for (int v = 0; v < nv; v++) s += e[m * nv + v]; hp[m] = s; }
        ull D = hp[full], E = e[full * nv + root], d = hp[full ^ ((size_t)1 << root)];
        free(e);
        ull *st = path_table(nv, out, root);
        ull s1 = 0, C = 0, U = 0;
        for (int v = 0; v < nv; v++) {
            s1 += st[full * nv + v];
            if ((out[v] >> root) & 1) C += st[full * nv + v];
        }
        for (size_t Y = 1; Y < S; Y++) {
            if (!((Y >> root) & 1)) continue;
            ull hs = 0;
            for (int v = 0; v < nv; v++) hs += st[Y * nv + v];
            U += hs * hp[full ^ Y];
        }
        if (nv < 3) C = 0;
        printf("%llu %llu %llu %llu %llu %llu\n", D, E, U, d, s1, C);
        fflush(stdout);
        free(st); free(hp);
    }
    return 0;
}

/* wl mode: minimal O/Z word length (from the two-vertex seed, d = 1) of a strong source with a U - b D = r, |r| <= R.
 * Explicit-stack depth-first enumeration of all words of length <= DEPTH. */
static int wl_mode(ll a, ll b, int depth, ll R) {
    unsigned char *best = malloc((size_t)(2 * R + 1));
    memset(best, 255, (size_t)(2 * R + 1));
    struct fr { ll U, D; int E, len; };
    struct fr *st = malloc(sizeof(struct fr) * (size_t)(2 * depth + 8));
    int sp = 0;
    st[sp++] = (struct fr){2, 1, 0, 0};
    while (sp) {
        struct fr f = st[--sp];
        if (f.E == 1) {
            ll F = a * f.U - b * f.D;
            if (F >= -R && F <= R && best[F + R] > f.len) best[F + R] = (unsigned char)f.len;
        }
        if (f.len == depth) continue;
        st[sp++] = (struct fr){2 * f.U - 1, f.D - f.E + f.U, 1, f.len + 1};   /* Z */
        st[sp++] = (struct fr){2 * f.U, f.D, 0, f.len + 1};                   /* O */
    }
    for (ll r = -R; r <= R; r++) if (best[r + R] != 255) printf("W %lld %d\n", r, best[r + R]);
    free(best); free(st);
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 3 && !strcmp(argv[1], "census")) return census(atoi(argv[2]));
    if (argc >= 2 && !strcmp(argv[1], "graph")) return graph_mode();
    if (argc >= 6 && !strcmp(argv[1], "wl")) return wl_mode(atoll(argv[2]), atoll(argv[3]), atoi(argv[4]), atoll(argv[5]));
    fprintf(stderr, "usage: audit_census census LOG2N | graph | wl a b DEPTH R\n");
    return 1;
}

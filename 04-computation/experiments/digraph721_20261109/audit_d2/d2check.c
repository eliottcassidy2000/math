/* d2check.c - audit D2 (THM-5056), independent brute-force checker.  Written from scratch for the audit; shares no
 * code with wcheck.c or dom_even.c.
 *
 * Input: lines "tag set n k | u-v u-v ..." (the wcheck format; "set" is ignored here).
 * For each graph:
 *   (1) d(G) by a Gray-code walk over all 2^n subsets S, maintaining cnt[v] = |N[v] cap S| and the number of
 *       undominated vertices (a different algorithm from the prefix-OR DP of wcheck.c / oddbuild.py);
 *   (2) the order: n must equal the bit length of k (= ceil(log2(k+1)));
 *   (3) oddness: an explicit certificate, i.e. a permutation p that is verified to be an automorphism and whose
 *       number of reversed edges (edges a<b with p(a)>p(b)) is counted directly and found odd.  First every
 *       transposition is tried; if none is an odd automorphism, Aut(G) is enumerated completely by backtracking
 *       (node cap) and every element's sign is computed.
 * Mode "full" (argv[1] = full): always enumerate Aut(G) completely and print |Aut| and the number of odd elements.
 * Output: one line per failure (or per graph in full/verbose mode) and a summary line.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 30
static int n;
static uint32_t adj[MAXN];          /* open neighbourhoods */

static int edge(int a, int b) { return (adj[a] >> b) & 1; }

static uint64_t dcount_gray(void) {
    uint32_t cl[MAXN]; int cnt[MAXN]; int und = n; uint64_t c = (n == 0);
    for (int v = 0; v < n; v++) { cl[v] = adj[v] | (1u << v); cnt[v] = 0; }
    uint32_t S = 0;
    uint64_t N = 1ull << n;
    for (uint64_t i = 1; i < N; i++) {
        int b = __builtin_ctzll(i);
        uint32_t r = cl[b];
        if ((S >> b) & 1) {             /* remove b */
            S &= ~(1u << b);
            while (r) { int v = __builtin_ctz(r); r &= r - 1; if (--cnt[v] == 0) und++; }
        } else {                        /* add b */
            S |= 1u << b;
            while (r) { int v = __builtin_ctz(r); r &= r - 1; if (cnt[v]++ == 0) und--; }
        }
        if (und == 0) c++;
    }
    return c;
}

/* reversed-edge count of a permutation p (assumed to be an automorphism) */
static int reversed_edges(const int *p) {
    int rev = 0;
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if (edge(a, b) && p[a] > p[b]) rev++;
    return rev;
}
static int is_automorphism(const int *p) {
    int seen[MAXN] = {0};
    for (int a = 0; a < n; a++) { if (p[a] < 0 || p[a] >= n || seen[p[a]]) return 0; seen[p[a]] = 1; }
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) if (edge(a, b) != edge(p[a], p[b])) return 0;
    return 1;
}

/* full backtracking enumeration of Aut(G) */
static int deg[MAXN], img[MAXN], used[MAXN];
static long long naut, nodd, nodes, cap;
static int found_odd, stop_at_first_odd, capped;
static int cert[MAXN];
static void bt(int v) {
    if (capped) return;
    if (++nodes > cap) { capped = 1; return; }
    if (v == n) {
        naut++;
        if (!is_automorphism(img)) { fprintf(stderr, "internal error: non-automorphism\n"); exit(2); }
        if (reversed_edges(img) & 1) { nodd++; if (!found_odd) { found_odd = 1; memcpy(cert, img, sizeof(int) * n); } }
        return;
    }
    for (int c = 0; c < n; c++) {
        if (used[c] || deg[c] != deg[v]) continue;
        int ok = 1;
        for (int j = 0; j < v; j++) if (edge(v, j) != edge(c, img[j])) { ok = 0; break; }
        if (!ok) continue;
        img[v] = c; used[c] = 1;
        bt(v + 1);
        used[c] = 0; img[v] = -1;
        if (stop_at_first_odd && found_odd) return;
        if (capped) return;
    }
}
static void enumerate_aut(int first_odd_only, long long nodecap) {
    for (int v = 0; v < n; v++) { deg[v] = __builtin_popcount(adj[v]); img[v] = -1; used[v] = 0; }
    naut = nodd = nodes = 0; found_odd = 0; capped = 0; cap = nodecap; stop_at_first_odd = first_odd_only;
    bt(0);
}

/* returns 1 = odd (certificate in cert[], kind in *kind), 0 = even (complete enumeration), -1 = undecided (cap) */
static int odd_certificate(char *kind) {
    int p[MAXN];
    for (int u = 0; u < n; u++) for (int v = u + 1; v < n; v++) {
        for (int a = 0; a < n; a++) p[a] = a;
        p[u] = v; p[v] = u;
        if (!is_automorphism(p)) continue;
        if (reversed_edges(p) & 1) { memcpy(cert, p, sizeof(int) * n); *kind = 'T'; return 1; }
    }
    enumerate_aut(1, 200000000LL);
    if (found_odd) { *kind = 'A'; return 1; }
    if (capped) { *kind = '?'; return -1; }
    *kind = 'E'; return 0;
}

static int bitlen(unsigned long long k) { int L = 0; while (k) { L++; k >>= 1; } return L; }

int main(int argc, char **argv) {
    int full = (argc > 1 && !strcmp(argv[1], "full"));
    int verbose = (argc > 1 && !strcmp(argv[1], "verbose"));
    static char line[1 << 17];
    long total = 0, bad = 0, nT = 0, nA = 0, maxn = 0;
    while (fgets(line, sizeof line, stdin)) {
        char tag[128], set[16]; unsigned long long k; int pos = 0;
        if (sscanf(line, "%127s %15s %d %llu |%n", tag, set, &n, &k, &pos) < 4 || pos == 0) continue;
        if (n < 0 || n > 24) { printf("SKIP %s n=%d (too large)\n", tag, n); continue; }
        memset(adj, 0, sizeof adj);
        char *s = line + pos; int a, b, m;
        while (sscanf(s, " %d-%d%n", &a, &b, &m) == 2) {
            if (a < 0 || b < 0 || a >= n || b >= n || a == b) { printf("BADEDGE %s %d-%d\n", tag, a, b); bad++; }
            adj[a] |= 1u << b; adj[b] |= 1u << a; s += m;
        }
        total++; if (n > maxn) maxn = n;
        uint64_t d = dcount_gray();
        int ordok = (n == bitlen(k));
        if (full) {
            enumerate_aut(0, 2000000000LL);
            printf("%s n=%d k=%llu d=%llu |Aut|=%lld odd_elements=%lld parity=%s%s\n", tag, n, k,
                   (unsigned long long)d, naut, nodd, nodd ? "ODD" : "EVEN", capped ? " (CAPPED)" : "");
            continue;
        }
        char kind; int odd = odd_certificate(&kind);
        if (kind == 'T') nT++; else if (kind == 'A') nA++;
        int ok = (d == k) && ordok && (odd == 1);
        if (!ok) bad++;
        if (!ok || verbose)
            printf("%s %s n=%d k=%llu d=%llu order_ok=%d odd=%d cert=%c\n", ok ? "OK" : "FAIL", tag, n, k,
                   (unsigned long long)d, ordok, odd, kind);
    }
    if (!full)
        printf("d2check: graphs=%ld failures=%ld (odd certificates: transposition %ld, other automorphism %ld; max n %ld)\n",
               total, bad, nT, nA, maxn);
    return bad ? 1 : 0;
}

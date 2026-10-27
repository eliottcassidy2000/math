// tcbatch.c -- read tournaments "n hex(out_0) ... hex(out_{n-1})" (one per line), and for each strong one compute
// H (pull Held-Karp, n <= 20), |Aut|, tc, the Aut-orbits, tc of every orbit-induced subtournament, and test
//   (1) tc not in {7, 21};  (2) orbit product lemma  tc(T) >= prod_O tc(T[O]);  (3) every T[O] regular (VT).
// Prints one summary line per input and a final tally.
#include "../common/tourlib.h"

static uint64_t *F = NULL; static int Fn = 0;
static uint64_t H_of(int n, const uint32_t *out) {
    if (n > Fn) { free(F); F = malloc(sizeof(uint64_t) * ((size_t)1 << n) * n); Fn = n; }
    uint32_t in[32]; make_in(n, out, in); return ham_pull(n, in, F);
}
static int same_orbit(int n, const uint32_t *out, const uint32_t *in, const uint32_t *cand, int v, int x) {
    uint32_t un = (1u << n) - 1, c2[32];
    return assign_prop(n, out, in, cand, un, v, x, c2) && ext_exists(n, out, in, c2, un & ~(1u << v));
}
int main(void) {
    char line[4096]; long long cnt = 0, strongc = 0, hits = 0, lemmafail = 0, vtfail = 0, divfail = 0;
    long long tchist[64] = {0};
    while (fgets(line, sizeof line, stdin)) {
        int n; char *p = line; if (sscanf(p, "%d", &n) != 1) continue;
        while (*p && *p != ' ') p++;
        uint32_t out[32] = {0};
        for (int v = 0; v < n; v++) { while (*p == ' ') p++; out[v] = (uint32_t)strtoul(p, &p, 16); }
        if (!is_tournament(n, out)) { printf("bad input\n"); continue; }
        cnt++;
        if (!strong_bfs(n, out)) continue;
        strongc++;
        uint64_t H = H_of(n, out); long long A = aut_order(n, out);
        if (H % (uint64_t)A) { divfail++; continue; }
        uint64_t tc = H / (uint64_t)A;
        if (tc < 64) tchist[tc]++;
        if (tc == 7 || tc == 21) { hits++; printf("HIT tc=%llu n=%d\n", (unsigned long long)tc, n); }
        // Aut orbits
        uint32_t in[32]; make_in(n, out, in); uint64_t col[32]; refine_colours(n, out, col);
        uint32_t cand[32]; for (int a = 0; a < n; a++) { cand[a] = 0; for (int x = 0; x < n; x++) if (col[x] == col[a]) cand[a] |= 1u << x; }
        int orb[32]; for (int v = 0; v < n; v++) orb[v] = -1;
        int no = 0; double prod = 1.0; char desc[512]; int dl = 0; dl += snprintf(desc + dl, sizeof desc - dl, "orbits:");
        for (int v = 0; v < n; v++) {
            if (orb[v] >= 0) continue;
            orb[v] = no; uint32_t O = 1u << v;
            for (int x = v + 1; x < n; x++) if (orb[x] < 0 && (cand[v] >> x & 1u) && same_orbit(n, out, in, cand, v, x)) { orb[x] = no; O |= 1u << x; }
            // induced subtournament
            int m = pc32(O), map[32], k = 0; for (uint32_t t = O; t; t &= t - 1) map[k++] = ctz32(t);
            uint32_t sub[32] = {0};
            for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) if (i != j && ((out[map[i]] >> map[j]) & 1u)) sub[i] |= 1u << j;
            int reg = 1; for (int i = 0; i < m; i++) if (pc32(sub[i]) * 2 != m - 1) reg = 0;
            if (!reg) vtfail++;
            uint64_t h = H_of(m, sub); long long a = aut_order(m, sub);
            prod *= (double)(h / (uint64_t)a);
            dl += snprintf(desc + dl, sizeof desc - dl, " %d(tc%llu)", m, (unsigned long long)(h / (uint64_t)a));
            no++;
        }
        if ((double)tc < prod - 0.5) { lemmafail++; printf("LEMMA FAIL n=%d tc=%llu prod=%.0f\n", n, (unsigned long long)tc, prod); }
        printf("n=%d H=%llu Aut=%lld tc=%llu prod_orbit_tc=%.0f %s\n", n, (unsigned long long)H, A, (unsigned long long)tc, prod, desc);
    }
    printf("TALLY inputs=%lld strong=%lld hits=%lld lemma-failures=%lld nonregular-orbits=%lld H%%Aut-fail=%lld small tc:", cnt, strongc, hits, lemmafail, vtfail, divfail);
    for (int t = 1; t < 64; t++) if (tchist[t]) printf(" %d:%lld", t, tchist[t]);
    printf("\n");
    return 0;
}

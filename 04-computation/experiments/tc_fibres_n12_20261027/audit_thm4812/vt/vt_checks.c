// vt_checks.c -- independent audit of the THM-4812 (c) vertex-transitive data.
//   ./vt_checks census N          read `gentourng -q -dK -DK N` (regular tournaments) on stdin; report H, |Aut|,
//                                 vertex-transitivity (orbit of every vertex), tc
//   ./vt_checks cayley m1 m2      all Cayley tournaments of Z_m1 x Z_m2; distinct (H,|Aut|,tc) values, min tc, tc <= 21
//   ./vt_checks affine15 A|B      G_A = {(i,x) -> (i+s, x+v_i): v in F_3^5, sum v = 0} on Z_5 x Z_3, or
//                                 G_B = same with Z_3 x Z_5 and v in F_5^3, sum v = 0; every invariant tournament
#include "../common/tourlib.h"

static int orbit_size(int n, const uint32_t *out, int v) {
    uint32_t in[32]; make_in(n, out, in);
    uint64_t col[32]; refine_colours(n, out, col);
    uint32_t cand[32]; for (int a = 0; a < n; a++) { cand[a] = 0; for (int x = 0; x < n; x++) if (col[x] == col[a]) cand[a] |= 1u << x; }
    uint32_t un = (1u << n) - 1; int cnt = 0;
    for (int x = 0; x < n; x++) {
        if (!(cand[v] >> x & 1u)) continue;
        uint32_t c2[32];
        if (assign_prop(n, out, in, cand, un, v, x, c2) && ext_exists(n, out, in, c2, un & ~(1u << v))) cnt++;
    }
    return cnt;
}

static uint64_t *F;
static void report(int n, const uint32_t *out, uint64_t *H, long long *A, int *vt) {
    uint32_t in[32]; make_in(n, out, in);
    *H = ham_pull(n, in, F); *A = aut_order(n, out);
    *vt = 1; for (int v = 0; v < n; v++) if (orbit_size(n, out, v) != n) { *vt = 0; break; }
}

int main(int argc, char **argv) {
    if (argc < 2) return 1;
    F = malloc(sizeof(uint64_t) * ((size_t)1 << 17) * 17);
    if (!strcmp(argv[1], "census")) {
        int n = atoi(argv[2]); char line[1024]; int idx = 0, nvt = 0;
        while (fgets(line, sizeof line, stdin)) {
            if (line[0] != '0' && line[0] != '1') continue;
            uint32_t out[32] = {0}; int k = 0;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (line[k] == '1') out[i] |= 1u << j; else out[j] |= 1u << i; k++; }
            int reg = 1; for (int v = 0; v < n; v++) if (pc32(out[v]) != (n - 1) / 2) reg = 0;
            uint64_t H; long long A; int vt; report(n, out, &H, &A, &vt); idx++;
            if (vt) nvt++;
            printf("class %2d: regular=%d H=%llu |Aut|=%lld vertex-transitive=%d tc=%llu%s\n", idx, reg, (unsigned long long)H, A, vt,
                   (unsigned long long)(H / A), H % A ? " (H not divisible!)" : "");
        }
        printf("classes %d, vertex-transitive %d\n", idx, nvt);
        return 0;
    }
    if (!strcmp(argv[1], "cayley")) {
        int m1 = atoi(argv[2]), m2 = atoi(argv[3]), n = m1 * m2;
        // elements g = (a,b) -> index a*m2 + b; choose for each pair {g,-g} one element
        int rep[32], nr = 0, used[32] = {0}; used[0] = 1;
        for (int g = 1; g < n; g++) { if (used[g]) continue; int a = g / m2, b = g % m2; int h = ((m1 - a) % m1) * m2 + (m2 - b) % m2; used[g] = used[h] = 1; rep[nr++] = g; }
        uint64_t minTc = ~0ull; int small = 0;
        uint64_t seenH[512]; long long seenA[512]; int seenN[512], ns = 0;
        for (int mask = 0; mask < (1 << nr); mask++) {
            int S[32];
            for (int k = 0; k < nr; k++) { int g = rep[k], a = g / m2, b = g % m2; S[k] = (mask >> k & 1) ? ((m1 - a) % m1) * m2 + (m2 - b) % m2 : g; }
            uint32_t out[32] = {0};
            for (int u = 0; u < n; u++) { int ua = u / m2, ub = u % m2; for (int k = 0; k < nr; k++) { int sa = S[k] / m2, sb = S[k] % m2; out[u] |= 1u << (((ua + sa) % m1) * m2 + (ub + sb) % m2); } }
            if (!is_tournament(n, out)) { printf("not a tournament\n"); return 1; }
            uint32_t in[32]; make_in(n, out, in);
            uint64_t H = ham_pull(n, in, F); long long A = aut_order(n, out);
            if (H % A) { printf("H not divisible by |Aut|\n"); return 1; }
            uint64_t tc = H / A; if (tc < minTc) minTc = tc; if (tc <= 21) small++;
            int f = -1; for (int i = 0; i < ns; i++) if (seenH[i] == H && seenA[i] == A) f = i;
            if (f < 0) { seenH[ns] = H; seenA[ns] = A; seenN[ns] = 1; ns++; } else seenN[f]++;
        }
        printf("Cayley Z_%d x Z_%d: %d connection sets, min tc %llu, sets with tc <= 21: %d; distinct (H,|Aut|,tc):", m1, m2, 1 << nr, (unsigned long long)minTc, small);
        for (int i = 0; i < ns; i++) printf(" (%llu,%lld,%llu)x%d", (unsigned long long)seenH[i], seenA[i], (unsigned long long)(seenH[i] / seenA[i]), seenN[i]);
        printf("\n");
        return 0;
    }
    if (!strcmp(argv[1], "affine15")) {
        int nb, bs; if (argv[2][0] == 'A') { nb = 5; bs = 3; } else { nb = 3; bs = 5; }
        int n = nb * bs;   // point (i,x) -> i*bs + x
        // group elements: s in Z_nb, v in F_bs^nb with sum 0
        static int G[1300][15]; int ng = 0;
        int v[8];
        long long total = 1; for (int i = 0; i < nb; i++) total *= bs;
        for (int s = 0; s < nb; s++) for (long long code = 0; code < total; code++) {
            long long c = code; int sum = 0; for (int i = 0; i < nb; i++) { v[i] = (int)(c % bs); c /= bs; sum += v[i]; }
            if (sum % bs) continue;
            for (int i = 0; i < nb; i++) for (int x = 0; x < bs; x++) G[ng][i * bs + x] = ((i + s) % nb) * bs + (x + v[i]) % bs;
            ng++;
        }
        // element orders: is there an element of order n (a regular cyclic subgroup)?
        int maxord = 0;
        for (int g = 0; g < ng; g++) { int p[15], q[15]; for (int a = 0; a < n; a++) p[a] = a; int o = 0; do { for (int a = 0; a < n; a++) q[a] = G[g][p[a]]; memcpy(p, q, sizeof p); o++; int id = 1; for (int a = 0; a < n; a++) if (p[a] != a) id = 0; if (id) break; } while (o < 1000); if (o > maxord) maxord = o; }
        // transitivity
        int reach[15] = {0}; for (int g = 0; g < ng; g++) reach[G[g][0]] = 1; int tr = 1; for (int a = 0; a < n; a++) if (!reach[a]) tr = 0;
        printf("G_%c: |G| = %d, transitive=%d, max element order %d\n", argv[2][0], ng, tr, maxord);
        // orbitals
        static int orb[15][15]; for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) orb[a][b] = -1;
        int no = 0;
        for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if (a != b && orb[a][b] < 0) { for (int g = 0; g < ng; g++) orb[G[g][a]][G[g][b]] = no; no++; }
        int pairA[200], pairB[200], np = 0, done[400] = {0};
        for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if (a != b && !done[orb[a][b]]) { int t = orb[a][b], ti = orb[b][a]; if (t == ti) { printf("self-paired orbital\n"); return 1; } done[t] = done[ti] = 1; pairA[np] = t; pairB[np] = ti; np++; }
        printf("  orbital pairs %d -> %d invariant tournaments\n", np, 1 << np);
        for (int m = 0; m < (1 << np); m++) {
            uint32_t out[32] = {0};
            for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if (a != b) {
                for (int k = 0; k < np; k++) { int want = (m >> k & 1) ? pairB[k] : pairA[k]; if (orb[a][b] == want) out[a] |= 1u << b; }
            }
            if (!is_tournament(n, out)) { printf("not tournament\n"); return 1; }
            // blocks as modules?
            int modok = 1;
            for (int i = 0; i < nb; i++) for (int j = 0; j < nb; j++) if (i != j) { int d = (out[i * bs] >> (j * bs)) & 1u; for (int x = 0; x < bs; x++) for (int y = 0; y < bs; y++) if ((int)((out[i * bs + x] >> (j * bs + y)) & 1u) != d) modok = 0; }
            uint64_t H; long long A; int vt; report(n, out, &H, &A, &vt);
            printf("  tournament %d: blocks are modules=%d H=%llu |Aut|=%lld VT=%d tc=%llu\n", m, modok, (unsigned long long)H, A, vt, (unsigned long long)(H / A));
        }
        return 0;
    }
    return 1;
}

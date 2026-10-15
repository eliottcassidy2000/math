/* Audit E part 2 (2026-10-15): brute-force validation of THM-4655 (7) / HYP-9279 (Z1) on random cut-vertex tournaments,
   with NO use of the transfer chain.  T = z + (D_1 => D_2 => ... => D_k), each D_j a random strong tournament (rejection
   sampling), all arcs between blocks forward, z attached by a random pattern conditioned on T strong (z -> some vertex of
   D_1, some vertex of D_k -> z).  Exact counts by forward/backward subset DP: N(v) = #Hamiltonian paths with v before z,
     N(v) = sum over S (z not in S, v in S) of [sum_{x in S, x -> z} f[S][x]] * g[V \ S][z],
   f[S][x] = paths covering S ending at x, g[R][y] = paths covering R starting at y.  Checks 1/3 <= N/H <= 2/3 for every v
   in an interior block (2 <= i <= k-1), separately recording interior blocks of size <= 7 (covered by the theorem) and
   the minimum balance seen.  Also independently checks that T - z has exactly the blocks as strong components.
   Usage: ./z1_bruteforce_random SEED TRIALS MODE   (MODE 0: random sizes; MODE 1: large end blocks, small interior)   */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int n; static unsigned out[20];
static u64 *f, *g;
static unsigned long long rs;
static unsigned rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (unsigned)(rs >> 11); }
static unsigned reach(unsigned S, int s) { unsigned r = 1u << s, fr = r; while (fr) { int v = __builtin_ctz(fr); fr &= fr - 1; unsigned nw = out[v] & S & ~r; r |= nw; fr |= nw; } return r; }
static int strong_on(unsigned S) { for (int v = 0; v < n; v++) if (S >> v & 1) if (reach(S, v) != S) return 0; return 1; }
int main(int argc, char **argv) {
    rs = 88172645463325252ULL ^ (u64)atoll(argv[1]) * 0x9E3779B97F4A7C15ULL; if (!rs) rs = 1;
    long trials = atol(argv[2]); int mode = atoi(argv[3]);
    f = malloc(sizeof(u64) * (1 << 17) * 17); g = malloc(sizeof(u64) * (1 << 17) * 17);
    long tested = 0, vchecked = 0, vsmall = 0, fails = 0, fails_small = 0, compfail = 0;
    u64 bn = 1, bd = 1; char bdesc[256] = "";
    int strongsizes[] = {1, 3, 4, 5, 6, 7};
    for (long t = 0; t < trials; t++) {
        int k, sz[12];
        for (;;) {
            if (mode == 0) { k = 3 + rnd() % 3; int tot = 1; for (int j = 0; j < k; j++) { sz[j] = strongsizes[rnd() % 6]; tot += sz[j]; } if (tot <= 16) break; }
            else { k = 3 + rnd() % 2; int tot = 1; sz[0] = 4 + rnd() % 4; sz[k - 1] = 4 + rnd() % 4;
                   for (int j = 1; j < k - 1; j++) sz[j] = strongsizes[rnd() % 3]; for (int j = 0; j < k; j++) tot += sz[j]; if (tot <= 16) break; }
        }
        n = 1; for (int j = 0; j < k; j++) n += sz[j];
        int blk[20]; int start[12]; int pos = 1; for (int j = 0; j < k; j++) { start[j] = pos; for (int q = 0; q < sz[j]; q++) blk[pos + q] = j; pos += sz[j]; }
        memset(out, 0, sizeof out);
        /* blocks: random strong tournaments */
        for (int j = 0; j < k; j++) {
            unsigned M = 0; for (int q = 0; q < sz[j]; q++) M |= 1u << (start[j] + q);
            for (;;) {
                for (int a = start[j]; a < start[j] + sz[j]; a++) out[a] &= ~M;
                for (int a = start[j]; a < start[j] + sz[j]; a++) for (int b = a + 1; b < start[j] + sz[j]; b++) { if (rnd() & 1) out[a] |= 1u << b; else out[b] |= 1u << a; }
                if (sz[j] == 1 || strong_on(M)) break;
            }
        }
        for (int a = 1; a < n; a++) for (int b = 1; b < n; b++) if (blk[a] < blk[b]) out[a] |= 1u << b;
        /* z = vertex 0 */
        for (;;) {
            out[0] = 0; for (int a = 1; a < n; a++) out[a] &= ~1u;
            for (int a = 1; a < n; a++) { if (rnd() & 1) out[0] |= 1u << a; else out[a] |= 1u; }
            int ok1 = 0, okk = 0;
            for (int a = start[0]; a < start[0] + sz[0]; a++) if (out[0] >> a & 1) ok1 = 1;
            for (int a = start[k - 1]; a < start[k - 1] + sz[k - 1]; a++) if (out[a] & 1u) okk = 1;
            if (ok1 && okk) break;
        }
        unsigned full = (1u << n) - 1;
        if (!strong_on(full)) { compfail++; continue; }
        /* components of T - z must be exactly the blocks */
        for (int a = 1; a < n; a++) { unsigned R = full ^ 1u; unsigned ra = reach(R, a); unsigned C = 0;
            for (int b = 1; b < n; b++) if ((ra >> b & 1) && (reach(R, b) >> a & 1)) C |= 1u << b;
            unsigned Bm = 0; for (int q = 0; q < sz[blk[a]]; q++) Bm |= 1u << (start[blk[a]] + q);
            if (C != Bm) compfail++; }
        /* subset DP */
        int N = 1 << n;
        memset(f, 0, sizeof(u64) * N * n); memset(g, 0, sizeof(u64) * N * n);
        for (int v = 0; v < n; v++) { f[(1 << v) * n + v] = 1; g[(1 << v) * n + v] = 1; }
        for (int S = 1; S < N; S++) {
            if (!(S & (S - 1))) continue;
            for (int v = 0; v < n; v++) if (S >> v & 1) {
                int R = S ^ (1 << v); u64 a = 0, b = 0;
                for (int u = 0; u < n; u++) if ((R >> u & 1) && (out[u] >> v & 1)) a += f[R * n + u];
                for (int u = 0; u < n; u++) if ((R >> u & 1) && (out[v] >> u & 1)) b += g[R * n + u];
                f[S * n + v] = a; g[S * n + v] = b;
            }
        }
        u64 H = 0; for (int v = 0; v < n; v++) H += f[(N - 1) * n + v];
        u64 Nv[20] = {0};
        for (int S = 2; S < N; S += 2) {      /* S without z (bit 0) */
            u64 w = 0; for (int x = 1; x < n; x++) if ((S >> x & 1) && (out[x] & 1u)) w += f[S * n + x];
            if (!w) continue;
            u64 tail = g[((N - 1) ^ S) * n + 0];
            if (!tail) continue;
            for (int v = 1; v < n; v++) if (S >> v & 1) Nv[v] += w * tail;
        }
        tested++;
        for (int v = 1; v < n; v++) {
            int j = blk[v]; if (j == 0 || j == k - 1) continue;
            vchecked++; int small = sz[j] <= 7; vsmall += small;
            u64 c = Nv[v], m = c < H - c ? c : H - c;
            if (3 * c < H || 3 * c > 2 * H) { fails++; if (small) fails_small++; }
            if ((unsigned __int128)m * bd < (unsigned __int128)bn * H) { bn = m; bd = H;
                int p = 0; p += snprintf(bdesc + p, 256 - p, "n=%d k=%d sizes", n, k); for (int q = 0; q < k; q++) p += snprintf(bdesc + p, 256 - p, " %d", sz[q]);
                p += snprintf(bdesc + p, 256 - p, " v in block %d", j + 1); }
        }
    }
    u64 a = bn, b = bd; while (b) { u64 t = a % b; a = b; b = t; }
    printf("mode %d: tournaments tested %ld (construction failures %ld); interior vertices checked %ld (in blocks of size <= 7: %ld); "
           "violations of 1/3 <= P(v before z) <= 2/3: %ld (in blocks <= 7: %ld); minimum balance seen %llu/%llu = %.6f (%s)\n",
           mode, tested, compfail, vchecked, vsmall, fails, fails_small, bn / a, bd / a, (double)bn / bd, bdesc);
    return 0;
}

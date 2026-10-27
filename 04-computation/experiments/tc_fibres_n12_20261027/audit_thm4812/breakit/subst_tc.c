// subst_tc.c -- structured attack on THM-4812: exact tc of every substitution S = Q[B_1..B_s] with Q a strong
// tournament class on s vertices (read from `gentourng -c -q s` on stdin) and blocks B_i in {K1, C3, R5, P7},
// restricted to orders nmin..nmax and at most `maxbig` blocks of type R5/P7.  H by the subset-run DP (C3 compressed,
// R5/P7 by subsets), |Aut| by tourlib's counter; cross-check H by Held-Karp when n <= 18 (every case).
// Usage: gentourng -c -q s | ./subst_tc s nmin nmax maxbig
#include "../common/tourlib.h"

static const int SZ[4] = {1, 3, 5, 7};
static int s, n, bkind[16], bsize[16], bverts[16][8];
static uint64_t hloc[16][128];
static int rad[16]; static long long mulc[16];
typedef struct { int newc; uint64_t w; } Tr;
static Tr trl[16][128][128]; static int trn[16][128];
static uint64_t *dp = NULL; static size_t dpcap = 0;

static void block_arcs(int kind, int m, uint32_t *lo) {   // local out-masks of the standard copy
    for (int x = 0; x < m; x++) lo[x] = 0;
    if (kind == 1) for (int x = 0; x < 3; x++) lo[x] |= 1u << ((x + 1) % 3);
    if (kind == 2) for (int x = 0; x < 5; x++) { lo[x] |= 1u << ((x + 1) % 5); lo[x] |= 1u << ((x + 2) % 5); }
    if (kind == 3) for (int x = 0; x < 7; x++) { lo[x] |= 1u << ((x + 1) % 7); lo[x] |= 1u << ((x + 2) % 7); lo[x] |= 1u << ((x + 4) % 7); }
}

static uint64_t dp_H(const int q[16][16]) {
    long long S = 1;
    for (int k = 0; k < s; k++) { rad[k] = bkind[k] == 0 ? 2 : bkind[k] == 1 ? 4 : (1 << bsize[k]); mulc[k] = S; S *= rad[k]; }
    size_t need = (size_t)S * s; if (need > dpcap) { free(dp); dp = malloc(sizeof(uint64_t) * need); dpcap = need; }
    memset(dp, 0, sizeof(uint64_t) * need);
    for (int k = 0; k < s; k++) for (int t = 0; t < trn[k][0]; t++) dp[(size_t)(trl[k][0][t].newc * mulc[k]) * s + k] += trl[k][0][t].w;
    long long full = S - 1; int comp[16];
    for (long long code = 1; code < full; code++) {
        uint64_t *row = dp + (size_t)code * s; int any = 0; for (int l = 0; l < s; l++) if (row[l]) { any = 1; break; } if (!any) continue;
        long long c = code; for (int k = 0; k < s; k++) { comp[k] = (int)(c % rad[k]); c /= rad[k]; }
        for (int l = 0; l < s; l++) { uint64_t v = row[l]; if (!v) continue;
            for (int k = 0; k < s; k++) { if (k == l || !q[l][k]) continue; int cc = comp[k]; if (cc == rad[k] - 1) continue;
                long long base = code - (long long)cc * mulc[k];
                for (int t = 0; t < trn[k][cc]; t++) dp[(size_t)(base + (long long)trl[k][cc][t].newc * mulc[k]) * s + k] += v * trl[k][cc][t].w; } }
    }
    uint64_t H = 0; for (int l = 0; l < s; l++) H += dp[(size_t)full * s + l]; return H;
}

static void setup_blocks(void) {
    for (int k = 0; k < s; k++) {
        int m = bsize[k]; uint32_t lo[8]; block_arcs(bkind[k], m, lo);
        for (int R = 1; R < (1 << m); R++) {
            uint64_t f[128][7]; memset(f, 0, sizeof f); uint64_t tot = 0;
            for (int a = 0; a < m; a++) if (R >> a & 1) f[1 << a][a] = 1;
            for (int msk = 1; msk < (1 << m); msk++) { if ((msk & R) != msk) continue;
                for (int a = 0; a < m; a++) if (f[msk][a]) for (int b = 0; b < m; b++) if ((R >> b & 1) && !(msk >> b & 1) && (lo[a] >> b & 1u)) f[msk | 1 << b][b] += f[msk][a]; }
            for (int a = 0; a < m; a++) tot += f[R][a]; hloc[k][R] = tot;
        }
        for (int c = 0; c < 128; c++) trn[k][c] = 0;
        if (bkind[k] == 0) { trl[k][0][0].newc = 1; trl[k][0][0].w = 1; trn[k][0] = 1; }
        else if (bkind[k] == 1) { static const int bn[4][4] = {{1,0,0,0},{1,1,0,0},{1,2,1,0},{1,3,3,1}}; uint64_t hs[4] = {0, 1, 1, 3};
            for (int u = 0; u < 3; u++) { int rem = 3 - u; for (int j = 1; j <= rem; j++) { trl[k][u][trn[k][u]].newc = u + j; trl[k][u][trn[k][u]].w = (uint64_t)bn[rem][j] * hs[j]; trn[k][u]++; } } }
        else { int full = (1 << m) - 1; for (int U = 0; U < full; U++) { int cp = full & ~U; for (int R = cp; R; R = (R - 1) & cp) { trl[k][U][trn[k][U]].newc = U | R; trl[k][U][trn[k][U]].w = hloc[k][R]; trn[k][U]++; } } }
    }
}

int main(int argc, char **argv) {
    s = atoi(argv[1]); int nmin = atoi(argv[2]), nmax = atoi(argv[3]), maxbig = atoi(argv[4]);
    char line[256]; long long combos = 0, hits = 0, hkchecks = 0, hkfail = 0, divfail = 0;
    long long mintc[64]; for (int i = 0; i < 64; i++) mintc[i] = -1;
    uint64_t *F = malloc(sizeof(uint64_t) * ((size_t)1 << 18) * 18);
    while (fgets(line, sizeof line, stdin)) {
        if (s > 1 && line[0] != '0' && line[0] != '1') continue;
        int q[16][16] = {{0}}; int k = 0;
        for (int i = 0; i < s; i++) for (int j = i + 1; j < s; j++) { if (line[k] == '1') q[i][j] = 1; else q[j][i] = 1; k++; }
        long long nas = 1; for (int i = 0; i < s; i++) nas *= 4;
        for (long long a = 0; a < nas; a++) {
            long long c = a; int big = 0; n = 0;
            for (int i = 0; i < s; i++) { bkind[i] = (int)(c % 4); c /= 4; bsize[i] = SZ[bkind[i]]; n += bsize[i]; if (bkind[i] >= 2) big++; }
            if (n < nmin || n > nmax || big > maxbig) continue;
            if (s == 1 && bkind[0] == 0) continue;
            int pos = 0; for (int i = 0; i < s; i++) { for (int x = 0; x < bsize[i]; x++) bverts[i][x] = pos + x; pos += bsize[i]; }
            uint32_t out[32] = {0};
            for (int i = 0; i < s; i++) { uint32_t lo[8]; block_arcs(bkind[i], bsize[i], lo);
                for (int x = 0; x < bsize[i]; x++) { for (uint32_t t = lo[x]; t; t &= t - 1) out[bverts[i][x]] |= 1u << bverts[i][ctz32(t)];
                    for (int j = 0; j < s; j++) if (j != i && q[i][j]) for (int y = 0; y < bsize[j]; y++) out[bverts[i][x]] |= 1u << bverts[j][y]; } }
            if (!is_tournament(n, out)) { printf("not tournament\n"); return 1; }
            if (!strong_bfs(n, out)) continue;
            setup_blocks();
            uint64_t H = dp_H(q);
            if (n <= 18) { uint32_t in[32]; make_in(n, out, in); uint64_t h2 = ham_pull(n, in, F); hkchecks++; if (h2 != H) { hkfail++; printf("HK MISMATCH\n"); } }
            long long A = aut_order(n, out);
            combos++;
            if (H % (uint64_t)A) { divfail++; continue; }
            long long tc = (long long)(H / (uint64_t)A);
            if (tc == 7 || tc == 21) { hits++; printf("HIT tc=%lld n=%d\n", tc, n); }
            if (mintc[n] < 0 || tc < mintc[n]) mintc[n] = tc;
        }
    }
    printf("s=%d orders %d..%d maxbig=%d: substitutions %lld, hits %lld, HK checks %lld (fail %lld), H%%Aut fail %lld; min tc by order:", s, nmin, nmax, maxbig, combos, hits, hkchecks, hkfail, divfail);
    for (int i = nmin; i <= nmax; i++) if (mintc[i] >= 0) printf(" %d:%lld", i, mintc[i]);
    printf("\n");
    return 0;
}

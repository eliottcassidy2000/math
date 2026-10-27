// window_enum.c -- independent audit enumerator for THM-4812 (d) window jobs written by window_groups.py.
// For every G-invariant tournament of a job (plain binary enumeration over orbital pairs):
//   strong by Landau's criterion; H by
//     (i)  pull Held-Karp when n <= 16 (every job, module or not),
//     (ii) for module jobs (every orbital a full rectangle, i.e. a substitution Q[B_1..B_s]) a subset-run DP:
//          a Hamiltonian path is a sequence of maximal runs inside blocks; state = (used part of every block, block
//          of the last run); a run in block k covering R adds weight h_k(R) = #Hamiltonian paths of T[R]; C3 blocks
//          are compressed to their used count (h of a C3-subset depends only on its size; asserted), K1 to a flag,
//          R5/P7 blocks keep the used subset.  (Different from the repository's run-count x ordered-path-cover DP.)
//     (iii) "spot" mode: on sampled strong tournaments of module jobs with n >= 17 compare (ii) with Held-Karp
//          (n <= 20) or with Karp's inclusion-exclusion formula H = sum_S (-1)^(n-|S|) walks_(n-1)(S) (n = 21, 22).
//   |Aut| by the orbit-stabilizer counter of tourlib.h for EVERY strong tournament (no filter); checks |G| | |Aut|,
//   |Aut| | H; reports tc hits in {7, 21}, min tc, and replays the repository's candidate counter.
// Build: cc -O3 -march=native -o window_enum window_enum.c
// Run:   ./window_enum job [spot NS]
#include "../common/tourlib.h"
#include <time.h>

static int n, s, K; static long long Gsize; static int jobindex; static char kindstr[16];
static int bkind[32];      // 0 K1, 1 C3, 2 R5, 3 P7
static int bsize[32], bverts[32][8];
static int fixu[1024], fixv[1024], nfix = 0;
static int *orbu[1024], *orbv[1024], olen[1024];

static uint64_t rng = 0x243F6A8885A308D3ULL;
static uint64_t xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

// ---------------- subset-run DP ----------------
static int rad[32]; static long long mulc[32]; static long long Scodes;
static uint64_t hloc[32][128];          // h_k(R) for local subsets R of block k
typedef struct { int newc; uint64_t w; } Tr;
static Tr *trl[32][128]; static int trn[32][128];   // transitions per block per component value
static uint64_t *dpbuf = NULL; static size_t dpcap = 0;

static void setup_dp(const uint32_t *out) {
    // block-local Hamiltonian path counts (independent tiny Held-Karp per block)
    for (int k = 0; k < s; k++) {
        int m = bsize[k];
        for (int R = 1; R < (1 << m); R++) {
            // count Hamiltonian paths of T[R] by brute recursion over local vertices
            uint64_t f[128][7]; memset(f, 0, sizeof f);
            for (int a = 0; a < m; a++) if (R >> a & 1) f[1 << a][a] = 1;
            uint64_t tot = 0;
            for (int msk = 1; msk < (1 << m); msk++) {
                if ((msk & R) != msk) continue;
                for (int a = 0; a < m; a++) if (f[msk][a]) {
                    for (int b = 0; b < m; b++) if ((R >> b & 1) && !(msk >> b & 1) && ((out[bverts[k][a]] >> bverts[k][b]) & 1u)) f[msk | 1 << b][b] += f[msk][a];
                }
            }
            for (int a = 0; a < m; a++) tot += f[R][a];
            hloc[k][R] = tot;
        }
    }
    Scodes = 1;
    for (int k = 0; k < s; k++) {
        if (bkind[k] == 0) rad[k] = 2; else if (bkind[k] == 1) rad[k] = 4; else rad[k] = 1 << bsize[k];
        mulc[k] = Scodes; Scodes *= rad[k];
    }
    for (int k = 0; k < s; k++) {
        for (int c = 0; c < rad[k]; c++) { trn[k][c] = 0; free(trl[k][c]); trl[k][c] = NULL; }
        if (bkind[k] == 0) {
            trl[k][0] = malloc(sizeof(Tr)); trl[k][0][0].newc = 1; trl[k][0][0].w = 1; trn[k][0] = 1;
        } else if (bkind[k] == 1) {
            // assert C3 subset counts depend only on size
            uint64_t hs[4] = {0, 0, 0, 0}; int seen[4] = {0};
            for (int R = 1; R < 8; R++) { int z = pc32(R); if (!seen[z]) { hs[z] = hloc[k][R]; seen[z] = 1; } else if (hs[z] != hloc[k][R]) { fprintf(stderr, "C3 compression invalid\n"); exit(3); } }
            static const int binom[4][4] = {{1,0,0,0},{1,1,0,0},{1,2,1,0},{1,3,3,1}};
            for (int u = 0; u < 3; u++) {
                int rem = 3 - u; trl[k][u] = malloc(sizeof(Tr) * 3);
                for (int j = 1; j <= rem; j++) { trl[k][u][trn[k][u]].newc = u + j; trl[k][u][trn[k][u]].w = (uint64_t)binom[rem][j] * hs[j]; trn[k][u]++; }
            }
        } else {
            int m = bsize[k], full = (1 << m) - 1;
            for (int U = 0; U < full; U++) {
                int comp = full & ~U; trl[k][U] = malloc(sizeof(Tr) * (1 << m));
                for (int R = comp; R; R = (R - 1) & comp) { trl[k][U][trn[k][U]].newc = U | R; trl[k][U][trn[k][U]].w = hloc[k][R]; trn[k][U]++; }
            }
        }
    }
    size_t need = (size_t)Scodes * s;
    if (need > dpcap) { free(dpbuf); dpbuf = malloc(sizeof(uint64_t) * need); dpcap = need; }
}

static uint64_t subset_run_dp(const uint32_t *out) {
    int q[32][32];
    for (int i = 0; i < s; i++) for (int j = 0; j < s; j++) q[i][j] = (i != j) && ((out[bverts[i][0]] >> bverts[j][0]) & 1u);
    memset(dpbuf, 0, sizeof(uint64_t) * (size_t)Scodes * s);
    for (int k = 0; k < s; k++) for (int t = 0; t < trn[k][0]; t++) dpbuf[(size_t)(trl[k][0][t].newc * mulc[k]) * s + k] += trl[k][0][t].w;
    int comp[32];
    long long fullcode = 0; for (int k = 0; k < s; k++) fullcode += (long long)(rad[k] - 1) * mulc[k];
    for (long long code = 1; code < fullcode; code++) {
        uint64_t *row = dpbuf + (size_t)code * s;
        int any = 0; for (int l = 0; l < s; l++) if (row[l]) { any = 1; break; }
        if (!any) continue;
        long long c = code; for (int k = 0; k < s; k++) { comp[k] = (int)(c % rad[k]); c /= rad[k]; }
        for (int l = 0; l < s; l++) {
            uint64_t v = row[l]; if (!v) continue;
            for (int k = 0; k < s; k++) {
                if (k == l || !q[l][k]) continue;
                int cc = comp[k]; if (cc == rad[k] - 1) continue;
                long long base = code - (long long)cc * mulc[k];
                for (int t = 0; t < trn[k][cc]; t++) dpbuf[(size_t)(base + (long long)trl[k][cc][t].newc * mulc[k]) * s + k] += v * trl[k][cc][t].w;
            }
        }
    }
    uint64_t H = 0; for (int l = 0; l < s; l++) H += dpbuf[(size_t)fullcode * s + l];
    return H;
}

// Karp inclusion-exclusion, arithmetic mod 2^64 (exact because the true H < 2^64)
static uint64_t ham_karp(const uint32_t *out) {
    uint32_t in[32]; make_in(n, out, in);
    uint64_t total = 0; uint32_t full = (1u << n) - 1;
    for (uint32_t S = 1; S <= full; S++) {
        uint64_t w[32], w2[32]; int verts[32], m = 0;
        for (uint32_t t = S; t; t &= t - 1) verts[m++] = ctz32(t);
        for (int i = 0; i < m; i++) w[verts[i]] = 1;
        for (int step = 1; step < n; step++) {
            for (int i = 0; i < m; i++) { int v = verts[i]; uint64_t acc = 0; for (uint32_t us = in[v] & S; us; us &= us - 1) acc += w[ctz32(us)]; w2[v] = acc; }
            for (int i = 0; i < m; i++) w[verts[i]] = w2[verts[i]];
        }
        uint64_t W = 0; for (int i = 0; i < m; i++) W += w[verts[i]];
        if ((n - m) & 1) total -= W; else total += W;
        if (S == full) break;
    }
    return total;
}

int main(int argc, char **argv) {
    FILE *fp = fopen(argv[1], "r"); if (!fp) { perror("job"); return 1; }
    int spot = 0, NS = 0; if (argc > 3 && !strcmp(argv[2], "spot")) { spot = 1; NS = atoi(argv[3]); }
    char buf[1 << 16];
    if (!fgets(buf, sizeof buf, fp)) return 1;
    sscanf(buf, "%d %d %d %lld %d %15s", &n, &s, &K, &Gsize, &jobindex, kindstr);
    int module = !strcmp(kindstr, "module");
    if (!fgets(buf, sizeof buf, fp)) return 1;
    { char *p = strtok(buf, " \n"); int k = 0;
      while (p) { char kind[4] = {p[0], p[1], 0}; bkind[k] = !strcmp(kind, "K1") ? 0 : !strcmp(kind, "C3") ? 1 : !strcmp(kind, "R5") ? 2 : 3;
                  char *q = p + 3; bsize[k] = 0; char *tok = q; while (*tok) { bverts[k][bsize[k]++] = atoi(tok); while (*tok && *tok != '.') tok++; if (*tok == '.') tok++; }
                  k++; p = strtok(NULL, " \n"); }
      if (k != s) { fprintf(stderr, "orbit parse\n"); return 1; } }
    if (!fgets(buf, sizeof buf, fp)) return 1;
    { char *p = strtok(buf, " \n"); while (p) { sscanf(p, "%d,%d", &fixu[nfix], &fixv[nfix]); nfix++; p = strtok(NULL, " \n"); } }
    for (int k = 0; k < K; k++) {
        if (!fgets(buf, sizeof buf, fp)) { fprintf(stderr, "orbital parse\n"); return 1; }
        int cnt = 0; for (char *c = buf; *c; c++) if (*c == ',') cnt++;
        orbu[k] = malloc(sizeof(int) * cnt); orbv[k] = malloc(sizeof(int) * cnt); olen[k] = 0;
        char *p = strtok(buf, " \n"); while (p) { sscanf(p, "%d,%d", &orbu[k][olen[k]], &orbv[k][olen[k]]); olen[k]++; p = strtok(NULL, " \n"); }
    }
    fclose(fp);
    uint64_t *f = NULL; if (n <= 16) f = malloc(sizeof(uint64_t) * ((size_t)1 << n) * n);
    long long strongc = 0, cand = 0, hits = 0, divfail = 0, gdivfail = 0, cmpn = 0, cmpfail = 0, spotn = 0, spotfail = 0, nontour = 0;
    uint64_t minH = ~0ull; long long minTc = -1;
    long long autmax = 0;
    int dp_ready = 0;
    clock_t t0 = clock();
    // spot sampling: choose NS random codes among strong ones (reservoir-free: sample codes with probability ~ NS/2^K)
    for (long long code = 0; code < (1ll << K); code++) {
        uint32_t out[32]; memset(out, 0, sizeof out);
        for (int i = 0; i < nfix; i++) out[fixu[i]] |= 1u << fixv[i];
        for (int k = 0; k < K; k++) {
            int rev = (code >> k) & 1;
            for (int j = 0; j < olen[k]; j++) { if (rev) out[orbv[k][j]] |= 1u << orbu[k][j]; else out[orbu[k][j]] |= 1u << orbv[k][j]; }
        }
        if (code < 64 && !is_tournament(n, out)) { nontour++; continue; }
        if (!strong_landau(n, out)) continue;
        strongc++;
        uint64_t H = 0, Hhk = 0, Hdp = 0;
        if (module) { if (!dp_ready) { setup_dp(out); dp_ready = 1; } Hdp = subset_run_dp(out); }
        if (n <= 16) { uint32_t in[32]; make_in(n, out, in); Hhk = ham_pull(n, in, f); }
        if (module && n <= 16) { cmpn++; if (Hdp != Hhk) cmpfail++; }
        H = module ? Hdp : Hhk;
        if (!module && n > 16) { fprintf(stderr, "non-module job beyond 16\n"); return 1; }
        if (spot && module && n > 16 && (xr() % (uint64_t)(((1ll << K) / (NS > 0 ? NS : 1)) + 1)) == 0 && spotn < NS) {
            uint64_t Hx;
            if (n <= 20) { uint64_t *g = malloc(sizeof(uint64_t) * ((size_t)1 << n) * n); uint32_t in[32]; make_in(n, out, in); Hx = ham_pull(n, in, g); free(g); }
            else Hx = ham_karp(out);
            spotn++; if (Hx != H) { spotfail++; printf("SPOT MISMATCH code=%lld dp=%llu other=%llu\n", code, (unsigned long long)H, (unsigned long long)Hx); }
        }
        if (H < minH) minH = H;
        long long A = aut_order(n, out);
        if (A > autmax) autmax = A;
        if (A % Gsize) gdivfail++;
        if (H % (uint64_t)A) { divfail++; continue; }
        long long tc = (long long)(H / (uint64_t)A);
        if (minTc < 0 || tc < minTc) minTc = tc;
        if (tc == 7 || tc == 21) { hits++; printf("HIT tc=%lld H=%llu Aut=%lld code=%lld\n", tc, (unsigned long long)H, A, code); }
        for (int t = 7; t <= 21; t += 14) { if (H % t) continue; uint64_t a = H / t; if (a % (uint64_t)Gsize) continue; cand++; }
    }
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("%s: n=%d orbits=%d K=%d |G|=%lld index=%d %s strong=%lld minH=%llu repo-style-autcalls=%lld hits=%lld minTc=%lld max|Aut|=%lld "
           "dp-vs-HK=%lld/%lld-fail Aut%%G-fail=%lld H%%Aut-fail=%lld spot=%lld/%lld-fail nontour=%lld %.1fs\n",
           argv[1], n, s, K, Gsize, jobindex, module ? "module" : "nonmodule", strongc, (unsigned long long)minH, cand, hits, minTc, autmax,
           cmpn, cmpfail, gdivfail, divfail, spotn, spotfail, nontour, secs);
    return 0;
}

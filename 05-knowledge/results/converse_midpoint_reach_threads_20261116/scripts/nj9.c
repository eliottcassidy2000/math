/* nj9.c <L> <targets.txt>: for every tournament class of orders 1..9 (orders >= 3 from gentourng -q on stdin, orders 1, 2
   built in) and every root t: h = H(T), e_t = paths ending at t, B_t = root response.  Marks the near-join values
   h + 2^k B_t (k >= 1) <= L in bitmap nj9.bits, and the dead-end appends h + B_t (k = 0, only when e_t = 0) in
   njext9.bits (which also contains every k >= 1 value).  Writes pairs7.txt (order h B e for orders <= 7) and
   hp8.txt (order strong h P_D(2) for every class of orders <= 8).  Records the first representation of each target. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
static int n; static u64 adj[12], radj[12]; static i64 Hm[1<<10], dp[1<<10][10], dpr[1<<10][10];
static unsigned char *nj, *njx; static long L;
#define MAXT 64
static long targets[MAXT]; static int nt = 0; static char rep[MAXT][128];
static u64 reachfrom(u64 s, u64 *A) { u64 r = s, p = 0; while (r != p) { p = r; for (int v = 0; v < n; v++) if ((r>>v)&1) r |= A[v]; } return r; }
static void setbit(unsigned char *b, long v) { b[v >> 3] |= 1 << (v & 7); }
static FILE *fp7, *fhp;
static void process(void) {
    for (int v = 0; v < n; v++) { radj[v] = 0; for (int u = 0; u < n; u++) if ((adj[u]>>v)&1) radj[v] |= 1ULL<<u; }
    u64 full = (1ULL<<n)-1; int strong = (reachfrom(1, adj) == full && reachfrom(1, radj) == full);
    int N = 1 << n;
    for (int m = 0; m < N; m++) for (int l = 0; l < n; l++) { dp[m][l] = 0; dpr[m][l] = 0; }
    for (int v = 0; v < n; v++) { dp[1<<v][v] = 1; dpr[1<<v][v] = 1; }
    for (int m = 1; m < N; m++) for (int l = 0; l < n; l++) {
        i64 a = dp[m][l]; if (a) { u64 nx = adj[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; dp[m|(1<<v)][v] += a; } }
        i64 b = dpr[m][l]; if (b) { u64 nx = radj[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; dpr[m|(1<<v)][v] += b; } }
    }
    for (int m = 1; m < N; m++) { i64 s = 0; for (int l = 0; l < n; l++) s += dp[m][l]; Hm[m] = s; } Hm[0] = 1;
    i64 h = Hm[N-1];
    if (n <= 8) { i64 P = 0; for (int m = 0; m < N; m++) P += Hm[m] * Hm[(N-1) & ~m]; fprintf(fhp, "%d %d %lld %lld\n", n, strong, h, P); }
    for (int t = 0; t < n; t++) {
        i64 B = 0; for (int m = 1; m < N; m++) if ((m>>t)&1) B += dpr[m][t] * Hm[(N-1) & ~m];
        i64 e = dp[N-1][t];
        if (n <= 7) fprintf(fp7, "%d %lld %lld %lld\n", n, h, B, e);
        if (e == 0 && h + B <= L) { setbit(njx, h + B); for (int i = 0; i < nt; i++) if (targets[i] == h + B && !rep[i][0]) snprintf(rep[i], 128, "k=0 dead-end append: order %d h %lld B %lld", n, h, B); }
        for (int k = 1; k < 40; k++) { i64 v = h + (B << k); if (v > L) break; setbit(nj, v); setbit(njx, v);
            for (int i = 0; i < nt; i++) if (targets[i] == v && !rep[i][0]) snprintf(rep[i], 128, "order %d h %lld B %lld k %d", n, h, B, k); }
    }
}
int main(int argc, char **argv) {
    L = atol(argv[1]); nj = calloc(L/8 + 2, 1); njx = calloc(L/8 + 2, 1);
    FILE *ft = fopen(argv[2], "r"); long v; while (ft && fscanf(ft, "%ld", &v) == 1 && nt < MAXT) { targets[nt] = v; rep[nt][0] = 0; nt++; } if (ft) fclose(ft);
    fp7 = fopen("pairs7.txt", "w"); fhp = fopen("hp8.txt", "w");
    long classes = 0;
    n = 1; adj[0] = 0; process(); classes++;
    n = 2; adj[0] = 2; adj[1] = 0; process(); classes++;
    char line[256];
    for (n = 3; n <= 9; n++) {
        char cmd[64]; snprintf(cmd, 64, "gentourng -q %d", n); FILE *g = popen(cmd, "r");
        while (fgets(line, sizeof line, g)) {
            int len = strlen(line); while (len && (line[len-1]=='\n' || line[len-1]=='\r')) line[--len] = 0; if (!len) continue; classes++;
            for (int i = 0; i < n; i++) adj[i] = 0; int k = 0;
            for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++) { if (line[k++] == '1') adj[i] |= 1ULL<<j; else adj[j] |= 1ULL<<i; }
            process();
        }
        pclose(g);
    }
    fclose(fp7); fclose(fhp);
    FILE *o = fopen("nj9.bits", "wb"); fwrite(nj, 1, L/8 + 2, o); fclose(o);
    o = fopen("njext9.bits", "wb"); fwrite(njx, 1, L/8 + 2, o); fclose(o);
    printf("classes of orders 1-9: %ld; L = %ld\n", classes, L);
    for (int pass = 0; pass < 2; pass++) {
        unsigned char *b = pass ? njx : nj; const char *name = pass ? "near-join + dead-end appends (k >= 0)" : "near-join (k >= 1)";
        long first = -1, miss = 0; long firstlist[40]; int nf = 0;
        for (long x = 69; x <= L; x += 2) if (!(b[x>>3] & (1 << (x&7)))) { miss++; if (first < 0) first = x; if (nf < 40) firstlist[nf++] = x; }
        printf("%s: every odd <= %ld except the seven is realised (first miss above 67: %ld); odd misses in (67, %ld]: %ld; first misses:", name, first - 2, first, L, miss);
        for (int i = 0; i < nf; i++) printf(" %ld", firstlist[i]); printf("\n");
        printf("  gap density of odd numbers by dyadic window [2^j, 2^(j+1)):");
        for (int j = 14; (1L << j) < L; j++) { long lo = 1L << j, hi = (1L << (j+1)) < L ? (1L << (j+1)) : L; long m = 0, tot = 0; for (long x = lo | 1; x < hi; x += 2) { tot++; if (!(b[x>>3] & (1 << (x&7)))) m++; } printf(" j=%d: %ld/%ld", j, m, tot); }
        printf("\n");
    }
    for (int i = 0; i < nt; i++) printf("target %ld: %s\n", targets[i], rep[i][0] ? rep[i] : "NOT reached from bases <= 9");
    return 0;
}

/* odometer.c <n> <L>: generalised binary port-move (odometer) search, containing wave 45's w45_tworoot.py as the special
   case P_i = {x_i} (z beats only the root x_i) and its mirror.  Base: every tournament class of order n (gentourng, stdin).
   A letter i (i = 0, 1) appends a vertex z that beats exactly the original vertices in P_i (any subset) and, among the
   previously appended vertices, either beats all of them (r_i = 1) or loses to all of them (r_i = 0).
   Binary odometer: H(w i) = 2 H(w) + d_i for every word w of length < L, with d_0 != d_1.  Reports every surviving
   (base, P_0, r_0, P_1, r_1) with its digits, and the counts of candidates surviving each length. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64; typedef long long i64;
static int n; static u64 adj[20];
static i64 dpv[1<<12][12];
static i64 ham(void) {
    int N = 1 << n; for (int m = 0; m < N; m++) for (int l = 0; l < n; l++) dpv[m][l] = 0;
    for (int v = 0; v < n; v++) dpv[1<<v][v] = 1;
    for (int m = 1; m < N; m++) for (int l = 0; l < n; l++) { i64 a = dpv[m][l]; if (!a) continue; u64 nx = adj[l] & ~(u64)m; while (nx) { int v = __builtin_ctzll(nx); nx &= nx-1; dpv[m|(1<<v)][v] += a; } }
    i64 s = 0; for (int l = 0; l < n; l++) s += dpv[(N-1)][l]; return s;
}
static int n0; static u64 base[20];
static i64 Hword(int *w, int len, int P[2], int r[2]) {
    n = n0; for (int i = 0; i < n0; i++) adj[i] = base[i];
    for (int t = 0; t < len; t++) {
        int i = w[t]; int z = n;
        adj[z] = 0;
        for (int v = 0; v < n0; v++) { if ((P[i] >> v) & 1) adj[z] |= 1ULL << v; else adj[v] |= 1ULL << z; }
        for (int v = n0; v < n; v++) { if (r[i]) adj[z] |= 1ULL << v; else adj[v] |= 1ULL << z; }
        n++;
    }
    return ham();
}
int main(int argc, char **argv) {
    n0 = atoi(argv[1]); int L = atoi(argv[2]); char line[256]; long classes = 0, cand = 0, surv1 = 0, surv2 = 0, surv3 = 0, found = 0;
    while (fgets(line, sizeof line, stdin)) {
        int len = strlen(line); while (len && (line[len-1]=='\n' || line[len-1]=='\r')) line[--len] = 0; if (!len) continue; classes++;
        int k = 0; for (int i = 0; i < n0; i++) base[i] = 0;
        for (int i = 0; i < n0; i++) for (int j = i+1; j < n0; j++) { if (line[k++] == '1') base[i] |= 1ULL<<j; else base[j] |= 1ULL<<i; }
        int w[8]; i64 H0 = Hword(w, 0, (int[]){0,0}, (int[]){0,0});
        int N0 = 1 << n0;
        /* precompute single-letter digits for every (P, r) */
        i64 *d1 = malloc(sizeof(i64) * N0 * 2);
        for (int P = 0; P < N0; P++) for (int r = 0; r < 2; r++) { int Pa[2] = {P, P}, ra[2] = {r, r}; w[0] = 0; d1[P*2+r] = Hword(w, 1, Pa, ra) - 2*H0; }
        for (int P0 = 0; P0 < N0; P0++) for (int r0 = 0; r0 < 2; r0++) for (int P1 = 0; P1 < N0; P1++) for (int r1 = 0; r1 < 2; r1++) {
            if (P1*2+r1 <= P0*2+r0) continue; /* unordered letter pairs */
            i64 dd[2] = {d1[P0*2+r0], d1[P1*2+r1]}; if (dd[0] == dd[1]) continue;
            cand++; int P[2] = {P0, P1}, r[2] = {r0, r1}; int ok = 1;
            /* length-2 words: w = (a), then letter b */
            for (int a = 0; a < 2 && ok; a++) { w[0] = a; i64 hw = Hword(w, 1, P, r); for (int b = 0; b < 2 && ok; b++) { w[1] = b; if (Hword(w, 2, P, r) - 2*hw != dd[b]) ok = 0; } }
            if (!ok) continue; surv1++;
            for (int a = 0; a < 2 && ok; a++) for (int b = 0; b < 2 && ok; b++) { w[0] = a; w[1] = b; i64 hw = Hword(w, 2, P, r); for (int c = 0; c < 2 && ok; c++) { w[2] = c; if (Hword(w, 3, P, r) - 2*hw != dd[c]) ok = 0; } }
            if (!ok) continue; surv2++;
            if (L >= 4) { for (int a = 0; a < 2 && ok; a++) for (int b = 0; b < 2 && ok; b++) for (int c = 0; c < 2 && ok; c++) { w[0]=a; w[1]=b; w[2]=c; i64 hw = Hword(w, 3, P, r); for (int e = 0; e < 2 && ok; e++) { w[3] = e; if (Hword(w, 4, P, r) - 2*hw != dd[e]) ok = 0; } } if (!ok) continue; surv3++; }
            found++;
            if (found <= 40) printf("ODOMETER base=%s H0=%lld letters (P0=%d r0=%d d0=%lld) (P1=%d r1=%d d1=%lld)\n", line, H0, P0, r0, dd[0], P1, r1, dd[1]);
        }
        free(d1);
    }
    printf("order %d: %ld classes; letter pairs with distinct digits %ld; surviving words of length 2: %ld; length 3: %ld; length 4: %ld; odometers found: %ld\n", n0, classes, cand, surv1, surv2, surv3, found);
    return 0;
}

/* Audit2 (D): H of a vertex-transitive tournament by inclusion-exclusion over vertex subsets, one subset per orbit of a
   given permutation group G (orbit size |G|/|Stab(X)|), __int128 arithmetic.  Input on stdin:
     n  m
     upper-triangle code ('1' at (i,j), i<j, iff i -> j)
     m lines, each a permutation of 0..n-1 (the group elements, identity included)
   n <= 30.  Used for the Paley tournament on F_27 with G = translations (Z_3)^3. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef __int128 i128;
static int n, m, A[32][32];
static uint32_t tab[64][4][256];   /* image of byte b in byte position p under group element g */

static void print_i128(i128 x) {
    char buf[64]; int k = 0, neg = 0;
    if (x < 0) { neg = 1; x = -x; }
    if (x == 0) buf[k++] = '0';
    while (x > 0) { buf[k++] = '0' + (int)(x % 10); x /= 10; }
    if (neg) putchar('-');
    while (k > 0) putchar(buf[--k]);
}

static i128 walks(uint32_t X) {
    i128 u[32], w[32]; int vs[32], k = 0;
    for (int v = 0; v < n; v++) if ((X >> v) & 1) vs[k++] = v;
    for (int i = 0; i < k; i++) u[i] = 1;
    for (int step = 0; step < n - 1; step++) {
        for (int j = 0; j < k; j++) {
            i128 s = 0; int b = vs[j];
            for (int i = 0; i < k; i++) if (A[vs[i]][b]) s += u[i];
            w[j] = s;
        }
        memcpy(u, w, sizeof(i128) * k);
    }
    i128 t = 0;
    for (int i = 0; i < k; i++) t += u[i];
    return t;
}

int main(void) {
    char code[1024];
    if (scanf("%d %d", &n, &m) != 2) return 1;
    if (scanf("%1000s", code) != 1) return 1;
    int k = 0;
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (code[k] == '1') A[i][j] = 1; else A[j][i] = 1; k++; }
    static int perm[64][32];
    for (int g = 0; g < m; g++) for (int v = 0; v < n; v++) if (scanf("%d", &perm[g][v]) != 1) return 1;
    for (int g = 0; g < m; g++)
        for (int p = 0; p < 4; p++)
            for (int b = 0; b < 256; b++) {
                uint32_t img = 0;
                for (int t = 0; t < 8; t++) {
                    int v = 8 * p + t;
                    if (((b >> t) & 1) && v < n) img |= (uint32_t)1 << perm[g][v];
                }
                tab[g][p][b] = img;
            }
    i128 H = 0;
    uint32_t M = (uint32_t)1 << n;
    long orbits = 0;
    for (uint32_t X = 1; X < M; X++) {
        int stab = 0, minimal = 1;
        for (int g = 0; g < m; g++) {
            uint32_t y = tab[g][0][X & 255] | tab[g][1][(X >> 8) & 255] | tab[g][2][(X >> 16) & 255] | tab[g][3][(X >> 24) & 255];
            if (y < X) { minimal = 0; break; }
            if (y == X) stab++;
        }
        if (!minimal) continue;
        orbits++;
        int sz = __builtin_popcount(X);
        i128 W = walks(X) * (m / stab);
        if ((n - sz) & 1) H -= W; else H += W;
    }
    fprintf(stderr, "orbits: %ld\n", orbits);
    print_i128(H); printf("\n");
    return 0;
}

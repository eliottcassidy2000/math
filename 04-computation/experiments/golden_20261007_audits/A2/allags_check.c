/* A2 audit: q_inf(T) = N_T / 2^T.
   (1) N_T = #distinct (a, c mod 3^a) over all y mod 2^T, c = 2^T T^T(y) - 3^a y.
   (2) brute force for small T: fraction of y mod 2^T with no r in [1, 2^(T+2)] such that T^T(Y+r) = T^T(Y)
       for two independent large lifts Y = y + 2^T M (identity merge by time T). */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef __int128 i128;
static int cmpu(const void *a, const void *b) { uint64_t x = *(uint64_t *)a, y = *(uint64_t *)b; return x < y ? -1 : x > y; }
static i128 Tn(i128 x, int n, int *a) { int c = 0; for (int i = 0; i < n; i++) { if (x & 1) { x = (3 * x + 1) / 2; c++; } else x /= 2; } if (a) *a = c; return x; }
int main(int argc, char **argv) {
    int TMAX = atoi(argv[1]), TB = atoi(argv[2]);
    uint64_t p3[64]; p3[0] = 1; for (int i = 1; i < 40; i++) p3[i] = p3[i - 1] * 3;
    for (int T = 4; T <= TMAX; T += 2) {
        uint64_t n = 1ULL << T; uint64_t *key = malloc(n * sizeof(uint64_t));
        for (uint64_t y = 0; y < n; y++) {
            int a; i128 x = Tn((i128)y, T, &a);
            i128 c = ((i128)x << T) - (i128)p3[a] * (i128)y;
            if (c < 0) { fprintf(stderr, "neg c\n"); return 1; }
            uint64_t cm = (uint64_t)(c % (i128)p3[a]);
            key[y] = (uint64_t)a * p3[T] + cm;
        }
        qsort(key, n, sizeof(uint64_t), cmpu);
        uint64_t N = 1; for (uint64_t i = 1; i < n; i++) if (key[i] != key[i - 1]) N++;
        printf("T=%2d N_T=%10llu q_inf=N_T/2^T=%.4f\n", T, (unsigned long long)N, (double)N / n);
        free(key);
    }
    /* brute force */
    srand(12345);
    for (int T = 4; T <= TB; T += 2) {
        uint64_t n = 1ULL << T, surv = 0; i128 M1 = ((i128)rand() << 30) ^ rand(), M2 = ((i128)rand() << 31) ^ rand() ^ 77777;
        for (uint64_t y = 0; y < n; y++) {
            i128 Y1 = (i128)y + (M1 << T), Y2 = (i128)y + (M2 << T);
            i128 t1 = Tn(Y1, T, 0), t2 = Tn(Y2, T, 0);
            int merged = 0;
            for (uint64_t r = 1; r <= (n << 2) && !merged; r++) {
                if (Tn(Y1 + r, T, 0) == t1 && Tn(Y2 + r, T, 0) == t2) merged = 1;
            }
            if (!merged) surv++;
        }
        printf("brute T=%2d: fraction never merging with any y+r (r<=2^(T+2)) = %.4f\n", T, (double)surv / n);
    }
    return 0;
}

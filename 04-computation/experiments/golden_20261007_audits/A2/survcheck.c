/* HYP-9231: survivors (3^(a_j) > 2^j for all j<=s) of equal length s and weight a never collide mod 3^a. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static uint64_t P3[64]; static int SMAX;
static uint64_t *keys[64]; static uint64_t cnt[64], cap[64];
static void push(int s, uint64_t k) { if (cnt[s] == cap[s]) { cap[s] = cap[s] ? 2 * cap[s] : 1024; keys[s] = realloc(keys[s], cap[s] * 8); } keys[s][cnt[s]++] = k; }
static void dfs(int s, int a, u128 c) {
    if (s > 0) push(s, (uint64_t)a * P3[SMAX] + (uint64_t)(c % P3[a]));
    if (s == SMAX) return;
    /* even: a, c unchanged; survivor iff 3^a > 2^(s+1) */
    if ((u128)P3[a] > ((u128)1 << (s + 1))) dfs(s + 1, a, c);
    /* odd: always survives (3^(a+1) > 2^(s+1) since 3^a > 2^s) */
    dfs(s + 1, a + 1, 3 * c + ((u128)1 << s));
}
static int cmp(const void *x, const void *y) { uint64_t a = *(uint64_t *)x, b = *(uint64_t *)y; return a < b ? -1 : a > b; }
int main(int argc, char **argv) {
    SMAX = atoi(argv[1]); P3[0] = 1; for (int i = 1; i < 40; i++) P3[i] = P3[i - 1] * 3;
    dfs(0, 0, 0);
    for (int s = 1; s <= SMAX; s++) {
        qsort(keys[s], cnt[s], 8, cmp); uint64_t col = 0;
        for (uint64_t i = 1; i < cnt[s]; i++) if (keys[s][i] == keys[s][i - 1]) col++;
        printf("s=%2d survivors=%10llu collisions=%llu\n", s, (unsigned long long)cnt[s], (unsigned long long)col);
    }
    return 0;
}

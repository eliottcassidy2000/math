/* In-order caterpillar automaton on (h,g,r,s) at the current spine end:
   L (leaf at the end):        (h+2r, g, r+s-g, r+s)
   S (extend spine, port moves): (h+2r, h, s, h+r)
   start (1,1,0,1).  BFS over states with h <= B; report which odd h <= B are reached. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef long long ll;
typedef struct { ll h, g, r, s; } St;
#define HS (1 << 26)
static St *tab; static unsigned char *used;
static unsigned long hsh(St x) { unsigned long k = (unsigned long)x.h * 1000003UL ^ (unsigned long)x.g * 998244353UL ^ (unsigned long)x.r * 2654435761UL ^ (unsigned long)x.s * 40503UL; return (k ^ (k >> 29)) & (HS - 1); }
static int insert(St x) { unsigned long i = hsh(x); while (used[i]) { if (!memcmp(&tab[i], &x, sizeof x)) return 0; i = (i + 1) & (HS - 1); } used[i] = 1; tab[i] = x; return 1; }
int main(int argc, char **argv) {
    ll B = atoll(argv[1]);
    tab = malloc(sizeof(St) * HS); used = calloc(HS, 1);
    unsigned char *hit = calloc(B + 1, 1);
    St *q = malloc(sizeof(St) * 40000000); long qh = 0, qt = 0;
    St s0 = {1, 1, 0, 1}; insert(s0); q[qt++] = s0; hit[1] = 1;
    while (qh < qt) {
        St x = q[qh++];
        St a = {x.h + 2 * x.r, x.g, x.r + x.s - x.g, x.r + x.s};
        St b = {x.h + 2 * x.r, x.h, x.s, x.h + x.r};
        St nx[2] = {a, b};
        for (int k = 0; k < 2; k++) { St y = nx[k]; if (y.h > B) continue;
            if (y.h == x.h && y.g == x.g && y.r == x.r && y.s == x.s) continue;
            if (insert(y)) { if (qt >= 40000000) { fprintf(stderr, "queue full\n"); goto done; } q[qt++] = y; hit[y.h] = 1; } }
    }
done:;
    ll bound = 0; int nm = 0; printf("states %ld; misses (odd, <= B):", qt);
    for (ll h = 1; h <= B; h += 2) if (!hit[h]) { if (nm < 25) printf(" %lld", h); nm++; }
    printf("  [total misses %d]\n", nm);
    for (ll h = 1; h <= B; h += 2) if (!hit[h] && h != 7 && h != 21) { bound = h - 2; break; }
    printf("complete except 7, 21 up to %lld (B = %lld)\n", bound, B);
}

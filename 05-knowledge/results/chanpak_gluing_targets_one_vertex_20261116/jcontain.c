/* jcontain.c -- which prime pairs (p,q) can sit, as induced subtournaments,
 * inside a tournament whose count is exactly N = p + q - 1 ?
 *
 * For every tournament T of order n (gentourng input) compute H of every
 * induced subtournament.  For every pair of odd primes p <= q with
 * p + q - 1 = H(T), report whether
 *   C: T has induced subtournaments X, Y with H(X)=p, H(Y)=q (any position);
 *   K: moreover X u Y = V(T) and neither contains the other (cover form);
 *   D: moreover X = V-a, Y = V-b for two vertices (vertex-exchange form).
 * Output one line per (p,q,N,form) the first time it is seen, with an example.
 *
 * Completeness (see JGLUE.md): if p,q >= 3 are primes with p + q - 1 = N and
 * X, Y are induced witnesses in any tournament with H = N, their strong cores
 * lie in one strong component S with H(S) = N, and f(|S|) <= N where f is the
 * Moon-Busch strong floor.  For N < f(10) = 125, |S| <= 9, so the census of
 * orders <= 9 decides C for all such N.
 * usage: gentourng n | ./jcontain n
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#define MAXN 10
static int n, N;
static uint32_t out_[MAXN], H[1 << MAXN], dp[1 << MAXN][MAXN];
static char isprime[1 << 16];

static void sieve(void) { memset(isprime, 1, sizeof isprime); isprime[0] = isprime[1] = 0;
    for (int i = 2; i * i < (1 << 16); i++) if (isprime[i]) for (int j = i * i; j < (1 << 16); j += i) isprime[j] = 0; }
static void allH(void) {
    H[0] = 1;
    for (int S = 1; S < N; S++) { uint32_t tot = 0;
        for (int v = 0; v < n; v++) { if (!(S >> v & 1)) continue; int R = S & ~(1 << v); uint32_t c;
            if (!R) c = 1; else { c = 0; for (int u = 0; u < n; u++) if ((R >> u & 1) && (out_[u] >> v & 1)) c += dp[R][u]; }
            dp[S][v] = c; tot += c; }
        H[S] = tot; }
}
static unsigned char flag[2048][2048]; /* flag[p/2][q/2] bit0=C bit1=K bit2=D */
int main(int argc, char **argv) {
    n = atoi(argv[1]); N = 1 << n; sieve();
    char line[256]; int V = N - 1;
    static unsigned char present[1 << 13];
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == '>') continue;
        int len = strlen(line); while (len && (line[len-1]=='\n'||line[len-1]=='\r')) line[--len]=0;
        if (len != n*(n-1)/2) continue;
        memset(out_, 0, sizeof out_);
        int k = 0; for (int i = 0; i < n; i++) for (int j = i+1; j < n; j++, k++) { if (line[k]=='1') out_[i] |= 1u<<j; else out_[j] |= 1u<<i; }
        allH();
        uint32_t HV = H[V];
        if (getenv("ROWS")) {
            if (HV < 9 || HV > 4000) continue;
            static unsigned char pres[1 << 13];
            memset(pres, 0, HV + 4);
            for (int S = 1; S < V; S++) if (H[S] <= HV + 3) pres[H[S]] = 1;
            static int c_of[3] = {1, -1, -3};
            for (int ci = 0; ci < 3; ci++) {
                int s = (int)HV + c_of[ci];  /* p + q = N + c */
                for (int p = 3; 2 * p <= s; p += 2) {
                    int q = s - p;
                    if (p == 7 || q == 7 || !isprime[p] || !isprime[q]) continue;
                    if (pres[p] && pres[q]) printf("R %u %d %d %d %d\n", HV, p, q, c_of[ci], n);
                }
            }
            continue;
        }
        if (HV < 21) continue;
        memset(present, 0, HV + 1);
        for (int S = 1; S < V; S++) if (H[S] <= HV) present[H[S]] = 1;
        for (uint32_t p = 3; 2 * p <= HV + 1; p += 2) {
            uint32_t q = HV + 1 - p;
            if (!isprime[p] || !isprime[q]) continue;
            int f = 0;
            if (present[p] && present[q]) f |= 1;
            if (f) {
                /* cover and deletion forms */
                for (int A = 1; A < V && !(f & 2); A++) {
                    if (H[A] != p && H[A] != q) continue;
                    uint32_t need = (H[A] == p) ? q : p;
                    uint32_t compA = V & ~A;
                    for (uint32_t X = A;; X = (X - 1) & A) {
                        uint32_t B = compA | X;
                        if (X != A && (A & ~B) && H[B] == need) { f |= 2; break; }
                        if (X == 0) break;
                    }
                }
                for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) if (a != b) {
                    uint32_t A = V & ~(1u << a), B = V & ~(1u << b);
                    if (H[A] == p && H[B] == q) f |= 4;
                }
            }
            int old = flag[p/2][q/2];
            if ((f | old) != old) {
                for (int bit = 0; bit < 3; bit++) if ((f >> bit & 1) && !(old >> bit & 1))
                    printf("%c p=%u q=%u N=%u n=%d T=%s\n", "CKD"[bit], p, q, HV, n, line);
                flag[p/2][q/2] = old | f;
            }
        }
    }
    return 0;
}

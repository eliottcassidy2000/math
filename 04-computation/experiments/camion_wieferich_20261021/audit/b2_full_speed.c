/* Audit b2 (independent code): THM-4750 statement 2 (full speed at Wieferich primes, claimed threshold m >= k)
 * and the general-base claim.  s_a(n) mod p^K is computed from the Hurwitz INVERSE of T_a(n) = a^C(n,2):
 *     W(0) = 1,  W(n) = -sum_{k=1}^{n} C(n,k) T_a(k) W(n-k),   s_a(n) = -W(n)   (n >= 1)
 * (a different recurrence from the session's Wright recurrence).  For a = 2 the circulant count is subtracted (D = s - c);
 * for a != 2 the plain s_a is used.  eps = (a/p) by Euler's criterion.
 * For every base m in [1, MTOP] the Mahler coefficients a_i(m) = Delta^i A_m(0), A_m(l) = eps^(-l) D(m + l(p-1)),
 * i = 0..2K-1, are computed mod p^K and compared with
 *     FULL  speed bound min(K, i)          (claimed for m >= K at a base-a Wieferich prime),
 *     HALF  speed bound min(K, ceil(i/2))  (claimed for m >= 2K-1 always).
 * Prints the bases m violating each bound, and the minimum over the p-1 classes in the window m in [K, K+p-2].
 * Compile into scratch: cc -O2 -o b2 b2_full_speed.c ; usage: ./b2 p K [a] */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static u64 Q;
static inline u64 mm(u64 x, u64 y){ return (u64)(((u128)x * y) % Q); }
static u64 pwq(u64 b, u64 e){ u64 r = 1 % Q; b %= Q; while(e){ if(e & 1) r = mm(r, b); b = mm(b, b); e >>= 1; } return r; }
static int val(u64 x, u64 p, int K){ if(x == 0) return K; int v = 0; while(x % p == 0){ x /= p; v++; } return v; }
int main(int argc, char **argv){
  if(argc < 3){ fprintf(stderr, "usage: %s p K [a]\n", argv[0]); return 1; }
  u64 p = strtoull(argv[1], 0, 10); int K = atoi(argv[2]); u64 a = argc > 3 ? strtoull(argv[3], 0, 10) : 2;
  Q = 1; for(int i = 0; i < K; i++) Q *= p;
  /* eps = (a/p): a^((p-1)/2) mod p */
  u64 e = 1, b = a % p, ex = (p - 1) / 2; while(ex){ if(ex & 1) e = e * b % p; b = b * b % p; ex >>= 1; }
  int eps = (e == 1) ? 1 : -1;
  int R = 2 * K;                                     /* i = 0..2K-1 */
  long MTOP = (long)K + (long)p + 2;                 /* bases m = 1..MTOP (covers the window [K, K+p-2] and more) */
  long N = MTOP + (long)(R - 1) * (long)(p - 1);
  u64 *T = malloc((N + 1) * 8), *W = malloc((N + 1) * 8), *row = calloc(N + 2, 8), *Dv = malloc((N + 1) * 8);
  for(long n = 0; n <= N; n++){ u128 c2 = (u128)n * (u128)(n > 0 ? n - 1 : 0) / 2; T[n] = pwq(a, (u64)c2); }
  W[0] = 1; row[0] = 1;
  for(long n = 1; n <= N; n++){
    for(long i = n; i >= 1; i--){ u64 t = row[i] + row[i - 1]; row[i] = (t >= Q) ? t - Q : t; }   /* row = C(n, .) */
    u128 acc = 0;
    for(long k = 1; k <= n; k++){ acc += (u128)mm(row[k], T[k]) * W[n - k]; if((k & 31) == 0) acc %= Q; }
    acc %= Q; W[n] = (Q - (u64)acc) % Q;
  }
  for(long n = 1; n <= N; n++){
    u64 sn = (Q - W[n]) % Q;
    u64 cn = (a == 2 && (n & 1)) ? pwq(2, (u64)(n - 1) / 2) : 0;
    Dv[n] = (sn + Q - cn) % Q;
  }
  u64 qa = (pwq(a, p - 1) + Q - 1) % Q; int wq = val(qa, p, K) - 1;
  printf("p=%llu K=%d base=%llu eps=%+d v_p(q_a)=%s%d  (s to n=%ld via Hurwitz inverse)\n", (unsigned long long)p, K,
         (unsigned long long)a, eps, (wq >= K - 1) ? ">=" : "", wq, N);
  int minwin[64]; for(int i = 0; i < R; i++) minwin[i] = K;
  long nfull_bad_ge = 0, nhalf_bad_ge = 0; char fb[4096]; int fl = 0; fb[0] = 0; char lb[4096]; int ll = 0; lb[0] = 0;
  for(long m = 1; m <= MTOP; m++){
    u64 t[64];
    for(int l = 0; l < R; l++){ u64 d = Dv[m + (long)l * (long)(p - 1)]; t[l] = (eps == 1 || (l % 2 == 0)) ? d : (Q - d) % Q; }
    int v[64];
    for(int i = 0; i < R; i++){ v[i] = val(t[0], p, K); for(int l = 0; l < R - 1 - i; l++) t[l] = (t[l + 1] + Q - t[l]) % Q; }
    int fullbad = 0, halfbad = 0;
    for(int i = 0; i < R; i++){ int fbd = i < K ? i : K; int hbd = (i + 1) / 2 < K ? (i + 1) / 2 : K;
      if(v[i] < fbd) fullbad = 1; if(v[i] < hbd) halfbad = 1; }
    if(m >= K && m <= K + (long)p - 2) for(int i = 0; i < R; i++) if(v[i] < minwin[i]) minwin[i] = v[i];
    if(fullbad){ if(m >= K) nfull_bad_ge++; if(fl < 3900) fl += snprintf(fb + fl, 4096 - fl, " %ld", m); }
    if(halfbad){ if(m >= 2 * K - 1) nhalf_bad_ge++; if(ll < 3900) ll += snprintf(lb + ll, 4096 - ll, " %ld", m); }
  }
  printf("  bases m in [1,%ld] violating FULL bound min(K,i): %s  (count at m >= K: %ld)\n", MTOP, fl ? fb : " none", nfull_bad_ge);
  printf("  bases m in [1,%ld] violating HALF bound min(K,ceil(i/2)): %s  (count at m >= 2K-1: %ld)\n", MTOP, ll ? lb : " none", nhalf_bad_ge);
  printf("  min over the p-1 classes, window m in [%d,%ld], of v_p(a_i) (cap K): ", K, (long)K + (long)p - 2);
  for(int i = 0; i < R; i++) printf("%d ", minwin[i]); printf("\n");
  return 0;
}

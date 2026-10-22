/* Audit b8 (independent code): THM-4750 statement 7 (connected 3-uniform hypergraphs b(n)).
 * b(n) mod p^K is computed through the rooted Hurwitz identity  Delta b = (Delta T) * T^(*-1)  (T(n) = 2^C(n,3)):
 *     W = T^(*-1):  W(0) = 1, W(n) = -sum_{k=1}^n C(n,k) T(k) W(n-k);   b(n+1) = sum_{k=0}^n C(n,k) T(k+1) W(n-k),
 * a different route from the session's direct recurrence 2^C(n,3) = sum_k C(n-1,k-1) b(k) T(n-k).
 * For every prime p <= PMAX with gcd(ord_p 2, 6) = 1 (the setting of statement 7):
 *   (1) mod p^2 (K=2): for every base m in [1, p+4], the "cubic law" p | a_1, a_2, a_3 and p^2 | a_4..a_6 (multiplier 1);
 *       the bases m that FAIL are listed (claimed sufficient threshold: m >= r(k-1)+2 = 5);
 *   (2) class 1: b(p) = 2^((p-1)(p-2)/6) - 1 (mod p) and = 0 mod p;
 *   (3) every zero class (base m in [6, p+4], as in the session): number of roots in F_p of the cubic
 *       g(l) = c0 + c1 l + c2 C(l,2) + c3 C(l,3), c_i = a_i/p, histogrammed by class type (1, 2, sporadic);
 *   (4) class 2: does the cubic of base m = p+1 vanish at l = -1 (the point n = 2, where b(2) = 0 exactly)?
 *       and does the base m = 2 itself already satisfy the cubic law (i.e. is n = 2 inside the polynomial range)?
 *   (5) for p <= P3 also mod p^3 (K=3): Mahler speed 1/3, v(a_i) >= min(3, ceil(i/3)), i <= 9, listing failing bases
 *       (claimed sufficient threshold m >= r(k-1)+2 = 8).
 * Compile into scratch: cc -O2 -o b8 b8_hyper3.c ; usage: ./b8 PMAX P3 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static u64 Q;
static inline u64 mm(u64 a, u64 b){ return (u64)((u128)a * b % Q); }
static u64 pwm(u64 b, u64 e, u64 m){ u64 r = 1 % m; b %= m; while(e){ if(e & 1) r = (u64)((u128)r * b % m); b = (u64)((u128)b * b % m); e >>= 1; } return r; }
static int isprime(u64 n){ if(n < 2) return 0; for(u64 d = 2; d * d <= n; d++) if(n % d == 0) return 0; return 1; }
static u64 ord2(u64 p){ u64 o = 1, x = 2 % p; while(x != 1){ x = x * 2 % p; o++; } return o; }
static u64 gcdu(u64 a, u64 b){ while(b){ u64 t = a % b; a = b; b = t; } return a; }
static int vpK(u64 x, u64 p, int K){ if(x == 0) return K; int v = 0; while(x % p == 0){ x /= p; v++; } return v; }
/* b(0..N) mod q */
static u64 *b_upto(long N, u64 q, u64 p){
  Q = q; u64 phi = q / p * (p - 1);
  u64 *T = malloc((N + 2) * 8), *W = malloc((N + 2) * 8), *row = calloc(N + 3, 8), *b = calloc(N + 2, 8);
  for(long n = 0; n <= N + 1; n++){ u128 c3 = (u128)n * (u128)(n >= 1 ? n - 1 : 0) * (u128)(n >= 2 ? n - 2 : 0) / 6; T[n] = pwm(2, (u64)(c3 % phi), q); }
  W[0] = 1; row[0] = 1;
  /* W and b together; row = C(n, .) */
  b[1] = T[1] % q;                       /* Delta b(0) = T(1) W(0) */
  for(long n = 1; n <= N; n++){
    for(long i = n; i >= 1; i--){ u64 t = row[i] + row[i - 1]; row[i] = t >= q ? t - q : t; }
    u128 accW = 0, accB = 0;
    for(long k = 1; k <= n; k++){ u64 c = mm(row[k], T[k]); accW += (u128)c * W[n - k]; if((k & 31) == 0) accW %= q; }
    accW %= q; W[n] = (q - (u64)accW) % q;
    for(long k = 0; k <= n; k++){ u64 c = mm(row[k], T[k + 1]); accB += (u128)c * W[n - k]; if((k & 31) == 0) accB %= q; }
    accB %= q; if(n + 1 <= N) b[n + 1] = (u64)accB;
  }
  free(T); free(W); free(row); return b;
}
int main(int argc, char **argv){
  u64 PMAX = argc > 1 ? strtoull(argv[1], 0, 10) : 2200, P3 = argc > 2 ? strtoull(argv[2], 0, 10) : 200;
  long nprimes = 0, law_fail_ge5 = 0, cls1_bad = 0, nzc = 0; long hist[3][5] = {{0}};
  long c2_forced = 0, c2_zero = 0, base2_inside = 0, nclasses = 0, notcubic = 0;
  char failset[64] = {0}; /* failing bases m <= 63 seen at K=2 */
  for(u64 p = 5; p <= PMAX; p++){
    if(!isprime(p)) continue; u64 t = ord2(p); if(gcdu(t, 6) != 1) continue;
    nprimes++;
    u64 q = p * p; long N = (long)p + 4 + 6 * (long)(p - 1) + 2;
    u64 *b = b_upto(N, q, p);
    /* (2) class 1 */
    u64 e = ((p - 1) * (p - 2) / 6); u64 pred = (pwm(2, e, p) + p - 1) % p;
    if(b[p] % p != pred || b[p] % p != 0) cls1_bad++;
    char fl[512]; int fll = 0; fl[0] = 0;
    for(long m = 1; m <= (long)p + 4; m++){
      u64 A[7]; for(int l = 0; l < 7; l++) A[l] = b[m + (long)l * (long)(p - 1)];
      u64 a[7], tt[7]; for(int l = 0; l < 7; l++) tt[l] = A[l];
      for(int i = 0; i < 7; i++){ a[i] = tt[0]; for(int l = 0; l < 6 - i; l++) tt[l] = (tt[l + 1] + q - tt[l]) % q; }
      int ok = (a[1] % p == 0) && (a[2] % p == 0) && (a[3] % p == 0) && a[4] == 0 && a[5] == 0 && a[6] == 0;
      if(!ok){ if(m >= 5) law_fail_ge5++; if(m < 64) failset[m] = 1; if(fll < 480) fll += snprintf(fl + fll, 512 - fll, " %ld", m); }
      if(m == 2 && ok) base2_inside++;
      if(m >= 6){
        nclasses++;
        if(!ok) continue;
        if(a[3] == 0) notcubic++;
        if(A[0] % p) continue;
        nzc++;
        u64 c0 = A[0] / p, c1 = a[1] / p, c2 = a[2] / p, c3 = a[3] / p;
        u64 i2 = (p + 1) / 2, i3 = pwm(3, p - 2, p); int roots = 0;
        for(u64 l = 0; l < p; l++){ u64 C2 = l * ((l + p - 1) % p) % p * i2 % p; u64 C3 = C2 * ((l + p - 2) % p) % p * i3 % p;
          if((c0 + c1 * l + c2 * C2 + c3 * C3) % p == 0) roots++; }
        long cls = m % (long)(p - 1); int typ = (cls == 1) ? 0 : (cls == 2) ? 1 : 2;
        hist[typ][roots < 4 ? roots : 4]++;
        if(cls == 2){ c2_zero++;
          u64 l = p - 1; /* l = -1 mod p */ u64 C2 = l * ((l + p - 1) % p) % p * i2 % p; u64 C3 = C2 * ((l + p - 2) % p) % p * i3 % p;
          if((c0 + c1 * l + c2 * C2 + c3 * C3) % p == 0) c2_forced++; }
      }
    }
    printf("p=%llu ord2=%llu: K=2 bases failing the cubic law:%s\n", (unsigned long long)p, (unsigned long long)t, fll ? fl : " none");
    free(b);
    if(p <= P3){ /* (5) mod p^3, Mahler speed 1/3 */
      u64 q3 = q * p; long N3 = (long)p + 8 + 10 * (long)(p - 1) + 2; u64 *b3 = b_upto(N3, q3, p);
      char f3[512]; int f3l = 0; f3[0] = 0; int minv[10]; for(int i = 0; i < 10; i++) minv[i] = 3;
      for(long m = 1; m <= (long)p + 8; m++){
        u64 tt[10]; for(int l = 0; l < 10; l++) tt[l] = b3[m + (long)l * (long)(p - 1)];
        int bad = 0;
        for(int i = 0; i < 10; i++){ int v = vpK(tt[0] % q3, p, 3); int bd = (i + 2) / 3 < 3 ? (i + 2) / 3 : 3; if(v < bd) bad = 1;
          if(m >= 8 && m <= 8 + (long)p - 2 && v < minv[i]) minv[i] = v;
          for(int l = 0; l < 9 - i; l++) tt[l] = (tt[l + 1] + q3 - tt[l]) % q3; }
        if(bad && f3l < 480) f3l += snprintf(f3 + f3l, 512 - f3l, " %ld", m);
      }
      printf("   K=3: bases failing v(a_i) >= min(3, ceil(i/3)):%s ; min over classes (m in [8, p+6]) of v(a_i), i=0..9:", f3l ? f3 : " none");
      for(int i = 0; i < 10; i++) printf(" %d", minv[i]); printf("\n");
      free(b3);
    }
    fflush(stdout);
  }
  printf("SUMMARY primes=%ld (gcd(ord_p 2, 6) = 1, 5 <= p <= %llu); classes (m in [6,p+4]) = %ld, not cubic (p^2 | a_3) = %ld\n", nprimes, (unsigned long long)PMAX, nclasses, notcubic);
  printf("  cubic-law failures at bases m >= 5: %ld; failing bases m < 64 seen anywhere:", law_fail_ge5); for(int m = 1; m < 64; m++) if(failset[m]) printf(" %d", m); printf("\n");
  printf("  class 1: b(p) = 2^((p-1)(p-2)/6) - 1 = 0 mod p failures: %ld\n", cls1_bad);
  printf("  zero classes %ld; root histogram [0,1,2,3,>3]: class1 %ld %ld %ld %ld %ld | class2 %ld %ld %ld %ld %ld | sporadic %ld %ld %ld %ld %ld\n", nzc,
         hist[0][0], hist[0][1], hist[0][2], hist[0][3], hist[0][4], hist[1][0], hist[1][1], hist[1][2], hist[1][3], hist[1][4], hist[2][0], hist[2][1], hist[2][2], hist[2][3], hist[2][4]);
  printf("  class 2: %ld zero classes, cubic vanishes at l = -1 (n = 2) in %ld; base m = 2 already obeys the cubic law at %ld of %ld primes\n",
         c2_zero, c2_forced, base2_inside, nprimes);
  return 0;
}

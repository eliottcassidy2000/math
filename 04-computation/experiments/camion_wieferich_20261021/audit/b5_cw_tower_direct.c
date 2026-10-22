/* Audit b5 (independent code): is p Camion-Wieferich, and what is its tower height h(p) = v_p(Lambda_p)?
 *
 * Direct computation, no interpolation polynomial of degree 2K-2:
 *   (1) s(n) mod p^3 for n <= p from the Hurwitz inverse W = T^(*-1) (u128 arithmetic)  ->  D(p) mod p^3, v_p(D(p));
 *   (2) s(n) mod p^2 for n <= 3p-2 (u64 arithmetic, lazy reduction)  ->  D(2p-1), D(3p-2) mod p^2;
 *   (3) the class-1 Camion function at base m = p evaluated at z = -p/(p-1) (the point n = 0) mod p^3:
 *         F_p(z) = sum_i alpha_i C(z,i),  v(alpha_i) >= ceil(i/2) (THM-4750 (1)/(4e)),  v(C(z,i)) >= 1 for 1 <= i < p,
 *       so mod p^3 only i <= 2 survive:  F_p(z) = D(p) + a_1 z + a_2 C(z,2)  (mod p^3),  a_1, a_2 needed only mod p^2.
 *       Lambda_p = eps * F_p(z) (sign irrelevant for valuations).  h(p) = v_p(Lambda_p) if < 3.
 * Also prints D(p) mod p^2 / p (the Camion quotient) for comparison with cw_fast.c.
 * Compile into scratch: cc -O2 -o b5 b5_cw_tower_direct.c ; usage: ./b5 p [p ...] */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128; typedef __int128 i128;
static u64 pw(u64 b, u64 e, u64 m){ u64 r = 1 % m; b %= m; while(e){ if(e & 1) r = (u64)((u128)r * b % m); b = (u64)((u128)b * b % m); e >>= 1; } return r; }
static u64 inv(u64 a, u64 m){ i128 t = 0, nt = 1, r = m, nr = a % m; while(nr){ i128 q = r / nr, x; x = t - q * nt; t = nt; nt = x; x = r - q * nr; r = nr; nr = x; } if(t < 0) t += m; return (u64)t; }
static int vp(u64 x, u64 p){ if(x == 0) return 99; int v = 0; while(x % p == 0){ x /= p; v++; } return v; }
/* s(n) mod q, n = 0..N, via the Hurwitz inverse; q < 2^62 */
static u64 *s_inverse(long N, u64 q){
  u64 *T = malloc((N + 1) * 8), *W = malloc((N + 1) * 8), *row = calloc(N + 2, 8), *ct = malloc((N + 1) * 8);
  for(long n = 0; n <= N; n++) T[n] = pw(2, (u64)n * (u64)(n ? n - 1 : 0) / 2, q);
  W[0] = 1; row[0] = 1;
  int small = (q < (1ULL << 30));
  for(long n = 1; n <= N; n++){
    for(long i = n; i >= 1; i--){ u64 t = row[i] + row[i - 1]; row[i] = t >= q ? t - q : t; }
    if(small){
      u64 acc = 0; /* terms < 2^60 after one reduction; reduce every 8 terms */
      for(long k = 1; k <= n; k++){ u64 c = row[k] * T[k] % q; acc += c * W[n - k]; if((k & 7) == 0) acc %= q; }
      acc %= q; W[n] = (q - acc) % q;
    } else {
      u128 acc = 0;
      for(long k = 1; k <= n; k++){ u64 c = (u64)((u128)row[k] * T[k] % q); acc += (u128)c * W[n - k]; if((k & 31) == 0) acc %= q; }
      acc %= q; W[n] = (q - (u64)acc) % q;
    }
  }
  u64 *s = malloc((N + 1) * 8); s[0] = 0; for(long n = 1; n <= N; n++) s[n] = (q - W[n]) % q;
  free(T); free(W); free(row); free(ct); return s;
}
static u64 cmod(long n, u64 q){ return (n & 1) ? pw(2, (u64)(n - 1) / 2, q) : 0; }
int main(int argc, char **argv){
  for(int ai = 1; ai < argc; ai++){
    u64 p = strtoull(argv[ai], 0, 10); u64 p2 = p * p, p3 = p2 * p;
    int eps = (p % 8 == 1 || p % 8 == 7) ? 1 : -1;
    u64 *s3 = s_inverse((long)p, p3);
    u64 Dp3 = (s3[p] + p3 - cmod((long)p, p3)) % p3;
    u64 *s2 = s_inverse(3 * (long)p - 2, p2);
    u64 A0 = Dp3 % p2;
    { long bad = 0; for(long n = 1; n <= (long)p; n++) if(s3[n] % p2 != s2[n]) bad++;
      printf("p=%llu: consistency of the u128 (mod p^3) and u64 (mod p^2) recurrences on n <= p: %s\n", (unsigned long long)p, bad ? "MISMATCH" : "ok"); }
    u64 D1 = (s2[2 * p - 1] + p2 - cmod((long)(2 * p - 1), p2)) % p2, D2 = (s2[3 * p - 2] + p2 - cmod((long)(3 * p - 2), p2)) % p2;
    u64 A1 = (eps == 1) ? D1 : (p2 - D1) % p2, A2 = D2;              /* A(l) = eps^(-l) D(p + l(p-1)) */
    u64 a1 = (A1 + p2 - A0) % p2, a2 = (A2 + 2 * p2 - (2 * A1) % p2 + A0) % p2;
    printf("p=%llu eps=%+d: D(p) mod p^3 = %llu, v_p(D(p)) = %d%s; Camion quotient D(p)/p mod p = %llu; p|a1: %s, p|a2: %s\n",
           (unsigned long long)p, eps, (unsigned long long)Dp3, vp(Dp3, p), vp(Dp3, p) >= 3 ? "+ (>= 3)" : "",
           (unsigned long long)((Dp3 % p2) / p % p), a1 % p ? "NO" : "yes", a2 % p ? "NO" : "yes");
    /* z = -p/(p-1) mod p^3 ; C(z,2) = z(z-1)/2 */
    u64 z = (p3 - (u64)((u128)p * inv(p - 1, p3) % p3)) % p3;
    u64 Cz2 = (u64)((u128)z * ((z + p3 - 1) % p3) % p3 * inv(2, p3) % p3);
    /* a1, a2 are known mod p^2 and multiply elements of p Z_p: lift them to residues mod p^3 (any lift is fine) */
    u64 F = (Dp3 + (u64)((u128)a1 * z % p3) + (u64)((u128)a2 * Cz2 % p3)) % p3;
    int h = vp(F, p);
    printf("   F_p(-p/(p-1)) mod p^3 = %llu  ->  v_p(Lambda_p) = %s%d  => h(p) %s %d\n", (unsigned long long)F,
           h >= 3 ? ">=" : "", h >= 3 ? 3 : h, h >= 3 ? ">=" : "=", h >= 3 ? 3 : h);
    if(vp(Dp3, p) >= 2){ u64 c1 = a1 / p % p, c2 = a2 / p % p, i2 = (p + 1) / 2; u64 lin = (c1 + p - c2 * i2 % p) % p;
      printf("   g_p(l) = c0 + c1 l + c2 C(l,2) with c0 = 0: c1 = %llu, c2 = %llu, c1 - c2/2 = %llu (root l=0 %s)\n",
             (unsigned long long)c1, (unsigned long long)c2, (unsigned long long)lin, lin ? "simple" : "DOUBLE"); }
    free(s3); free(s2); fflush(stdout);
  }
  return 0;
}

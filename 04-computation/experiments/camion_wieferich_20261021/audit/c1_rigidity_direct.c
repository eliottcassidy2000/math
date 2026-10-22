/* c1_rigidity_direct.c -- second-pass audit of THM-4750 Statement 5' (Teichmueller rigidity), DIRECT values, no interpolation.
 *
 * For a base a (any integer, reduced mod M = p^(kmax+1)) computes s_a(n) mod M for n <= p^kmax with the recurrence coming
 * from T' = S' * T^2 (EGF identity for T = 1/(1-S)); this is a different identity from the Wright / Pascal-row form
 * s(n) = T(n) - sum C(n,k) s(k) T(n-k) used by the session:
 *      T2(m)   = sum_j C(m,j) T(j) T(m-j)                     (Hurwitz square of T)
 *      s(m+1)  = T(m+1) - sum_{k=0}^{m-1} C(m,k) s(k+1) T2(m-k)
 * with T(n) = a^C(n,2).  It prints, for k = 1..kmax,  eps^k s_a(p^k) mod p^(k+1)  (which by the tower relation for s equals
 * sigma_p(a) mod p^(k+1)) and eps^k D_a(p^k) mod p^(k+1) (= Lambda_p(a) mod p^(k+1)), with D_a = s_a - a^((n-1)/2).
 *
 * Usage:  ./c1 p kmax lifts a0          -> bases a0, a0+p, a0+2p, a0+5p, a0-p+p^2, omega(a0), a0+p*R (R pseudo-random)
 *         ./c1 p kmax explicit b1 b2 ... -> the given bases (decimal, reduced mod M)
 *         ./c1 p kmax selftest           -> compare the recurrence with the Wright recurrence for n <= 60 (base 2 and 3)
 * Second-pass auditor, 2026-10-21.  Compile: cc -O2 -o c1 c1_rigidity_direct.c                                         */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static u64 M;                                  /* modulus p^(kmax+1) < 2^32 is enforced, so products of residues fit in u64 */
static u64 mulm(u64 x, u64 y){ return (x * y) % M; }
static u64 pw(u64 b, u64 e){ u64 r = 1 % M; b %= M; while(e){ if(e & 1) r = mulm(r, b); b = mulm(b, b); e >>= 1; } return r; }
/* s[1..N] for base a, by T' = S' T^2 */
static void strong_TSprime(u64 a, long N, u64 *s){
  u64 *T = malloc((N + 2) * 8), *T2 = malloc((N + 2) * 8), *row = calloc(N + 2, 8);
  u64 phi_hint = 0; (void)phi_hint;
  /* T(n) = a^C(n,2): iterate T(n+1) = T(n) * a^n */
  u64 an = 1; T[0] = 1 % M;
  for(long n = 0; n <= N; n++){ T[n + 1] = mulm(T[n], an); an = mulm(an, a % M); }
  s[0] = 0;
  row[0] = 1 % M;                              /* row = C(m, .) */
  for(long m = 0; m < N; m++){
    if(m >= 1){ for(long i = m; i >= 1; i--){ row[i] += row[i - 1]; if(row[i] >= M) row[i] -= M; } }
    u128 acc = 0;
    for(long j = 0; j <= m; j++) acc += (u128)mulm(row[j], T[j]) * T[m - j];
    T2[m] = (u64)(acc % M);
    acc = 0;
    for(long k = 0; k <= m - 1; k++) acc += (u128)mulm(row[k], s[k + 1]) * T2[m - k];
    s[m + 1] = (T[m + 1] + M - (u64)(acc % M)) % M;
  }
  free(T); free(T2); free(row);
}
/* plain Wright recurrence, for the self-test */
static void strong_wright(u64 a, long N, u64 *s){
  u64 *T = malloc((N + 2) * 8), *row = calloc(N + 2, 8); u64 an = 1; T[0] = 1 % M;
  for(long n = 0; n <= N; n++){ T[n + 1] = mulm(T[n], an); an = mulm(an, a % M); }
  row[0] = 1 % M; s[0] = 0;
  for(long n = 1; n <= N; n++){
    for(long i = n; i >= 1; i--){ row[i] += row[i - 1]; if(row[i] >= M) row[i] -= M; }   /* row = C(n, .) */
    u128 acc = 0; for(long k = 1; k < n; k++) acc += (u128)mulm(row[k], s[k]) * T[n - k];
    s[n] = (T[n] + M - (u64)(acc % M)) % M;
  }
  free(T); free(row);
}
static int vp(u64 x, u64 p, int cap){ int v = 0; if(x == 0) return cap; while(x % p == 0 && v < cap){ x /= p; v++; } return v; }
int main(int argc, char **argv){
  if(argc < 4){ fprintf(stderr, "usage\n"); return 1; }
  u64 p = strtoull(argv[1], 0, 10); int kmax = atoi(argv[2]);
  u64 pk[40]; pk[0] = 1; for(int i = 1; i <= kmax + 1; i++) pk[i] = pk[i - 1] * p;
  M = pk[kmax + 1]; if(M >= (1ULL << 32)){ fprintf(stderr, "modulus too large\n"); return 1; }
  long N = (long)pk[kmax];
  u64 bases[64]; int nb = 0;
  if(!strcmp(argv[3], "selftest")){
    long n0 = 60; u64 *s1 = calloc(n0 + 2, 8), *s2 = calloc(n0 + 2, 8); int bad = 0;
    for(u64 a = 2; a <= 5; a++){ strong_TSprime(a, n0, s1); strong_wright(a, n0, s2);
      for(long n = 1; n <= n0; n++) if(s1[n] != s2[n]) bad++; }
    printf("selftest p=%llu M=%llu: T'=S'T^2 recurrence vs Wright recurrence, bases 2..5, n<=60: %d mismatches\n",
           (unsigned long long)p, (unsigned long long)M, bad);
    printf("   s_2(1..8) mod M: "); strong_TSprime(2, 8, s1); for(int n = 1; n <= 8; n++) printf("%llu ", (unsigned long long)s1[n]); printf("\n");
    return 0;
  }
  if(!strcmp(argv[3], "lifts")){
    u64 a0 = strtoull(argv[4], 0, 10) % p;
    u64 om = a0; for(int i = 0; i < kmax + 4; i++) om = pw(om, p);           /* Teichmueller lift mod M */
    u64 R = 0x9E3779B97F4A7C15ULL % pk[kmax];
    u64 L[7] = { a0, a0 + p, a0 + 2 * p, a0 + 5 * p, (a0 + M - p + p * p) % M, om, (a0 + p * R) % M };
    for(int i = 0; i < 7; i++) bases[nb++] = L[i] % M;
  } else {
    for(int i = 4; i < argc && nb < 64; i++) bases[nb++] = strtoull(argv[i], 0, 10) % M;
  }
  u64 *s = calloc(N + 2, 8);
  u64 sig[64][40], lam[64][40];
  for(int b = 0; b < nb; b++){
    u64 a = bases[b];
    int eps = (pw(a, (p - 1) / 2) % p == 1) ? 1 : -1;       /* (a/p), per base */
    strong_TSprime(a, N, s);
    for(int k = 1; k <= kmax; k++){
      u64 mod = pk[k + 1], n = pk[k];
      u64 sv = s[n] % mod; u64 c = pw(a, (n - 1) / 2) % mod;
      u64 D = (sv + mod - c) % mod;
      int sg = (k % 2 == 0 || eps == 1) ? 1 : -1;          /* eps^k */
      sig[b][k] = sg == 1 ? sv : (mod - sv) % mod;
      lam[b][k] = sg == 1 ? D : (mod - D) % mod;
    }
    printf("p=%llu eps=%+d base a=%llu (mod p^%d):", (unsigned long long)p, eps, (unsigned long long)a, kmax + 1);
    for(int k = 1; k <= kmax; k++) printf("  k=%d: eps^k s(p^k)=%llu, eps^k D(p^k)=%llu [v=%d]", k,
           (unsigned long long)sig[b][k], (unsigned long long)lam[b][k], vp(lam[b][k], p, k + 1));
    printf("\n");
    /* tower relation for s and for D: eps^k X(p^k) = eps^(k-1) X(p^(k-1)) mod p^k */
    int tw = 1; for(int k = 2; k <= kmax; k++){ if((sig[b][k] - sig[b][k - 1] + pk[k + 1]) % pk[k] != 0) tw = 0; if((lam[b][k] + pk[k + 1] - lam[b][k - 1]) % pk[k] != 0) tw = 0; }
    printf("   tower relations (s and D) for k=2..%d: %s\n", kmax, tw ? "OK" : "FAIL");
  }
  /* rigidity: compare eps^k s(p^k) mod p^(k+1) across the bases (only meaningful when all bases share the residue mod p) */
  int same_res = 1; for(int b = 1; b < nb; b++) if(bases[b] % p != bases[0] % p) same_res = 0;
  if(same_res){
    for(int k = 1; k <= kmax; k++){
      int ok = 1; for(int b = 1; b < nb; b++) if(sig[b][k] != sig[0][k]) ok = 0;
      printf("RIGIDITY p=%llu residue %llu level k=%d (mod p^%d): sigma equal for all %d lifts: %s  (sigma = %llu)\n",
             (unsigned long long)p, (unsigned long long)(bases[0] % p), k, k + 1, nb, ok ? "YES" : "NO", (unsigned long long)sig[0][k]);
    }
    /* consistency of the levels: sigma at level k+1 reduces to sigma at level k */
    int cons = 1; for(int k = 2; k <= kmax; k++) if(sig[0][k] % pk[k] != sig[0][k - 1]) cons = 0;
    printf("   level consistency sigma_(k+1) = sigma_k mod p^(k+1): %s\n", cons ? "OK" : "FAIL");
  }
  free(s); return 0;
}

"""S(k) = sum_n C(2n,n)C(4n,2n)/((kn+1)64^n).

KEY REDUCTION (exact):  C(2n,n)C(4n,2n)/64^n = (1/4)_n (3/4)_n / (n!)^2.
Proof: C(2n,n)=4^n (1/2)_n/n!; C(4n,2n)=16^n (1/2)_{2n}/(2n)!; and duplication
(1/2)_{2n} = 4^n (1/4)_n(3/4)_n, (2n)! = 4^n n! (1/2)_n.
Hence  sum_n C(2n,n)C(4n,2n) x^n/64^n = 2F1(1/4,3/4;1;x)   (Ramanujan signature 4),
and with 1/(kn+1) = int_0^1 t^{kn} dt,
      S(k) = int_0^1 2F1(1/4,3/4;1;t^k) dt = (1/k) int_0^1 s^{1/k-1} 2F1(1/4,3/4;1;s) ds
           = 3F2(1/4, 3/4, 1/k;  1, 1+1/k;  1).
"""
from mpmath import mp, mpf, binomial, hyp2f1, hyp3f2, quad, pi, sqrt, log, atan, nstr

mp.dps = 60

# verify the binomial identity exactly for many n
from fractions import Fraction as Fr
from math import comb
def poch(a, n):
    r = Fr(1)
    for i in range(n):
        r *= (a + i)
    return r
ok = all(Fr(comb(2*n, n)*comb(4*n, 2*n), 64**n)
         == poch(Fr(1, 4), n)*poch(Fr(3, 4), n)/Fr(comb(n, n))**1 / Fr(1)
         * Fr(1, 1) / (Fr(1)*Fr(1))**1 / Fr(1)
         if False else
         Fr(comb(2*n, n)*comb(4*n, 2*n), 64**n) ==
         poch(Fr(1, 4), n)*poch(Fr(3, 4), n)/(Fr(1)*Fr(int(__import__('math').factorial(n))))**2
         for n in range(0, 40))
print("identity C(2n,n)C(4n,2n)/64^n == (1/4)_n(3/4)_n/(n!)^2 for n<40:", ok)

def S_series(k, N=200000):
    """direct series with exact rational terms -> mpf (slow but honest)"""
    tot = mpf(0)
    a = mpf(1)
    for n in range(N):
        tot += a/(k*n+1)
        a *= mpf((n+mpf(1)/4)*(n+mpf(3)/4))/mpf((n+1)**2)
    return tot

def S(k):
    return hyp3f2(mpf(1)/4, mpf(3)/4, mpf(1)/k, 1, 1+mpf(1)/k, 1)

print("\nk   S(k) [3F2]                                  check vs claimed")
claims = {
 1: 8*sqrt(2)/(3*pi),
 2: 4/pi*log(1+sqrt(2)),
 3: -(sqrt(3)*log(5-2*sqrt(6)) + 2*atan(sqrt(2)/5))/pi,
}
for k in range(1, 6):
    v = S(k)
    c = claims.get(k)
    tag = "" if c is None else ("MATCH" if abs(v-c) < mpf(10)**(-40) else f"MISMATCH {nstr(v-c,5)}")
    print(f"{k}   {nstr(v, 40)}   {tag}")
print("\nseries cross-check S(2):", nstr(S_series(2, 60000), 12), " vs 3F2:", nstr(S(2), 12))

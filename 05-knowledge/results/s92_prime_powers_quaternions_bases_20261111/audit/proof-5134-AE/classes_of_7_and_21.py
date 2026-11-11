# classes_of_7_and_21.py -- auditor proof-5134-AE.
# THM-5134 (E), first sentence, on concrete classes: the class of N = 7 and of N = 21 modulo m = 441 = 9 * 49 and
# m = 2205 = 9 * 5 * 49 (every prime of m to a power higher than v_p(n) for n = 1, 7, 11, 21), both sheets.
# Reads census_28_big.json (my census at X = 2^28 for the moduli 441, 2205).
# Printed: measured class factor of the single class, the first-order value Phi(g0) from the measured types,
# the limit L_m(N - c_sigma) and the product formula of (E).
import json, math
from fractions import Fraction as Fr

HERE = "C:/tmp/wt-s92/05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/proof-5134-AE/"


def primes_of(n):
    ps = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        ps.append(n)
    return ps


def c_d(m, d):
    v = Fr(1, d * d)
    for p in primes_of(m // d):
        v *= (1 - Fr(1, p * p))
    return v


def vp(n, p):
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n


def main():
    res = json.load(open(HERE + "census_28_big.json"))
    for ms, dat in res["data"].items():
        m = int(ms)
        divs = dat["divisors"]
        for N in (7, 21):
            for sh in range(2):
                e = N + 1 if sh == 0 else N          # M = N - c_sigma
                n = odd_part(e)
                off = dat["off"][sh]
                dia = dat["diag"][sh]
                T = 2 * sum(sum(r) for r in off) + sum(sum(r) for r in dia)
                r = e % m
                cnt = 2 * sum(off[t][r] for t in range(len(divs))) + sum(dia[t][r] for t in range(len(divs)))
                g0 = math.gcd(m, r) if r else m
                D = {d: 2 * sum(off[t]) for t, d in enumerate(divs)}
                phi = float(sum(d * Fr(D[d], T) for d in divs if g0 % d == 0))
                lim = sum(d * c_d(m, d) for d in divs if g0 % d == 0)
                prod = Fr(1)
                for p in primes_of(m):
                    v = vp(n, p)
                    assert v < vp(m, p)
                    prod *= (1 - Fr(1, p * p)) * sum(Fr(1, p ** i) for i in range(v + 1))
                assert prod == lim
                print("m = %4d  N = %2d  sheet %d  (n = %2d, g0 = %3d): class count %7d  measured factor %.4f (noise %.4f)  first-order %.4f  limit L_m = product of (E) = %.4f  [full series S = %.4f]"
                      % (m, N, sh, n, g0, cnt, m * cnt / T, math.sqrt(m / T), phi, float(lim), 8 / math.pi ** 2 * sum(1.0 / d for d in range(1, n + 1) if n % d == 0)))


if __name__ == "__main__":
    main()

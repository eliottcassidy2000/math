"""Exact digit arithmetic, ternary clocks, and Collatz hostile controls.

Only Python's standard library is used.  Run from the repository root:
  python 04-computation/experiments/ternary_digits_20260925.py
The output is deterministic; this file makes no primality claims about giants.
"""

from itertools import product
from math import isqrt
from pathlib import Path


def valuation(n, p):
    assert n != 0 and p > 1
    n = abs(n)
    a = 0
    while n % p == 0:
        n //= p
        a += 1
    return a


def trial_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


def odd_step(n):
    assert n > 0 and n % 2 == 1
    a = valuation(3 * n + 1, 2)
    return (3 * n + 1) >> a, a


def N(k):
    assert k >= 1
    return (10**k - 7) // 3


def rotate(n):
    s = str(n)
    return int(s[1:] + s[0])


def sieve(limit):
    a = bytearray(b"\x01") * limit
    a[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit - 1) + 1):
        if a[p]:
            a[p * p : limit : p] = b"\x00" * (((limit - 1 - p * p) // p) + 1)
    return a


def main():
    lines = []

    def report(s):
        lines.append(str(s))

    report("TERNARY DIGIT AUDIT 2026-09-25: exact integers, no random samples")
    report("Universe: seven prime candidates by exhaustive trial division; circular primes below 10^6 by sieve and independent alphabet enumeration; ternary quotient levels 1..8; valuation k,t=1..60; Collatz decimal k=2..300; composite rising controls m=2..100.")
    initial = [(k, N(k), trial_prime(N(k))) for k in range(2, 10)]
    assert all(p for _, _, p in initial[:-1]) and not initial[-1][2]
    assert N(9) == 17 * 19607843
    assert trial_prime(17) and trial_prime(19607843)
    report(f"N_k initial run: {initial}")
    report("Seven consecutive primes k=2..8; first composite k=9: 333333331=17*19607843.")
    assert pow(10, 16, 17) == 1
    assert all(pow(10, j, 17) != 1 for j in range(1, 16))
    assert pow(10, 9, 17) == 7
    assert all((N(k) % 17 == 0) == (k % 16 == 9) for k in range(1, 301))
    assert all(pow(2, j, 641) != 1 for j in range(1, 64))
    assert pow(2, 32, 641) == 640 and pow(2, 64, 641) == 1
    assert 2**32 + 1 == 641 * 6700417
    assert trial_prime(641) and trial_prime(6700417)
    assert 641 == 5 * 2**7 + 1 == 5**4 + 2**4
    report("ord_17(10)=16, 17|N_k iff k=9 mod16; ord_641(2)=64, F5=641*6700417.")
    a = [2*n**3 + 4*n**2 + n for n in range(1, 6)]
    assert a == [7, 34, 93, 196, 355]
    assert all(2*n**3+4*n**2+n == n*(2*(n+1)**2-1) for n in range(1, 301))
    primes12 = [p for p in range(2, 38) if trial_prime(p)]
    assert len(primes12) == 12 and sum(primes12) == 197
    report(f"Cubic values={a}; a_n=n*(2(n+1)^2-1), composite for n>=2.")
    report(f"First twelve primes={primes12}; sum=197, whereas a4=196=14^2.")

    prime = sieve(10**6)
    circular = []
    for n in range(2, 10**6):
        if not prime[n]:
            continue
        s = str(n)
        good = True
        for j in range(len(s)):
            x = int(s[j:] + s[:j])
            if not prime[x]:
                good = False
                break
        if good:
            circular.append(n)
    independent = {2, 3, 5, 7}
    for k in range(2, 7):
        for digits in product("1379", repeat=k):
            s = "".join(digits)
            rotations = [int(s[j:]+s[:j]) for j in range(k)]
            if all(trial_prime(x) for x in rotations):
                independent.add(int(s))
    assert set(circular) == independent and len(circular) == 55
    assert all(rotate(p) == p for p in (2, 3, 5, 7))
    assert all(str(p) != "1" * len(str(p)) for p in (2, 3, 5, 7))
    reps = sorted({min(int(str(n)[j:]+str(n)[:j]) for j in range(len(str(n)))) for n in circular})
    assert len(reps) == 19
    counts = {k: sum(len(str(n)) == k for n in circular) for k in range(1, 7)}
    report(f"Circular count below 10^6={len(circular)}; counts by length={counts}; orbit representatives={reps}")
    report("One-digit exception: 2,3,5,7 are prime rotation fixed points and are not repunits; the prime-fixed-point repunit theorem requires k>=2.")
    assert trial_prime(331) and not trial_prime(133) and 133 == 7*19
    report("Hostile: 331 is prime but rotation133=7*19; the 333...31 prime run is not a circular-prime run.")

    for k in range(1, 61):
        for t in range(1, 61):
            assert valuation(N(k+t)-N(k), 3) == 1 + valuation(t, 3)
            assert valuation(10**t-1, 3) == 2 + valuation(t, 3)
            for n in (0, 1, 3, 5, 17):
                S_t = (4**t*(3*n+1)-1)//3
                assert valuation(S_t-n, 3) == valuation(t, 3)
                A_t = 10**t*n + 10*(10**t-1)//9
                assert valuation(A_t-n, 3) == valuation(t, 3)
    report("Valuations verified: N_{k+t}-N_k has v3=1+v3(t); S^t(n)-n and A^t(n)-n have v3=v3(t), S=4n+1, A=10n+10.")
    maps = {}
    for r in range(1, 9):
        q = 3**r
        h = {}
        s, a = 1, 0
        for _ in range(q):
            assert s not in h
            h[s] = a
            s, a = (4*s+1) % q, (10*a+10) % q
        assert (s, a) == (1, 0) and len(set(h.values())) == q
        assert all(h[(4*n+1) % q] == (10*h[n]+10) % q for n in range(q))
        if r > 1:
            previous = maps[r-1]
            assert all(h[n] % (q//3) == previous[n % (q//3)] for n in range(q))
        maps[r] = h
    report("Compatible bijective conjugacies H_r(S^j(1))=A^j(0) verified for every residue at r=1..8.")
    report(f"H_2 on inputs0..8: {[maps[2][n] for n in range(9)]}; H_1(n)=n-1 mod3, so raw unit guard is not preserved.")

    for k in range(2, 7):
        M = 10**k-1
        v = 2+valuation(k, 3)
        for n in range(10**(k-1), min(10**k, 10**(k-1)+1000)):
            d = n//10**(k-1)
            rho = rotate(n)
            assert rho == 10*n-d*M
            assert (rho-10*n) % 3**v == 0
            assert (rho-n) % 9 == 0
        n = 10**(k-1)
        assert (rotate(n)-10*n) % 3**(v+1) != 0
    assert 337 % 27 != 373 % 27 and 373 % 27 == (10*337) % 27
    assert rotate(13) % 27 != (10*13) % 27
    report("Rotation: rho=10n-d(10^k-1); carry-free multiplication mod3^r iff r<=2+v3(k), for all k-digit strings. Always identity mod9; not identity mod27 (337->373). At length2 mod27 even multiplication needs carry (13->31).")

    decimal_samples = []
    for k in range(2, 301):
        n = N(k)
        y, a1 = odd_step(n)
        z, a2 = odd_step(y)
        assert a1 == 1 and y > n
        if k >= 4:
            assert z < n
        if k == 4:
            assert (y, z, a2) == (4997, 937, 4)
        if k >= 5:
            assert a2 == 3 and z == (3*10**k)//16-1
        if k <= 10:
            decimal_samples.append((k, n, y, z, a1, a2))
    report(f"Decimal two-step samples (k,n,U(n),U²(n),a1,a2): {decimal_samples}")
    report("For every k>=4, U²(N_k)<N_k; for k>=5, gap=(7*10^k-64)/48. Prime factorization is unused.")
    for m in range(2, 101):
        k = 2*m
        n = 2**k-1
        assert n % 3 == 0 and n > 3
        x = n
        for j in range(1, k):
            x, a = odd_step(x)
            assert a == 1 and x == 3**j*2**(k-j)-1 and x > n
    report("Composite hostile: all n=2^(2m)-1, m=2..100, have first2m-1 accelerated odd steps strictly growing with exponent1; formula proves arbitrary lengths.")
    report("ALL ASSERTIONS PASSED. Infinite claims are proved in the companion note; external giant primality certificates are CITED, not rerun here.")
    output = "\n".join(lines) + "\n"
    destination = Path(__file__).resolve().parents[2] / "05-knowledge/results/ternary_digits_20260925.out"
    destination.write_text(output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()

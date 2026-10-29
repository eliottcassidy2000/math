"""Audit of the post-launch THM-4845 (B) addition: self-converse tournaments (A002785) as a twisted Gaussian average.

(1) Direct twisted Burnside over all sigma in S_n (n <= 8): #{T : sigma.T = T^op} via orbits of
    tau(a,b) = (sigma b, sigma a) on ordered pairs (0 if an orbit contains a pair and its reverse, else 2^(#orbits/2)).
(2) The claimed formula: coefficient of x^n in E[(1 + W_1 x) exp(sum_{j = 2 mod 4} e^{X_j} x^j/j)], i.e.
    sum over cycle types with at most one fixed point and all other parts = 2 mod 4 of 2^((m^T G m - m_1)/2)/z_m.
"""
import itertools
from fractions import Fraction
from math import factorial, gcd
from audit_gf import partitions, zee


def direct(n):
    tot = 0
    pairs = [(a, b) for a in range(n) for b in range(n) if a != b]
    for s in itertools.permutations(range(n)):
        seen = set()
        norb = 0
        ok = True
        for p in pairs:
            if p in seen:
                continue
            orb = []
            q = p
            while q not in seen:
                seen.add(q)
                orb.append(q)
                q = (s[q[1]], s[q[0]])
            if any((b, a) in set(orb) for a, b in orb):
                ok = False
                break
            norb += 1
        if ok:
            assert norb % 2 == 0
            tot += 2 ** (norb // 2)
    assert tot % factorial(n) == 0
    return tot // factorial(n)


def formula(n):
    tot = Fraction(0)
    for m in partitions(n, [1] + [j for j in range(2, n + 1) if j % 4 == 2]):
        if m.get(1, 0) > 1:
            continue
        q = sum(gcd(a, b) * ka * kb for a, ka in m.items() for b, kb in m.items())
        e = Fraction(q - m.get(1, 0), 2)
        assert e.denominator == 1
        tot += Fraction(2 ** int(e), zee(m))
    assert tot.denominator == 1
    return int(tot)


if __name__ == '__main__':
    d = [direct(n) for n in range(1, 9)]
    f = [formula(n) for n in range(1, 17)]
    print('direct twisted Burnside n=1..8:', d)
    print('claimed Gaussian formula n=1..16:', f)
    print('agree on n<=8:', d == f[:8])

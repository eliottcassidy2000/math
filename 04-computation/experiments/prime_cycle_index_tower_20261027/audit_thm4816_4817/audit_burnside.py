"""Auditor's independent A(n) (tournament classes) and SC(n) (self-converse classes) by Burnside's lemma.

For every cycle type lambda of S_n a representative permutation sigma is built and the number of labelled
tournaments T with sigma(T) = T (resp. sigma(T) = T^op) is obtained by explicit constraint propagation on
ordered pairs (union-find with parity), not by a closed formula:
    automorphism:      x(sigma a, sigma b) = x(a, b)
    anti-automorphism: x(sigma a, sigma b) = x(b, a) = 1 - x(a, b)
with x(b, a) = 1 - x(a, b).  A(n) = sum_lambda fix(lambda)/z_lambda, SC(n) = sum_lambda antifix(lambda)/z_lambda.
"""
from fractions import Fraction
from math import factorial
import sys


def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield []
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield [k] + rest


def zee(part):
    z = 1
    from collections import Counter
    for l, m in Counter(part).items():
        z *= l ** m * factorial(m)
    return z


def perm_of(part):
    sigma, v = [], 0
    for l in part:
        for s in range(l):
            sigma.append(v + (s + 1) % l)
        v += l
    return sigma


def fixcount(sigma, anti):
    """number of labelled tournaments fixed by sigma (anti=False) or mapped to their converse (anti=True)"""
    n = len(sigma)
    # variables: unordered pairs {a<b}, value x(a,b) for a<b; x(b,a) = 1 - x(a,b)
    idx = {}
    for a in range(n):
        for b in range(a + 1, n):
            idx[(a, b)] = len(idx)
    parent = list(range(len(idx)))
    parity = [0] * len(idx)   # parity relative to parent

    def find(i):
        if parent[i] == i:
            return i, 0
        r, p = find(parent[i])
        parent[i] = r
        parity[i] ^= p
        return r, parity[i]

    def lit(a, b):  # returns (var, flip) with x(a,b) = var xor flip
        return (idx[(a, b)], 0) if a < b else (idx[(b, a)], 1)

    ok = True
    for a in range(n):
        for b in range(a + 1, n):
            v1, f1 = lit(a, b)
            v2, f2 = lit(sigma[a], sigma[b])
            # constraint: x(sigma a, sigma b) = x(a,b) xor anti
            want = f1 ^ f2 ^ (1 if anti else 0)
            r1, p1 = find(v1)
            r2, p2 = find(v2)
            if r1 == r2:
                if (p1 ^ p2) != want:
                    ok = False
            else:
                parent[r1] = r2
                parity[r1] = p1 ^ p2 ^ want
    if not ok:
        return 0
    comps = len({find(i)[0] for i in range(len(idx))})
    return 2 ** comps


def counts(n):
    A = Fraction(0)
    SC = Fraction(0)
    for part in partitions(n):
        sigma = perm_of(part)
        z = zee(part)
        A += Fraction(fixcount(sigma, False), z)
        SC += Fraction(fixcount(sigma, True), z)
    assert A.denominator == 1 and SC.denominator == 1
    return int(A), int(SC)


if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    for n in range(1, N + 1):
        a, sc = counts(n)
        print(f"n={n} A={a} SC={sc}", flush=True)

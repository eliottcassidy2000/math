#!/usr/bin/env python3
"""JC(2) golden-corner census, reproduced, and its metallic/Markov stratification.

The mac-mini-S137 reflection
(07-reflections/the-hurwitz-principle-jc2-golden-corner-lrc-ghost-convergents-
gbonacci-macmini-S137.md) reports: under the conservative citable battery
[Magnus gcd = 1 => invertible; m | n divisibility reduction; Moh max <= 100],
among all surviving degree pairs <= 300 the maximal Euclid-proxy chain
(length 9) is achieved by EXACTLY ONE pair, (178,288) = 2*(F_11,F_12), the
best phi-approximant in range -- "the reduction-resistant corner of the JC_2
landscape is golden, on the nose."  The census OUTPUT FILE referenced there
(jc2_golden_corner_macmini_S137.out) is absent from 05-knowledge/results/ and
no script exists, so this reproduces it from scratch.

It then tests the natural extension: does the corner stratify by the METALLIC
ratios (x^2 = nx+1, i.e. continued fraction [n;n,n,...]) / the Markov-Lagrange
spectrum sqrt5, sqrt8, sqrt(221)/5, ...?

SCOPE WARNING.  This is the Euclid PROXY, not the true Abhyankar-Moh polygon
calculus, and MISTAKE-237 explicitly quarantines the promotion of
"Lame-for-polygons" to an equivalent form of the JC(2) obstruction.  Nothing
here is a reduction; it is a frame statistic.
"""
from math import gcd, isqrt

LIMIT = 300
MOH = 100


def euclid_chain(a, b):
    """Number of division steps of the Euclidean algorithm (the proxy)."""
    n = 0
    while b:
        a, b = b, a % b
        n += 1
    return n


def cf(a, b):
    """Continued-fraction partial quotients of a/b."""
    q = []
    while b:
        q.append(a // b)
        a, b = b, a % b
    return q


def survives(m, n):
    """Conservative citable battery."""
    if gcd(m, n) == 1:
        return False                      # Magnus: gcd 1 => invertible
    if n % m == 0 or m % n == 0:
        return False                      # divisibility reduction
    if max(m, n) <= MOH:
        return False                      # Moh: degrees <= 100 settled
    return True


def fib(limit):
    f = [1, 1]
    while f[-1] < limit:
        f.append(f[-1] + f[-2])
    return f


def pell(limit):
    p = [1, 2]
    while p[-1] < limit:
        p.append(2 * p[-1] + p[-2])
    return p


def metallic_seq(nn, limit):
    """x_{k+1} = nn*x_k + x_{k-1};  ratio -> (nn+sqrt(nn^2+4))/2."""
    s = [1, nn]
    while s[-1] < limit:
        s.append(nn * s[-1] + s[-2])
    return s


if __name__ == "__main__":
    pairs = [(m, n) for m in range(1, LIMIT + 1) for n in range(m, LIMIT + 1)
             if survives(m, n)]
    print(f"battery: Magnus gcd=1, divisibility, Moh max<={MOH}; degrees <= {LIMIT}")
    print(f"surviving pairs: {len(pairs)}   (S137 reports 14,661)")

    scored = sorted(((euclid_chain(n, m), m, n) for m, n in pairs), reverse=True)
    best = scored[0][0]
    top = [(m, n) for c, m, n in scored if c == best]
    print(f"\nmaximal Euclid-proxy chain length: {best}")
    print(f"attained by {len(top)} pair(s): {top}")
    for m, n in top:
        g = gcd(m, n)
        print(f"   ({m},{n}) = {g}*({m//g},{n//g})   cf(n/m) = {cf(n, m)}"
              f"   n/m = {n/m:.9f}")
    print(f"   phi = {(1+5**0.5)/2:.9f}")

    print("\nchain-length distribution (top of the census):")
    from collections import Counter
    dist = Counter(c for c, _, _ in scored)
    for c in sorted(dist, reverse=True)[:6]:
        ex = [(m, n) for cc, m, n in scored if cc == c][:4]
        print(f"   length {c}: {dist[c]:6d} pair(s)   e.g. {ex}")

    # --- the metallic question -------------------------------------------
    print("\n--- metallic / Markov stratification test ---")
    F, P = fib(LIMIT), pell(LIMIT)
    print(f"Fibonacci <= {LIMIT}: {F}")
    print(f"Pell      <= {LIMIT}: {P}   (ratio -> 1+sqrt2, the silver ratio)")
    print("\nlongest chain achievable by each metallic family (any scaling):")
    for nn in (1, 2, 3, 4):
        seq = metallic_seq(nn, LIMIT)
        bestc, bestp = 0, None
        for i in range(len(seq) - 1):
            a, b = seq[i], seq[i + 1]
            if a == 0:
                continue
            for s in range(1, LIMIT // b + 1 if b else 1):
                m, n = s * a, s * b
                if n > LIMIT or not survives(m, n):
                    continue
                c = euclid_chain(n, m)
                if c > bestc:
                    bestc, bestp = c, (m, n)
        x = (nn + (nn * nn + 4) ** 0.5) / 2
        print(f"   metallic n={nn} (x={x:.6f}, cf=[{nn};{nn},{nn},...]): "
              f"best surviving chain {bestc} at {bestp}")
    print("\nLagrange spectrum (Hurwitz constants): sqrt5=%.6f, sqrt8=%.6f, "
          "sqrt(221)/5=%.6f" % (5 ** 0.5, 8 ** 0.5, (221 ** 0.5) / 5))
    print("""
READING.  Euclid chain length and Hurwitz approximability are DIFFERENT
extremal problems.  Chain length is maximised by all-ones partial quotients,
so the golden family is uniquely extremal and every other metallic family is
strictly WORSE (larger partial quotients terminate the algorithm faster).
The Lagrange/Markov spectrum sqrt5 < sqrt8 < ... instead stratifies
approximation quality, where the silver ratio IS the runner-up.  So the two
readings of the "golden corner" diverge below the top, and which secondary
corner is meaningful depends on which statistic the true Abhyankar-Moh
polygon calculus actually tracks -- the census proxy tracks the first.""")

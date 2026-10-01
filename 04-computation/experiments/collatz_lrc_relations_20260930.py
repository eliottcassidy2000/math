#!/usr/bin/env python3
"""collatz_lrc_relations_20260930.py -- Beck-Everett's lonely runner relations checked on small tight instances, and the
Collatz cycle half read as a per-period finite check (session collatz-posets-zeta5-20260927, opus, 2026-09-30, sixteenth note).

 (1) lon(n) = max_t min_j ||t n_j|| computed exactly for integer speeds (the maximum is attained at a breakpoint: a tent peak
     t = (2a+1)/(2 n_j) or a crossing t = a/(n_j +- n_l)); tight instances (lon = 1/(k+1)) among k-subsets of {1..N} for
     k = 2..5; for each, the shortest harmful relation (m . n = 0, sum m odd) in 1-norm and 2-norm, against Beck-Everett's
     bounds ||m||_1 <= 2k+3 and ||m||_2 <= 2(k+1)/sqrt(k-1) (Theorem 2 of arXiv:2609.06259).
 (2) The converse direction is not a theorem: sets with a short harmful relation that are nevertheless strictly lonely.
 (3) The Collatz cycle half as a per-period finite check: for shape (A, p) every cycle element is bounded by the carry
     maximum over the clock, so the cycle question at fixed period is a finite computation (as LRC at fixed k is, by Tao's
     finite-checking theorem); the bound tabulated for p <= 12 against the eleventh note's census.
Usage: python3 collatz_lrc_relations_20260930.py
"""
import itertools, math, time
from fractions import Fraction
from math import comb

T0 = time.time()


def dist_to_int(x):
    x = x - math.floor(x)
    return min(x, 1 - x)


def lon(n):
    """exact loneliness for integer speeds: evaluate min_j ||t n_j|| at all breakpoints of the piecewise-linear function."""
    cands = set()
    k = len(n)
    for j in range(k):
        for a in range(0, 2 * n[j]):
            cands.add(Fraction(2 * a + 1, 2 * n[j]))
        for l in range(j + 1, k):
            for s in (n[j] + n[l], abs(n[j] - n[l])):
                if s:
                    for a in range(0, s + 1):
                        cands.add(Fraction(a, s))
    best = Fraction(0)
    for t in cands:
        t = t - math.floor(t)
        m = min(min(Fraction(x * t - math.floor(x * t)), 1 - (x * t - math.floor(x * t))) for x in n)
        if m > best:
            best = m
    return best


def l1_ball(k, r):
    """integer vectors of dimension k with 1-norm <= r (recursive)."""
    if k == 0:
        yield ()
        return
    for x in range(-r, r + 1):
        for rest in l1_ball(k - 1, r - abs(x)):
            yield (x,) + rest


def shortest_harmful(n, bound1):
    """shortest (in 1-norm, then 2-norm) relation m . n = 0 with odd coordinate sum and ||m||_1 <= bound1; (m, l1, l2) or None."""
    best = None
    for m in l1_ball(len(n), bound1):
        l1 = sum(abs(x) for x in m)
        if l1 == 0 or sum(m) % 2 == 0 or sum(x * y for x, y in zip(m, n)) != 0:
            continue
        l2 = math.sqrt(sum(x * x for x in m))
        if best is None or (l1, l2) < (best[1], best[2]):
            best = (m, l1, l2)
    return best


def part0():
    print("== (0) Beck-Everett Lemma 3 identities and the LRC(14) numbers ==")
    import cmath
    for k in (3, 13):
        h = (k - 1) / (2 * (k + 1))
        # a_r = |psi_hat(r)|^2 with psi_h(x) = sqrt(2/h) cos(pi x/h) on |x| <= h/2 (numerical quadrature)
        M = 20000
        xs = [-h / 2 + h * (i + 0.5) / M for i in range(M)]
        def psi_hat(r):
            return sum(math.sqrt(2 / h) * math.cos(math.pi * x / h) * cmath.exp(-2j * math.pi * r * x) for x in xs) * (h / M)
        a = {r: abs(psi_hat(r)) ** 2 for r in range(-60, 61)}
        print(" k = %2d, h = %.4f: sum a_r = %.5f (should be 1), sum r^2 a_r = %.4f (should be 1/(4h^2) = %.4f)" % (k, h, sum(a.values()), sum(r * r * v for r, v in a.items()), 1 / (4 * h * h)))
    k = 13
    print(" LRC(14) (k = 13): a counterexample or tight instance has a harmful relation with ||m||_1 <= %d and ||m||_2 <= %.3f (sum m_i^2 <= %d), against THM-4009's ||a||_2 < 14 (sum <= 195), ||a||_1 <= 50 without the parity" % (2 * k + 3, 2 * (k + 1) / math.sqrt(k - 1), math.floor((2 * (k + 1)) ** 2 / (k - 1))))
    pairs_repo = [(a, b) for a in range(1, 14) for b in range(a + 1, 14) if math.gcd(a, b) == 1 and a * a + b * b < 196]
    pairs_be = [(a, b) for a in range(1, 9) for b in range(a + 1, 9) if math.gcd(a, b) == 1 and a * a + b * b <= 65 and (a + b) % 2 == 1]
    print(" support-two branch: coprime ratios a:b with a^2 + b^2 < 196: %d; with a^2 + b^2 <= 65 and a - b odd (harmful): %d %s" % (len(pairs_repo), len(pairs_be), pairs_be))


def part1():
    print("== (1) tight instances and their harmful relations (Beck-Everett Theorem 2) ==")
    for k, N in ((2, 24), (3, 18), (4, 14), (5, 11)):
        tight = []; viol = []
        for n in itertools.combinations(range(1, N + 1), k):
            L = lon(list(n))
            if L == Fraction(1, k + 1):
                rel = shortest_harmful(list(n), 2 * k + 3)
                tight.append((n, rel))
                if rel is None or rel[2] > 2 * (k + 1) / math.sqrt(k - 1) + 1e-9:
                    viol.append((n, rel))
            elif L < Fraction(1, k + 1):
                print("  COUNTEREXAMPLE?? n = %s lon = %s" % (n, L))
        print(" k = %d, speeds <= %d: tight instances %d; Beck-Everett bounds ||m||_1 <= %d, ||m||_2 <= %.3f; violations: %d (%.0fs)" % (
            k, N, len(tight), 2 * k + 3, 2 * (k + 1) / math.sqrt(k - 1), len(viol), time.time() - T0))
        for n, rel in tight[:5]:
            print("   n = %s: shortest harmful relation m = %s, ||m||_1 = %d, ||m||_2 = %.2f" % (n, rel[0] if rel else None, rel[1] if rel else -1, rel[2] if rel else -1))


def part2():
    print("== (2) the converse fails: strictly lonely sets with short harmful relations ==")
    k = 3; N = 12; have = 0; strict = 0
    for n in itertools.combinations(range(1, N + 1), k):
        rel = shortest_harmful(list(n), 2 * k + 3)
        if rel:
            have += 1
            if lon(list(n)) > Fraction(1, k + 1):
                strict += 1
    print(" k = 3, speeds <= 12: sets with a harmful relation of 1-norm <= 9: %d, of which strictly lonely (lon > 1/4): %d -> the relation is necessary for tightness, not sufficient" % (have, strict))


def part3():
    print("== (3) the Collatz cycle half is a per-period finite check ==")
    LOG23 = math.log2(3)
    rows = []
    for p in range(1, 13):
        A = math.floor(p * LOG23) + 1          # the least A with 2^A > 3^p: the positive clocks of period p
        D = 2 ** A - 3 ** p
        # the carry of any word of shape (A, p) is < 3^(p-1) * sum_i (2/3)^i * 2^(d_i)/... crude bound: S < 2^A * p
        Smax = 0
        # exact maximum carry over words of shape (A, p): the carry is maximised by putting the halvings late? enumerate small
        if comb(A - 1, p - 1) <= 200000:
            for cuts in itertools.combinations(range(1, A), p - 1):
                d = (0,) + cuts; S = sum(3 ** (p - 1 - i) * 2 ** d[i] for i in range(p))
                Smax = max(Smax, S)
        rows.append((p, A, D, Smax, (Smax // D) if Smax else None))
    print(" (period p, least A with positive clock, clock, max carry, max possible cycle element S_max/D):")
    for r in rows:
        print("  ", r)
    print(" every positive cycle of period p has elements <= S_max/D for its shape, so at fixed p the cycle question is a finite check (Steiner for one run, Simons-de Weger and Hercher for m runs); the divergence half has no such parameter")


if __name__ == "__main__":
    part0(); part1(); part2(); part3()
    print("total %.0fs" % (time.time() - T0))

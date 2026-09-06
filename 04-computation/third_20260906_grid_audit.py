"""Independent exact audit of translated-grid excess and inert pair bounds.

No producer imports; all requirements remain active under python -O.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm
import hashlib
import json
import sys

sys.stdout.reconfigure(newline="\n")
GATES = 0


def need(value, message):
    global GATES
    GATES += 1
    if not value:
        raise ArithmeticError(message)


def ceilq(x):
    return -(-x.numerator // x.denominator)


def frac(x):
    return x - x.numerator // x.denominator


def bad(x, u):
    a = (x.numerator * u) % x.denominator
    return 14 * min(a, x.denominator - a) < x.denominator


def phase_bank(U, t):
    # lambda=t*theta in [0,1). Every change is retained, including endpoints.
    cuts = {F(0), F(1)}
    for u in U:
        for k in range(u):
            for s in (-1, 1):
                cuts.add(frac(F(t * (14 * k + s), 14 * u)))
    cuts = sorted(cuts)
    return cuts[:-1] + [(a + b) / 2 for a, b in zip(cuts, cuts[1:])]


def grid_controls():
    rows = 0
    phase_rows = 0
    Ubank = [U for b in range(1, 5) for U in combinations(range(1, 8), b)]
    Ubank += [(1, 2, 3, 4, 5, 6, 7), (1, 1, 2, 2), (2, 4, 6, 8)]
    for U in Ubank:
        for t in range(1, 19):
            ds = [gcd(t, u) for u in U]
            B = sum(d * ceilq(F(t, 7 * d)) for d in ds)
            R0 = sum(ceilq(F(t, 7 * u)) - 1 for u in U[1:])
            for lam in phase_bank(U, t):
                points = [(lam + j) / t for j in range(t)]
                masks = [[bad(x, u) for x in points] for u in U]
                safe = sum(not any(mask[j] for mask in masks) for j in range(t))
                actual_sizes = sum(sum(mask) for mask in masks)
                actual_excess = actual_sizes - (t - safe)
                origin_layers = sum(
                    sum(bad(x, 1) and min(frac(x), 1 - frac(x)) < F(1, 14 * u)
                        for x in points)
                    for u in U[1:]
                )
                need(actual_sizes <= B, "strict single-mask upper bound")
                need(actual_excess >= origin_layers >= R0, "nested origin excess with strict endpoints")
                need(safe >= t - B + R0, "full translated-grid weak-safety count")
                phase_rows += 1
            rows += 1
    # Touching endpoints belong to the weak-safe set, not the strict bad set.
    need(not bad(F(1, 14), 1), "positive boundary is weak-safe")
    need(not bad(F(13, 14), 1), "negative boundary is weak-safe")
    need(ceilq(F(1)) - 1 == 0, "an open interval of one grid spacing may contain no point")
    return rows, phase_rows


def cap_controls():
    caps = (90, 30, 9, 4, 2, 1, 1)
    rows = []
    for tau in range(1, 7):
        weights = []
        for d in range(1, 91):
            if d % 7 == 0:
                continue
            n = sum(d <= c for c in caps)
            weights += [(d * ((-tau * pow(d, -1, 7)) % 7), d)] * n
        top = sorted(weights, reverse=True)[:7]
        ds = [d for _, d in top]
        E = F(sum(w for w, _ in top), 7)
        for k, c in enumerate(caps, 1):
            need(all(gcd(*S) <= c for S in combinations(ds, k)), "top divisors obey every subset cap")
        T = lcm(*ds)
        T *= (tau * pow(T, -1, 7)) % 7
        need(T % 7 == tau, "declared residue realization")
        actual = sum(d * ceilq(F(T, 7 * d)) for d in ds) - T
        need(actual == E, "abstract gcd-cap maximum is realized")
        rows.append((tau, int(E), ds))
    need([r[1] for r in rows] == [429, 426, 438, 435, 447, 445], "six exact E ceilings")
    seven_rows = []
    for tau in range(1, 7):
        weights = []
        for h in range(1, 13):
            if h % 7:
                weights += [h * ((-tau * pow(h, -1, 7)) % 7)] * sum(7 * h <= c for c in caps)
        seven_rows.append(sum(sorted(weights, reverse=True)[:3]))
    need(seven_rows == [151, 134, 152, 135, 146, 122], "full-v7=1 relaxed top-three maxima")
    need(7 * 6 == 42 and 7 ** 3 > 90, "higher-v7 residual ceilings")
    return rows, seven_rows


def allowed_sum(n):
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            if p % 3 != 2 or e > 2:
                return False
        p += 1
    return n == 1 or n % 3 == 2


def bernoulli2(x):
    x = frac(x)
    return x * x - x + F(1, 6)


def overlap(p, q):
    # Independent center-resonance computation, including contained intervals.
    s = p + q
    R = (s - 1) // 14
    return F(p + sum(min(2 * p, s - 14 * r) for r in range(1, R + 1)), 7 * p * q)


def literal_overlap(p, q):
    cuts = {F(0), F(1)}
    for u in (p, q):
        for k in range(u):
            cuts.add(frac(F(14 * k + 1, 14 * u)))
            cuts.add(frac(F(14 * k - 1, 14 * u)))
    cuts = sorted(cuts)
    flags = [bad((a + b) / 2, p) and bad((a + b) / 2, q)
             for a, b in zip(cuts, cuts[1:])]
    mu = sum((b - a for a, b, yes in zip(cuts, cuts[1:], flags) if yes), F(0))
    J = sum(yes and not flags[i - 1] for i, yes in enumerate(flags))
    return mu, J


def pair_controls():
    atlas = []
    literal_rows = 0
    for p in range(1, 178):
        for q in range(p + 1, 357 - p):
            if gcd(p, q) != 1 or not allowed_sum(p + q):
                continue
            mu = overlap(p, q)
            B = F(1, 49) + (bernoulli2(F(q - p, 14)) - bernoulli2(F(q + p, 14))) / (p * q)
            J = 2 * ((p + q - 1) // 14) + 1
            need(mu == B, "clipped center-resonance equals independent Bernoulli formula")
            need(mu >= F(1, 91), "parent's coarse uniform measure bound")
            need(J <= 51, "parent's coarse component bound")
            if q <= 50 or (p, q) in [(5, 348), (4, 349), (6, 347), (177, 179)]:
                need(literal_overlap(p, q) == (mu, J), "literal arc sweep checks measure and cyclic component count")
                literal_rows += 1
            atlas.append((F(447 + 30 * J) / mu, p, q, J, mu))
    need(len(atlas) == 5855, "complete actual inert ratio atlas")
    need(len({p + q for _, p, q, _, _ in atlas}) == 94, "independent 94 admissible-sum count")
    best = max(atlas)
    need(best == (F(6019965, 62), 5, 348, 51, F(62, 3045)), "exact universal pair cutoff optimizer")
    need(sum(row[0] == best[0] for row in atlas) == 1, "unique maximizing pair")
    need(best[0] < 97097 and 97096 < best[0], "correct integer strict threshold")
    need(min((mu, p, q) for _, p, q, _, mu in atlas) == (F(1, 70), 1, 10), "actual atlas minimum overlap")
    # Explicit correction control for the inherited microscopic addendum.
    need(overlap(1, 2) == F(1, 14) != F(3, 28), "containment clipping cannot be omitted")
    small = [(p, q) for p in range(1, 27) for q in range(p + 1, 27)
             if gcd(p, q) == 1 and p * q <= 26]
    need(len(small) == 36, "complete coarse-bound small-product universe")
    need(min(overlap(p, q) for p, q in small) == F(1, 91), "coarse bound closes every small-product exception")
    need(F(1, 49) - F(1, 4 * 27) > F(1, 91), "coarse Fourier envelope closes every larger product")
    window_cuts = sorted({F(0), F(1, 14)} | {
        F(14 * k + s, 196) for k in range(15) for s in (-1, 1)
        if 0 < F(14 * k + s, 196) < F(1, 14)
    })
    window = sum((b - a for a, b in zip(window_cuts, window_cuts[1:])
                  if bad((a + b) / 2, 1) and bad((a + b) / 2, 14)), F(0))
    need(window == F(1, 98) != overlap(1, 14) / 14,
         "windowed resonance pieces must be clipped; far-ratio bulk statement fails")
    phase_rows = 0
    for p in range(1, 8):
        for q in range(p + 1, 13):
            if gcd(p, q) != 1:
                continue
            mu = overlap(p, q)
            J = 2 * ((p + q - 1) // 14) + 1
            for h in range(1, 4):
                u, v = h * p, h * q
                for t in range(1, 13):
                    e = gcd(t, h)
                    for lam in phase_bank((u, v), t):
                        actual = sum(bad((lam + j) / t, u) and bad((lam + j) / t, v)
                                     for j in range(t))
                        need(actual >= t * mu - e * J, "translated pair count with exact gcd multiplicity")
                        phase_rows += 1
    # Endpoints prevent replacing the error J by zero even for positive measure.
    need(sum(bad(F(j, 7), 1) and bad(F(j, 7), 2) for j in range(7)) == 1,
         "literal positive overlap control")
    return atlas, literal_rows, phase_rows


def main():
    grid_rows, grid_phases = grid_controls()
    caps, seven = cap_controls()
    atlas, literal_rows, pair_phases = pair_controls()
    serial = [[str(a), p, q, j, str(mu)] for a, p, q, j, mu in atlas]
    digest = hashlib.sha256(json.dumps(serial, separators=(",", ":")).encode()).hexdigest()
    print("AUDIT SCOPE: sufficient translated-grid weak-safety count and balanced actual-entry scale bound; no LRC(14) closure.")
    print("Grid universe: all positive distinct U subsets of1..7 of size1..4, plus three repeated/seven-label controls; t1..18; every phase wall and open cell.")
    print("Grid rows:", grid_rows, "exact phase rows:", grid_phases)
    print("Residue-wise abstract gcd-cap maxima:", caps)
    print("Full-v7=1 relaxed ceilings:", seven, "higher ceilings42,0")
    print("Actual inert atlas:", len(atlas), "ratios;94 sums; literal arc sweeps:", literal_rows)
    print("Pair grid phase rows:", pair_phases)
    print("Exact optimizer:", max(atlas))
    print("Uniform sufficient integer cutoff: t>=97097. Actual atlas minimum overlap1/70 at(1,10).")
    print("Inherited microscopic-overlap correction: p1,q2 has measure1/14, not3/28; main Bernoulli formula passes all5855 rows.")
    print("Windowed correction: p1,q14 on[0,1/14) has overlap1/98, while full-circle overlap/14 is1/686.")
    print("atlas_semantic_sha256:", digest)
    print("Exact gates:", GATES)
    print("PASS")


if __name__ == "__main__":
    main()

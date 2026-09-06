"""Exact signed-pair virtual walls and a constructive primitive LRC(14) family.

Standalone standard-library verifier. No giant-height grid scan: family
controls use the virtual denominator 14d and literal physical inequalities.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce
from collections import Counter
import hashlib
import json
import sys

sys.stdout.reconfigure(newline="\n")
GATES = 0
Q = 91 ** 6


def need(test, label):
    global GATES
    GATES += 1
    if not test:
        raise RuntimeError(label)


def norm(x):
    x %= 1
    return min(x, 1 - x)


def mask(w, y):
    return tuple(j for j in range(3) if norm(w * (y + j) / 3) < F(1, 14))


def full_bad(T, y):
    return set().union(*(mask(w, y) for w in T)) == {0, 1, 2}


def signed_difference(u, v):
    sign = -1 if (u - v) % 3 == 0 else 1
    d = abs(u + sign * v) // 3
    need(d > 0 and (u + sign * v) % 3 == 0, "legal positive signed pair frequency")
    return sign, d


def wall_bad(c, d, numerator=1):
    modulus = 14 * d
    return {k for k in range(d)
            if min(c * (14 * k + numerator) % modulus,
                   -c * (14 * k + numerator) % modulus) < d}


def wall_count(c, d, numerator=1):
    g = gcd(c, d)
    q, h = d // g, c // g
    return g * ((q - 1 - numerator * h) // 14
                + (q - 1 + numerator * h) // 14 + 1)


def family(d, m):
    b, c = 3 * d + 1, 3 * d + 2
    L = 42 * d * b * c
    h = 1 + m * L
    C = tuple(sorted((d, 2 * d, 3 * d, 4 * d, 14, 14 * b, 14 * c, h, 2 * h, 1)))
    T = (1, b, c)
    k = (35 * d + 18) // 42
    y = F(14 * k + 1, 14 * d)
    return C, T, L, h, k, y


def check_family(d, m, scan=False, clocks=False):
    C, T, L, h, k, y = family(d, m)
    need(d >= 31 and gcd(d, 42) == 1 and m >= 1 and h > Q, "family domain")
    need(len(set(C)) == 10 and reduce(gcd, C) == 1, "primitive ten-pack")
    S = tuple(sorted(set([3 * z for z in C] + list(T))))
    need(len(S) == 13 and reduce(gcd, S) == 1, "actual primitive thirteen-speed row")
    need(all(w % 3 for w in T) and T[0] + T[1] == T[2], "additive ternary-unit tail")
    need(all((h - 1) % z == 0 for z in C if z not in (h, 2 * h)),
         "h preserves every base body residue")
    cross = min(max(x // gcd(x, z), z // gcd(x, z))
                for x in (h, 2 * h) for z in C if z not in (h, 2 * h))
    need(cross == h > Q, "exact necessary cross-height filter")
    need(max(T) < 11 * max(C), "necessary width-cone filter")
    need((h // gcd(h, 2 * h), 2 * h // gcd(h, 2 * h)) == (1, 2),
         "small primitive pair ratio")
    want_shift = F(5, 21 * d) if d % 6 == 1 else -F(2, 21 * d)
    need(y == F(5, 6) + want_shift and 0 <= k < d, "explicit virtual grid address")
    need(norm(d * y) == F(1, 14), "d is an actual body wall")
    need(all(norm(z * y) > F(1, 14) for z in C if z != d), "d is unique endpoint owner")
    need(all(norm(w * y) < F(3, 14) for w in T), "all three tails genuinely active")
    need(all(mask(w, y) == (2,) for w in T), "owner word is (2,2,2), not a permutation")
    x = y / 3
    need(min(norm(z * x) for z in S) == F(1, 14), "literal full physical witness")
    need(all(norm(z * (y + 1) / 3) >= F(1, 14) for z in S), "second free sheet")
    eta = min([F(1, 14 * d)]
              + [(norm(z * y) - F(1, 14)) / (2 * z) for z in C if z != d]
              + [(F(3, 14) - norm(w * y)) / (2 * w) for w in T])
    need(eta > 0, "positive addressed component cell")
    for point in (y, y + eta / 2, y + eta):
        need(all(norm(z * point) >= F(1, 14) for z in C), "literal body-safe cell")
        need(all(mask(w, point) == (2,) for w in T), "constant owner word on cell")
        need(all(norm(z * point / 3) >= F(1, 14) for z in S), "physical cell lift")
    need(min(norm(z * (y + eta / 2) / 3) for z in S) > F(1, 14),
         "strict full-row loneliness inside the addressed cell")
    for w in T:
        need(14 * w in C and wall_count(14 * w, w, 3) == w,
             "each actual positive tail event coset is wholly blocked")
        need(wall_count(14 * w, w, -3) == w, "negative tail coset also blocked")
    multiples = [z for z in C if z % d == 0]
    others = [z for z in C if z % d]
    need(multiples == [d, 2 * d, 3 * d, 4 * d], "four zero-blocking body multiples")
    need(len(others) == 6 and all(gcd(z, d) == 1 for z in others), "six coprime residual frequencies")
    need(6 * ((d + 6) // 7) < d, "phase-free sufficient virtual capacity")
    if scan:
        bad_sets = [wall_bad(z, d) for z in C]
        survivors = set(range(d)) - set().union(*bad_sets)
        need(k in survivors and len(survivors) > 0, "explicit witness versus entire virtual code")
        need(all(len(bad) == wall_count(z, d) for z, bad in zip(C, bad_sets)),
             "huge-height body residues retain exact virtual counts")
        # This is still a finite endpoint code, not a scan to max(C).
        for w in T:
            need(len(wall_bad(14 * w, w, 3)) == w, "literal actual-tail coset blocker")
    if clocks:
        for q in range(2, 15):
            need(any(z % q == 0 for z in S), "every denominator 2..14 is blocked")
        ten_gcds = [(P, reduce(gcd, P)) for P in combinations(S, 10)]
        nonprimitive = [(P, g) for P, g in ten_gcds if g > 1]
        need(nonprimitive == [(tuple(3 * z for z in C), 3)],
             "only original ten-pack has a nontrivial clock, exactly three")
        need(all(reduce(gcd, P) == 1 for P in combinations(S, 11)), "all eleven-subsets primitive")
        need(all(reduce(gcd, P) == 1 for P in combinations(S, 12)), "all twelve-subsets primitive")
    return {"d": d, "m": m, "h": h, "kappa": cross, "virtual_k": k,
            "y": str(y), "x": str(x), "eta": str(eta), "C": C, "T": T,
            "owners": [mask(w, y) for w in T]}


def main():
    # Exact point/cell checks of all pair bands in the finite tail universe.
    unit = [w for w in range(1, 15) if w % 3]
    tail_rows = 0
    for T in combinations(unit, 3):
        pairs = [signed_difference(u, v) for u, v in combinations(T, 2)]
        breaks = {F(0), F(1)}
        for w in T:
            breaks.update(F(14 * k + sign * 3, 14 * w) % 1
                          for k in range(w) for sign in (-1, 1))
        for _, d in pairs:
            breaks.update(F(21 * k + sign * 4, 21 * d) % 1
                          for k in range(d) for sign in (-1, 1))
        points = sorted(breaks)
        cells = points + [(a + b) / 2 for a, b in zip(points, points[1:])]
        for y in cells:
            masks = [mask(w, y) for w in T]
            need(all(len(B) <= 1 for B in masks), "literal ternary-unit single ownership")
            if set().union(*masks) == {0, 1, 2}:
                need(all(norm(d * y) > F(4, 21) for _, d in pairs),
                     "all signed virtual pair bands, including strict boundaries")
            for (u, v), (_, d) in zip(combinations(T, 2), pairs):
                if norm(d * y) <= F(4, 21):
                    A, B = mask(u, y), mask(v, y)
                    need(not A or not B or A == B, "pair collision or inactivity at virtual band")
        tail_rows += 1

    # Sharpness through a primitive additive sequence, without a limiting census.
    for N in range(2, 102, 2):
        T = (9 * N - 1, 33 * N - 1, 42 * N - 2)
        y = F(1, T[2])
        d = 8 * N
        need(reduce(gcd, T) == 1 and all(w % 3 for w in T), "primitive sharp additive sequence")
        need(full_bad(T, y), "actual owner permutation in sharpness control")
        need(norm(d * y) == F(4 * N, 21 * N - 1), "sharp pair-band approach")
        need(norm(d * y) - F(4, 21) == F(4, 21 * (21 * N - 1)), "exact sharpness error")
    need(all(norm(d * F(1, 4)) > F(4, 21) for d in (1, 2, 3))
         and not full_bad((1, 4, 5), F(1, 4)), "all pair bands are necessary, not sufficient")

    # New numerator-one code independently checked against exact modular distance.
    for d in range(1, 81):
        for c in list(range(1, 101)) + [14 * d, 28 * d, 13 * d]:
            actual = wall_bad(c, d)
            g, q = gcd(c, d), d // gcd(c, d)
            need(len(actual) == wall_count(c, d), "exact numerator-one body-wall count")
            need(len(actual) <= g * ((q + 6) // 7), "gcd-sensitive virtual count cap")
            need((len(actual) == d) == (c % (14 * d) == 0), "single-owner virtual cover boundary")
    for r in range(37, 44):
        need(6 * ((r + 6) // 7) < r, "seven-residue induction base for d>=37")
    need(6 * ((31 + 6) // 7) < 31, "coprime-to-42 exception d=31")
    need(all(gcd(d, 42) > 1 for d in range(32, 37)), "omitted small residues are outside domain")

    # Deleted whole components refute a false d-no-owner=>redundancy bridge.
    G14 = [(F(14 * k + 1, 196), F(14 * k + 13, 196)) for k in range(14)]
    surviving = [(a, b) for a, b in G14 if a >= F(1, 14) and b <= F(13, 14)]
    need(len(surviving) == 12 and surviving == G14[1:13], "two whole safe components deleted")
    need(all(norm(x) != F(1, 14) for I in surviving for x in I), "no surviving d=1-owned endpoint")
    need(norm(14 * F(1, 14)) < F(1, 14), "both d walls blocked by frequency14")
    padded = tuple(range(1, 10)) + (14,)
    need(all(norm(c * F(1, 28)) >= F(1, 14) for c in padded if c != 1)
         and norm(F(1, 28)) < F(1, 14), "primitive ten-pack deleted-component embedding")
    need(all(norm(c * F(1, 11)) >= F(1, 14) for c in padded), "same ten-pack remains body-safe")

    controls = []
    for d in range(31, 181):
        if gcd(d, 42) == 1:
            L = 42 * d * (3 * d + 1) * (3 * d + 2)
            m0 = Q // L + 1
            for m in (m0, m0 + 1):
                controls.append(check_family(d, m, scan=(d <= 43)))
    strong = []
    for r in range(3):
        d = 715 * (11 + 168 * r)
        L = 42 * d * (3 * d + 1) * (3 * d + 2)
        m = Q // L + 1
        strong.append(check_family(d, m, scan=(r == 0), clocks=True))
    payload = {"finite_tail_rows": tail_rows, "family_controls": controls, "strong": strong}
    print("STATUS: PASS; unconditional constructive family, LRC(14) remains open")
    print("FINITE TAIL UNIVERSE:", tail_rows, "ternary-unit triples, max<=14, every critical point and cell")
    print("PAIR BAND: full spoil implies ||dy||>4/21; primitive sharpness sequence verified")
    print("CODE UNIVERSE: d1..80,c1..100 plus full-cover boundaries; exact numerator-one counts")
    print("DELETION HOSTILE: C=(1,14) loses two whole components; no d-owned surviving endpoint")
    print("FAMILY CONTROLS:", len(controls), "plus", len(strong), "all-small-clock hostile controls")
    for row in (controls[0], strong[0]):
        print(json.dumps(row, sort_keys=True))
    print("ALL THREE ACTUAL TAIL EVENT COSETS BLOCKED; virtual endpoint owners=(2,2,2); two sheets free")
    print("ENTRY SCOPE: necessary numeric filters only; no THM3818 decoder-entry or W=Vdec claim")
    print("SEMANTIC SHA256:", hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest())
    print("ACTIVE GATES:", GATES)


if __name__ == "__main__":
    main()

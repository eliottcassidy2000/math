#!/usr/bin/env python3
"""Exact arithmetic audit of the disconnected-low affine-ray quotient.

This companion proves only the finite ray bookkeeping and the strict carrier
cutoff arithmetic.  The carrier inequalities and literal physical census are
separate parts of the THM-3355 candidate.
"""

from fractions import Fraction as F
from hashlib import sha256
from math import gcd


DMAX = F(186_636_088_362, 11_773_143_757_375)
TARGET = DMAX / 5
P_MIN = 264
P_MAX = 14_913
MAX_C = 46
PERTURBATION = 888
EXPECTED_PRIMITIVE = "cc7cc412a79baf2bd72fb3a58df0e3e0ef306de1372acc4fbfa945e9279cbb2d"
EXPECTED_SEMANTIC = "7a4ff17393137c5e3afcf24411d3d6cf2d0078255534c0b22b82110f167abcec"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def original_rays():
    rows = []
    for d in range(1, 9):
        for a in range(7 * d + 1):
            for c in tuple(range(-MAX_C, 0)) + tuple(range(1, MAX_C + 1)):
                if (a == 0 and c < 0) or (a == 7 * d and c > 0):
                    continue
                for p0 in range(1, d + 1):
                    if (a * p0 + c) % d == 0:
                        q0 = p0 + (a * p0 + c) // d
                        rows.append((d, a, c, p0, q0))
    require(len(rows) == 22_890, len(rows))
    return tuple(rows)


def primitive_row(row):
    d, a, c, p0, q0 = row
    scale = gcd(gcd(d, a), abs(c))
    dd, aa, cc = d // scale, a // scale, c // scale
    pp = next(x for x in range(1, dd + 1) if (aa * x + cc) % dd == 0)
    qq = pp + (aa * pp + cc) // dd
    require((p0 - pp) % dd == 0, (row, dd, pp))
    require(d * q0 - (d + a) * p0 == c, row)
    require(dd * qq - (dd + aa) * pp == cc, (dd, aa, cc, pp, qq))
    # The old ray is a residue subsequence of the divided ray.  Since a
    # solvable congruence has gcd(d,a)|c, division makes gcd(dd,aa)=1.
    require(gcd(dd, aa) == 1, (row, scale, dd, aa, cc))
    require(d // gcd(d, d + a) == dd, (row, dd))
    require((d + a) // gcd(d, d + a) == dd + aa, (row, dd + aa))
    require(c // gcd(d, d + a) == cc, (row, cc))
    return dd, aa, cc, pp, qq


def carrier_data(row):
    d, a, c, _, _ = row
    D = d + a
    require(gcd(d, D) == 1, row)
    if d + D >= 8:
        continuum_floor = F(1, 105)
        coefficient = F(12 * d, 35) + 2 * abs(c)
        chamber = "high"
    elif abs(c) >= 2:
        continuum_floor = F(120, 6811)
        coefficient = F(89 * d, 84) + 3 * abs(c)
        chamber = "low-nonresonant"
    else:
        continuum_floor = F(25, 1596)
        coefficient = F(89 * d, 84) + 3 * abs(c)
        chamber = "low-unit"
    require(continuum_floor > TARGET, (row, continuum_floor, TARGET))
    ratio = coefficient / (continuum_floor - TARGET)
    cutoff = ratio.numerator // ratio.denominator + 1
    require(continuum_floor - coefficient / cutoff > TARGET, (row, cutoff))
    if cutoff > 1:
        require(continuum_floor - coefficient / (cutoff - 1) <= TARGET, (row, cutoff))
    return chamber, continuum_floor, coefficient, cutoff


def main():
    original = original_rays()
    primitive = tuple(sorted({primitive_row(row) for row in original}))
    require(len(primitive) == 14_168, len(primitive))
    by_d = tuple(sum(row[0] == d for row in primitive) for d in range(1, 9))
    require(by_d == (644, 644, 1288, 1288, 2576, 1288, 3864, 2576), by_d)

    primitive_digest = sha256()
    maximum_cutoff = 0
    maximum_rows = []
    carrier_counts = {}
    relevant_occurrences = 0
    residual_occurrences = 0
    residual_by_d = [0] * 9
    for row in primitive:
        primitive_digest.update((repr(row) + "\n").encode())
        chamber, _, _, cutoff = carrier_data(row)
        carrier_counts[chamber] = carrier_counts.get(chamber, 0) + 1
        if cutoff > maximum_cutoff:
            maximum_cutoff, maximum_rows = cutoff, [row]
        elif cutoff == maximum_cutoff:
            maximum_rows.append(row)

        d, a, c, p0, q0 = row
        D = d + a
        n = max(0, (P_MIN - p0 + d - 1) // d)
        while True:
            p = p0 + d * n
            if p >= cutoff:
                break
            q = q0 + D * n
            if p < q < 8 * p and gcd(p, q) <= 3:
                relevant_occurrences += 1
                universal_A_lower = max(0, 168 * abs(c) - PERTURBATION)
                # The Dirichlet many-turn lemma is available only in its
                # short-rotation regime 9|c| <= p.  A large-c affine witness
                # may share this physical pair, but it cannot justify the
                # analytic skip.
                many_turn = (
                    9 * abs(c) <= p
                    and (p // d) * universal_A_lower >= 5 * 168 * p
                )
                if not many_turn:
                    residual_occurrences += 1
                    residual_by_d[d] += 1
            n += 1

    require(maximum_cutoff == P_MAX, maximum_cutoff)
    require(len(maximum_rows) == 56, len(maximum_rows))
    require(residual_occurrences == 8_079_264, residual_occurrences)
    require(F(1, 294) > TARGET, (F(1, 294), TARGET))
    require(5 * F(1, 294) > DMAX, (5 * F(1, 294), DMAX))

    primitive_hash = primitive_digest.hexdigest()
    semantic = sha256(
        repr(
            (
                len(original), len(primitive), by_d, carrier_counts,
                maximum_cutoff, tuple(maximum_rows), relevant_occurrences,
                residual_occurrences, tuple(residual_by_d[1:]),
            )
        ).encode()
    ).hexdigest()
    require(primitive_hash == EXPECTED_PRIMITIVE, primitive_hash)
    require(semantic == EXPECTED_SEMANTIC, semantic)
    print("LRC14 DISCONNECTED AFFINE-RAY QUOTIENT AUDIT")
    print("original_rays", len(original), "primitive_rays", len(primitive), "by_d", by_d)
    print("carrier_chambers", tuple(sorted(carrier_counts.items())))
    print("maximum_strict_cutoff", maximum_cutoff, "maximizing_rows", len(maximum_rows))
    print("carrier_relevant_occurrences", relevant_occurrences)
    print("post_universal_many_turn_occurrences", residual_occurrences, "by_d", tuple(residual_by_d[1:]))
    print("primitive_ray_sha256", primitive_hash)
    print("semantic_sha256", semantic)
    print("status=FINITE-EXACT ray quotient and cutoff arithmetic; carrier proof and physical census separate")


if __name__ == "__main__":
    main()

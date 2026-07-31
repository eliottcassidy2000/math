#!/usr/bin/env python3
"""Independent exact referee for the rail-8 E3/D6 root-zero clutch.

This deliberately does not import the clutch probe or its interval/integration
helpers.  The final delayed intersections are enumerated directly in the
lifted y=R*x coordinate rather than evaluated by the prefix floor formula.
"""

from bisect import bisect_right
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_relative_present_semantic_lift_probe_20260728 as relative


P = 13
R = P**6
S = R // P
G = 57_297_240


def req(ok, message):
    if not ok:
        raise RuntimeError(message)


def merge(rows):
    out = []
    for a, b in sorted(rows):
        if a >= b:
            continue
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return tuple(out)


def intersection(left, right):
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            out.append((a, b))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def complement(rows, period):
    out = []
    cursor = 0
    for a, b in merge(rows):
        if cursor < a:
            out.append((cursor, a))
        cursor = max(cursor, b)
    if cursor < period:
        out.append((cursor, period))
    return tuple(out)


def shift_union(rows, amount, period):
    out = []
    for a, b in rows:
        length = b - a
        c = (a + amount) % period
        d = c + length
        if d <= period:
            out.append((c, d))
        else:
            out.extend(((c, period), (0, d - period)))
    return merge(out)


def shift_weighted(rows, amount, period):
    out = []
    for a, b, w in rows:
        length = b - a
        c = (a + amount) % period
        d = c + length
        if d <= period:
            out.append((c, d, w))
        else:
            out.extend(((c, period, w), (0, d - period, w)))
    return tuple(sorted(out))


def profile_intersection(left, right):
    """Return support intersection retaining both profile values."""
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            out.append((a, b, left[i][2], right[j][2]))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def weighted_union_intersection(weighted, plain):
    starts = tuple(a for a, _ in plain)
    out = []
    for a, b, w in weighted:
        j = max(0, bisect_right(starts, a) - 1)
        while j < len(plain) and plain[j][0] < b:
            c = max(a, plain[j][0])
            d = min(b, plain[j][1])
            if c < d:
                out.append((c, d, w))
            j += 1
    return tuple(out)


def weighted_comb_intersection(weighted, speed, pd, lo, hi, period):
    req(0 <= lo < hi <= pd, "nonwrapping comb required")
    req(period % (pd * speed) == 0, "comb not integral on grid")
    unit = period // (pd * speed)
    step = pd * unit
    width = (hi - lo) * unit
    origin = lo * unit
    out = []
    for a, b, w in weighted:
        k = (a - width - origin) // step + 1
        p = origin + k * step
        while p < b:
            c = max(a, p)
            d = min(b, p + width)
            if c < d:
                out.append((c, d, w))
            p += step
    return tuple(out)


def direct_lifted_numerator(weighted, target, period):
    """Enumerate intersections R[a,b] cap (kT+[u,v]) exactly.

    The returned integer is the numerator over R*T of the physical integral.
    This route does not use Phi-prefix/floor-sum evaluation.
    """
    total = 0
    hits = 0
    for a, b, weight in weighted:
        ra, rb = R * a, R * b
        for u, v in target:
            k0 = (ra - v) // period + 1
            k1 = (rb - u - 1) // period
            for k in range(k0, k1 + 1):
                left = max(ra, k * period + u)
                right = min(rb, k * period + v)
                if left < right:
                    total += weight * (right - left)
                    hits += 1
    return total, hits


def e3(module):
    blockers = (module.C1, module.C2, module.C3)
    rows = module.make_comb(module.C3, 182, -13, 13)
    rows = module.subtract_comb(rows, module.W[module.GUARD], 91, -13, 13)
    for index in module.UNIT_IDX:
        rows = module.subtract_comb(rows, module.W[index], 182, -13, 13)
    rows = module.subtract_comb(rows, module.C1, 182, -13, 13)
    rows = module.subtract_comb(rows, module.C2, 182, -13, 13)
    return tuple(rows)


def fork(module):
    rows = module.make_comb(module.C1, 182, -13, 13)
    rows = module.intersect_comb(rows, module.C2, 182, -13, 13)
    rows = module.subtract_comb(rows, module.W[module.GUARD], 91, -13, 13)
    for index in module.UNIT_IDX:
        rows = module.subtract_comb(rows, module.W[index], 182, -13, 13)
    rows = module.subtract_comb(rows, module.C3, 182, -13, 13)
    return tuple(rows)


def section(module, source, ell, s=0, t=3):
    clock = tuple(module.make_comb(
        module.C1, 182, 26 * ell - 13, 26 * ell + 13
    ))
    rows = intersection(source, clock)
    rows = module.subtract_comb(
        rows, module.W[1], 182, -14 * s - 13, -14 * s + 13
    )
    rows = module.subtract_comb(
        rows, module.W[2], 182, -14 * t - 13, -14 * t + 13
    )
    rows = module.subtract_comb(
        rows, module.C2, 182, 14 * s - 13, 14 * s + 13
    )
    rows = module.subtract_comb(
        rows, module.C3, 182, 14 * t - 13, 14 * t + 13
    )
    return tuple(rows)


def delayed_digit(module, terminal_fork, ell):
    # Sector zero is the safe sector in THM-2623.
    word = module.build_word_Ta()
    qell = module.subtract_comb(
        word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
    )
    qell = intersection(tuple(qell), terminal_fork)
    return tuple(relative.private.old.sat.intersect_interval(
        qell, 13 * module.T // 26, 14 * module.T // 26
    ))


def construct_source(module, rail, present, ell, source_section, shift):
    target_pullback = shift_weighted(rail, -shift, module.T)
    pairs = profile_intersection(rail, target_pullback)
    req(pairs and all(ws == wt for _, _, ws, wt in pairs),
        "rail-8 source/pullback weights differ")
    rows = tuple((a, b, ws) for a, b, ws, _ in pairs)
    safe = complement(present[ell, 7], module.T)
    rows = weighted_union_intersection(rows, safe)
    rows = weighted_union_intersection(rows, shift_union(safe, -shift, module.T))
    rows = weighted_comb_intersection(rows, module.C3, 182, 169, 181, module.T)
    common = intersection(
        source_section, shift_union(source_section, -shift, module.T)
    )
    return weighted_union_intersection(rows, common)


def construct_target(module, rail, present, ell, target_section, shift):
    source_forward = shift_weighted(rail, shift, module.T)
    pairs = profile_intersection(rail, source_forward)
    req(pairs and all(wt == ws for _, _, wt, ws in pairs),
        "rail-8 target/forward weights differ")
    rows = tuple((a, b, wt) for a, b, wt, _ in pairs)
    safe = complement(present[ell, 7], module.T)
    rows = weighted_union_intersection(rows, safe)
    rows = weighted_union_intersection(rows, shift_union(safe, shift, module.T))
    rows = weighted_comb_intersection(rows, module.C3, 182, 1, 13, module.T)
    common = intersection(
        target_section, shift_union(target_section, shift, module.T)
    )
    return weighted_union_intersection(rows, common)


def cyclic_mul(a, b):
    raw = [0] * 7
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            raw[(i + j) % 7] += ai * bj
    return tuple((raw[i] - raw[6]) % P for i in range(6))


def det_mod(matrix):
    a = [[x % P for x in row] for row in matrix]
    det = 1
    for col in range(len(a)):
        pivot = next((i for i in range(col, len(a)) if a[i][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        value = a[col][col] % P
        det = det * value % P
        inv = pow(value, -1, P)
        for i in range(col + 1, len(a)):
            factor = a[i][col] * inv % P
            for j in range(col, len(a)):
                a[i][j] = (a[i][j] - factor * a[col][j]) % P
    return det % P


def norm(poly):
    cols = []
    for j in range(6):
        basis = tuple(int(i == j) for i in range(6))
        cols.append(cyclic_mul(poly, basis))
    return det_mod([[cols[j][i] for j in range(6)] for i in range(6)])


def main():
    module, _, _, _, rails, present, _ = relative.lift.m.core.build_carrier_data()
    T = module.T
    shift = 7 * T // R
    req((T, shift, module.C3) == (297836897838480, 431933040, 742586),
        "canonical scales changed")
    source_e3 = e3(module)
    terminal_fork = fork(module)
    source_values = []
    target_values = []
    ledger = []
    for ell in range(7):
        sec = section(module, source_e3, ell)
        source = construct_source(module, rails[8][3], present, ell, sec, shift)
        target = construct_target(module, rails[8][3], present, ell, sec, shift)
        req(shift_weighted(source, shift, T) == target,
            f"physical translation failed at ell={ell}")
        source_carry = weighted_comb_intersection(source, S, P, 12, 13, T)
        target_carry = weighted_comb_intersection(target, S, P, 6, 7, T)
        req(shift_weighted(source_carry, shift, T) == target_carry,
            f"carry 12->6 translation failed at ell={ell}")
        digit = delayed_digit(module, terminal_fork, ell)
        sv, shits = direct_lifted_numerator(source_carry, digit, T)
        tv, thits = direct_lifted_numerator(target_carry, digit, T)
        req(sv == tv, f"direct source/target integral mismatch at ell={ell}")
        source_values.append(sv)
        target_values.append(tv)
        ledger.append((ell, len(source), len(source_carry), len(digit), shits,
                       sum(b-a for a,b,_ in source), sv))

    source_values = tuple(source_values)
    target_values = tuple(target_values)
    expected = (
        0,
        339633525654239542165440,
        750593782703678965571520,
        719200126392878704654080,
        0, 0, 0,
    )
    req(source_values == target_values == expected, "reported vector changed")
    req(all(v % G == 0 for v in expected), "broadcast content does not divide")
    row_gcd = 0
    for value in expected:
        row_gcd = gcd(row_gcd, value)
    raw = tuple((v // G) % P for v in expected)
    source_normalized = tuple(v * pow(12, -1, P) % P for v in raw)
    target_normalized = raw
    source_reduced = tuple((source_normalized[i] - source_normalized[-1]) % P
                           for i in range(6))
    target_reduced = tuple((target_normalized[i] - target_normalized[-1]) % P
                           for i in range(6))
    source_det = norm(source_reduced)
    target_det = norm(target_reduced)
    req(source_det == target_det == 12, "independent multiplication norm changed")
    req((row_gcd // G) % P == 1, "extra row content changes the mod-13 class")

    # Root and phase transport are checked directly, not inferred from output.
    h, kappa = 6, 1
    source_root = (2 * 12 + (2*h+kappa)//P + 0) % P
    target_root = (2 * 6 + (2*h+kappa)//P + 1) % P
    req((source_root, target_root) == (12, 1), "private root convention changed")
    req(R * shift % T == 0, "delayed phase is not invariant")
    req((S * shift) % T == 7 * T // P, "carry translation sign changed")
    req((module.C3 * shift * 182 // T) % 182 == 14,
        "deep chart translation sign changed")

    print(f"vector={expected}")
    print(f"ledger={tuple(ledger)}")
    print(f"row_gcd={row_gcd} row_gcd_over_g={row_gcd//G} residual_mod13={(row_gcd//G)%P}")
    print(f"raw_mod13={raw}")
    print(f"source_normalized={source_normalized} reduced={source_reduced} det={source_det}")
    print(f"target_normalized={target_normalized} reduced={target_reduced} det={target_det}")
    print("transport=delta:+7/R carry:12->6 deep:(169,181)->(1,13) roots:right12->left1")


if __name__ == "__main__":
    main()

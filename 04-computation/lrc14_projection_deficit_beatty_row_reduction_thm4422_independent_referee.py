#!/usr/bin/env python3
"""Independent clean-room exact referee for THM-4422.

This standard-library program does not import the canonical verifier or any
projection scout.  It rebuilds generic carriers from the integer kernel,
reconstructs all three one-dimensional row compilers, and separately checks
the signed norm-four rays and their layer-cake formula.
"""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from json import dumps
from math import gcd
import argparse
import sys


TARGET = Q(6, 77)
AP_LEADER = Q(24, 343)
CHECKS = 0


class RefereeFailure(RuntimeError):
    pass


def need(condition, payload):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RefereeFailure(payload)


def strict_floor(numerator, denominator):
    """Largest integer strictly below numerator / denominator."""
    return (numerator - 1) // denominator


def floor_q(value):
    return value.numerator // value.denominator


def ceil_q(value):
    return -floor_q(-value)


def ceil_ratio(numerator, denominator):
    return (numerator + denominator - 1) // denominator


def ff(value):
    return f"{value.numerator}/{value.denominator}"


def ternary_units(height):
    return tuple(x for x in range(1, height + 1, 2) if x % 3)


def primitive(w):
    return gcd(gcd(w[0], w[1]), w[2]) == 1


def roofs(w, carrier):
    """Integer roof numerators p_i, indexed by the omitted coordinate i."""
    a, b, c = w
    x, y, z = carrier
    return (
        3 * (b + c) - 14 * abs(x),
        3 * (a + c) - 14 * abs(y),
        3 * (a + b) - 14 * abs(z),
    )


def exhaustive_carriers(w):
    """Generic exact kernel-box enumeration, independent of every short ray."""
    a, b, c = w
    bx = strict_floor(3 * (b + c), 14)
    by = strict_floor(3 * (a + c), 14)
    bz = strict_floor(3 * (a + b), 14)
    result = set()
    for x in range(-bx, bx + 1):
        for y in range(-by, by + 1):
            numerator = -(a * x + b * y)
            if numerator % c:
                continue
            z = numerator // c
            carrier = (x, y, z)
            if abs(z) <= bz and all(v % 3 for v in carrier) and all(p > 0 for p in roofs(w, carrier)):
                result.add(carrier)
    return result


def projection_sums(w, carriers):
    q = Q(3, 7 * w[2])
    result = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        result.append(sum(
            min(q, Q(roofs(w, carrier)[i], 14 * w[j] * w[k]))
            for carrier in carriers
        ))
    return tuple(result)


def normalized_deficits(w, carriers):
    """Return direct cap deficits and the three boundary-strip formulas."""
    a, b, c = w
    direct = [Q(0), Q(0), Q(0)]
    boundary = [Q(0), Q(0), Q(0)]
    for carrier in carriers:
        x, y, z = carrier
        p = roofs(w, carrier)
        direct[0] += max(Q(0), 1 - Q(c * p[0], 6 * b * c))
        direct[1] += max(Q(0), 1 - Q(c * p[1], 6 * a * c))
        direct[2] += max(Q(0), 1 - Q(c * p[2], 6 * a * b))
        boundary[0] += Q(max(0, 14 * abs(x) - 3 * (c - b)), 6 * b)
        boundary[1] += Q(max(0, 14 * abs(y) - 3 * (c - a)), 6 * a)
        boundary[2] += Q(max(0, 14 * c * abs(z) - 3 * c * (a + b) + 6 * a * b), 6 * a * b)
    return tuple(direct), tuple(boundary)


def extended_gcd(a, b):
    if b == 0:
        return a, 1, 0
    g, x, y = extended_gcd(b, a % b)
    return g, y, x - (a // b) * y


def progression_open(lower, upper, residue):
    """All integers t == residue (mod 3) in the strict interval."""
    if lower >= upper:
        return ()
    first_integer = floor_q(lower) + 1
    last_integer = ceil_q(upper) - 1
    first = first_integer + ((residue - first_integer) % 3)
    if first > last_integer:
        return ()
    return tuple(range(first, last_integer + 1, 3))


def compile_projection_rows(w, fixed_index):
    """Independently reconstruct every carrier in fixed-coordinate rows."""
    i = fixed_index
    j, k = (i + 1) % 3, (i + 2) % 3
    A, B, C = w[i], w[j], w[k]
    d = gcd(B, C)
    B0, C0 = B // d, C // d
    g, u, v = extended_gcd(B0, C0)
    need(g == 1 and B0 * u + C0 * v == 1,
         ("bezout", w, i, B0, C0, g, u, v))
    need(gcd(A, d) == 1 and d % 3 and B0 % 3 and C0 % 3,
         ("primitive-row-units", w, i, A, d, B0, C0))

    fixed_bound = strict_floor(3 * (B + C), 14 * d)
    compiled = set()
    weighted_sum = Q(0)
    max_multiplicity = 0
    max_cap = 0
    nonempty_rows = 0
    candidate_rows = 0
    for n in range(-fixed_bound, fixed_bound + 1):
        if n % 3 == 0:
            continue
        candidate_rows += 1
        residue_candidates = []
        for residue in range(3):
            trial = [0, 0, 0]
            trial[i] = d * n
            trial[j] = -A * n * u + C0 * residue
            trial[k] = -A * n * v - B0 * residue
            if all(value % 3 for value in trial):
                residue_candidates.append(residue)
        need(len(residue_candidates) == 1,
             ("unique-row-residue", w, i, n, residue_candidates))
        residue = residue_candidates[0]

        Aj = Q(3 * (A + C), 14)
        Ak = Q(3 * (A + B), 14)
        lower = max(Q(-Aj + A * n * u, C0), Q(-A * n * v - Ak, B0))
        upper = min(Q(Aj + A * n * u, C0), Q(-A * n * v + Ak, B0))
        ts = progression_open(lower, upper, residue)
        M = len(ts)
        if M:
            nonempty_rows += 1
        row_cap = ceil_ratio(d * (A + C), 7 * C)
        need(M <= row_cap, ("row-cap", w, i, n, M, row_cap, lower, upper))
        max_multiplicity = max(max_multiplicity, M)
        max_cap = max(max_cap, row_cap)

        for t in ts:
            carrier = [0, 0, 0]
            carrier[i] = d * n
            carrier[j] = -A * n * u + C0 * t
            carrier[k] = -A * n * v - B0 * t
            carrier = tuple(carrier)
            need(sum(w[r] * carrier[r] for r in range(3)) == 0,
                 ("compiled-relation", w, i, n, t, carrier))
            need(all(value % 3 for value in carrier) and all(p > 0 for p in roofs(w, carrier)),
                 ("compiled-live", w, i, n, t, carrier))
            need((B0 * carrier[j] - A * n) % 3 == 0
                 and (C0 * carrier[k] - A * n) % 3 == 0,
                 ("row-F3-identity", w, i, n, t, carrier))
            compiled.add(carrier)

        roof = Q(3 * (B + C) - 14 * d * abs(n), 14 * B * C)
        weighted_sum += M * min(Q(3, 7 * w[2]), roof)

    return (compiled, weighted_sum, max_multiplicity, max_cap,
            nonempty_rows, candidate_rows)


def norm_family(w):
    a, b, c = w
    labels = []
    if c == 2 * a + b:
        labels.append("2a+b")
    if c == 2 * b - a:
        labels.append("2b-a")
    if c == a + 2 * b:
        labels.append("a+2b")
    need(len(labels) <= 1, ("overlapping-norm-four-families", w, labels))
    return labels[0] if labels else None


def family_geometry(w, label):
    a, b, c = w
    if label == "2a+b":
        need(c == 2 * a + b, ("family", w, label))
        v = (2, 1, -1)
        selected = 0
        A, B = Q(3 * a, 14), Q(3 * (a + b), 14)
        radius_sum = Q(3 * (b + c), 28) + Q(3 * (a + c), 14)
        need(radius_sum == Q(3 * c, 7) < Q(c, 2),
             ("transverse-separation", w, label, radius_sum))
    elif label == "2b-a":
        need(c == 2 * b - a, ("family", w, label))
        v = (1, -2, 1)
        selected = 1
        A, B = Q(3 * (b - a), 14), Q(3 * b, 14)
        transverse_width = Q(3 * (b + c), 14) + Q(3 * (a + b), 14)
        need(transverse_width == Q(6 * b, 7) < b,
             ("transverse-width", w, label, transverse_width))
    elif label == "a+2b":
        need(c == a + 2 * b, ("family", w, label))
        v = (1, 2, -1)
        selected = 1
        A, B = Q(3 * b, 14), Q(3 * (a + b), 14)
        radius_sum = Q(3 * (b + c), 14) + Q(3 * (a + c), 28)
        need(radius_sum == Q(3 * c, 7) < Q(c, 2),
             ("transverse-separation", w, label, radius_sum))
    else:
        raise RefereeFailure(("unknown-family", w, label))
    need(A >= 0 and B > A and A + B == Q(3 * c, 14),
         ("AB-geometry", w, label, A, B))
    return v, selected, A, B


def predicted_ray(w, label):
    v, selected, A, B = family_geometry(w, label)
    K = strict_floor(B.numerator, B.denominator)
    carriers = {
        tuple(sign * k * coordinate for coordinate in v)
        for k in range(1, K + 1)
        if k % 3
        for sign in (-1, 1)
    }
    return carriers, v, selected, A, B, K


def R(t):
    """Positive integers <=t that are nonzero mod 3, for t>=0."""
    n = floor_q(t)
    return n - n // 3


def integrate_R(A, B):
    """Exact integration by independent constant-step partition."""
    points = {A, B}
    for integer in range(floor_q(A), ceil_q(B) + 1):
        point = Q(integer)
        if A < point < B:
            points.add(point)
    ordered = sorted(points)
    total = Q(0)
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        total += R(midpoint) * (right - left)
    return total


def family_rows(height):
    values = ternary_units(height)
    seen = set()
    for a in values:
        for b in values:
            if b <= a or gcd(a, b) != 1:
                continue
            for label, c in (
                ("2a+b", 2 * a + b),
                ("2b-a", 2 * b - a),
                ("a+2b", a + 2 * b),
            ):
                w = (a, b, c)
                if c <= height and c > b and c % 2 and c % 3:
                    need(w not in seen, ("duplicate-family-row", w, label))
                    seen.add(w)
                    yield w, label


def referee_generic(height):
    values = ternary_units(height)
    triples = 0
    sparse = 0
    dense = []
    dense_pairwise = 0
    dense_norm = []
    row_compilers = 0
    compiled_nonempty_rows = 0
    max_row_multiplicity = 0
    max_row_cap = 0
    first_binary_checks = 0
    candidate_rows_by_index = [0, 0, 0]
    max_m_by_index = [0, 0, 0]
    semantic = []
    for w in combinations(values, 3):
        if not primitive(w):
            continue
        triples += 1
        carriers = exhaustive_carriers(w)
        sums = projection_sums(w, carriers)
        direct_D, boundary_D = normalized_deficits(w, carriers)
        need(direct_D == boundary_D, ("boundary-deficit", w, direct_D, boundary_D))
        q = Q(3, 7 * w[2])
        need(all(sums[i] == q * (len(carriers) - direct_D[i]) for i in range(3)),
             ("deficit-identity", w, sums, direct_D, len(carriers)))
        need((min(sums) <= TARGET)
             == (max(direct_D) >= Q(len(carriers)) - Q(2 * w[2], 11)),
             ("deficit-equivalence", w, sums, direct_D, len(carriers)))

        for i in range(3):
            (compiled, compiled_sum, max_m, max_cap,
             nonempty, candidate_rows) = compile_projection_rows(w, i)
            need(compiled == carriers, ("row-compiler-set", w, i, compiled ^ carriers))
            need(compiled_sum == sums[i], ("row-compiler-sum", w, i, compiled_sum, sums[i]))
            row_compilers += 1
            compiled_nonempty_rows += nonempty
            max_row_multiplicity = max(max_row_multiplicity, max_m)
            max_row_cap = max(max_row_cap, max_cap)
            candidate_rows_by_index[i] += candidate_rows
            max_m_by_index[i] = max(max_m_by_index[i], max_m)
            if i == 0 and gcd(w[1], w[2]) == 1:
                need(max_m <= 1, ("binary-first-row", w, max_m))
                first_binary_checks += 1

        if Q(len(carriers)) <= Q(2 * w[2], 11):
            sparse += 1
        else:
            dense.append(w)
            if gcd(w[0], w[1]) == gcd(w[0], w[2]) == gcd(w[1], w[2]) == 1:
                dense_pairwise += 1
            label = norm_family(w)
            if label:
                dense_norm.append((w, label))
        semantic.append((w, tuple(sorted(carriers)), tuple(ff(x) for x in sums),
                         tuple(ff(x) for x in direct_D)))

    need(triples == 2910 and sparse == 2796 and len(dense) == 114,
         ("H79-counts", triples, sparse, len(dense)))
    need(candidate_rows_by_index[0] == 81792 and max_m_by_index[0] == 3,
         ("first-row-canonical-counts", candidate_rows_by_index, max_m_by_index))
    need(dense_pairwise == 114 and len(dense_norm) == 113,
         ("dense-classification-counts", dense_pairwise, len(dense_norm)))
    exceptions = [w for w in dense if norm_family(w) is None]
    need(exceptions == [(19, 23, 29)], ("dense-exception", exceptions))
    hostile = (19, 23, 29)
    hostile_carriers = exhaustive_carriers(hostile)
    expected_hostile = {
        tuple(sign * x for x in vector)
        for vector in ((1, 8, -7), (10, -7, -1), (11, 1, -8))
        for sign in (-1, 1)
    }
    need(hostile_carriers == expected_hostile,
         ("hostile-carriers", hostile_carriers ^ expected_hostile))
    hostile_sums = projection_sums(hostile, hostile_carriers)
    need(hostile_sums == (Q(156, 4669), Q(192, 3857), Q(3840, 88711)),
         ("hostile-sums", hostile_sums))
    digest = sha256(dumps(semantic, sort_keys=True).encode()).hexdigest()
    return {
        "triples": triples,
        "sparse": sparse,
        "dense": len(dense),
        "dense_pairwise": dense_pairwise,
        "dense_norm": len(dense_norm),
        "row_compilers": row_compilers,
        "nonempty_rows": compiled_nonempty_rows,
        "max_m": max_row_multiplicity,
        "max_cap": max_row_cap,
        "first_binary": first_binary_checks,
        "candidate_rows_by_index": tuple(candidate_rows_by_index),
        "max_m_by_index": tuple(max_m_by_index),
        "hostile_sums": hostile_sums,
        "sha": digest,
    }


def referee_norm_four(height):
    rows = 0
    family_counts = Counter()
    equality = []
    base35 = []
    base67 = []
    exhaustive_base_checks = 0
    large_generic_controls = []
    semantic = []
    for w, label in sorted(family_rows(height), key=lambda item: (item[0][2], item[0], item[1])):
        rows += 1
        family_counts[label] += 1
        carriers, vector, selected, A, B, K = predicted_ray(w, label)
        need(all(sum(w[i] * carrier[i] for i in range(3)) == 0
                 and all(x % 3 for x in carrier)
                 and all(p > 0 for p in roofs(w, carrier))
                 for carrier in carriers),
             ("predicted-ray-live", w, label))
        need(len(carriers) == 2 * (K - K // 3),
             ("ray-count", w, label, len(carriers), K))

        if w[2] < 67:
            direct = exhaustive_carriers(w)
            need(direct == carriers, ("base-exhaustion", w, label, direct ^ carriers))
            exhaustive_base_checks += 1
        elif len(large_generic_controls) < 12 and (sum(w) + 5 * w[0] + 7 * w[1]) % 97 == 0:
            direct = exhaustive_carriers(w)
            need(direct == carriers, ("large-exhaustion", w, label, direct ^ carriers))
            large_generic_controls.append((w, label))

        sums = projection_sums(w, carriers)
        selected_sum = sums[selected]
        q = Q(3, 7 * w[2])
        T = sum(
            min(Q(1), (B - k) / (B - A))
            for k in range(1, K + 1) if k % 3
        )
        need(selected_sum == 2 * q * T,
             ("selected-trapezoid", w, label, selected_sum, T, A, B))
        layer_integral = integrate_R(A, B)
        need(T == layer_integral / (B - A),
             ("layer-cake", w, label, T, layer_integral, A, B))
        upper_T = Q(w[2], 14) + Q(2, 3)
        need(T <= upper_T, ("R-row-upper", w, label, T, upper_T))
        selected_upper = Q(3, 49) + Q(4, 7 * w[2])
        need(selected_sum <= selected_upper,
             ("selected-upper", w, label, selected_sum, selected_upper))
        if w[2] >= 35:
            need(selected_upper <= TARGET,
                 ("large-height-target", w, label, selected_upper))
        else:
            base35.append((selected_sum, w, label))
            need(selected_sum <= TARGET,
                 ("small-height-target", w, label, selected_sum))
        if w[2] < 67:
            base67.append((selected_sum, w, label))
        else:
            need(selected_upper < AP_LEADER,
                 ("leader-tail", w, label, selected_upper, AP_LEADER))
        if selected_sum == TARGET:
            equality.append((w, label))
        semantic.append((w, label, vector, K, tuple(sorted(carriers)),
                         tuple(ff(x) for x in sums), ff(T), ff(selected_upper)))

    need(len(base35) == 25, ("base35-count", len(base35)))
    need(len(base67) == 105 and exhaustive_base_checks == 105,
         ("base67-count", len(base67), exhaustive_base_checks))
    leaders = sorted(base67, reverse=True)
    expected_leaders = [
        (TARGET, (1, 5, 11)),
        (Q(12, 161), (1, 11, 23)),
        (AP_LEADER, (1, 25, 49)),
    ]
    need([(value, w) for value, w, _ in leaders[:3]] == expected_leaders,
         ("base67-leaders", leaders[:3]))
    need(equality == [((1, 5, 11), "a+2b")],
         ("unique-equality", equality))
    need(len(large_generic_controls) == 12,
         ("large-control-count", len(large_generic_controls), large_generic_controls))
    digest = sha256(dumps(semantic, sort_keys=True).encode()).hexdigest()
    return {
        "height": height,
        "rows": rows,
        "family_counts": family_counts,
        "base35": len(base35),
        "base67": len(base67),
        "leaders": leaders[:3],
        "equality": equality,
        "large_controls": large_generic_controls,
        "sha": digest,
    }


def referee_no_fixed_average():
    w1, w2 = (1, 5, 7), (5, 11, 17)
    c1, c2 = exhaustive_carriers(w1), exhaustive_carriers(w2)
    expected1 = {tuple(sign * x for x in (2, 1, -1)) for sign in (-1, 1)}
    expected2 = {
        tuple(sign * k * x for x in (1, -2, 1))
        for k in (1, 2) for sign in (-1, 1)
    }
    need(c1 == expected1 and c2 == expected2,
         ("average-hostile-carriers", c1, c2))
    s1, s2 = projection_sums(w1, c1), projection_sums(w2, c2)
    need(s1 == (Q(8, 245), Q(6, 49), Q(4, 35)),
         ("average-hostile-sums-1", s1))
    need(s2 == (Q(122, 1309), Q(8, 119), Q(12, 119)),
         ("average-hostile-sums-2", s2))
    delta1 = tuple(value - TARGET for value in s1)
    delta2 = tuple(value - TARGET for value in s2)
    need(delta1 == (Q(-122, 2695), Q(24, 539), Q(2, 55)),
         ("average-delta-1", delta1))
    need(delta2 == (Q(20, 1309), Q(-14, 1309), Q(30, 1309)),
         ("average-delta-2", delta2))
    alpha1_floor = min(delta1[1], delta1[2]) / (
        -delta1[0] + min(delta1[1], delta1[2])
    )
    alpha2_floor = min(delta2[0], delta2[2]) / (
        -delta2[1] + min(delta2[0], delta2[2])
    )
    need(alpha1_floor == Q(49, 110) and alpha2_floor == Q(10, 17),
         ("average-alpha-floors", alpha1_floor, alpha2_floor))
    need(alpha1_floor + alpha2_floor == Q(1933, 1870) > 1,
         ("average-contradiction", alpha1_floor + alpha2_floor))
    return s1, s2, alpha1_floor, alpha2_floor


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--generic-height", type=int, default=79)
    parser.add_argument("--family-height", type=int, default=999)
    args = parser.parse_args()
    need(args.generic_height == 79, ("frozen-generic-height", args.generic_height))
    need(args.family_height >= 67, ("family-height", args.family_height))
    sys.stdout.reconfigure(newline="\n")

    generic = referee_generic(args.generic_height)
    norm = referee_norm_four(args.family_height)
    s1, s2, alpha1, alpha2 = referee_no_fixed_average()

    print("THM-4422 INDEPENDENT CLEAN-ROOM REFEREE")
    print("status=PASS; signed_norm_four=INDEPENDENTLY_AUDITED; LRC14=OPEN")
    print("generic_H79=triples:%d sparse:%d dense:%d pairwise_dense:%d norm4_dense:%d"
          % (generic["triples"], generic["sparse"], generic["dense"],
             generic["dense_pairwise"], generic["dense_norm"]))
    print("row_compilers=%d candidate_rows_by_index=%s nonempty_rows=%d max_m_by_index=%s max_cap=%d first_binary_checks=%d"
          % (generic["row_compilers"], generic["candidate_rows_by_index"],
             generic["nonempty_rows"], generic["max_m_by_index"],
             generic["max_cap"], generic["first_binary"]))
    print("dense_exception=(19,23,29); carriers=+/-{(1,8,-7),(10,-7,-1),(11,1,-8)}")
    print("dense_exception_sums=%s" % (tuple(ff(x) for x in generic["hostile_sums"]),))
    print("generic_semantic_sha256=" + generic["sha"])
    print("norm4_height=%d rows=%d counts=%s base_c_lt35=%d base_c_lt67=%d"
          % (norm["height"], norm["rows"], tuple(sorted(norm["family_counts"].items())),
             norm["base35"], norm["base67"]))
    print("carrier_rays=2a+b:(2,1,-1); 2b-a:(1,-2,1); a+2b:(1,2,-1)")
    print("layer_cake=T=(B-A)^-1*integral_A^B R(t)dt; R(t)<=2(t+1)/3")
    print("selected_bound=3/49+4/(7c)<=6/77_for_c>=35")
    print("base67_leaders=%s" % (tuple((w, ff(value), label) for value, w, label in norm["leaders"]),))
    print("unique_equality=%s" % (norm["equality"],))
    print("large_generic_controls=%s" % (tuple(norm["large_controls"]),))
    print("norm4_semantic_sha256=" + norm["sha"])
    print("fixed_average_hostile_sums=%s;%s"
          % (tuple(ff(x) for x in s1), tuple(ff(x) for x in s2)))
    print("fixed_average_forces=alpha1>=%s alpha2>=%s sum=%s>1"
          % (ff(alpha1), ff(alpha2), ff(alpha1 + alpha2)))
    print("checks=%d" % CHECKS)
    print("verdict=PASS")


if __name__ == "__main__":
    main()


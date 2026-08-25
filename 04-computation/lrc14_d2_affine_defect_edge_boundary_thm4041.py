#!/usr/bin/env python3
"""Independent exact audit of the THM-4024 d=2,c_2=9 boundary.

The companion compares three genuinely different views:

* literal rational wall-cell subdivision in the lift coordinate z;
* direct affine-centre intersection; and
* the one-edge defect lattice and its closed-form reduced-sum test.

All truth gates are explicit (no bare ``assert``), so ``python -O`` retains
the semantic audit.
"""

from fractions import Fraction as F
from itertools import combinations
from math import gcd

ONE = F(1)
THRESHOLD = F(1, 14)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def frac(x):
    return x % ONE


def norm(x):
    r = frac(x)
    return min(r, ONE - r)


def danger(w, x):
    return norm(w * x) < THRESHOLD


def lift_masks(alpha, beta, z):
    """Bad lift labels for the two odd exception speeds."""
    return (
        tuple(j for j in range(2) if danger(alpha, z + F(j, 2))),
        tuple(j for j in range(2) if danger(beta, z + F(j, 2))),
    )


def fully_spoiled(alpha, beta, z):
    masks = lift_masks(alpha, beta, z)
    return all(any(j in mask for mask in masks) for j in range(2))


def wall_points(alpha, beta):
    points = {F(0)}
    for w in (alpha, beta):
        for j in range(2):
            for n in range(w):
                centre = F(n, w) - F(j, 2)
                radius = F(1, 14 * w)
                points.add(frac(centre - radius))
                points.add(frac(centre + radius))
    return sorted(points)


def direct_wall_witness(alpha, beta):
    """Literal strict-wall-cell decision on z in R/Z."""
    walls = wall_points(alpha, beta)
    probes = list(walls)
    for i, left in enumerate(walls):
        right = walls[(i + 1) % len(walls)]
        if i + 1 == len(walls):
            right += ONE
        probes.append(frac((left + right) / 2))
    for z in probes:
        if fully_spoiled(alpha, beta, z):
            return z
    return None


def direct_wall_measure_and_components(alpha, beta):
    """Lebesgue measure and positive-length component count from wall cells."""
    walls = wall_points(alpha, beta)
    states = []
    lengths = []
    for i, left in enumerate(walls):
        right = walls[(i + 1) % len(walls)]
        if i + 1 == len(walls):
            right += ONE
        lengths.append(right - left)
        states.append(fully_spoiled(alpha, beta, frac((left + right) / 2)))
    measure = sum((length for length, state in zip(lengths, states) if state), F(0))
    if all(states):
        components = 1
    else:
        components = sum(
            state and not states[(i - 1) % len(states)]
            for i, state in enumerate(states)
        )
    return measure, components


def affine_witness(alpha, beta):
    """Find intersecting lifted components after fixing alpha on label zero."""
    radius = F(1, 14 * alpha) + F(1, 14 * beta)
    for A in range(alpha):
        c_alpha = F(A, alpha)
        # Simultaneous translation fixes A modulo alpha.  This deliberately
        # generous B box contains every centre within radius < 1 of c_alpha.
        for B in range(-beta, 2 * beta + 1):
            c_beta = F(B, beta) - F(1, 2)
            if abs(c_alpha - c_beta) < radius:
                N = 2 * beta * A - 2 * alpha * B + alpha * beta
                return A, B, N
    return None


def defect_witness(alpha, beta):
    """Enumerate the exact one-edge defect sidecars."""
    g = gcd(alpha, beta)
    bound = (alpha + beta - 1) // 7
    for N in range(-bound, bound + 1):
        if N == 0:
            continue
        if N % g == 0 and (N // g) % 2 != 0:
            return N
    return None


def reduced_sum_test(alpha, beta):
    g = gcd(alpha, beta)
    return alpha // g + beta // g > 7


def defect_phase_measure_and_components(alpha, beta):
    """Closed formula for the fully-spoiled phase set.

    Odd common dilation is a circle covering and preserves measure.  For the
    reduced coprime pair a<b, positive odd defects n index the two signs.
    Each fixed-orientation z interval has length

      min(1/(7b), (a+b-7n)/(14ab)).

    Adding the opposite label orientation gives the factor four after summing
    only positive defects.  The z-set has four components per positive defect;
    the quotient y=2z identifies half-turn partners and has two.  Pullback by
    the common odd gcd multiplies either component count by that gcd.  Haar
    measure is unchanged by the quotient.
    """
    g = gcd(alpha, beta)
    a, b = sorted((alpha // g, beta // g))
    total = F(0)
    positive_defects = 0
    for n in range(1, a + b, 2):
        if 7 * n >= a + b:
            break
        positive_defects += 1
        total += 4 * min(F(1, 7 * b), F(a + b - 7 * n, 14 * a * b))
    return total, 4 * g * positive_defects, 2 * g * positive_defects


def threshold_safe_phase(speeds):
    """Exact cell search for a phase with all clearances >= 1/14."""
    walls = {F(0)}
    for w in speeds:
        for n in range(w):
            centre = F(n, w)
            radius = F(1, 14 * w)
            walls.add(frac(centre - radius))
            walls.add(frac(centre + radius))
    walls = sorted(walls)
    probes = list(walls)
    for i, left in enumerate(walls):
        right = walls[(i + 1) % len(walls)]
        if i + 1 == len(walls):
            right += ONE
        probes.append(frac((left + right) / 2))
    for x in probes:
        clearance = min(norm(w * x) for w in speeds)
        if clearance >= THRESHOLD:
            return x, clearance
    return None


def audit_pair_universe(limit=79):
    odds = list(range(1, limit + 1, 2))
    profiles = 0
    positives = 0
    primitive = 0
    semantic_gates = 0
    first_positive = None
    first_nonprimitive_positive = None
    for alpha, beta in combinations(odds, 2):
        profiles += 1
        wall = direct_wall_witness(alpha, beta)
        affine = affine_witness(alpha, beta)
        defect = defect_witness(alpha, beta)
        closed = reduced_sum_test(alpha, beta)
        verdicts = (wall is not None, affine is not None, defect is not None, closed)
        semantic_gates += 4
        require(len(set(verdicts)) == 1,
                f"view mismatch {(alpha, beta)}: {verdicts}, {wall}, {affine}, {defect}")

        direct_measure = direct_wall_measure_and_components(alpha, beta)
        defect_measure = defect_phase_measure_and_components(alpha, beta)
        require(direct_measure == defect_measure[:2],
                f"phase-set formula mismatch {(alpha,beta)}: {direct_measure} != {defect_measure}")
        semantic_gates += 1

        g = gcd(alpha, beta)
        if g == 1:
            primitive += 1
        if closed:
            positives += 1
            if first_positive is None:
                first_positive = (alpha, beta, wall, affine, defect)
            if g > 1 and first_nonprimitive_positive is None:
                first_nonprimitive_positive = (alpha, beta, wall, affine, defect)

        # Every actual affine defect has exactly the claimed residue/gcd type.
        if affine is not None:
            A, B, N = affine
            require(N == 2 * beta * A - 2 * alpha * B + alpha * beta,
                    "defect definition drift")
            require(N % g == 0 and (N // g) % 2 != 0,
                    f"defect sidecar failed {(alpha, beta,A,B,N)}")
            require(7 * abs(N) < alpha + beta,
                    f"strict window failed {(alpha,beta,A,B,N)}")
            semantic_gates += 3

        # An attainable wall equality is impossible: after division by g,
        # even (odd+odd) cannot equal 7 times an odd integer.
        for N in range(-(alpha + beta), alpha + beta + 1):
            if N and N % g == 0 and (N // g) % 2 and 7 * abs(N) == alpha + beta:
                raise RuntimeError(f"unexpected attainable endpoint {(alpha,beta,N)}")
        semantic_gates += 1

    return {
        "limit": limit,
        "profiles": profiles,
        "primitive": primitive,
        "positives": positives,
        "semantic_gates": semantic_gates,
        "first_positive": first_positive,
        "first_nonprimitive_positive": first_nonprimitive_positive,
    }


def control_packet():
    # Exact geometry controls.
    require(direct_wall_witness(1, 5) is None, "(1,5) must be selector-safe")
    require(direct_wall_witness(3, 9) is None, "(3,9) must reduce to safe (1,3)")
    z35 = F(3, 16)
    require(fully_spoiled(3, 5, z35), "minimal-max hostile (3,5) failed")
    require(lift_masks(3, 5, z35) == ((1,), (0,)), "(3,5) labels drifted")
    z915 = F(1, 16)
    require(fully_spoiled(9, 15, z915), "odd dilation of (3,5) failed")
    require(lift_masks(9, 15, z915) == ((1,), (0,)), "(9,15) labels drifted")
    require(defect_phase_measure_and_components(3, 5) == (F(2, 105), 4, 2),
            "(3,5) phase-set formula drifted")
    require(defect_phase_measure_and_components(9, 15) == (F(2, 105), 12, 6),
            "odd dilation must preserve measure and triple components")

    # The requested H10 is a core.  The exact c_2=9 boundary needs one more
    # divided owner: nine divided body speeds plus two divided pair speeds.
    H10 = tuple(range(1, 11))
    K11 = H10 + (12,)
    y = F(1, 11)
    require(min(norm(w * y) for w in K11) == F(1, 11), "K11 phase is not safe")
    z = y / 2
    require(lift_masks(1, 11, z) == ((0,), (1,)), "typed selector masks drifted")
    require(fully_spoiled(1, 11, z), "typed selector phase was not spoiled")

    # Exact typed row: s=1,t=2,(p,q)=(1,12); nine even body owners and two
    # odd exceptions.  Dividing the pack gives K11=H10 union {12}.
    body = tuple(2 * h for h in range(2, 11)) + (1, 11)
    pair = (2, 24)
    row = tuple(sorted(body + pair))
    require(len(body) == 11 and len(row) == 13 and len(set(row)) == 13,
            "typed row cardinality/distinctness failed")
    require(gcd(*body) == 1, "typed body is not primitive")
    require(sum(u % 2 == 0 for u in body) == 9, "typed row is not c_2=9")
    safe = threshold_safe_phase(row)
    require(safe is not None, "typed selector hostile accidentally has no threshold-safe phase")
    strict_safe = (F(229, 560), F(5, 56))
    require(min(norm(w * strict_safe[0]) for w in row) == strict_safe[1],
            "typed strict-safe witness drifted")
    require(strict_safe[1] > THRESHOLD, "typed witness must be strictly safe")

    # A second literal H11 control uses the consecutive eleven-pack.
    H11 = tuple(range(1, 12))
    y11 = F(1, 12)
    require(min(norm(w * y11) for w in H11) == F(1, 12), "H11 phase is not safe")
    require(lift_masks(1, 11, y11 / 2) == ((0,), (1,)), "H11 masks drifted")

    return {
        "safe_pair": (1, 5),
        "safe_nonprimitive": (3, 9),
        "minimal_max_hostile": (3, 5, z35, lift_masks(3, 5, z35)),
        "dilated_hostile": (9, 15, z915, lift_masks(9, 15, z915)),
        "hostile_phase_set": defect_phase_measure_and_components(3, 5),
        "dilated_phase_set": defect_phase_measure_and_components(9, 15),
        "H10_core": H10,
        "K11": K11,
        "K11_y": y,
        "typed_exceptions": (1, 11),
        "typed_masks": lift_masks(1, 11, z),
        "typed_body": body,
        "typed_pair": pair,
        "typed_row_safe": safe,
        "typed_row_strict_safe": strict_safe,
        "H11_y": y11,
    }


def main():
    census = audit_pair_universe()
    controls = control_packet()
    print("THM4041 d=2,c2=9 affine defect-edge exact audit")
    print("status=PASS")
    print("universe=unordered distinct positive odd pairs <=", census["limit"])
    print("profiles=", census["profiles"])
    print("primitive_profiles=", census["primitive"])
    print("positive_profiles=", census["positives"])
    print("semantic_gates=", census["semantic_gates"])
    print("equivalence=literal_wall_cells=affine_centres=defect_sidecars=reduced_sum_gt_7")
    print("first_positive=", census["first_positive"])
    print("first_nonprimitive_positive=", census["first_nonprimitive_positive"])
    print("safe_pair=", controls["safe_pair"])
    print("safe_nonprimitive=", controls["safe_nonprimitive"])
    print("minimal_max_hostile=", controls["minimal_max_hostile"])
    print("dilated_hostile=", controls["dilated_hostile"])
    print("minimal_hostile_phase_measure_components=", controls["hostile_phase_set"])
    print("dilated_phase_measure_components=", controls["dilated_phase_set"])
    print("requested_H10_core=", controls["H10_core"])
    print("correct_typed_K11=", controls["K11"])
    print("K11_safe_phase=", controls["K11_y"])
    print("typed_exceptions=", controls["typed_exceptions"])
    print("typed_masks=", controls["typed_masks"])
    print("typed_body=", controls["typed_body"])
    print("typed_pair=", controls["typed_pair"])
    print("typed_row_safe_phase_clearance=", controls["typed_row_safe"])
    print("typed_row_strict_safe_phase_clearance=", controls["typed_row_strict_safe"])
    print("consecutive_H11_safe_phase=", controls["H11_y"])
    print("circuit_rank=0 (one defect edge; no nontrivial cycle equation)")
    print("scope=selector certificate only; pack-safe-phase intersection is forgotten")


if __name__ == "__main__":
    main()



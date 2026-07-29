#!/usr/bin/env python3
"""Exact body-projection screens for THM-2928's two-drift quotient.

For a literal six-body carrier

    G_F = T \\ union_{v in F} D_v,    F subset {1,...,14}, |F|=6,

let L=14*lcm(F), let J be the exact safe 1/L-cell word, and put
S_D=J mod D for D|L.  If the five-aligned residual is covered by two
distinct quotient combs D_{a_1},D_{a_2}, then

    (|S_D|/D) u_A <= mu(D_{a_1} union D_{a_2}) <= 25/91.

The proved five-comb floor u_A>=478/1365 therefore forces

    |S_D|/D <= 375/478.                                  (support screen)

If d_i=D/gcd(a_i,D)>1 and lcm(d_1,d_2)=D, then at fixed normalized
phase u the i-th comb covers at most

    (D/d_i) ceil(d_i/7)

cell residues.  This gives the independent denominator-cardinality screen

    |S_D|/D <= ceil(d_1/7)/d_1 + ceil(d_2/7)/d_2.

There is a sharper reflection screen when either denominator is 2.  The
involution r -> D-1-r preserves S_D and swaps parity, so each parity class
contains |S_D|/2 carrier residues.  If the other denominator is d, it must
therefore have capacity at least |S_D|/2:

    |S_D| <= 2(D/d)ceil(d/7).                           (parity screen)

Finally, when s=|S_D|/D>2535/3346, the pair overlap must be below 1/49.
Writing a_i=g*alpha_i with coprime alpha_1<alpha_2, LEM-043 gives

    alpha_2 <= 1/[7(s*478/1365-13/49)].

Together with a_1<(585/154)D this bounds g<=a_1 and hence a_2.  Thus these
high-support rows are finite in the actual quotient slopes, not merely in
projective pair shape.  The coprimality gcd(g,D)=1 follows separately from
gcd(a_i,D)=D/d_i and lcm(d_1,d_2)=D; it is not itself a size bound.

All arithmetic below is integral or Fraction-exact.  The selected-cell
projection is computed independently by a bitset union and by merging
cyclic integer intervals.  Explicit RuntimeError checks remain active under
python -O.
"""

from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from math import gcd, lcm


PAIR_UNION_CAP = Q(25, 91)
FIVE_SAFE_FLOOR = Q(478, 1365)
SUPPORT_CUTOFF = PAIR_UNION_CAP / FIVE_SAFE_FLOOR  # 375/478
SUBINDEPENDENCE_CUTOFF = Q(13, 49) / FIVE_SAFE_FLOOR  # 2535/3346


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def divisors(n):
    low = []
    high = []
    d = 1
    while d * d <= n:
        if n % d == 0:
            low.append(d)
            if d * d != n:
                high.append(n // d)
        d += 1
    return low + high[::-1]


def safe_cell_ranges(F):
    """Return L and half-open integer ranges comprising the exact cell word J."""
    L = lcm(*(14 * v for v in F))
    danger = []
    for v in F:
        half_tooth = L // (14 * v)
        period = L // v
        for k in range(v + 1):
            center = k * period
            danger.append(
                (max(0, center - half_tooth), min(L, center + half_tooth))
            )
    danger.sort()
    merged = []
    for left, right in danger:
        if merged and left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1] = (merged[-1][0], right)
        else:
            merged.append((left, right))
    safe = []
    cursor = 0
    for left, right in merged:
        if cursor < left:
            safe.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < L:
        safe.append((cursor, L))
    return L, safe


def support_size_bitset(D, ranges):
    """Compute |J mod D| by a cyclic bitset union."""
    mask = 0
    full = (1 << D) - 1
    for left, right in ranges:
        length = right - left
        if length >= D:
            return D
        residue = left % D
        if residue + length <= D:
            mask |= ((1 << length) - 1) << residue
        else:
            first = D - residue
            mask |= ((1 << first) - 1) << residue
            mask |= (1 << (length - first)) - 1
        if mask == full:
            return D
    return mask.bit_count()


def support_size_arcs(D, ranges):
    """Independent computation by splitting and merging cyclic integer arcs."""
    pieces = []
    for left, right in ranges:
        length = right - left
        if length >= D:
            return D
        residue = left % D
        endpoint = residue + length
        if endpoint <= D:
            pieces.append((residue, endpoint))
        else:
            pieces.append((residue, D))
            pieces.append((0, endpoint - D))
    pieces.sort()
    merged = []
    for left, right in pieces:
        if merged and left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1] = (merged[-1][0], right)
        else:
            merged.append((left, right))
    return sum(right - left for left, right in merged)


def denominator_pairs(D):
    ds = [d for d in divisors(D) if d > 1]
    return tuple(
        (d1, d2)
        for i, d1 in enumerate(ds)
        for d2 in ds[i:]
        if lcm(d1, d2) == D
    )


def numerator_floor(frac):
    return frac.numerator // frac.denominator


def main():
    support_hash = sha256()
    denominator_hash = sha256()
    parity_hash = sha256()
    high_support_hash = sha256()
    denominator_cache = {}

    body_count = 0
    total_rows = 0
    full_support_rows = 0
    support_killed_rows = 0
    support_hard_rows = 0
    support_hard_D = set()
    max_hard_divisors_per_body = 0
    minimum_hard_D = None
    minimum_hard_rows = []
    worst_support = (Q(1), None)

    denominator_raw_pairs = 0
    denominator_survivors = 0
    denominator_survivor_rows = 0
    denominator_survivor_D = set()
    max_denominator_survivors = (0, None)
    parity_candidates = 0
    parity_killed = 0
    parity_survivors = 0
    parity_survivor_rows = []
    survivors_after_parity = 0

    high_support_rows = 0
    low_support_rows = 0
    max_reduced_alpha2 = (0, None)
    max_actual_a2 = (0, None)

    control_support = None

    for F in combinations(range(1, 15), 6):
        body_count += 1
        L, ranges = safe_cell_ranges(F)
        safe_cells = sum(right - left for left, right in ranges)
        hard_this_body = 0
        control_rows = {}

        for D in divisors(L):
            total_rows += 1
            support_bits = support_size_bitset(D, ranges)
            support_arcs = support_size_arcs(D, ranges)
            require(
                support_bits == support_arcs,
                ("projection implementations disagree", F, L, D),
            )
            support = support_bits
            density = Q(support, D)
            support_hash.update(
                f"{','.join(map(str, F))}|{L}|{safe_cells}|{D}|{support}\n".encode()
            )
            if F == (1, 2, 3, 4, 5, 6):
                control_rows[D] = support
            if support == D:
                full_support_rows += 1
            if density < worst_support[0]:
                worst_support = (density, (F, L, D, support, Q(safe_cells, L)))

            if density > SUPPORT_CUTOFF:
                support_killed_rows += 1
                continue

            support_hard_rows += 1
            hard_this_body += 1
            support_hard_D.add(D)
            if minimum_hard_D is None or D < minimum_hard_D:
                minimum_hard_D = D
                minimum_hard_rows = [(F, L, D, support, density)]
            elif D == minimum_hard_D:
                minimum_hard_rows.append((F, L, D, support, density))

            if D not in denominator_cache:
                denominator_cache[D] = denominator_pairs(D)
            pairs = denominator_cache[D]
            denominator_raw_pairs += len(pairs)
            kept = []
            for d1, d2 in pairs:
                cap = Q((d1 + 6) // 7, d1) + Q((d2 + 6) // 7, d2)
                if density <= cap:
                    kept.append((d1, d2))
            denominator_survivors += len(kept)
            if kept:
                denominator_survivor_rows += 1
                denominator_survivor_D.add(D)
            if len(kept) > max_denominator_survivors[0]:
                max_denominator_survivors = (
                    len(kept),
                    (F, L, D, support, density),
                )
            denominator_hash.update(
                (
                    f"{','.join(map(str, F))}|{L}|{D}|{support}|"
                    + ";".join(f"{d1},{d2}" for d1, d2 in kept)
                    + "\n"
                ).encode()
            )
            kept_after_parity = []
            for d1, d2 in kept:
                if d1 == 2:
                    # Reflection r -> D-1-r swaps the two parity classes.
                    # Hence either parity contains exactly |S_D|/2 carrier
                    # residues, while the other comb has the displayed
                    # absolute residue capacity.
                    parity_candidates += 1
                    second_capacity = (D // d2) * ((d2 + 6) // 7)
                    if support > 2 * second_capacity:
                        parity_killed += 1
                        continue
                    parity_survivors += 1
                    parity_survivor_rows.append(
                        (F, L, D, support, density, d1, d2)
                    )
                kept_after_parity.append((d1, d2))
            survivors_after_parity += len(kept_after_parity)
            parity_hash.update(
                (
                    f"{','.join(map(str, F))}|{L}|{D}|{support}|"
                    + ";".join(f"{d1},{d2}" for d1, d2 in kept_after_parity)
                    + "\n"
                ).encode()
            )

            if density > SUBINDEPENDENCE_CUTOFF:
                high_support_rows += 1
                delta = density * FIVE_SAFE_FLOOR - Q(13, 49)
                require(delta > 0, ("nonpositive overlap defect", F, D))
                alpha2_cap = numerator_floor(Q(1, 1) / (7 * delta))
                a1_cap = (585 * D + 153) // 154 - 1
                a2_cap = a1_cap * alpha2_cap
                high_support_hash.update(
                    (
                        f"{','.join(map(str, F))}|{L}|{D}|{support}|"
                        f"{delta.numerator}/{delta.denominator}|"
                        f"{alpha2_cap}|{a1_cap}|{a2_cap}\n"
                    ).encode()
                )
                if alpha2_cap > max_reduced_alpha2[0]:
                    max_reduced_alpha2 = (
                        alpha2_cap,
                        (F, L, D, support, density, delta),
                    )
                if a2_cap > max_actual_a2[0]:
                    max_actual_a2 = (
                        a2_cap,
                        (F, L, D, support, density, delta, alpha2_cap, a1_cap),
                    )
            else:
                low_support_rows += 1

        max_hard_divisors_per_body = max(
            max_hard_divisors_per_body, hard_this_body
        )
        if F == (1, 2, 3, 4, 5, 6):
            require(Q(safe_cells, L) == Q(16, 35), "control body mass changed")
            control_support = control_rows

    require(body_count == 3003, "body universe changed")
    require(total_rows == 251536, "divisor-row universe changed")
    require(support_killed_rows == 240560, "support kill count changed")
    require(support_hard_rows == 10976, "support frontier count changed")
    require(full_support_rows == 223282, "full-support count changed")
    require(max_hard_divisors_per_body == 6, "per-body frontier changed")
    require(
        minimum_hard_rows
        == [
            (
                (1, 2, 3, 4, 6, 12),
                168,
                42,
                32,
                Q(16, 21),
            )
        ],
        "minimum support-hard row changed",
    )
    require(
        denominator_raw_pairs == 3066274,
        "denominator candidate count changed",
    )
    require(
        denominator_survivors == 23755,
        "denominator survivor count changed",
    )
    require(
        denominator_survivor_rows == 6292,
        "denominator survivor row count changed",
    )
    require(parity_candidates == 6756, "parity candidate count changed")
    require(parity_killed == 6754, "parity kill count changed")
    require(parity_survivors == 2, "parity survivor count changed")
    require(survivors_after_parity == 17001, "post-parity count changed")
    require(
        parity_survivor_rows
        == [
            (
                (1, 4, 5, 7, 9, 11),
                194040,
                194040,
                55392,
                Q(2308, 8085),
                2,
                194040,
            ),
            (
                (1, 5, 7, 8, 9, 11),
                388080,
                388080,
                109044,
                Q(3029, 10780),
                2,
                388080,
            ),
        ],
        "parity survivor rows changed",
    )
    require(high_support_rows == 1150, "high-support row count changed")
    require(low_support_rows == 9826, "low-support row count changed")
    require(
        control_support[280] == 210
        and control_support[420] == 282
        and control_support[840] == 384,
        "F=(1,...,6) support control changed",
    )

    print(f"support_cutoff={SUPPORT_CUTOFF}")
    print(f"subindependence_cutoff={SUBINDEPENDENCE_CUTOFF}")
    print(f"body_count={body_count}")
    print(f"divisor_rows={total_rows}")
    print(f"full_support_rows={full_support_rows}")
    print(f"support_killed_rows={support_killed_rows}")
    print(f"support_hard_rows={support_hard_rows}")
    print(f"unique_support_hard_D={len(support_hard_D)}")
    print(f"max_hard_divisors_per_body={max_hard_divisors_per_body}")
    print(f"minimum_support_hard_D={minimum_hard_D}")
    print(f"minimum_support_hard_rows={minimum_hard_rows}")
    print(f"worst_support={worst_support}")
    print(f"denominator_raw_pairs={denominator_raw_pairs}")
    print(f"denominator_survivors={denominator_survivors}")
    print(f"denominator_survivor_rows={denominator_survivor_rows}")
    print(f"unique_denominator_survivor_D={len(denominator_survivor_D)}")
    print(f"max_denominator_survivors={max_denominator_survivors}")
    print(f"parity_candidates={parity_candidates}")
    print(f"parity_killed={parity_killed}")
    print(f"parity_survivors={parity_survivors}")
    print(f"parity_survivor_rows={parity_survivor_rows}")
    print(f"survivors_after_parity={survivors_after_parity}")
    print(f"high_support_rows={high_support_rows}")
    print(f"low_support_rows={low_support_rows}")
    print(f"max_reduced_alpha2={max_reduced_alpha2}")
    print(f"max_actual_a2={max_actual_a2}")
    print(f"control_F_1_to_6_D280={control_support[280]}/280")
    print(f"control_F_1_to_6_D420={control_support[420]}/420")
    print(f"control_F_1_to_6_D840={control_support[840]}/840")
    print(f"support_ledger_sha256={support_hash.hexdigest()}")
    print(f"denominator_ledger_sha256={denominator_hash.hexdigest()}")
    print(f"parity_ledger_sha256={parity_hash.hexdigest()}")
    print(f"high_support_ledger_sha256={high_support_hash.hexdigest()}")


if __name__ == "__main__":
    main()

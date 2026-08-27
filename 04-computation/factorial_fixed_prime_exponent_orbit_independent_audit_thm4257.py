#!/usr/bin/env python3
"""Independent exact fixed-prime exponent-orbit atlas for odd-core suffix closure.

The central representation is the signed-five chart

    U(2^ell) = {(-1)^sigma 5^t : sigma in C2, t in C_(2^(ell-2))}

for ell >= 3 (with the evident C2 chart at ell=2).  Closure is evaluated
from explicit occupied output-position sets of a schoolbook digit multiplier;
the script never uses an integer product-and-mask test to decide membership.

Every load-bearing check uses require(), so python -O removes no gate.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from math import gcd


B_VALUES = (3, 5, 7, 9, 11, 13)
ELL_MIN = 1
ELL_MAX = 12


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fmt(values) -> str:
    return "{" + ",".join(str(x) for x in sorted(values)) + "}"


def divisors_power_two(n: int) -> tuple[int, ...]:
    require(n >= 1 and n & (n - 1) == 0, f"not a positive power of two: {n}")
    out = []
    d = 1
    while d <= n:
        out.append(d)
        d *= 2
    return tuple(out)


def signed_five_chart(ell: int):
    """Return n, encoding and decoding for C2 x C_n -> U(2^ell)."""
    require(ell >= 2, f"signed-five chart requires ell>=2, got {ell}")
    modulus = 1 << ell
    n = 1 if ell == 2 else 1 << (ell - 2)
    encode = {}
    decode = {}
    value = 1
    for t in range(n):
        for sigma, residue in ((0, value), (1, modulus - value)):
            coordinate = (sigma, t)
            require(coordinate not in encode, f"duplicate coordinate at ell={ell}")
            require(residue not in decode, f"signed-five collision ell={ell}, r={residue}")
            encode[coordinate] = residue
            decode[residue] = coordinate
        value = (5 * value) % modulus
    expected_units = tuple(range(1, modulus, 2))
    require(tuple(sorted(decode)) == expected_units, f"chart does not cover U(2^{ell})")
    return n, encode, decode


def coordinate_sum(left, right, n: int):
    return ((left[0] + right[0]) % 2, (left[1] + right[1]) % n)


def coordinate_multiple(coordinate, k: int, n: int):
    return ((coordinate[0] * k) % 2, (coordinate[1] * k) % n)


def coordinate_order(coordinate, n: int) -> int:
    state = (0, 0)
    for order in range(1, 2 * n + 1):
        state = coordinate_sum(state, coordinate, n)
        if state == (0, 0):
            return order
    raise RuntimeError(f"coordinate failed to return: coordinate={coordinate}, n={n}")


def binary_positions(value: int, ell: int) -> frozenset[int]:
    require(0 <= value < (1 << ell), f"word outside ell-bit range: {value}, ell={ell}")
    positions = set()
    remaining = value
    for position in range(ell):
        remaining, digit = divmod(remaining, 2)
        if digit:
            positions.add(position)
    require(remaining == 0, f"binary extraction left a tail for {value}")
    return frozenset(positions)


def occupied_product_positions(value: int, multiplier: int, ell: int) -> frozenset[int]:
    """Low ell output positions of multiplier*value via schoolbook carry."""
    input_positions = binary_positions(value, ell)
    output_positions = set()
    carry = 0
    for position in range(ell):
        input_digit = 1 if position in input_positions else 0
        carry, output_digit = divmod(carry + multiplier * input_digit, 2)
        if output_digit:
            output_positions.add(position)
    return frozenset(output_positions)


def closes_odd_core(value: int, b: int, ell: int) -> bool:
    require(b >= 3 and b % 2 == 1, f"invalid odd core {b}")
    product_positions = {
        multiplier: occupied_product_positions(value, multiplier, ell)
        for multiplier in range(1, b)
    }
    for s in range(1, (b - 1) // 2 + 1):
        if product_positions[s].isdisjoint(product_positions[b - s]):
            return False
    return True


def least_period(classes: frozenset[int], order: int) -> int:
    for period in divisors_power_two(order):
        if all(((k + period) % order in classes) == (k in classes) for k in range(order)):
            return period
    raise RuntimeError(f"no period found for order={order}, classes={sorted(classes)}")


def is_single_affine_coset(classes: frozenset[int], order: int) -> bool:
    if not classes or order % len(classes) != 0:
        return False
    step = order // len(classes)
    anchor = min(classes)
    candidate = frozenset((anchor + j * step) % order for j in range(len(classes)))
    return candidate == classes


def v2(value: int) -> int:
    require(value > 0, f"v2 requires a positive integer, got {value}")
    answer = 0
    while value % 2 == 0:
        answer += 1
        value //= 2
    return answer


def target_in_generated_subgroup_formula(generator, target, n: int) -> bool:
    """Solve k*(sigma,t)=(tau,a) in C2 x C_n without orbit enumeration."""
    sigma, t = generator
    tau, a = target
    if t == 0:
        if a != 0:
            return False
        if sigma == 0:
            return tau == 0
        return True
    common = gcd(t, n)
    if a % common != 0:
        return False
    if sigma == 0:
        return tau == 0
    quotient = n // common
    require(quotient >= 2 and quotient % 2 == 0, "unexpected nontrivial 2-group quotient")
    # After division by common, t/common is odd.  Its inverse modulo the even
    # quotient is odd, so the solution parity is exactly parity(a/common).
    return tau == ((a // common) % 2)


def universal_level(b: int) -> int:
    ell = 1
    while (1 << (ell - 1)) < b:
        ell += 1
    return ell


def collar_residues(b: int, ell: int) -> frozenset[int]:
    modulus = 1 << ell
    x = modulus // (2 * (b - 1))
    return frozenset(modulus - w for w in range(1, x + 1, 2))


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return value == divisor
        divisor += 1
    return True


def first_prime_in_classes(lower: int, modulus: int, classes: frozenset[int]) -> int:
    candidate = lower + 1
    while True:
        if is_prime(candidate) and candidate % modulus in classes:
            return candidate
        candidate += 1


CHARTS = {}
for chart_ell in range(2, ELL_MAX + 1):
    CHARTS[chart_ell] = signed_five_chart(chart_ell)


ORBIT_DATA = {}
coordinate_power_cells = 0
for ell in range(2, ELL_MAX + 1):
    modulus = 1 << ell
    n, encode, decode = CHARTS[ell]
    for base in sorted(decode):
        coordinate = decode[base]
        order = coordinate_order(coordinate, n)
        residues = tuple(encode[coordinate_multiple(coordinate, k, n)] for k in range(order))
        require(len(set(residues)) == order, f"orbit repeats early ell={ell}, p={base}")
        require(residues[0] == 1, f"orbit does not start at identity ell={ell}, p={base}")
        for k, residue in enumerate(residues):
            require(
                residue == pow(base, k, modulus),
                f"coordinate/direct power mismatch ell={ell}, p={base}, k={k}",
            )
            coordinate_power_cells += 1
        ORBIT_DATA[(ell, base)] = (order, residues, coordinate)


R_ATLAS = {}
closure_cells = 0
for b in B_VALUES:
    for ell in range(ELL_MIN, ELL_MAX + 1):
        units = (1,) if ell == 1 else tuple(sorted(CHARTS[ell][2]))
        accepted = set()
        for residue in units:
            if closes_odd_core(residue, b, ell):
                accepted.add(residue)
            closure_cells += 1
        R_ATLAS[(b, ell)] = frozenset(accepted)

require(
    all(1 not in R_ATLAS[(b, ell)] for b in B_VALUES for ell in range(ELL_MIN, ELL_MAX + 1)),
    "the identity residue unexpectedly closes in the declared universe",
)


FIRST_EXPECTED = {
    3: (2, frozenset({3})),
    5: (3, frozenset({5, 7})),
    7: (3, frozenset({7})),
    9: (4, frozenset({9, 15})),
    11: (4, frozenset({3, 15})),
    13: (4, frozenset({5, 15})),
}

UNIVERSAL_EXPECTED = {
    3: (3, frozenset({3, 7})),
    5: (4, frozenset({5, 7, 13, 15})),
    7: (4, frozenset({3, 5, 7, 15})),
    9: (5, frozenset({9, 15, 19, 25, 29, 31})),
    11: (5, frozenset({3, 7, 9, 15, 19, 21, 27, 31})),
    13: (5, frozenset({5, 15, 21, 31})),
}

for b in B_VALUES:
    first_ell, first_set = FIRST_EXPECTED[b]
    require(R_ATLAS[(b, first_ell)] == first_set, f"first set mismatch b={b}")
    require(
        all(not R_ATLAS[(b, ell)] for ell in range(1, first_ell)),
        f"claimed first level is not first for b={b}",
    )
    universal_ell, universal_set = UNIVERSAL_EXPECTED[b]
    require(universal_ell == universal_level(b), f"universal level formula mismatch b={b}")
    require(R_ATLAS[(b, universal_ell)] == universal_set, f"universal set mismatch b={b}")


K_CACHE = {}


def orbit_classes(b: int, ell: int, base: int):
    key = (b, ell, base)
    if key not in K_CACHE:
        order, residues, _ = ORBIT_DATA[(ell, base)]
        classes = frozenset(k for k, residue in enumerate(residues) if residue in R_ATLAS[(b, ell)])
        period = least_period(classes, order)
        reduced = frozenset(k for k in range(period) if k in classes)
        K_CACHE[key] = (order, classes, period, reduced)
    return K_CACHE[key]


def orbit_signature_groups(b: int, ell: int):
    groups = defaultdict(list)
    for base in sorted(CHARTS[ell][2]):
        order, classes, period, reduced = orbit_classes(b, ell, base)
        signature = (period, tuple(sorted(reduced))) if classes else (0, ())
        groups[signature].append(base)
    return groups


def cyclic_subgroups(ell: int):
    groups = defaultdict(list)
    for base in sorted(CHARTS[ell][2]):
        groups[frozenset(ORBIT_DATA[(ell, base)][1])].append(base)
    return groups


def maximal_empty_subgroups(b: int, ell: int):
    accepted = R_ATLAS[(b, ell)]
    empty = [(subgroup, tuple(generators)) for subgroup, generators in cyclic_subgroups(ell).items()
             if subgroup.isdisjoint(accepted)]
    return tuple(
        sorted(
            ((subgroup, generators) for subgroup, generators in empty
             if not any(subgroup < other for other, _ in empty)),
            key=lambda item: (len(item[0]), tuple(sorted(item[0]))),
        )
    )


def coordinate_cyclic_subgroup(element, n: int) -> frozenset[tuple[int, int]]:
    subgroup = {(0, 0)}
    state = (0, 0)
    while True:
        state = coordinate_sum(state, element, n)
        if state in subgroup:
            require(state == (0, 0), "coordinate cycle returned away from identity")
            return frozenset(subgroup)
        subgroup.add(state)


def all_coordinate_subgroups(n: int) -> frozenset[frozenset[tuple[int, int]]]:
    """Every subgroup of the rank-two group C2 x C_n, from at most two generators."""
    elements = tuple((sigma, t) for sigma in (0, 1) for t in range(n))
    cyclic = {element: coordinate_cyclic_subgroup(element, n) for element in elements}
    subgroups = set()
    for left in elements:
        for right in elements:
            generated = frozenset(
                coordinate_sum(x, y, n) for x in cyclic[left] for y in cyclic[right]
            )
            subgroups.add(generated)
    return frozenset(subgroups)


def maximal_group_cosets(b: int, ell: int):
    accepted = R_ATLAS[(b, ell)]
    n, encode, _ = CHARTS[ell]
    elements = tuple((sigma, t) for sigma in (0, 1) for t in range(n))
    candidates = set()
    for subgroup in all_coordinate_subgroups(n):
        for translate in elements:
            coordinate_coset = frozenset(
                coordinate_sum(translate, element, n) for element in subgroup
            )
            values = frozenset(encode[coordinate] for coordinate in coordinate_coset)
            if values and values.issubset(accepted):
                candidates.add((coordinate_coset, values))
    maximal = [
        candidate for candidate in candidates
        if not any(candidate[1] < other[1] for other in candidates)
    ]
    return tuple(sorted(maximal, key=lambda item: (-len(item[1]), tuple(sorted(item[0])))))


# Exact special lanes requested for hostile audit.
require(
    frozenset(base for base in CHARTS[5][2] if not orbit_classes(9, 5, base)[1])
    == frozenset({1, 7, 17, 23}),
    "b=9, ell=5 empty-orbit classification failed",
)
require(
    frozenset(base for base in CHARTS[5][2] if not orbit_classes(11, 5, base)[1])
    == frozenset({1, 17}),
    "b=11, ell=5 empty-orbit classification failed",
)
require(
    all(bool(orbit_classes(11, 5, base)[1]) == (base % 16 != 1) for base in CHARTS[5][2]),
    "b=11 orbit nonemptiness is not exactly p!=1 mod16",
)
EXPECTED_B13_MINIMAL_SIGNATURES = {
    (0, ()): [1, 3, 7, 9, 11],
    (2, (1,)): [15],
    (4, (1,)): [5],
    (4, (3,)): [13],
}
require(
    {key: value for key, value in orbit_signature_groups(13, 4).items()}
    == EXPECTED_B13_MINIMAL_SIGNATURES,
    "b=13, ell=4 fixed-prime family classification failed",
)


# Negative-collar inclusion and the coordinate membership criterion.
collar_formula_cells = 0
for b in B_VALUES:
    for ell in range(2, ELL_MAX + 1):
        collar = collar_residues(b, ell)
        require(collar.issubset(R_ATLAS[(b, ell)]), f"negative collar failed b={b}, ell={ell}")
        n, _, decode = CHARTS[ell]
        for base in sorted(decode):
            subgroup = frozenset(ORBIT_DATA[(ell, base)][1])
            generator_coordinate = decode[base]
            for target_residue in sorted(collar):
                target_coordinate = decode[target_residue]
                formula = target_in_generated_subgroup_formula(
                    generator_coordinate, target_coordinate, n
                )
                require(
                    formula == (target_residue in subgroup),
                    f"collar subgroup formula mismatch b={b}, ell={ell}, p={base}, target={target_residue}",
                )
                collar_formula_cells += 1


# Exact exponent-class lift law over every high residue class in the universe.
lift_base_cells = 0
lift_parent_exponent_cells = 0
ratio_hist = Counter()
transition_stats = {
    b: {"A": 0, "B": 0, "strict": 0, "awaken": 0, "first_awaken": None}
    for b in B_VALUES
}

for b in B_VALUES:
    for ell in range(2, ELL_MAX):
        low_modulus = 1 << ell
        high_ell = ell + 1
        for high_base in sorted(CHARTS[high_ell][2]):
            low_base = high_base % low_modulus
            low_order, low_classes, _, _ = orbit_classes(b, ell, low_base)
            high_order, high_classes, _, _ = orbit_classes(b, high_ell, high_base)
            require(high_order % low_order == 0, "order divisibility failed")
            ratio = high_order // low_order
            require(ratio in (1, 2), f"order jump outside 1/2: ell={ell}, p={high_base}")
            ratio_hist[ratio] += 1
            lift_base_cells += 1
            low_residues = ORBIT_DATA[(ell, low_base)][1]
            high_residues = ORBIT_DATA[(high_ell, high_base)][1]
            if ratio == 2:
                new_failed_children = 0
                for j in range(low_order):
                    low_residue = low_residues[j]
                    first_child = high_residues[j]
                    second_child = high_residues[j + low_order]
                    require(
                        frozenset({first_child, second_child})
                        == frozenset({low_residue, low_residue + low_modulus}),
                        f"two exponent lifts are not residue children ell={ell}, p={high_base}, j={j}",
                    )
                    child_hits = int(j in high_classes) + int(j + low_order in high_classes)
                    if j in low_classes:
                        require(child_hits == 2, "closed orbit parent did not have two closed children")
                    else:
                        require(child_hits in (0, 1), "failed orbit parent had two closed children")
                        new_failed_children += child_hits
                    lift_parent_exponent_cells += 1
                require(
                    len(high_classes) == 2 * len(low_classes) + new_failed_children,
                    "ratio-two orbit recurrence failed",
                )
                transition_stats[b]["A"] += new_failed_children
            else:
                new_observed_children = 0
                for j in range(low_order):
                    low_residue = low_residues[j]
                    high_residue = high_residues[j]
                    reduced_high = high_residue if high_residue < low_modulus else high_residue - low_modulus
                    require(reduced_high == low_residue, "unique observed child has wrong parent")
                    if j in low_classes:
                        require(j in high_classes, "closed orbit parent lost its observed child")
                    elif j in high_classes:
                        new_observed_children += 1
                    lift_parent_exponent_cells += 1
                require(
                    len(high_classes) == len(low_classes) + new_observed_children,
                    "ratio-one orbit recurrence failed",
                )
                transition_stats[b]["B"] += new_observed_children
            if len(high_classes) > ratio * len(low_classes):
                transition_stats[b]["strict"] += 1
            if not low_classes and high_classes:
                transition_stats[b]["awaken"] += 1
                if transition_stats[b]["first_awaken"] is None:
                    transition_stats[b]["first_awaken"] = (
                        ell,
                        low_base,
                        high_base,
                        low_order,
                        high_order,
                        tuple(sorted(high_classes)),
                    )


def level_summary(b: int, ell: int) -> str:
    bases = sorted(CHARTS[ell][2])
    records = [orbit_classes(b, ell, base) for base in bases]
    nonempty_bases = sum(bool(record[1]) for record in records)
    compressed_bases = sum(bool(record[1]) and record[2] < record[0] for record in records)
    coset_bases = sum(is_single_affine_coset(record[1], record[0]) for record in records)
    period_ratios = Counter(
        record[0] // record[2] for record in records if record[1]
    )
    subgroups = cyclic_subgroups(ell)
    empty_subgroups = sum(subgroup.isdisjoint(R_ATLAS[(b, ell)]) for subgroup in subgroups)
    collar = collar_residues(b, ell)
    collar_hit_subgroups = sum(not subgroup.isdisjoint(collar) for subgroup in subgroups)
    ratios = ",".join(f"{key}:{period_ratios[key]}" for key in sorted(period_ratios)) or "none"
    return (
        f"b={b} ell={ell} |R|={len(R_ATLAS[(b, ell)])}/{1 << (ell - 1)} "
        f"nonempty_bases={nonempty_bases}/{len(bases)} subgroups={len(subgroups)} "
        f"empty_subgroups={empty_subgroups} collar={len(collar)} "
        f"collar_hit_subgroups={collar_hit_subgroups} compressed_bases={compressed_bases} "
        f"single_coset_bases={coset_bases} period_compression={{{ratios}}}"
    )


print("FACTORIAL FIXED-PRIME EXPONENT ORBIT INDEPENDENT ATLAS")
print(f"universe odd_cores={B_VALUES} ell={ELL_MIN}..{ELL_MAX}")
print("representation=U(2^ell) signed-five coordinates; closure=schoolbook occupied-position sets")
print("K includes exponent class 0 when it closes; here 0 never closes because residue 1 is hostile")
print()

print("LEAST NONEMPTY AND LEAST UNIVERSAL SUFFIX LEVELS")
for b in B_VALUES:
    first_ell, first_set = FIRST_EXPECTED[b]
    universal_ell, universal_set = UNIVERSAL_EXPECTED[b]
    print(
        f"b={b} first_nonempty_ell={first_ell} R={fmt(first_set)} "
        f"least_universal_ell={universal_ell} R_universal={fmt(universal_set)}"
    )
print()

print("MINIMAL FIXED-PRIME CONGRUENCE FAMILIES (least exponent periods)")
for b in B_VALUES:
    ell = FIRST_EXPECTED[b][0]
    modulus = 1 << ell
    groups = orbit_signature_groups(b, ell)
    for signature in sorted(groups, key=lambda item: (item[0] == 0, item[0], item[1])):
        period, classes = signature
        bases = groups[signature]
        if period == 0:
            print(f"b={b} mod={modulus} empty_p_classes={fmt(bases)}")
        else:
            print(
                f"b={b} mod={modulus} p_classes={fmt(bases)} "
                f"k_mod={period} K={fmt(classes)}"
            )
print()

print("LEAST-UNIVERSAL ORBIT SIGNATURE ATLAS")
for b in B_VALUES:
    ell = universal_level(b)
    modulus = 1 << ell
    print(f"b={b} ell={ell} mod={modulus} R={fmt(R_ATLAS[(b, ell)])}")
    groups = orbit_signature_groups(b, ell)
    for signature in sorted(groups, key=lambda item: (item[0] == 0, item[0], item[1])):
        period, classes = signature
        bases = groups[signature]
        if period == 0:
            print(f"  empty p={fmt(bases)}")
        else:
            print(f"  p={fmt(bases)} k_mod={period} K={fmt(classes)}")
    hostile_classes = frozenset(base for base in CHARTS[ell][2] if not orbit_classes(b, ell, base)[1])
    hostile_prime = first_prime_in_classes(b, modulus, hostile_classes)
    print(f"  fixed-prime empty hostile: p={hostile_prime}, p_mod={hostile_prime % modulus}, K=empty")
print()

print("MAXIMAL EMPTY CYCLIC SUBGROUPS AT FIRST/UNIVERSAL LEVELS")
for b in B_VALUES:
    levels = sorted({FIRST_EXPECTED[b][0], universal_level(b)})
    for ell in levels:
        entries = maximal_empty_subgroups(b, ell)
        rendered = "; ".join(
            f"H={fmt(subgroup)} generators={fmt(generators)}" for subgroup, generators in entries
        )
        print(f"b={b} ell={ell}: {rendered}")
print()

print("MAXIMAL AFFINE GROUP COSETS INSIDE R AT LEAST-UNIVERSAL LEVEL")
for b in B_VALUES:
    ell = universal_level(b)
    pieces = []
    for coordinate_coset, values in maximal_group_cosets(b, ell):
        coordinates = "{" + ",".join(
            f"({sigma},{t})" for sigma, t in sorted(coordinate_coset)
        ) + "}"
        pieces.append(f"coordinates={coordinates},values={fmt(values)}")
    print(f"b={b} ell={ell}: " + "; ".join(pieces))
print()

print("SELECTED ALL-BASE LEVEL SUMMARIES")
for b in B_VALUES:
    selected_levels = sorted({
        FIRST_EXPECTED[b][0],
        universal_level(b),
        min(universal_level(b) + 1, ELL_MAX),
        8,
        10,
        ELL_MAX,
    })
    for ell in selected_levels:
        print(level_summary(b, ell))
print()

print("EXPONENT-CLASS LIFT TRANSITIONS")
print(
    f"base_transition_cells={lift_base_cells} parent_exponent_cells={lift_parent_exponent_cells} "
    f"order_ratio_hist={dict(sorted(ratio_hist.items()))}"
)
for b in B_VALUES:
    stats = transition_stats[b]
    first = stats["first_awaken"]
    require(first is not None, f"no empty-orbit awakening found for b={b}")
    ell, low_base, high_base, low_order, high_order, high_classes = first
    print(
        f"b={b} A_ratio2={stats['A']} B_ratio1={stats['B']} "
        f"strict_base_transitions={stats['strict']} empty_to_nonempty={stats['awaken']} "
        f"first_awaken=(ell {ell}->{ell + 1},p {low_base}->{high_base},"
        f"ord {low_order}->{high_order},K_high={fmt(high_classes)})"
    )
print()

print("REQUESTED EXACT HOSTILES")
print("b=9 ell=5 empty iff p mod32 in {1,7,17,23}")
print("b=11 ell=5 orbit nonempty iff p mod16 !=1")
print("b=13 ell=4 families: p=5,k=1 mod4; p=13,k=3 mod4; p=15,k=1 mod2; no others")
print("empty orbit is level-local: every b has an explicit empty-to-nonempty lift above")
print()

print("ALGEBRAIC COLLAR-HIT GATE")
print("coordinates generator=(sigma,t), target=(tau,a), n=2^(ell-2)")
print("t=0: require a=0 and either sigma=1 or tau=0")
print("t!=0: require gcd(t,n)|a; if sigma=0 also tau=0; if sigma=1 also tau=(a/gcd(t,n)) mod2")
print(f"collar_formula_cells={collar_formula_cells}")
print()

print("VERIFICATION COUNTS")
print(
    f"closure_cells={closure_cells} coordinate_power_cells={coordinate_power_cells} "
    f"K_records={len(K_CACHE)}"
)
print("status=PASS")

#!/usr/bin/env python3
"""Exact unrestricted positive-transverse zero-mode-cochain rank probe.

THM-3405 reduces every active zero-cochain family U=dV to two canonical
centres.  If g=gcd(q,d), Q=q/g, and d=g*d0, then multiplication by d0
permutes the Q-sheet quotient.  Every certificate is therefore a g-fold
pullback of a primitive gcd-one cover on Q sheets at centre 0 or 1/(2Q).

The primitive gcd-one condition is encoded exactly by adjoining one breaker
bit for every prime dividing the finite owner modulus.  A union-state BFS and
an independent exhaustive-combination solver classify Q=2..28.  Divisor
minimization then gives the unrestricted ranks for q=15..28.  Literal
positive-owner witnesses are replayed at c=1/(2q).  Exact quotient-order
capacity proves the universal rank floor four; the primitive Q=8,9 covers
are the only primitive rank-four covers.  Hence global rank four occurs
exactly when q is divisible by 8 or 9, with harmonic, Berggren-spine, and
Fibonacci pullbacks.

This is a zero-mode-cochain classification, not the larger synchronized
half-grid physical rank.  It gives no LRC(14) ledger decrement.
Runtime truth gates survive python -O.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINNED = (
    (
        "THM-3405-zero-mode-gauge",
        ROOT / "01-canon/theorems/THM-3405-common-centre-gcd-gauge-and-boolean-half-twist.md",
        "d3e7dbeeb85c6f897bd9e31270bd0b6602ae4feac3b46a45eb5ce23ae5d24fe0",
    ),
    (
        "THM-3401-fixed-zero-rank",
        ROOT / "01-canon/theorems/THM-3401-centered-transverse-sheet-cover-rank-fifteen-through-twenty-eight.md",
        "dae088cbda12fb64d24f84ab26a6879e94939e04cb03601d8fb996a48c077716",
    ),
    (
        "MISTAKE-389-halfgrid-repair",
        ROOT / "04-computation/lrc_all_owner_halfgrid_physical_chart_probe_20260815.py",
        "318705d41da80a69e1aa00f40c52d1677eaadaa00d57ea92bdafe108bc7087ab",
    ),
)

EXPECTED_PRIMITIVE = (
    (2, None, None),
    (3, None, None),
    (4, None, None),
    (5, None, None),
    (6, None, None),
    (7, None, None),
    (8, None, 4),
    (9, None, 4),
    (10, None, 5),
    (11, None, 6),
    (12, None, 5),
    (13, None, 7),
    (14, None, 7),
    (15, 6, 6),
    (16, 5, 5),
    (17, 8, 8),
    (18, 5, 5),
    (19, 9, 9),
    (20, 6, 6),
    (21, 8, 8),
    (22, 7, 6),
    (23, 11, 6),
    (24, 6, 5),
    (25, 11, 6),
    (26, 8, 7),
    (27, 10, 5),
    (28, 8, 7),
)

EXPECTED_GLOBAL = (
    (15, 6, ((15, 0), (15, 1))),
    (16, 4, ((8, 1),)),
    (17, 8, ((17, 0), (17, 1))),
    (18, 4, ((9, 1),)),
    (19, 9, ((19, 0), (19, 1))),
    (20, 5, ((10, 1),)),
    (21, 8, ((21, 0), (21, 1))),
    (22, 6, ((11, 1), (22, 1))),
    (23, 6, ((23, 1),)),
    (24, 4, ((8, 1),)),
    (25, 6, ((25, 1),)),
    (26, 7, ((13, 1), (26, 1))),
    (27, 4, ((9, 1),)),
    (28, 7, ((14, 1), (28, 1))),
)

EXPLICIT_HALF_WITNESSES = (
    (15, 15, (1, 4, 6, 7, 8, 10)),
    (16, 8, (2, 6, 10, 14)),
    (17, 17, (1, 3, 4, 5, 7, 8, 9, 11)),
    (18, 9, (2, 10, 12, 14)),
    (19, 19, (1, 3, 4, 5, 7, 8, 9, 11, 12)),
    (20, 10, (2, 6, 8, 14, 18)),
    (21, 21, (1, 4, 5, 6, 8, 11, 13, 14)),
    (22, 11, (2, 4, 6, 10, 14, 18)),
    (23, 23, (1, 4, 5, 7, 9, 11)),
    (24, 8, (3, 9, 15, 21)),
    (25, 25, (1, 9, 10, 11, 19, 21)),
    (26, 13, (2, 4, 6, 10, 14, 18, 22)),
    (27, 9, (3, 15, 18, 21)),
    (28, 14, (2, 6, 8, 10, 18, 22, 26)),
)

EXPECTED_FIXED_ZERO = (6, 5, 8, 5, 9, 6, 8, 7, 11, 6, 11, 8, 10, 8)
EXPECTED_HALFGRID = (3, 2, 8, 2, 9, 2, 3, 2, 6, 2, 5, 2, 3, 2)
EXPECTED_ZERO_CRITICAL_ORDERS = (3, 4, 5, 6, 15, 16, 17, 18, 29, 30)
EXPECTED_ZERO_CRITICAL_PAIRS = (
    (3, 3, 6, 1),
    (3, 4, 12, 1),
    (3, 5, 30, 1),
    (3, 6, 6, 0),
    (3, 15, 30, 1),
    (3, 16, 48, 1),
    (3, 17, 102, 1),
    (3, 18, 18, 0),
    (3, 29, 174, 1),
    (3, 30, 30, 0),
    (4, 4, 4, 0),
)
EXPECTED_HALF_CRITICAL_ORDERS = (5, 8, 9, 10, 11, 12, 15, 17, 22, 23, 24, 29, 36)
EXPECTED_HALF_RANK4_TAILS = (
    (5, 8, 8),
    (5, 8, 9),
    (8, 8, 8),
    (8, 8, 9),
    (8, 8, 10),
    (8, 8, 11),
    (8, 8, 12),
    (8, 8, 15),
    (8, 8, 17),
    (8, 8, 22),
    (8, 8, 23),
    (8, 8, 24),
    (8, 8, 29),
    (8, 8, 36),
    (8, 9, 9),
    (8, 9, 10),
    (8, 9, 15),
    (9, 9, 9),
)
EXPECTED_SEMANTIC_DIGEST = "233c092a9b73dcf8a40b9c21b52b99e322b059f003554bd79721bae317c30c7e"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def divisors(value):
    return tuple(divisor for divisor in range(1, value + 1) if value % divisor == 0)


def prime_factors(value):
    factors = []
    trial = 2
    remaining = value
    while trial * trial <= remaining:
        if remaining % trial == 0:
            factors.append(trial)
            while remaining % trial == 0:
                remaining //= trial
        trial += 1
    if remaining > 1:
        factors.append(remaining)
    return tuple(factors)


def type_danger_mask(q, residue, epsilon):
    modulus = 2 * q
    mask = 0
    for sheet in range(q):
        phase_word = residue * (2 * sheet + epsilon) % modulus
        if 14 * min(phase_word, modulus - phase_word) < modulus:
            mask |= 1 << sheet
    return mask


def direct_danger_mask(q, owner, centre):
    mask = 0
    for sheet in range(q):
        phase = (owner * (centre + Fraction(sheet, q))) % 1
        if 14 * min(phase, 1 - phase) < 1:
            mask |= 1 << sheet
    return mask


def augmented_type_bank(q, epsilon):
    require(epsilon in (0, 1), (q, epsilon))
    modulus = q if epsilon == 0 else 2 * q
    primes = prime_factors(modulus)
    grouped = {}
    raw_count = 0
    for residue in range(1, modulus):
        if residue % q == 0:
            continue
        sheet_mask = type_danger_mask(q, residue, epsilon)
        if not sheet_mask:
            continue
        raw_count += 1
        augmented = sheet_mask
        for offset, prime in enumerate(primes):
            if residue % prime:
                augmented |= 1 << (q + offset)
        if augmented not in grouped or residue < grouped[augmented]:
            grouped[augmented] = residue

    items = tuple(sorted(((mask, residue) for mask, residue in grouped.items()), key=lambda item: item[1]))
    maximal = tuple(
        item
        for item in items
        if not any(item[0] != other[0] and item[0] | other[0] == other[0] for other in items)
    )
    full = (1 << (q + len(primes))) - 1
    return modulus, primes, raw_count, items, maximal, full


def union_state_bfs(full, items):
    frontier = {0: ()}
    visited = {0}
    transitions = 0
    rank = 0
    while frontier:
        rank += 1
        following = {}
        for state, witness in frontier.items():
            for mask, residue in items:
                transitions += 1
                joined = state | mask
                if joined == state or joined in visited or joined in following:
                    continue
                candidate = witness + (residue,)
                following[joined] = candidate
        if full in following:
            return rank, tuple(sorted(following[full])), len(visited) + len(following), transitions
        visited.update(following)
        frontier = following
    return None, (), len(visited), transitions


def exhaustive_combinations(q, epsilon, full, items, event_hasher):
    tests = 0
    for rank in range(1, len(items) + 1):
        for chosen in combinations(items, rank):
            tests += 1
            joined = 0
            residues = []
            for mask, residue in chosen:
                joined |= mask
                residues.append(residue)
            event_hasher.update(
                repr((q, epsilon, rank, tuple(residues), joined)).encode("ascii") + bytes((10,))
            )
            if joined == full:
                return rank, tuple(residues), tests
    return None, (), tests


def lift_residues_to_gcd_one(modulus, residues):
    require(len(residues) >= 2, (modulus, residues))
    rest = gcd(*residues[1:])
    period = 1
    for prime in prime_factors(rest):
        period *= prime
    for shift in range(period + 1):
        first = residues[0] + shift * modulus
        lifted = (first,) + tuple(residues[1:])
        if gcd(*lifted) == 1:
            require(
                tuple(value % modulus for value in lifted)
                == tuple(value % modulus for value in residues),
                (modulus, residues, lifted),
            )
            return lifted
    raise RuntimeError(("gcd-one lift missing", modulus, residues, period))


def primitive_record(q, epsilon, event_hasher):
    modulus, primes, raw_count, items, maximal, full = augmented_type_bank(q, epsilon)
    bfs_rank, bfs_witness, bfs_states, bfs_transitions = union_state_bfs(full, maximal)
    combo_rank, combo_witness, combo_tests = exhaustive_combinations(
        q, epsilon, full, items, event_hasher
    )
    require(bfs_rank == combo_rank, (q, epsilon, bfs_rank, combo_rank))

    lifted = ()
    if bfs_rank is not None:
        require(bfs_rank >= 2, (q, epsilon, bfs_rank))
        lifted = lift_residues_to_gcd_one(modulus, bfs_witness)
        require(gcd(*lifted) == 1, (q, epsilon, lifted))
        union = 0
        augmented_union = 0
        grouped = dict((residue, mask) for mask, residue in items)
        for residue in bfs_witness:
            union |= type_danger_mask(q, residue, epsilon)
            augmented_union |= grouped[residue]
        require(union == (1 << q) - 1, (q, epsilon, bfs_witness, union))
        require(augmented_union == full, (q, epsilon, bfs_witness, augmented_union))

    return (
        q,
        epsilon,
        modulus,
        primes,
        raw_count,
        len(items),
        len(maximal),
        bfs_rank,
        bfs_witness,
        combo_witness,
        lifted,
        bfs_states,
        bfs_transitions,
        combo_tests,
    )


def verify_zero_mode_pullback(q, quotient, epsilon, primitive_lifts):
    scale = q // quotient
    owners = tuple(scale * value for value in primitive_lifts)
    require(len(owners) == len(set(owners)), (q, quotient, epsilon, owners))
    require(all(owner > 0 and owner % q for owner in owners), (q, quotient, epsilon, owners))
    owner_gcd = gcd(*owners)
    gauge_gcd = gcd(q, owner_gcd)
    require(owner_gcd == scale and gauge_gcd == scale, (q, quotient, owners))
    centre = Fraction(0) if epsilon == 0 else Fraction(1, 2 * q)
    scalar = Fraction(2 * q * owner_gcd) * centre
    require(scalar.denominator == 1, (q, quotient, epsilon, scalar))
    require(scalar.numerator == epsilon * gauge_gcd, (q, quotient, epsilon, scalar))

    masks = tuple(direct_danger_mask(q, owner, centre) for owner in owners)
    require(all(masks), (q, quotient, epsilon, owners, masks))
    union = 0
    for mask in masks:
        union |= mask
    require(union == (1 << q) - 1, (q, quotient, epsilon, owners, masks))

    for owner, value, mask in zip(owners, primitive_lifts, masks):
        centre_word = Fraction(2 * q * owner) * centre
        require(centre_word.denominator == 1, (q, owner, centre_word))
        require(gcd(q, owner) and centre_word.numerator % gcd(q, owner) == 0, (q, owner))
        base_mask = type_danger_mask(quotient, value, epsilon)
        expected = sum(1 << sheet for sheet in range(q) if base_mask & (1 << (sheet % quotient)))
        require(mask == expected, (q, quotient, epsilon, owner, mask, expected))

    return (
        quotient,
        epsilon,
        centre,
        owners,
        tuple(mask.bit_count() for mask in masks),
        owner_gcd,
        gauge_gcd,
        scalar.numerator,
    )


def verify_explicit_half_witness(q, quotient, owners):
    centre = Fraction(1, 2 * q)
    scale = q // quotient
    require(gcd(*owners) == scale, (q, quotient, owners, gcd(*owners)))
    primitive = tuple(owner // scale for owner in owners)
    require(gcd(*primitive) == 1, (q, quotient, primitive))
    record = verify_zero_mode_pullback(q, quotient, 1, primitive)
    require(record[3] == owners, (q, quotient, owners, record))
    return record


def capacity_profile(q, prime):
    records = []
    for epsilon in (0, 1):
        modulus = q if epsilon == 0 else 2 * q
        sizes_breaking = []
        sizes_divisible = []
        for residue in range(1, modulus):
            if residue % q == 0:
                continue
            size = type_danger_mask(q, residue, epsilon).bit_count()
            if not size:
                continue
            if residue % prime:
                sizes_breaking.append(size)
            else:
                sizes_divisible.append(size)
        records.append(
            (
                epsilon,
                max(sizes_breaking, default=0),
                max(sizes_divisible, default=0),
            )
        )
    return q, prime, tuple(records)


def zero_layer_order_count(order):
    return 1 + 2 * ((order - 1) // 14)


def half_odd_residue_count(order):
    largest = (order - 1) // 7
    return 2 * ((largest + 1) // 2)


def half_layer_order_count(order):
    odd_count = half_odd_residue_count(order)
    if order % 2 == 0:
        return odd_count
    return max(odd_count, zero_layer_order_count(order))


def lcm(left, right):
    return left // gcd(left, right) * right


def pisano_period(modulus):
    values = [0]
    previous, current = 0, 1
    while True:
        previous, current = current, (previous + current) % modulus
        if (previous, current) == (0, 1):
            return tuple(values)
        values.append(previous)


def rank_four_floor_audit(primitive_by_key):
    # At zero twist, a quotient-order-m owner sees every point j/m and hence
    # exactly z(m)=1+2 floor((m-1)/14) dangerous phases.
    zero_formula = tuple(
        (
            order,
            zero_layer_order_count(order),
            type_danger_mask(order, 1, 0).bit_count(),
        )
        for order in range(2, 501)
    )
    require(all(expected == actual for _, expected, actual in zero_formula), zero_formula)
    require(
        all(3 * zero_layer_order_count(order) <= order for order in range(3, 501)),
        "zero density exceeds one third",
    )
    require(
        tuple(
            order
            for order in range(3, 501)
            if 3 * zero_layer_order_count(order) == order
        )
        == (3,),
        "zero one-third equality",
    )

    # If one zero-layer owner has order two, the other two can meet its
    # remaining half-mass only in this finite list.  The interval formula
    # z(m)=2k+1 on 14k+1<=m<=14k+14 makes k<=2 exact.
    critical_orders = tuple(
        order
        for order in range(3, 501)
        if 6 * zero_layer_order_count(order) >= order
    )
    require(critical_orders == EXPECTED_ZERO_CRITICAL_ORDERS, critical_orders)
    require(14 * 3 + 1 > 12 * 3 + 6 and 14 > 12, "critical interval tail")

    critical_pairs = []
    for left in critical_orders:
        for right in critical_orders:
            if right < left:
                continue
            if (
                2 * zero_layer_order_count(left) * right
                + 2 * zero_layer_order_count(right) * left
                < left * right
            ):
                continue
            quotient = 2 * left // gcd(2, left)
            quotient = quotient * right // gcd(quotient, right)
            capacity = (
                quotient // 2
                + (quotient // left) * zero_layer_order_count(left)
                + (quotient // right) * zero_layer_order_count(right)
            )
            critical_pairs.append((left, right, quotient, capacity - quotient))
    critical_pairs = tuple(critical_pairs)
    require(critical_pairs == EXPECTED_ZERO_CRITICAL_PAIRS, critical_pairs)
    require(all(excess <= 1 for _, _, _, excess in critical_pairs), critical_pairs)
    require(
        all(
            quotient + excess - 2 < quotient
            for _, _, quotient, excess in critical_pairs
        ),
        critical_pairs,
    )

    # At half twist an order-two block is empty.  Every order m>=3 has at
    # most ceil(m/7) dangerous quotient phases, at most one third of m with
    # equality only at m=3.  All nonempty order-three types are the same
    # singleton, so the equality case cannot cover.
    require(type_danger_mask(2, 1, 1) == 0, "half order-two block")
    half_caps = tuple((order, (order + 6) // 7) for order in range(3, 501))
    require(all(3 * cap <= order for order, cap in half_caps), half_caps)
    require(tuple(order for order, cap in half_caps if 3 * cap == order) == (3,), half_caps)
    zero_order_three = {
        type_danger_mask(3, residue, 0)
        for residue in range(1, 3)
        if type_danger_mask(3, residue, 0)
    }
    half_order_three = {
        type_danger_mask(3, residue, 1)
        for residue in range(1, 6)
        if residue % 3 and type_danger_mask(3, residue, 1)
    }
    require(zero_order_three == {1 << 0}, zero_order_three)
    require(half_order_three == {1 << 1}, half_order_three)

    # Primitive Q=8 and Q=9 half-twist covers attain the universal floor.
    ray_records = []
    support = []
    for q in range(2, 501):
        ancestors = tuple(base for base in (8, 9) if q % base == 0)
        if not ancestors:
            continue
        support.append(q)
        pulls = []
        for quotient in ancestors:
            primitive = primitive_by_key[(quotient, 1)]
            require(primitive[7] == 4, (q, quotient, primitive[7]))
            pulls.append(verify_zero_mode_pullback(q, quotient, 1, primitive[10]))
        ray_records.append((q, ancestors, tuple(pulls)))
    ray_records = tuple(ray_records)
    require(len(ray_records) == 111, len(ray_records))
    require(sum(len(record[1]) for record in ray_records) == 117, ray_records)

    support_residues = tuple(
        residue for residue in range(72) if residue % 8 == 0 or residue % 9 == 0
    )
    require(len(support_residues) == 16, support_residues)

    berggren_residues = tuple(
        residue
        for residue in range(9)
        if (4 * residue * residue + 12 * residue + 11) % 9 == 0
    )
    require(berggren_residues == (1, 5), berggren_residues)

    period_eight = pisano_period(8)
    period_nine = pisano_period(9)
    require(len(period_eight) == 12, len(period_eight))
    require(len(period_nine) == 24, len(period_nine))
    require(
        tuple(index for index, value in enumerate(period_eight) if value == 0)
        == (0, 6),
        period_eight,
    )
    require(
        tuple(index for index, value in enumerate(period_nine) if value == 0)
        == (0, 12),
        period_nine,
    )

    return (
        sha256(repr(zero_formula).encode("ascii")).hexdigest(),
        critical_orders,
        critical_pairs,
        zero_order_three,
        half_order_three,
        len(ray_records),
        sha256(repr(ray_records).encode("ascii")).hexdigest(),
        support_residues,
        berggren_residues,
        period_eight,
        period_nine,
    )


def exact_order_augmented_bank(q, order):
    modulus = 2 * q
    primes = prime_factors(modulus)
    grouped = {}
    for residue in range(1, modulus):
        if q // gcd(q, residue) != order:
            continue
        sheet_mask = type_danger_mask(q, residue, 1)
        if not sheet_mask:
            continue
        augmented = sheet_mask
        for offset, prime in enumerate(primes):
            if residue % prime:
                augmented |= 1 << (q + offset)
        if augmented not in grouped or residue < grouped[augmented]:
            grouped[augmented] = residue
    return tuple(sorted(((mask, residue) for mask, residue in grouped.items()), key=lambda item: item[1]))


def rank_four_profile_record(orders):
    q = 1
    for order in orders:
        q = lcm(q, order)
    primes = prime_factors(2 * q)
    full = (1 << (q + len(primes))) - 1
    banks = {order: exact_order_augmented_bank(q, order) for order in set(orders)}
    grouped_choices = []
    for order in sorted(banks):
        multiplicity = orders.count(order)
        grouped_choices.append(tuple(combinations(banks[order], multiplicity)))

    tests = 0
    witnesses = []
    for choice_groups in product(*grouped_choices):
        selected = tuple(item for group in choice_groups for item in group)
        tests += 1
        union = 0
        for mask, _ in selected:
            union |= mask
        if union == full:
            witnesses.append(tuple(sorted(residue for _, residue in selected)))
    return (
        orders,
        q,
        tuple((order, len(banks[order])) for order in sorted(banks)),
        tests,
        tuple(witnesses),
    )


def primitive_rank_four_classification_audit():
    # Zero twist: without order two, useful order-three and order-four blocks
    # are each unique.  The density hierarchy forces either duplicated
    # order-three/four blocks or total mass below one.  With one order-two
    # block, restrict the other owners to the odd sheet coset.
    require(
        all(4 * zero_layer_order_count(order) <= order for order in range(4, 501)),
        "zero quarter density",
    )
    require(
        tuple(
            order
            for order in range(4, 501)
            if 4 * zero_layer_order_count(order) == order
        )
        == (4,),
        "zero quarter equality",
    )
    require(
        all(5 * zero_layer_order_count(order) <= order for order in range(5, 501)),
        "zero fifth density",
    )
    zero_order_four = {
        type_danger_mask(4, residue, 0)
        for residue in range(1, 4)
        if gcd(4, residue) == 1 and type_danger_mask(4, residue, 0)
    }
    require(zero_order_four == {1 << 0}, zero_order_four)

    odd_coset_density = []
    for order in range(3, 501):
        if order % 2:
            numerator = zero_layer_order_count(order)
            denominator = order
        else:
            local_odd = 2 * ((((order - 1) // 14) + 1) // 2)
            numerator = 2 * local_odd
            denominator = order
        odd_coset_density.append((order, numerator, denominator))
    require(
        all(3 * numerator <= denominator for _, numerator, denominator in odd_coset_density),
        odd_coset_density,
    )
    require(
        tuple(
            order
            for order, numerator, denominator in odd_coset_density
            if 3 * numerator == denominator
        )
        == (3,),
        odd_coset_density,
    )

    # Half twist: h(m) is the exact maximum block count on the order-m
    # quotient.  Odd owner words permute odd residues modulo 2m; for odd m,
    # even owner words additionally realize the zero-grid count z(m).
    half_formula = []
    for order in range(2, 501):
        direct_counts = [type_danger_mask(order, 1, 1).bit_count()]
        if order % 2:
            direct_counts.append(type_danger_mask(order, 2, 1).bit_count())
        actual = max(direct_counts)
        expected = half_layer_order_count(order)
        require(actual == expected, (order, expected, actual))
        half_formula.append((order, expected))
    half_formula = tuple(half_formula)
    require(
        all(
            4 * count <= order
            for order, count in half_formula
            if order != 3
        ),
        "half quarter density",
    )
    require(
        tuple(
            order
            for order, count in half_formula
            if order != 3 and 4 * count == order
        )
        == (8,),
        half_formula,
    )
    require(
        all(
            9 * count <= 2 * order
            for order, count in half_formula
            if order not in (3, 8)
        ),
        "half two-ninth density",
    )
    require(
        tuple(
            order
            for order, count in half_formula
            if order not in (3, 8) and 9 * count == 2 * order
        )
        == (9,),
        half_formula,
    )

    critical_orders = tuple(
        order
        for order in range(4, 501)
        if 6 * half_layer_order_count(order) >= order
    )
    require(critical_orders == EXPECTED_HALF_CRITICAL_ORDERS, critical_orders)
    one_eight_pairs = tuple(
        (left, right)
        for left in critical_orders
        for right in critical_orders
        if left <= right and left != 8 and right != 8
        and 12
        * (
            half_layer_order_count(left) * right
            + half_layer_order_count(right) * left
        )
        >= 5 * left * right
    )
    require(one_eight_pairs == ((5, 9), (9, 9), (9, 10), (9, 15)), one_eight_pairs)

    tails = []
    for third in critical_orders:
        tails.append(tuple(sorted((8, 8, third))))
    for left, right in one_eight_pairs:
        tails.append(tuple(sorted((8, left, right))))
    tails.append((9, 9, 9))
    tails = tuple(sorted(set(tails)))
    require(tails == EXPECTED_HALF_RANK4_TAILS, tails)

    profiles = ((8, 8, 8, 8),) + tuple((3,) + tail for tail in tails)
    profile_records = tuple(rank_four_profile_record(profile) for profile in profiles)
    positive_profiles = tuple(
        (orders, q, witnesses)
        for orders, q, _, _, witnesses in profile_records
        if witnesses
    )
    require(
        positive_profiles
        == (
            ((8, 8, 8, 8), 8, ((1, 3, 5, 7),)),
            ((3, 9, 9, 9), 9, ((1, 5, 6, 7),)),
        ),
        positive_profiles,
    )

    return (
        sha256(repr(odd_coset_density).encode("ascii")).hexdigest(),
        sha256(repr(half_formula).encode("ascii")).hexdigest(),
        critical_orders,
        one_eight_pairs,
        tails,
        profile_records,
        positive_profiles,
        sum(record[3] for record in profile_records),
    )


def main():
    for name, path, expected in PINNED:
        require(lf_hash(path) == expected, ("dependency changed", name, path))

    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)
        ),
        "floating literal",
    )

    event_hasher = sha256()
    primitive_records = tuple(
        primitive_record(q, epsilon, event_hasher)
        for q in range(2, 29)
        for epsilon in (0, 1)
    )
    primitive_by_key = {(record[0], record[1]): record for record in primitive_records}
    primitive_table = tuple(
        (
            q,
            primitive_by_key[(q, 0)][7],
            primitive_by_key[(q, 1)][7],
        )
        for q in range(2, 29)
    )
    require(primitive_table == EXPECTED_PRIMITIVE, ("primitive table", primitive_table))

    global_records = []
    for q in range(15, 29):
        candidates = []
        for quotient in divisors(q):
            if quotient < 2:
                continue
            for epsilon in (0, 1):
                primitive = primitive_by_key[(quotient, epsilon)]
                rank = primitive[7]
                if rank is None:
                    continue
                pullback = verify_zero_mode_pullback(q, quotient, epsilon, primitive[10])
                candidates.append((rank, quotient, epsilon, pullback))
        require(candidates, ("no global candidate", q))
        best = min(candidate[0] for candidate in candidates)
        minimizers = tuple(
            (quotient, epsilon)
            for rank, quotient, epsilon, _ in candidates
            if rank == best
        )
        global_records.append((q, best, minimizers, tuple(candidates)))
    global_records = tuple(global_records)
    global_table = tuple((q, rank, minimizers) for q, rank, minimizers, _ in global_records)
    require(global_table == EXPECTED_GLOBAL, ("global table", global_table))

    explicit_records = tuple(
        verify_explicit_half_witness(q, quotient, owners)
        for q, quotient, owners in EXPLICIT_HALF_WITNESSES
    )
    require(
        tuple(len(record[3]) for record in explicit_records)
        == tuple(rank for _, rank, _ in EXPECTED_GLOBAL),
        explicit_records,
    )

    zero_mode_sequence = tuple(rank for _, rank, _ in EXPECTED_GLOBAL)
    fixed_gap = tuple(
        fixed - mobile for fixed, mobile in zip(EXPECTED_FIXED_ZERO, zero_mode_sequence)
    )
    halfgrid_gap = tuple(
        mobile - physical for mobile, physical in zip(zero_mode_sequence, EXPECTED_HALFGRID)
    )
    require(fixed_gap == (0, 1, 0, 1, 0, 1, 0, 1, 5, 2, 5, 1, 6, 1), fixed_gap)
    require(halfgrid_gap == (3, 2, 0, 2, 0, 3, 5, 4, 0, 2, 1, 5, 1, 5), halfgrid_gap)

    capacity_controls = (
        capacity_profile(25, 5),
        capacity_profile(9, 3),
        capacity_profile(27, 3),
    )
    rank_four_floor = rank_four_floor_audit(primitive_by_key)
    primitive_rank_four_classification = primitive_rank_four_classification_audit()

    solver_audit = (
        sum(record[11] for record in primitive_records),
        sum(record[12] for record in primitive_records),
        sum(record[13] for record in primitive_records),
        event_hasher.hexdigest(),
    )

    semantic = sha256(
        repr(
            (
                primitive_records,
                global_records,
                explicit_records,
                fixed_gap,
                halfgrid_gap,
                capacity_controls,
                rank_four_floor,
                primitive_rank_four_classification,
                solver_audit,
            )
        ).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_DIGEST:
        require(semantic == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", semantic))

    print("LRC UNRESTRICTED ZERO-MODE-COCHAIN DIVISOR-RANK PROBE")
    print(f"source_sha256_lf={lf_hash(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINNED)}")
    print("status=PROVED-ANALYTIC divisor_pullback_reduction_and_universal_rank_floor_4_and_primitive_rank4_iff_Q=8_or_9_and_global_rank4_iff_8|q_or_9|q;FINITE-EXACT q15..28 unrestricted_positive_transverse_zero_mode_cochain_ranks;INDEPENDENT_BFS_AND_COMBINATION_SOLVERS;NOT_halfgrid_physical_rank;no_LRC14_decrement")
    print("reduction=U=dV;g=gcd(q,d);Q=q/g;epsilon_in_{0,1};every_certificate_is_g_fold_pullback_of_primitive_gcd_one_Q_cover_at_0_or_1/(2Q)")
    print("primitive_gcd_gate=gcd(M_epsilon,r_1,...,r_s)=1;implemented_as_one_prime_breaker_bit_per_prime_dividing_M_epsilon")
    print(f"primitive_rank_table_Q2_Q28=(Q,zero,half)={primitive_table}")
    print(f"global_rank_table_q15_q28=(q,rank,minimizing_(Q,epsilon))={global_table}")
    print(f"explicit_half_centre_witnesses=(Q,epsilon,c,U,block_sizes,d,g,a)={explicit_records}")
    print(f"three_layer_sequences=(fixed_zero,genuine_zero_mode,halfgrid_physical)={(EXPECTED_FIXED_ZERO, zero_mode_sequence, EXPECTED_HALFGRID)}")
    print(f"fixed_minus_zero_mode={fixed_gap};zero_mode_minus_halfgrid={halfgrid_gap}")
    print(f"capacity_controls_(Q,p,((epsilon,max_breaker,max_divisible),...))={capacity_controls}")
    print("universal_floor=every_positive_transverse_zero_mode_cochain_cover_has_at_least_4_owners;zero_twist_forced_common_sheet_and_finite_order_capacity;half_twist_order2_empty_and_order3_equality_collapses_to_one_block")
    print(f"rank4_floor_audit=(zero_formula_sha256,critical_orders,critical_pairs,zero_m3_masks,half_m3_masks,ray_count,ray_sha256,residues_mod72,berggren_n_mod9,pisano8,pisano9)={rank_four_floor}")
    print(f"primitive_rank4_classification=(zero_odd_coset_sha256,half_order_formula_sha256,critical_orders,one8_pairs,critical_tails,profile_records,positive_profiles,total_profile_tests)={primitive_rank_four_classification}")
    print("exact_rank4_classification=rho_ZMC(q)=4_iff_8|q_or_9|q;full_support_density=2/9;reciprocal_mass=(2/9)log_N+O(1);Berggren_Q_n_rank4_for_n=1,5_mod9;Fibonacci_F_n_rank4_for_6|n")
    print(f"solver_audit=(BFS_states,BFS_transitions,combination_tests,event_sha256)={solver_audit}")
    print("scope=distinct_positive_transverse_owners_of_unrestricted_size;allowing_q_divisible_owners_would_trivialize_fixed_zero_to_rank_one")
    print(f"semantic_sha256={semantic}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()

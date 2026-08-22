#!/usr/bin/env python3
"""Exact referee for THM-3429's prime-fibre activity descent.

The unbounded proof is elementary.  Its cover universe is transverse:
Q does not divide any selected residue, so residues 0 and Q modulo 2Q are
excluded.  On a p-point fibre of Z/Q -> Z/(Q/p),
an owner with p|r is a full pullback, while an owner with p not dividing r
hits at most ceil(p/7) points.  This companion reconstructs those statements
directly for every transverse residue on every odd composite Q<=315, freezes
the activity/defect arithmetic, and audits the sharp Q=51 positive and Q=39
nonprimitive hostile.

No floating-point decision occurs.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from math import gcd, isqrt, lcm


TARGET_BASES = (8, 9, 10, 11, 12, 15, 23, 25)
Q51_RESIDUES = (1, 11, 12, 18, 23, 34, 35)
Q13_RESIDUES = (1, 2, 3, 5, 7, 9, 11)
EXPECTED_SEMANTIC_SHA256 = "c3b8596233d20f55e5c0bf0419c5f3f627cbf371e3004a5f40ab8a582a09b9f8"
TRANSVERSE_UNIVERSE = "Q_does_not_divide_each_selected_residue"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def prime_factors(n: int) -> tuple[int, ...]:
    factors = []
    candidate = 2
    while candidate * candidate <= n:
        if n % candidate == 0:
            factors.append(candidate)
            while n % candidate == 0:
                n //= candidate
        candidate += 1 if candidate == 2 else 2
    if n > 1:
        factors.append(n)
    return tuple(factors)


def primes_through(limit: int) -> tuple[int, ...]:
    return tuple(
        n for n in range(2, limit + 1)
        if all(n % p for p in range(2, isqrt(n) + 1))
    )


def danger(q: int, residue: int, sheet: int) -> bool:
    modulus = 2 * q
    phase = (residue * (2 * sheet + 1)) % modulus
    distance = min(phase, modulus - phase)
    return 14 * distance < modulus


def danger_mask(q: int, residue: int) -> int:
    return sum(1 << sheet for sheet in range(q) if danger(q, residue, sheet))


def quotient_order(q: int, residue: int) -> int:
    return q // gcd(q, residue)


def joint_period(q: int, residues: tuple[int, ...]) -> int:
    return lcm(*(quotient_order(q, residue) for residue in residues))


def covers(q: int, residues: tuple[int, ...]) -> bool:
    joined = 0
    for residue in residues:
        joined |= danger_mask(q, residue)
    return joined == (1 << q) - 1


def fibre_cap(p: int) -> int:
    return (p + 6) // 7


def activity_floor(p: int) -> int:
    cap = fibre_cap(p)
    return (p + cap - 1) // cap


def is_target_free(q: int) -> bool:
    return all(q % base for base in TARGET_BASES)


def exact_fibre_universe(limit: int = 315) -> tuple[int, int, int, int, int]:
    moduli = 0
    prime_fibres = 0
    pullback_cells = 0
    active_cells = 0
    largest_active_hit = 0
    for q in range(9, limit + 1, 2):
        factors = prime_factors(q)
        if len(factors) == 1 and factors[0] == q:
            continue
        moduli += 1
        for p in factors:
            prime_fibres += 1
            m = q // p
            for residue in range(2 * q):
                if residue % q == 0:
                    continue
                if residue % p == 0:
                    reduced = residue // p
                    for sheet in range(q):
                        require(
                            danger(q, residue, sheet)
                            == danger(m, reduced, sheet % m),
                            ("pullback", q, p, residue, sheet),
                        )
                        pullback_cells += 1
                else:
                    for base_sheet in range(m):
                        hits = sum(
                            danger(q, residue, base_sheet + step * m)
                            for step in range(p)
                        )
                        require(hits <= fibre_cap(p),
                                ("active cap", q, p, residue, base_sheet, hits))
                        active_cells += p
                        largest_active_hit = max(largest_active_hit, hits)
    return moduli, prime_fibres, pullback_cells, active_cells, largest_active_hit


def threshold_arithmetic() -> tuple[tuple[tuple[int, int, int], ...], tuple[int, ...], tuple[int, ...]]:
    primes = tuple(p for p in primes_through(997) if p % 2)
    rows = tuple((p, fibre_cap(p), activity_floor(p)) for p in primes)
    below_seven = tuple(p for p, _, floor in rows if floor < 7)
    equality_primes = tuple(p for p, cap, floor in rows if cap * floor == p)
    require(below_seven == (3, 5, 11, 17, 23, 29), below_seven)
    require(equality_primes == (3, 5, 7), equality_primes)
    require(all(floor == 7 for p, _, floor in rows if p in (7, 13, 19) or p >= 31), rows)
    selected = tuple(row for row in rows if row[0] <= 43)
    return selected, below_seven, equality_primes


def fibre_record(q: int, residues: tuple[int, ...], p: int) -> tuple[object, ...]:
    m = q // p
    active = tuple(index for index, residue in enumerate(residues) if residue % p)
    dormant = tuple(index for index, residue in enumerate(residues) if residue % p == 0)
    dormant_base_mask = 0
    for index in dormant:
        dormant_base_mask |= danger_mask(m, residues[index] // p)
    missed = tuple(sheet for sheet in range(m) if not (dormant_base_mask >> sheet) & 1)
    hit_profiles = []
    for base_sheet in missed:
        hits = tuple(
            sum(
                danger(q, residues[index], base_sheet + step * m)
                for step in range(p)
            )
            for index in active
        )
        union = 0
        for index in active:
            for step in range(p):
                sheet = base_sheet + step * m
                if danger(q, residues[index], sheet):
                    union |= 1 << step
        require(union == (1 << p) - 1, (q, p, base_sheet, hits, union))
        hit_profiles.append((base_sheet, hits, sum(hits) - p))
    return p, m, active, dormant, missed, tuple(hit_profiles)


def q51_positive() -> tuple[object, ...]:
    q = 51
    residues = Q51_RESIDUES
    require(is_target_free(q), (q, "not target free"))
    require(covers(q, residues), (q, residues, "not a cover"))
    orders = tuple(quotient_order(q, residue) for residue in residues)
    require(Counter(orders) == Counter({3: 1, 17: 2, 51: 4}), orders)
    require(joint_period(q, residues) == q, orders)
    activities = tuple((p, sum(residue % p != 0 for residue in residues)) for p in (3, 17))
    require(activities == ((3, 5), (17, 6)), activities)
    defects = tuple((p, tuple(i for i, r in enumerate(residues) if r % p == 0)) for p in (3, 17))
    full_order_indices = tuple(i for i, order in enumerate(orders) if order == q)
    require(full_order_indices == (0, 1, 4, 6), full_order_indices)
    return orders, activities, defects, full_order_indices, fibre_record(q, residues, 3), fibre_record(q, residues, 17)


def q39_joint_period_hostile() -> tuple[object, ...]:
    q = 39
    residues = tuple(3 * residue for residue in Q13_RESIDUES)
    require(is_target_free(q), (q, "not target free"))
    require(covers(q, residues), (q, residues, "not a cover"))
    orders = tuple(quotient_order(q, residue) for residue in residues)
    require(joint_period(q, residues) == 13, orders)
    require(all(residue % 3 == 0 for residue in residues), residues)
    return residues, orders, joint_period(q, residues), sum(residue % 3 != 0 for residue in residues)


def boolean_defect_budget() -> tuple[object, ...]:
    exceptional = (3, 5, 7, 17, 29)
    budgets = tuple((p, 7 - activity_floor(p)) for p in exceptional)
    require(budgets == ((3, 4), (5, 2), (7, 0), (17, 1), (29, 1)), budgets)
    pair_full_order_floors = tuple(
        (left, right, max(0, 7 - (7 - activity_floor(left)) - (7 - activity_floor(right))))
        for left, right in ((3, 17), (3, 29), (5, 17), (5, 29), (17, 29))
    )
    require(
        pair_full_order_floors
        == ((3, 17, 2), (3, 29, 2), (5, 17, 4), (5, 29, 4), (17, 29, 5)),
        pair_full_order_floors,
    )
    return budgets, pair_full_order_floors


def main() -> None:
    fibre_universe = exact_fibre_universe()
    thresholds = threshold_arithmetic()
    q51 = q51_positive()
    q39 = q39_joint_period_hostile()
    boolean_budget = boolean_defect_budget()
    transverse_controls = tuple(
        (q, danger_mask(q, 0).bit_count(), danger_mask(q, q).bit_count())
        for q in (13, 31, 51)
    )
    require(transverse_controls == ((13, 13, 0), (31, 31, 0), (51, 51, 0)),
            transverse_controls)
    semantic_surface = (
        TRANSVERSE_UNIVERSE,
        fibre_universe,
        thresholds,
        q51,
        q39,
        boolean_budget,
        transverse_controls,
    )
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256,
                (semantic_digest, EXPECTED_SEMANTIC_SHA256))

    print("THM-3429 prime-fibre activity descent -- exact referee")
    print("status=VERIFIED_EXACT_SUPPORT_FOR_PROVED_UNBOUNDED_ELEMENTARY_PROOF;odd_half_twist;cap7;joint_period_typed;no_LRC14_decrement")
    print(f"typed_universe={TRANSVERSE_UNIVERSE};r=0_universal_and_r=Q_empty_are_excluded")
    print(f"fibre_universe=(composite_moduli,prime_fibres,pullback_cells,active_cells,max_active_hit)={fibre_universe}")
    print(f"threshold_arithmetic=(selected_rows,activity_floor_lt7,equality_primes)={thresholds}")
    print(f"boolean_defect_budget=(single_prime,pair_full_order_floors)={boolean_budget}")
    print(f"Q51_sharp_positive=(orders,activities,defects,full_order_indices,p3,p17)={q51}")
    print(f"Q39_joint_period_hostile=(residues,orders,joint_period,a3)={q39}")
    print(f"transverse_boundary_controls=(Q,r0_size,rQ_size)={transverse_controls}")
    print(f"semantic_sha256={semantic_digest}")
    print("scope=prime_fibre_descent_and_five_lane_reduction_only;mixed_small_prime_lanes_open;fixed_zero_open;arbitrary_time_open;LRC14_open")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact residual 17-adic tower probe downstream of THM-3429.

The unbounded part is an elementary fibre-incidence proof.  For a p=17
active half-twist block, every 17-fibre contains two or three points.  The
number of third-point events is therefore its total mass minus twice the
number of base fibres.  THM-3429 supplies exactly six active owners and one
even inactive owner in the residual towers.  Comparing the required and
available third-point events excludes 17^a for a>=2 and 5*17^a for a>=2.

The sole coarse-boundary modulus Q=85 is closed twice: a cross-prime p=17/p=5
fixed-fibre incidence deficit and an unrestricted exact literal/joint-period
set-cover search.  The discarded coarse profile's overlap census is retained
as a hostile check, not as a needed third proof.
Every truth gate uses ``require`` and therefore survives ``python -O``.
"""

from __future__ import annotations

import ast
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
THM3429 = ROOT / "01-canon/theorems/THM-3429-prime-fibre-activity-descent-for-mixed-order-half-twist-seven-covers.md"
THM3416 = ROOT / "01-canon/theorems/THM-3416-zero-mode-cochain-global-rank-six-support.md"
THM3421 = ROOT / "01-canon/theorems/THM-3421-prime-half-twist-rank-seven-classification.md"
PINNED = (
    ("THM-3429", THM3429, "58ebf850fc79fc9afed57966b7599e7376a6684fa3bbc5a2aa2e1a8e6e0ca148"),
    ("THM-3416", THM3416, "42a9309145de51d1bb6fca0b7c1945302ff37a63a3183e1dfed838c07118e8bf"),
    ("THM-3421", THM3421, "2f577354a06628660d90f70aca34378a186b7e67a55795e6c4a5c32a255d9736"),
)
EXPECTED_SEMANTIC_SHA256 = "2939e6278cc64659599ec030d620cd52b12574cbfb2f531a3e34d8ebf8df1fcd"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def prime_factors(value: int) -> tuple[int, ...]:
    factors = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            factors.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1
    if value > 1:
        factors.append(value)
    return tuple(factors)


def quotient_order(q: int, residue: int) -> int:
    return q // gcd(q, residue)


def half_mask(q: int, residue: int) -> int:
    """Literal strict half-twist danger mask on the q odd sheets."""
    modulus = 2 * q
    result = 0
    for sheet in range(q):
        if is_danger(q, residue, sheet):
            result |= 1 << sheet
    return result


def is_danger(q: int, residue: int, sheet: int) -> bool:
    modulus = 2 * q
    word = residue * (2 * sheet + 1) % modulus
    return 14 * min(word, modulus - word) < modulus


def h(order: int) -> int:
    """THM-3416 equation (7), exact maximal half-twist block mass."""
    zero_branch = 1 + 2 * ((order - 1) // 14)
    half_word = (order - 1) // 7
    half_branch = 2 * ((half_word + 1) // 2)
    return half_branch if order % 2 == 0 else max(half_branch, zero_branch)


def divisors(value: int) -> tuple[int, ...]:
    result = [1]
    remaining = value
    prime = 2
    while prime * prime <= remaining:
        if remaining % prime:
            prime += 1
            continue
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        old = tuple(result)
        multiplier = 1
        for _ in range(exponent):
            multiplier *= prime
            result.extend(base * multiplier for base in old)
        prime += 1
    if remaining > 1:
        old = tuple(result)
        result.extend(base * remaining for base in old)
    return tuple(sorted(result))


def fibre_mask(q: int, prime: int, residue: int, base: int) -> int:
    require(q % prime == 0, (q, prime))
    cofactor = q // prime
    result = 0
    for lift in range(prime):
        sheet = base + cofactor * lift
        if is_danger(q, residue, sheet):
            result |= 1 << lift
    return result


def active_fibre_profile(q: int, prime: int, residue: int) -> tuple[int, ...]:
    require(q % prime == 0 and residue % prime, (q, prime, residue))
    cofactor = q // prime
    profile = tuple(fibre_mask(q, prime, residue, base).bit_count() for base in range(cofactor))
    require(sum(profile) == half_mask(q, residue).bit_count(), (q, prime, residue, profile))
    return profile


def divisor_density_audit():
    """Prove the two downstairs density maxima from THM-3416's h(m)."""
    require(h(5) == 1 and h(17) == 3, (h(5), h(17)))

    # Pure tower: m=17 is equality.  For every larger possible order m>=289,
    # h(m)/m <= (m+6)/(7m) <= 3/17 because m>=26.
    require(Fraction(289 + 6, 7 * 289) <= Fraction(3, 17), "pure tail")
    pure_samples = tuple((power, h(17**power), Fraction(h(17**power), 17**power)) for power in range(1, 9))
    require(max(row[2] for row in pure_samples) == Fraction(3, 17), pure_samples)

    # Mixed tower: order 5 is equality.  Every other nontrivial divisor order
    # is at least 17, and h(m)/m <= (m+6)/(7m) <= 1/5 for m>=15.
    require(Fraction(17 + 6, 7 * 17) <= Fraction(1, 5), "mixed tail")
    mixed_orders = tuple(sorted({5, *(17**power for power in range(1, 8)), *(5 * 17**power for power in range(1, 8))}))
    mixed_samples = tuple((order, h(order), Fraction(h(order), order)) for order in mixed_orders)
    require(max(row[2] for row in mixed_samples) == Fraction(1, 5), mixed_samples)
    return pure_samples, mixed_samples


def tower_invoice_audit():
    """Freeze the all-height third-point inequalities and exact h samples."""
    # Q=17N.  Six full-order active owners supply at most
    # (18N+36)/7 third events.  An inactive block covers b<=3N/17 base
    # fibres, so at least 70N/17 events are required.  Their gap is
    # (184N-612)/119, positive already at N=17 and increasing.
    pure_base_n = 17
    pure_supply_upper = Fraction(18 * pure_base_n + 36, 7)
    pure_demand_lower = Fraction(70 * pure_base_n, 17)
    pure_gap = pure_demand_lower - pure_supply_upper
    require(pure_gap == Fraction(2516, 119) > 0, pure_gap)
    require(184 > 0, "pure gap slope")

    # Q=17N with N=5*17^(a-1).  A 17-active block has coindex 1 or 5,
    # hence mass <=(Q+30)/7.  The inactive downstairs density is <=1/5.
    # Supply is <=(18N+180)/7 and demand >=4N; the gap is
    # (10N-180)/7, positive for N>=85.
    mixed_base_n = 85
    mixed_supply_upper = Fraction(18 * mixed_base_n + 180, 7)
    mixed_demand_lower = Fraction(4 * mixed_base_n, 1)
    mixed_gap = mixed_demand_lower - mixed_supply_upper
    require(mixed_gap == Fraction(670, 7) > 0, mixed_gap)
    require(10 > 0, "mixed gap slope")

    exact_pure = []
    exact_mixed = []
    for exponent in range(2, 13):
        q = 17**exponent
        n = q // 17
        downstairs_orders = tuple(17**power for power in range(1, exponent))
        b_max = max(n // order * h(order) for order in downstairs_orders)
        supply = 6 * h(q) - 12 * n
        demand = 5 * (n - b_max)
        require(supply < demand, ("pure", exponent, q, n, supply, demand, b_max))
        exact_pure.append((exponent, q % 14, n, h(q), b_max, supply, demand, demand - supply))

        q5 = 5 * q
        n5 = q5 // 17
        downstairs_orders5 = tuple(order for order in divisors(n5) if order > 1)
        b5_max = max(n5 // order * h(order) for order in downstairs_orders5)
        active_max = max(h(q5), 5 * h(q))
        supply5 = 6 * active_max - 12 * n5
        demand5 = 5 * (n5 - b5_max)
        require(supply5 < demand5, ("mixed", exponent, q5, n5, supply5, demand5, b5_max))
        exact_mixed.append((exponent, q5 % 14, n5, active_max, b5_max, supply5, demand5, demand5 - supply5))

    symbolic = (
        ("pure", "T<=6h(Q)-12N<=(18N+36)/7", "T>=5(N-b)>=70N/17", "gap=(184N-612)/119", "N>=17"),
        ("mixed", "T<=6(Q+30)/7-12N=(18N+180)/7", "T>=5(N-b)>=4N", "gap=(10N-180)/7", "N>=85"),
    )
    return symbolic, pure_gap, mixed_gap, tuple(exact_pure), tuple(exact_mixed)


def fibre_lemma_hostiles():
    """Directly test the 2/3 fibre law, pullbacks, and the Q=51 boundary."""
    active_cells = 0
    profile_histograms = []
    for q in (51, 85, 289, 1445):
        require(q % 17 == 0, q)
        profiles = Counter()
        for residue in range(1, 2 * q):
            if residue % q == 0 or residue % 17 == 0:
                continue
            profile = active_fibre_profile(q, 17, residue)
            require(set(profile) <= {2, 3}, (q, residue, profile))
            profiles[(profile.count(2), profile.count(3), half_mask(q, residue).bit_count())] += 1
            active_cells += len(profile)
        profile_histograms.append((q, tuple(sorted(profiles.items()))))

    pullback_cells = 0
    for q in (51, 85, 289, 1445):
        n = q // 17
        for residue in range(17, 2 * q, 17):
            if residue % q == 0:
                continue
            downstairs = half_mask(n, residue // 17)
            upstairs = half_mask(q, residue)
            reconstructed = 0
            for base in range(n):
                if downstairs >> base & 1:
                    for lift in range(17):
                        reconstructed |= 1 << (base + n * lift)
                pullback_cells += 17
            require(upstairs == reconstructed, (q, residue))

    q51_residues = (1, 11, 12, 18, 23, 34, 35)
    q51_inactive = 34
    q51_base = tuple(base for base in range(3) if half_mask(3, q51_inactive // 17) >> base & 1)
    q51_active = tuple(residue for residue in q51_residues if residue % 17)
    q51_fibres = tuple(
        (
            base,
            tuple(active_fibre_profile(51, 17, residue)[base] for residue in q51_active),
        )
        for base in range(3)
    )
    # On the two inactive-missed fibres all six active blocks have size three
    # and their union is all 17 points; this is the sharp affine-lift hostile.
    q51_union_rows = []
    for base in range(3):
        joined = 0
        incidence = 0
        for residue in q51_active:
            local = fibre_mask(51, 17, residue, base)
            joined |= local
            incidence += local.bit_count()
        q51_union_rows.append((base, incidence, joined.bit_count()))
    require(q51_base == (1,), q51_base)
    require(tuple(q51_union_rows) == ((0, 18, 17), (1, 14, 7), (2, 18, 17)), q51_union_rows)
    return active_cells, pullback_cells, tuple(profile_histograms), q51_residues, q51_base, q51_fibres, tuple(q51_union_rows)


def q85_profile_arithmetic():
    """Derive the unique Q=85 parity/order profile before any mask search."""
    # x/y: even/odd order-17 owners (p5-inactive); e/o: even/odd full
    # order-85 owners.  There are six p17-active owners and at most two are
    # p5-inactive.  Their third-event supplies are 5,0,3,2 respectively.
    raw_profiles = []
    after_covered_fibre = []
    for x in range(3):
        for y in range(3 - x):
            for e in range(7 - x - y):
                o = 6 - x - y - e
                if o < 0:
                    continue
                third_events = 5 * x + 3 * e + 2 * o
                row = (x, y, e, o, third_events)
                if third_events >= 20:
                    raw_profiles.append(row)
                    # Every even active owner, including every even full-order
                    # owner, has size three on the reflection-fixed p17 fibre.
                    # Hence x+e third events are necessarily spent on the sole
                    # fibre already covered by the inactive order-5 owner.
                    if third_events - x - e >= 20:
                        after_covered_fibre.append(row)
    require(tuple(raw_profiles) == ((1, 0, 5, 0, 20), (2, 0, 2, 2, 20), (2, 0, 3, 1, 21), (2, 0, 4, 0, 22)), raw_profiles)
    require(tuple(after_covered_fibre) == (), after_covered_fibre)
    return tuple(raw_profiles), tuple(after_covered_fibre)


def q85_type_banks():
    grouped = defaultdict(dict)
    for residue in range(1, 170):
        if residue % 85 == 0:
            continue
        if residue % 17 == 0:
            order_type = "p17_inactive_order5"
        elif residue % 5 == 0:
            order_type = "p17_active_order17"
        else:
            order_type = "full_order85"
        parity = "even" if residue % 2 == 0 else "odd"
        grouped[(order_type, parity)].setdefault(half_mask(85, residue), residue)
    return grouped


def q85_cross_prime_audit():
    """Hostile overlap bank for the coarse profile killed by the invoice."""
    grouped = q85_type_banks()
    counts = tuple(sorted((key, len(value)) for key, value in grouped.items()))
    expected_counts = (
        (("full_order85", "even"), 32),
        (("full_order85", "odd"), 32),
        (("p17_active_order17", "even"), 8),
        (("p17_active_order17", "odd"), 8),
        (("p17_inactive_order5", "even"), 1),
        (("p17_inactive_order5", "odd"), 1),
    )
    require(counts == expected_counts, counts)

    inactive_items = tuple(grouped[("p17_inactive_order5", "even")].items())
    lower_items = tuple(grouped[("p17_active_order17", "even")].items())
    full_items = tuple(grouped[("full_order85", "even")].items())
    require(len(inactive_items) == 1 and len(lower_items) == 8 and len(full_items) == 32, counts)
    inactive_mask, inactive_residue = inactive_items[0]
    require(inactive_residue == 34 and inactive_mask.bit_count() == 17, (inactive_residue, inactive_mask.bit_count()))

    covered_p17_bases = tuple(base for base in range(5) if fibre_mask(85, 17, inactive_residue, base))
    require(covered_p17_bases == (2,), covered_p17_bases)
    missed_p17_bases = tuple(base for base in range(5) if base not in covered_p17_bases)

    p5_fixed_fibre = tuple(sorted(8 + 17 * lift for lift in range(5)))
    require(tuple(sorted(sheet % 5 for sheet in p5_fixed_fibre)) == (0, 1, 2, 3, 4), p5_fixed_fibre)
    require(sum(sheet % 5 == 2 for sheet in p5_fixed_fibre) == 1, p5_fixed_fibre)

    lower_records = []
    for mask, residue in lower_items:
        require(residue % 10 == 0 and quotient_order(85, residue) == 17, residue)
        downstairs_residue = residue // 5
        require(downstairs_residue % 2 == 0, (residue, downstairs_residue))
        downstairs = half_mask(17, downstairs_residue)
        require(downstairs >> 8 & 1, (residue, downstairs))
        require(all(mask >> sheet & 1 for sheet in p5_fixed_fibre), (residue, p5_fixed_fibre))
        profile = active_fibre_profile(85, 17, residue)
        require(profile == (3, 3, 3, 3, 3), (residue, profile))
        lower_records.append((residue, downstairs.bit_count(), profile))

    pair_overlaps = []
    outside_inactive = ((1 << 85) - 1) ^ inactive_mask
    for (left_mask, left_residue), (right_mask, right_residue) in combinations(lower_items, 2):
        overlap = left_mask & right_mask
        outside_count = (overlap & outside_inactive).bit_count()
        missed_fibre_counts = tuple(
            (fibre_mask(85, 17, left_residue, base) & fibre_mask(85, 17, right_residue, base)).bit_count()
            for base in missed_p17_bases
        )
        require(outside_count >= 4 and missed_fibre_counts == (1, 1, 1, 1), (left_residue, right_residue, outside_count, missed_fibre_counts))
        pair_overlaps.append((left_residue, right_residue, overlap.bit_count(), outside_count, missed_fibre_counts))
    outside_histogram = tuple(sorted(Counter(row[3] for row in pair_overlaps).items()))
    missed_pattern_histogram = tuple(sorted(Counter(row[4] for row in pair_overlaps).items()))
    global_overlap_histogram = tuple(sorted(Counter(row[2] for row in pair_overlaps).items()))
    require(outside_histogram == ((4, 28),), pair_overlaps)
    require(missed_pattern_histogram == (((1, 1, 1, 1), 28),), pair_overlaps)
    require(global_overlap_histogram == ((5, 28),), pair_overlaps)
    pair_summary = (len(pair_overlaps), global_overlap_histogram, outside_histogram, missed_pattern_histogram)
    return counts, inactive_residue, covered_p17_bases, missed_p17_bases, p5_fixed_fibre, tuple(lower_records), pair_summary


def generic_joint_period_bank(q: int):
    """Sheet masks augmented only by prime-activity bits for L(R)=Q."""
    factors = prime_factors(q)
    grouped = {}
    raw = 0
    for residue in range(1, 2 * q):
        # The positive-transverse universe excludes the zero/order-one
        # classes 0 and Q modulo 2Q.  The latter is empty at half twist.
        if residue % q == 0:
            continue
        sheet_mask = half_mask(q, residue)
        if not sheet_mask:
            continue
        raw += 1
        augmented = sheet_mask
        for offset, prime in enumerate(factors):
            if residue % prime:
                augmented |= 1 << (q + offset)
        grouped.setdefault(augmented, residue)
    unique = tuple(sorted(((mask, residue) for mask, residue in grouped.items()), key=lambda row: row[1]))
    maximal = tuple(
        item
        for item in unique
        if not any(item[0] != other[0] and item[0] | other[0] == other[0] for other in unique)
    )
    full = (1 << (q + len(factors))) - 1
    return factors, raw, unique, maximal, full


def exact_joint_period_cover(q: int, cap: int):
    factors, raw, unique, maximal, full = generic_joint_period_bank(q)
    masks = tuple(mask for mask, _ in maximal)
    residues = tuple(residue for _, residue in maximal)
    width = q + len(factors)
    by_bit = tuple(tuple(index for index, mask in enumerate(masks) if mask >> bit & 1) for bit in range(width))
    require(all(by_bit), (q, tuple(bit for bit, entries in enumerate(by_bit) if not entries)))
    nodes = 0
    branches = 0

    @lru_cache(maxsize=None)
    def solve(state: int, slots: int):
        nonlocal nodes, branches
        nodes += 1
        if state == full:
            return ()
        if slots == 0:
            return None
        missing = full ^ state
        gains = sorted(((mask & missing).bit_count() for mask in masks if mask | state != state), reverse=True)
        if not gains or sum(gains[:slots]) < missing.bit_count():
            return None
        missing_bits = tuple(bit for bit in range(width) if missing >> bit & 1)
        pivot = min(
            missing_bits,
            key=lambda bit: (sum(masks[index] | state != state for index in by_bit[bit]), bit),
        )
        candidates = sorted(
            (index for index in by_bit[pivot] if masks[index] | state != state),
            key=lambda index: (-(masks[index] & missing).bit_count(), residues[index]),
        )
        for index in candidates:
            branches += 1
            suffix = solve(state | masks[index], slots - 1)
            if suffix is not None:
                return (index,) + suffix
        return None

    answer = solve(0, cap)
    witness = None if answer is None else tuple(sorted(residues[index] for index in answer))
    if witness is not None:
        joined = 0
        orders = []
        for residue in witness:
            joined |= half_mask(q, residue)
            orders.append(quotient_order(q, residue))
        require(joined == (1 << q) - 1, (q, witness, joined.bit_count()))
        joint = 1
        for order in orders:
            joint = joint * order // gcd(joint, order)
        require(joint == q and len(witness) <= cap, (q, witness, orders, joint))
    return (q, factors, cap, raw, len(unique), len(maximal), witness, nodes, branches, solve.cache_info().hits)


def q85_forced_profile_census():
    """Exhaust the discarded coarse profile as an independent hostile."""
    grouped = q85_type_banks()
    inactive = tuple(grouped[("p17_inactive_order5", "even")].keys())
    lower = tuple(grouped[("p17_active_order17", "even")].keys())
    full = tuple(grouped[("full_order85", "even")].keys())
    require((len(inactive), len(lower), len(full)) == (1, 8, 32), (len(inactive), len(lower), len(full)))
    four_unions = tuple(first | second | third | fourth for first, second, third, fourth in combinations(full, 4))
    histogram = Counter()
    covers = 0
    for left, right in combinations(lower, 2):
        base = inactive[0] | left | right
        for four_union in four_unions:
            covered = (base | four_union).bit_count()
            histogram[covered] += 1
            if covered == 85:
                covers += 1
    tests = len(tuple(combinations(lower, 2))) * len(four_unions)
    require(tests == 1_006_880 and sum(histogram.values()) == tests, (tests, sum(histogram.values())))
    require(covers == 0 and max(histogram) == 77, (covers, max(histogram)))
    expected_histogram = (
        (49, 64), (51, 800), (53, 5920), (55, 31200), (57, 92880),
        (59, 181696), (61, 254072), (63, 217440), (65, 139040),
        (67, 58624), (69, 20288), (71, 3776), (73, 928), (75, 128), (77, 24),
    )
    require(tuple(sorted(histogram.items())) == expected_histogram, histogram)
    return len(four_unions), tests, covers, max(histogram), tuple(sorted(histogram.items()))


def main() -> None:
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert found")

    dependencies = tuple((label, lf_sha256(path)) for label, path, _ in PINNED)
    for label, path, expected in PINNED:
        require(lf_sha256(path) == expected, (label, lf_sha256(path), expected))

    density = divisor_density_audit()
    towers = tower_invoice_audit()
    fibres = fibre_lemma_hostiles()
    q85_profiles = q85_profile_arithmetic()
    q85_cross = q85_cross_prime_audit()
    generic_searches = (
        exact_joint_period_cover(51, 7),
        exact_joint_period_cover(85, 7),
        exact_joint_period_cover(289, 7),
    )
    require(generic_searches[0][6] == (1, 11, 12, 18, 23, 34, 35), generic_searches[0])
    require(generic_searches[1][6] is None and generic_searches[2][6] is None, generic_searches)
    q85_census = q85_forced_profile_census()

    consequence = (
        "NO_TARGET_FREE_LITERAL_JOINT_PERIOD_CAP7_COVER_FOR_Q_17POW_A_GE_2",
        "NO_TARGET_FREE_LITERAL_JOINT_PERIOD_CAP7_COVER_FOR_Q_5_TIMES_17POW_A_GE_1",
    )
    semantic_surface = (
        dependencies,
        density,
        towers,
        fibres,
        q85_profiles,
        q85_cross,
        generic_searches,
        q85_census,
        consequence,
    )
    semantic_digest = sha256(repr(semantic_surface).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256, (semantic_digest, EXPECTED_SEMANTIC_SHA256))

    print("LRC rank-seven residual 17-adic tower probe")
    print("status=VERIFIED_EXACT_INDEPENDENT_COMPANION_FOR_PROVED_THM3434;closes_both_THM3429_residual_towers;no_LRC14_decrement")
    print(f"dependency_sha256_lf={dependencies}")
    print("typed_universe=Q_odd_target_free;one_to_seven_nontrivial_residues_mod_2Q;residue_not_0_mod_Q;literal_strict_half_twist_sheet_union;joint_quotient_period_Q;multisets_allowed_so_negative_is_stronger;no_positive_lift_parity_assumed")
    print("inherited_profile=THM3429_plus_THM3421_fixed_fibre_gate:for_p17_exactly_one_inactive_owner;it_is_even;exactly_six_active_owners;inactive_mask_is_full_fibre_pullback")
    print("fibre_lemma=each_p17_active_block_has_2_or_3_points_on_every_17_fibre;T=sum_active_sizes-12N_counts_third_point_events;every_inactive_missed_fibre_needs_at_least_5_events;equivalently_active_incidence_at_least_17(N-b)+12b=17N-5b")
    print(f"downstairs_density_samples=(pure_17powers,mixed_5x17powers)={density}")
    print(f"all_height_invoices=(symbolic,pure_base_gap,mixed_base_gap)={towers[:3]}")
    print(f"exact_h_height_samples_pure={towers[3]}")
    print(f"exact_h_height_samples_mixed={towers[4]}")
    print(f"fibre_and_pullback_hostiles=(active_cells,pullback_cells,profile_histograms,Q51_residues,Q51_inactive_base,Q51_active_profiles,Q51_union_rows)={fibres}")
    print(f"Q85_profile_arithmetic=(coarse_mass_feasible,after_fixed_fibre_invoice)={q85_profiles}")
    print(f"Q85_cross_prime=(type_counts,inactive_residue,covered_p17_bases,missed_p17_bases,p5_fixed_fibre,lower_records,pair_summary)={q85_cross}")
    print("Q85_cross_prime_mechanism=every_coarse_mass_feasible_profile_spends_one_fixed_fibre_third_event_per_even_active_owner_and_none_retains_the_20_events_needed_on_four_missed_fibres;the_separate_order17_pair_overlap_bank_is_a_hostile_sidecar_not_the_proof")
    print(f"generic_joint_period_searches=(Q,factors,cap,raw,unique,maximal,witness,nodes,branches,hits)={generic_searches}")
    print(f"Q85_forced_profile_census=(four_unions,tests,covers,max_union,histogram)={q85_census}")
    print(f"consequence={consequence}")
    print("positive_hostile=Q51_mixed_atom_survives_and_is_found_exactly;its_order3_inactive_density_and_affine_lift_character_show_why_prime_support_alone_does_not_imply_the_tower_invoice")
    print("scope=literal_half_twist_common_centre_only;uses_proved_THM3429_profile_and_THM3416_h_formula;does_not_classify_arbitrary_common_time;does_not_supply_physical_runner_transport;LRC14_remains_open")
    print(f"semantic_sha256={semantic_digest}")
    print(f"script_sha256_lf={lf_sha256(source)}")
    print("reproducibility=standard_library_only;no_elapsed_fields;all_truth_gates_survive_python_O")
    print("commands=python -B 04-computation/lrc_rank7_residual_17adic_tower_probe_20260815.py;python -B -O 04-computation/lrc_rank7_residual_17adic_tower_probe_20260815.py")


if __name__ == "__main__":
    main()

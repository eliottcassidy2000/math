#!/usr/bin/env python3
"""Exact all-owner synchronized half-grid physical rank by affine charts.

At a physical time c satisfying 2*q*u*c in Z for every active owner, write
theta=q*c=a/b in lowest terms.  Then b|2u.  If b is odd write u=bv; if
b=2d write u=dv.  Exact danger blocks depend on v(a+b*ell) modulo q and
v(a+2d*ell) modulo 2q.  Sheet-affine permutations and unit reindexing
normalize these to

    O(q,g): v(1+g*ell) mod q,   g=gcd(b,q), g odd,
    E(q,g): v(1+2g*ell) mod 2q, g=gcd(d,q),

for proper divisors g of q.  Thus every positive-owner synchronized
half-grid physical cover belongs to one of finitely many exact set-cover
charts.

MISTAKE-389 records the repaired type: the half-grid condition is necessary
but not sufficient for c to be a THM-3398 mode centre.  THM-3405 adds the
missing mode divisibility.  This artifact therefore computes half-grid
physical ranks, not zero-mode-cochain ranks.  It has no LRC(14) ledger
consequence.  It also checks divisor-pullback families proving half-grid rank
two for even q>=8, rank three for odd 3|q, and rank five for odd 5|q with
3 not dividing q (q>=25).  On the even rank-two branch it reconstructs the
unique containing modes and proves the exact quadratic leakage ladder
P=-a*q^2/2 for every odd 1<=a<q/7.
Runtime gates survive python -O.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINNED = (
    (
        "THM-3398",
        ROOT / "01-canon/theorems/THM-3398-general-finite-mode-sheet-cover-cochain.md",
        "01901da2bb382184cfe4466550afe79255598f580f00a761fc32731a52ec9378",
    ),
    (
        "THM-3402-constructive-locus",
        ROOT / "01-canon/theorems/THM-3402-atomized-sheet-covers-and-constructive-cochain-locus.md",
        "0aaff0ffe66042ccae8de3158b1cb7ece056264fe96753b0bb166167728d472f",
    ),
    (
        "THM-3405-zero-mode-gauge",
        ROOT / "01-canon/theorems/THM-3405-common-centre-gcd-gauge-and-boolean-half-twist.md",
        "d3e7dbeeb85c6f897bd9e31270bd0b6602ae4feac3b46a45eb5ce23ae5d24fe0",
    ),
)

EXPECTED_GLOBAL_MINIMA = (
    (15, 3),
    (16, 2),
    (17, 8),
    (18, 2),
    (19, 9),
    (20, 2),
    (21, 3),
    (22, 2),
    (23, 6),
    (24, 2),
    (25, 5),
    (26, 2),
    (27, 3),
    (28, 2),
)
EXPECTED_SEMANTIC_DIGEST = "97d17ab3b16bb595b4cf765e30365b0afe37d916a1eaf3291dd635212471e578"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def divisors(value):
    return tuple(divisor for divisor in range(1, value + 1) if value % divisor == 0)


def direct_mask(q, owner, centre):
    mask = 0
    for sheet in range(q):
        phase = (owner * (centre + Fraction(sheet, q))) % 1
        distance = min(phase, 1 - phase)
        if 14 * distance < 1:
            mask |= 1 << sheet
    return mask


def containing_single_atom_mode_centre(q, owner, block, physical_time):
    order = q // gcd(q, owner)
    residue_classes = {sheet % order for sheet in block}
    require(len(residue_classes) == 1, (q, owner, block, residue_classes))
    atom = next(iter(residue_classes))
    scaled = owner * (physical_time + Fraction(atom, q))
    floor_scaled = scaled.numerator // scaled.denominator
    candidates = []
    for tooth in range(floor_scaled - 2, floor_scaled + 3):
        centre = Fraction(tooth, owner) - Fraction(atom, q)
        if abs(physical_time - centre) < Fraction(1, 14 * owner):
            candidates.append(centre)
    require(len(candidates) == 1, (q, owner, block, physical_time, candidates))
    return candidates[0]


def chart_masks(q, kind, g):
    require(q % g == 0 and g < q, (q, kind, g))
    quotient = q // g
    modulus = q if kind == "O" else 2 * q
    owner_period = modulus
    representatives = {}
    for v in range(1, owner_period):
        if v % quotient == 0:
            continue
        mask = 0
        for sheet in range(q):
            affine = 1 + (g if kind == "O" else 2 * g) * sheet
            residue = v * affine % modulus
            if 14 * min(residue, modulus - residue) < modulus:
                mask |= 1 << sheet
        if mask:
            representatives[mask] = min(v, representatives.get(mask, v))

    maximal = {}
    masks = tuple(representatives)
    for mask in masks:
        if any(mask != other and mask | other == other for other in masks):
            continue
        maximal[mask] = representatives[mask]
    return tuple(sorted(maximal.items(), key=lambda item: (item[1], item[0])))


def greedy_cover(full, masks):
    covered = 0
    witness = []
    while covered != full:
        candidate = max(
            (item for item in masks if item[0] | covered != covered),
            key=lambda item: ((item[0] | covered).bit_count() - covered.bit_count(), -item[1]),
            default=None,
        )
        require(candidate is not None, (full, masks, covered))
        witness.append(candidate)
        covered |= candidate[0]
    return witness


def minimum_cover(q, mask_items):
    full = (1 << q) - 1
    require(mask_items, (q, "empty mask bank"))
    union = 0
    for mask, _ in mask_items:
        union |= mask
    require(union == full, (q, "chart noncover", mask_items))

    greedy = greedy_cover(full, mask_items)
    best = [len(greedy), tuple(greedy)]
    seen_depth = {}
    nodes = [0]

    coverers = tuple(
        tuple(item for item in mask_items if item[0] & (1 << sheet))
        for sheet in range(q)
    )

    def search(covered, chosen):
        nodes[0] += 1
        depth = len(chosen)
        if depth >= best[0]:
            return
        old = seen_depth.get(covered)
        if old is not None and old <= depth:
            return
        seen_depth[covered] = depth
        if covered == full:
            best[:] = [depth, tuple(chosen)]
            return

        uncovered = full ^ covered
        gains = tuple((mask & uncovered).bit_count() for mask, _ in mask_items)
        max_gain = max(gains)
        if depth + (uncovered.bit_count() + max_gain - 1) // max_gain >= best[0]:
            return

        uncovered_sheets = tuple(sheet for sheet in range(q) if uncovered & (1 << sheet))
        pivot = min(
            uncovered_sheets,
            key=lambda sheet: sum(1 for mask, _ in coverers[sheet] if mask | covered != covered),
        )
        candidates = tuple(
            sorted(
                (item for item in coverers[pivot] if item[0] | covered != covered),
                key=lambda item: (-(item[0] & uncovered).bit_count(), item[1]),
            )
        )
        for item in candidates:
            search(covered | item[0], chosen + (item,))

    # A solution of size best-1 must be found to improve the greedy bound.
    search(0, ())
    rank, witness = best
    require(rank == len(witness), (q, rank, witness))
    covered = 0
    for mask, _ in witness:
        covered |= mask
    require(covered == full, (q, "witness cover", witness))
    return rank, witness, nodes[0], len(seen_depth)


def combination_cover(q, mask_items):
    full = (1 << q) - 1
    tested = 0
    for rank in range(1, len(mask_items) + 1):
        for chosen in combinations(mask_items, rank):
            tested += 1
            covered = 0
            for mask, _ in chosen:
                covered |= mask
            if covered == full:
                return rank, chosen, tested
    raise RuntimeError((q, "combination route found no cover"))


def raw_affine_mask(q, modulus, offset, slope, v, sheet_map=None):
    mask = 0
    for sheet in range(q):
        source_sheet = sheet if sheet_map is None else sheet_map(sheet)
        residue = v * (offset + slope * source_sheet) % modulus
        if 14 * min(residue, modulus - residue) < modulus:
            mask |= 1 << sheet
    return mask


def find_normalization(q, modulus, offset, slope, g):
    for shift in range(q):
        constant = (offset + slope * shift) % modulus
        if gcd(constant, modulus) != 1:
            continue
        scalar = pow(constant, -1, modulus)
        for multiplier in range(q):
            if gcd(multiplier, q) != 1:
                continue
            if scalar * slope * multiplier % modulus == (g if modulus == q else 2 * g) % modulus:
                return shift, multiplier, scalar
    raise RuntimeError((q, modulus, offset, slope, g, "normalization"))


def normalization_audit(q):
    odd_charts = 0
    odd_owner_rows = 0
    even_charts = 0
    even_owner_rows = 0

    odd_slopes = tuple(
        slope
        for slope in range(q)
        if q % 2 or slope % 2 == 1
    )
    for slope in odd_slopes:
        g = gcd(slope, q)
        if g == q:
            continue
        quotient = q // g
        for offset in range(q):
            if gcd(offset, slope, q) != 1:
                continue
            shift, multiplier, scalar = find_normalization(q, q, offset, slope, g)
            inverse_scalar = pow(scalar, -1, q)
            odd_charts += 1
            for v in range(1, q):
                if slope * v % q == 0:
                    continue
                normalized_v = v * inverse_scalar % q
                require(normalized_v % quotient, (q, "odd transverse transport"))
                original = raw_affine_mask(
                    q,
                    q,
                    offset,
                    slope,
                    v,
                    lambda sheet, shift=shift, multiplier=multiplier: (shift + multiplier * sheet) % q,
                )
                canonical = raw_affine_mask(q, q, 1, g, normalized_v)
                require(original == canonical, (q, "odd normalization", offset, slope, v))
                odd_owner_rows += 1

    for d in range(q):
        g = gcd(d, q)
        if g == q:
            continue
        quotient = q // g
        slope = 2 * d
        for offset in range(1, 2 * q, 2):
            if gcd(offset, d, q) != 1:
                continue
            shift, multiplier, scalar = find_normalization(q, 2 * q, offset, slope, g)
            inverse_scalar = pow(scalar, -1, 2 * q)
            even_charts += 1
            for v in range(1, 2 * q):
                if d * v % q == 0:
                    continue
                normalized_v = v * inverse_scalar % (2 * q)
                require(normalized_v % quotient, (q, "even transverse transport"))
                original = raw_affine_mask(
                    q,
                    2 * q,
                    offset,
                    slope,
                    v,
                    lambda sheet, shift=shift, multiplier=multiplier: (shift + multiplier * sheet) % q,
                )
                canonical = raw_affine_mask(q, 2 * q, 1, 2 * g, normalized_v)
                require(original == canonical, (q, "even normalization", offset, d, v))
                even_owner_rows += 1

    return q, odd_charts, odd_owner_rows, even_charts, even_owner_rows


def universal_block_cap(q):
    return max(
        (q // order) * ((order + 6) // 7)
        for order in divisors(q)
        if order > 1
    )


def divisor_pullback_record(q, prime):
    require(q % prime == 0, (q, prime))
    d = q // prime
    centre = Fraction(1, 2 * d * q)
    if prime == 2:
        require(q >= 8 and q % 2 == 0, (q, prime))
        v_values = (1, q - 1)
    elif prime == 3:
        require(q >= 9 and q % 2 == 1, (q, prime))
        if q == 9:
            v_values = (1, 5, 7)
        elif d % 3 == 0:
            v_values = (1, 2 * d - 1, 2 * d - 2)
        elif d % 3 == 1:
            v_values = (1, 2 * d, 2 * d - 1)
        else:
            v_values = (1, 2 * d - 2, 2 * d)
    else:
        require(prime == 5 and q >= 25 and q % 2 == 1 and q % 3, (q, prime))
        v_values = (1,) + tuple(
            value for value in range(2 * d - 2, 2 * d + 3) if value % 5
        )

    require(len(v_values) == prime and len(set(v_values)) == prime, (q, prime, v_values))
    owners = tuple(d * value for value in v_values)
    require(all(owner % q for owner in owners), (q, prime, owners))
    masks = tuple(direct_mask(q, owner, centre) for owner in owners)
    full = (1 << q) - 1
    require(all(mask.bit_count() == d for mask in masks), (q, prime, "block size"))
    require(sum(mask.bit_count() for mask in masks) == q, (q, prime, "mass"))
    union = 0
    for mask in masks:
        require(not union & mask, (q, prime, "overlap"))
        union |= mask
    require(union == full, (q, prime, "partition"))
    require(universal_block_cap(q) == d, (q, prime, universal_block_cap(q), d))
    residue_classes = tuple(
        tuple(sorted({sheet % prime for sheet in range(q) if mask & (1 << sheet)}))
        for mask in masks
    )
    require(sorted(item[0] for item in residue_classes) == list(range(prime)), (q, prime, residue_classes))
    return q, prime, centre, owners, residue_classes


def even_parity_leakage_record(q, scalar):
    """Recover the exact mode sidecar on one rank-two parity partition."""
    require(q >= 8 and q % 2 == 0, (q, scalar))
    require(scalar >= 1 and scalar % 2 == 1 and 7 * scalar < q, (q, scalar))
    half = q // 2
    physical_time = Fraction(scalar, q * q)
    owners = (half, half * (q - 1))
    masks = tuple(direct_mask(q, owner, physical_time) for owner in owners)
    full = (1 << q) - 1
    require(tuple(mask.bit_count() for mask in masks) == (half, half), (q, scalar, masks))
    require(not masks[0] & masks[1] and masks[0] | masks[1] == full, (q, scalar, masks))
    blocks = tuple(
        tuple(sheet for sheet in range(q) if mask & (1 << sheet))
        for mask in masks
    )
    require(
        tuple(tuple(sorted({sheet % 2 for sheet in block})) for block in blocks)
        == ((0,), (1,)),
        (q, scalar, blocks),
    )
    mode_centres = tuple(
        containing_single_atom_mode_centre(q, owner, block, physical_time)
        for owner, block in zip(owners, blocks)
    )
    expected_centres = (Fraction(0), Fraction(scalar, q * (q - 1)))
    require(mode_centres == expected_centres, (q, scalar, mode_centres))
    cochain = 2 * q * owners[0] * owners[1] * (mode_centres[0] - mode_centres[1])
    require(cochain == Fraction(-scalar * q * q, 2), (q, scalar, cochain))
    require(cochain.denominator == 1, (q, scalar, cochain))

    owner_gcd = gcd(*owners)
    gauge_gcd = gcd(q, owner_gcd)
    gauge_scalar = Fraction(2 * q * owner_gcd) * physical_time
    require(gauge_scalar == scalar, (q, scalar, gauge_scalar))
    require(gauge_scalar.numerator % gauge_gcd, (q, scalar, gauge_gcd))

    # The next odd scalar crosses the strict source-radius wall for owner q/2.
    boundary_scalar = scalar + 2
    if 7 * boundary_scalar >= q:
        boundary_mask = direct_mask(q, owners[0], Fraction(boundary_scalar, q * q))
        require(boundary_mask == 0, (q, scalar, boundary_scalar, boundary_mask))

    return (
        q,
        scalar,
        physical_time,
        owners,
        mode_centres,
        cochain.numerator,
        owner_gcd,
        gauge_gcd,
        gauge_scalar.numerator,
    )


def chart_record(q, kind, g):
    mask_items = chart_masks(q, kind, g)
    rank, witness, nodes, states = minimum_cover(q, mask_items)
    combination_rank, combination_witness, combination_tests = combination_cover(q, mask_items)
    require(combination_rank == rank, (q, kind, g, rank, combination_rank))
    centre = Fraction(1, g * q) if kind == "O" else Fraction(1, 2 * g * q)
    owner_scale = g
    owners = tuple(owner_scale * v for _, v in witness)
    blocks = tuple(
        tuple(sheet for sheet in range(q) if mask & (1 << sheet))
        for mask, _ in witness
    )
    require(all(owner % q for owner in owners), (q, kind, g, owners))
    require(
        tuple(direct_mask(q, owner, centre) for owner in owners)
        == tuple(mask for mask, _ in witness),
        (q, kind, g, "direct block mismatch"),
    )
    require(set().union(*(set(block) for block in blocks)) == set(range(q)), (q, kind, g))
    return (
        kind,
        g,
        q // g,
        rank,
        len(mask_items),
        centre,
        Fraction(q) * centre % 1,
        owners,
        blocks,
        nodes,
        states,
        combination_tests,
        tuple(v for _, v in combination_witness),
    )


def all_charts(q):
    records = []
    for g in divisors(q):
        if g == q:
            continue
        if g % 2 == 1:
            records.append(chart_record(q, "O", g))
        records.append(chart_record(q, "E", g))
    return tuple(records)


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

    records = tuple((q, all_charts(q)) for q in range(15, 29))
    minima = tuple(
        (
            q,
            min(record[3] for record in charts),
            tuple((record[0], record[1]) for record in charts if record[3] == min(item[3] for item in charts)),
        )
        for q, charts in records
    )
    rank_table = tuple((q, rank) for q, rank, _ in minima)
    require(rank_table == EXPECTED_GLOBAL_MINIMA, ("global ranks", rank_table))

    normalization = tuple(normalization_audit(q) for q in range(15, 29))
    capacity = tuple(
        (
            q,
            universal_block_cap(q),
            (q + universal_block_cap(q) - 1) // universal_block_cap(q),
            rank,
        )
        for q, rank, _ in minima
    )
    capacity_sharp_support = tuple(q for q, _, lower, rank in capacity if lower == rank)
    require(
        capacity_sharp_support == (15, 16, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28),
        capacity_sharp_support,
    )

    # At q=17,19 the odd chart is THM-3401's fixed-zero sign-pair cover.
    # The even chart has one common central sheet plus disjoint noncentral
    # pairs, forcing and attaining (q-1)/2 owners.
    short_prime_pair_controls = []
    for q in (17, 19):
        charts = dict(((record[0], record[1]), record) for record in dict(records)[q])
        for key in (("O", 1), ("E", 1)):
            record = charts[key]
            masks = tuple(mask for mask, _ in chart_masks(q, *key))
            common = set(range(q))
            for mask in masks:
                common &= {sheet for sheet in range(q) if mask & (1 << sheet)}
            require(len(masks) == (q - 1) // 2, (q, key, len(masks)))
            require(len(common) == 1, (q, key, common))
            require(all(mask.bit_count() == 3 for mask in masks), (q, key, masks))
            require(
                all((left & right).bit_count() == 1 for left, right in combinations(masks, 2)),
                (q, key, "pair overlap"),
            )
            short_prime_pair_controls.append((q, key, tuple(sorted(common)), len(masks)))
    short_prime_pair_controls = tuple(short_prime_pair_controls)

    # MISTAKE-389 hostile: this exact physical half-grid partition is not a
    # zero-mode-cochain certificate.  THM-3405 requires gcd(q,d)|a, where
    # d is the active owner gcd and a=2*q*d*c.
    hostile_record = next(
        record
        for record in dict(records)[15]
        if (record[0], record[1]) == ("E", 5)
    )
    hostile_centre = hostile_record[5]
    hostile_owners = hostile_record[7]
    hostile_owner_gcd = gcd(*hostile_owners)
    hostile_gauge_gcd = gcd(15, hostile_owner_gcd)
    hostile_scalar = Fraction(2 * 15 * hostile_owner_gcd) * hostile_centre
    require(hostile_scalar.denominator == 1, ("half-grid scalar", hostile_scalar))
    require(
        hostile_scalar.numerator % hostile_gauge_gcd,
        ("hostile unexpectedly zero-mode", hostile_scalar, hostile_gauge_gcd),
    )
    mode_centre_divisibility_hostile = (
        15,
        ("E", 5),
        hostile_centre,
        hostile_owners,
        hostile_owner_gcd,
        hostile_gauge_gcd,
        hostile_scalar,
    )
    hostile_blocks = hostile_record[8]
    hostile_mode_centres = tuple(
        containing_single_atom_mode_centre(
            15, owner, block, hostile_centre
        )
        for owner, block in zip(hostile_owners, hostile_blocks)
    )
    hostile_cochain = tuple(
        (
            left,
            right,
            2
            * 15
            * hostile_owners[left]
            * hostile_owners[right]
            * (hostile_mode_centres[left] - hostile_mode_centres[right]),
        )
        for left, right in combinations(range(len(hostile_owners)), 2)
    )
    require(
        all(value.denominator == 1 for _, _, value in hostile_cochain),
        hostile_cochain,
    )
    hostile_cochain = tuple(
        (left, right, value.numerator) for left, right, value in hostile_cochain
    )
    hostile_cochain_values = {
        (left, right): value for left, right, value in hostile_cochain
    }
    require(
        hostile_owners[2] * hostile_cochain_values[(0, 1)]
        + hostile_owners[0] * hostile_cochain_values[(1, 2)]
        - hostile_owners[1] * hostile_cochain_values[(0, 2)]
        == 0,
        ("triangle closure", hostile_cochain),
    )
    hostile_l1 = sum(abs(value) for _, _, value in hostile_cochain)
    hostile_linf = max(abs(value) for _, _, value in hostile_cochain)
    require(
        (hostile_mode_centres, hostile_cochain, hostile_l1, hostile_linf)
        == (
            (Fraction(0), Fraction(1, 120), Fraction(1, 150)),
            ((0, 1, -50), (0, 2, -50), (1, 2, 100)),
            200,
            100,
        ),
        (hostile_mode_centres, hostile_cochain, hostile_l1, hostile_linf),
    )
    mode_centre_drift_hostile = (
        hostile_mode_centres,
        hostile_cochain,
        hostile_l1,
        hostile_linf,
    )

    divisor_pullback_records = (
        tuple(divisor_pullback_record(q, 2) for q in range(8, 501, 2)),
        tuple(divisor_pullback_record(q, 3) for q in range(9, 501, 6)),
        tuple(
            divisor_pullback_record(q, 5)
            for q in range(25, 501)
            if q % 2 == 1 and q % 5 == 0 and q % 3
        ),
    )
    divisor_pullback_audit = tuple(
        (
            prime,
            len(family),
            family[0][0],
            family[-1][0],
            sha256(repr(family).encode("ascii")).hexdigest(),
        )
        for prime, family in zip((2, 3, 5), divisor_pullback_records)
    )

    even_parity_leakage_records = tuple(
        even_parity_leakage_record(q, scalar)
        for q in range(8, 501, 2)
        for scalar in range(1, q, 2)
        if 7 * scalar < q
    )
    require(len(even_parity_leakage_records) == 4482, len(even_parity_leakage_records))
    canonical_even_leakage = tuple(
        record for record in even_parity_leakage_records if record[1] == 1
    )
    require(len(canonical_even_leakage) == 247, len(canonical_even_leakage))
    require(
        tuple((record[0], record[2], record[3]) for record in canonical_even_leakage)
        == tuple((record[0], record[2], record[3]) for record in divisor_pullback_records[0]),
        "canonical even leakage/pullback mismatch",
    )
    even_parity_leakage_audit = (
        len(even_parity_leakage_records),
        even_parity_leakage_records[0],
        even_parity_leakage_records[-1],
        max(abs(record[5]) for record in even_parity_leakage_records),
        sha256(repr(even_parity_leakage_records).encode("ascii")).hexdigest(),
    )

    semantic = sha256(
        repr(
            (
                records,
                minima,
                normalization,
                capacity,
                capacity_sharp_support,
                short_prime_pair_controls,
                mode_centre_divisibility_hostile,
                mode_centre_drift_hostile,
                divisor_pullback_records,
                even_parity_leakage_records,
            )
        ).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_DIGEST:
        require(semantic == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", semantic))

    print("LRC ALL-OWNER SYNCHRONIZED HALF-GRID PHYSICAL AFFINE-CHART PROBE")
    print(f"source_sha256_lf={lf_hash(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINNED)}")
    print("status=PROVED-ELEMENTARY all_q_half_grid_physical_affine_chart_reduction;exact_all_positive_owner_half_grid_ranks_q15..28;all_q_half_grid_rank_families_even=2_odd3multiple=3_odd5multiple_not3=5;exact_even_rank2_mode_leakage_P=-a*q^2/2;INDEPENDENT_COMBINATION_AND_NORMALIZATION_AUDITS;MISTAKE389_zero_mode_cochain_interpretation_retracted")
    print("odd_chart=theta=a/b,b_odd,u=bv;normalize_to_v(1+g*ell)_mod_q,g=gcd(b,q)_odd")
    print("even_chart=theta=a/(2d),u=dv;normalize_to_v(1+2g*ell)_mod_2q,g=gcd(d,q)")
    print("gauge=sheet_affine_permutation_plus_unit_owner_reindex;transverse_v_not_divisible_by_q/g")
    print(f"half_grid_minima_q15_q28=(q,rank,minimizing_charts)={minima}")
    print(f"universal_block_capacity=(q,max_block,rank_lower_bound,exact_rank)={capacity};sharp_support={capacity_sharp_support}")
    print(f"q17_q19_central_plus_pair_controls={short_prime_pair_controls}")
    print(f"mode_centre_divisibility_hostile=(q,chart,c,owners,d,g,a)={mode_centre_divisibility_hostile}")
    print(f"mode_centre_drift_hostile=(centres,pairs,L1,Linf)={mode_centre_drift_hostile}")
    print(f"divisor_pullback_family_audit_q_through_500=(prime,count,first_q,last_q,sha256)={divisor_pullback_audit}")
    print("even_rank2_mode_leakage_formula=q_even>=8,a_odd,1<=a<q/7;c=a/q^2;owners=(q/2,q(q-1)/2);centres=(0,a/(q(q-1)));P=-a*q^2/2;THM3405_g=q/2_not_dividing_a")
    print(f"even_rank2_mode_leakage_audit_q_through_500=(count,first,last,max_abs_P,sha256)={even_parity_leakage_audit}")
    print(f"normalization_audit=(q,odd_charts,odd_owner_rows,even_charts,even_owner_rows)={normalization}")
    for q, charts in records:
        print(f"q={q};charts={charts}")
    print("scope=all_positive_transverse_owners;maximal_actual_blocks_at_synchronized_physical_times_with_2quc_integral;NOT_zero_mode_cochain_rank;NOT_mobile_common_mode_centre_rank;minimum_rank_not_literal_certificate_count;no_LRC14_decrement")
    print(f"semantic_sha256={semantic}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()

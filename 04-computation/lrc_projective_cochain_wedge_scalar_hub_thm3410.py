#!/usr/bin/env python3
"""Exact companion for THM-3410.

The program has two independent jobs.

1.  It treats a realized THM-3398 cochain as the determinant system of the
    integral columns ``(u_i,A_i)``, contracts equal projective rays, and checks
    the exact minimum-spanning-tree and bottleneck-tree tariffs on a bounded
    hostile lattice universe.
2.  Starting from the original strict danger inequality, it reconstructs the
    scalar fibres over the canonical two-, three-, and five-colour half-grid
    partitions through q=40.  It does not feed predicted masks or centres into
    the direct physical test.

All arithmetic is integral or ``Fraction`` arithmetic.  There are no random
choices, floating literals, external solvers, network calls, or truth-bearing
``assert`` statements.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINS = (
    (
        "THM-3398",
        ROOT / "01-canon/theorems/THM-3398-general-finite-mode-sheet-cover-cochain.md",
        "01901da2bb382184cfe4466550afe79255598f580f00a761fc32731a52ec9378",
    ),
    (
        "THM-3405",
        ROOT / "01-canon/theorems/THM-3405-common-centre-gcd-gauge-and-boolean-half-twist.md",
        "d3e7dbeeb85c6f897bd9e31270bd0b6602ae4feac3b46a45eb5ce23ae5d24fe0",
    ),
    (
        "THM-3409",
        ROOT / "01-canon/theorems/THM-3409-q15-exceptional-edge-positive-cochain-rigidity-and-leakage-tariff.md",
        "d114beaf2eed6349588a0db0a135a9c70c5a215fcda46ee6555dbccb33566bff",
    ),
    (
        "THM-3409-output",
        ROOT / "05-knowledge/results/lrc_q15_exceptional_positive_cochain_tariff_thm3409.out",
        "ef1a0fa0f17411046f77c3c90fbe71009c6cddb14a91ae1df647e1a8ff2e1b20",
    ),
)
EXPECTED_SEMANTIC_DIGEST = "bc6cc2bd14fb566ee7e57adc27291e3dea1085cf542891128cc4feea93d44c51"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


def digest_repr(value):
    return sha256(repr(value).encode("ascii")).hexdigest()


def pair_order(size):
    return tuple(combinations(range(size), 2))


def wedge_values(speeds, lifts):
    return tuple(
        (left, right, lifts[left] * speeds[right] - speeds[left] * lifts[right])
        for left, right in pair_order(len(speeds))
    )


def edge_dictionary(values):
    return {(left, right): value for left, right, value in values}


def weighted_triangle_closure(speeds, values):
    edges = edge_dictionary(values)
    for left, middle, right in combinations(range(len(speeds)), 3):
        residue = (
            speeds[right] * edges[(left, middle)]
            - speeds[middle] * edges[(left, right)]
            + speeds[left] * edges[(middle, right)]
        )
        require(residue == 0, ("triangle closure", speeds, values, left, middle, right, residue))


def components_at_threshold(vertex_count, edges, threshold):
    parent = list(range(vertex_count))

    def find(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for left, right, weight in edges:
        if weight > threshold:
            continue
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[root_left] = root_right
    return len({find(vertex) for vertex in range(vertex_count)})


def kruskal_tariff(vertex_count, weighted_edges):
    if vertex_count <= 1:
        return (0, 0, ())
    parent = list(range(vertex_count))

    def find(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    chosen = []
    for left, right, weight in sorted(weighted_edges, key=lambda item: (item[2], item[0], item[1])):
        root_left = find(left)
        root_right = find(right)
        if root_left == root_right:
            continue
        parent[root_left] = root_right
        chosen.append((left, right, weight))
        if len(chosen) == vertex_count - 1:
            break
    require(len(chosen) == vertex_count - 1, ("disconnected graph", vertex_count, weighted_edges))
    thresholds = sorted({weight for _, _, weight in weighted_edges})
    bottleneck = next(
        threshold
        for threshold in thresholds
        if components_at_threshold(vertex_count, weighted_edges, threshold) == 1
    )
    require(max(weight for _, _, weight in chosen) == bottleneck, ("bottleneck mismatch", chosen, bottleneck))
    return (sum(weight for _, _, weight in chosen), bottleneck, tuple(chosen))


def projective_record(speeds, lifts):
    require(len(speeds) == len(lifts) and speeds and all(speed > 0 for speed in speeds), (speeds, lifts))
    values = wedge_values(speeds, lifts)
    weighted_triangle_closure(speeds, values)
    edge_map = edge_dictionary(values)
    contents = tuple(gcd(speed, abs(lift)) for speed, lift in zip(speeds, lifts))
    primitive_rays = tuple(
        (speed // content, lift // content)
        for speed, lift, content in zip(speeds, lifts, contents)
    )
    for left, right in pair_order(len(speeds)):
        determinant = (
            primitive_rays[left][1] * primitive_rays[right][0]
            - primitive_rays[left][0] * primitive_rays[right][1]
        )
        require(
            edge_map[(left, right)] == contents[left] * contents[right] * determinant,
            ("primitive factorization", speeds, lifts, left, right),
        )

    classes_by_ray = {}
    for index, ray in enumerate(primitive_rays):
        classes_by_ray.setdefault(ray, []).append(index)
    rays = tuple(sorted(classes_by_ray))
    classes = tuple(tuple(classes_by_ray[ray]) for ray in rays)
    class_contents = tuple(min(contents[index] for index in block) for block in classes)
    quotient_edges = []
    for left, right in pair_order(len(classes)):
        ray_left = rays[left]
        ray_right = rays[right]
        determinant = abs(ray_left[1] * ray_right[0] - ray_left[0] * ray_right[1])
        require(determinant > 0, ("distinct ray determinant", ray_left, ray_right))
        formula_weight = class_contents[left] * class_contents[right] * determinant
        direct_weight = min(
            abs(edge_map[tuple(sorted((i, j)))])
            for i in classes[left]
            for j in classes[right]
        )
        require(formula_weight == direct_weight, ("quotient edge", formula_weight, direct_weight))
        quotient_edges.append((left, right, formula_weight))

    raw_edges = tuple((left, right, abs(value)) for left, right, value in values)
    raw_tariff = kruskal_tariff(len(speeds), raw_edges)
    quotient_tariff = kruskal_tariff(len(classes), tuple(quotient_edges))
    require(raw_tariff[:2] == quotient_tariff[:2], ("ray contraction", speeds, lifts, raw_tariff, quotient_tariff))

    if len(classes) >= 2:
        ordered_contents = tuple(sorted(class_contents))
        lower_l1 = ordered_contents[0] * sum(ordered_contents[1:])
        lower_linf = ordered_contents[0] * ordered_contents[-1]
        require(raw_tariff[0] >= lower_l1, ("L1 content bound", raw_tariff, lower_l1))
        require(raw_tariff[1] >= lower_linf, ("Linf content bound", raw_tariff, lower_linf))
    else:
        lower_l1 = lower_linf = 0
    return (
        values,
        contents,
        primitive_rays,
        classes,
        class_contents,
        raw_tariff[:2],
        (lower_l1, lower_linf),
    )


def synthetic_projective_audit():
    cases = 0
    zero_ray_cases = 0
    multi_ray_cases = 0
    equality_l1 = 0
    equality_linf = 0
    summary = []
    for size in range(2, 6):
        for speeds in combinations(range(1, 8), size):
            for lifts in product(range(-2, 3), repeat=size):
                record = projective_record(speeds, lifts)
                cases += 1
                class_count = len(record[3])
                zero_ray_cases += class_count == 1
                multi_ray_cases += class_count > 1
                equality_l1 += class_count > 1 and record[5][0] == record[6][0]
                equality_linf += class_count > 1 and record[5][1] == record[6][1]
                if cases in (1, 525, 4900, 26775, 92400):
                    summary.append((cases, speeds, lifts, record[3], record[5], record[6]))
    require(cases == 92400, cases)
    equal_ray_hostile = projective_record((2, 4, 3, 6), (1, 2, -1, -2))
    require(len(equal_ray_hostile[3]) == 2, equal_ray_hostile)
    require(equal_ray_hostile[5] == (5, 5), equal_ray_hostile)
    return (
        cases,
        zero_ray_cases,
        multi_ray_cases,
        equality_l1,
        equality_linf,
        tuple(summary),
        equal_ray_hostile,
    )


def direct_danger_mask(q, owner, physical_time):
    mask = 0
    for sheet in range(q):
        phase = (owner * (physical_time + Fraction(sheet, q))) % 1
        if 14 * min(phase, 1 - phase) < 1:
            mask |= 1 << sheet
    return mask


def residue_mask(q, prime, residue):
    return sum(1 << sheet for sheet in range(q) if sheet % prime == residue)


def containing_single_atom_centre(q, owner, block, physical_time):
    order = q // gcd(q, owner)
    residues = {sheet % order for sheet in block}
    require(len(residues) == 1, ("not a single quotient atom", q, owner, block, residues))
    atom = next(iter(residues))
    scaled = owner * (physical_time + Fraction(atom, q))
    floor_scaled = scaled.numerator // scaled.denominator
    candidates = []
    for tooth in range(floor_scaled - 2, floor_scaled + 3):
        centre = Fraction(tooth, owner) - Fraction(atom, q)
        if abs(physical_time - centre) < Fraction(1, 14 * owner):
            candidates.append(centre)
    require(len(candidates) == 1, ("containing centre", q, owner, block, physical_time, candidates))
    return candidates[0]


def hub_owner_word(prime, degree):
    if prime == 2:
        require(degree >= 4, (prime, degree))
        return (1, 2 * degree - 1)
    if prime == 3:
        require(degree == 3 or (degree >= 5 and degree % 2 == 1), (prime, degree))
        if degree == 3:
            return (1, 5, 7)
        if degree % 3 == 0:
            return (1, 2 * degree - 1, 2 * degree - 2)
        if degree % 3 == 1:
            return (1, 2 * degree, 2 * degree - 1)
        return (1, 2 * degree - 2, 2 * degree)
    require(prime == 5 and degree >= 5 and degree % 2 == 1 and degree % 3 != 0, (prime, degree))
    return (1,) + tuple(
        value
        for value in range(2 * degree - 2, 2 * degree + 3)
        if value % 5 != 0
    )


def hub_delta(word, degree):
    return max((1,) + tuple(abs(2 * degree - value) for value in word[1:]))


def prufer_edges(word, vertex_count):
    if vertex_count == 2:
        return ((0, 1),)
    degree = [1] * vertex_count
    for vertex in word:
        degree[vertex] += 1
    edges = []
    for vertex in word:
        leaf = next(index for index in range(vertex_count) if degree[index] == 1)
        edges.append(tuple(sorted((leaf, vertex))))
        degree[leaf] -= 1
        degree[vertex] -= 1
    leaves = tuple(index for index in range(vertex_count) if degree[index] == 1)
    require(len(leaves) == 2, (word, degree, leaves))
    edges.append(tuple(sorted(leaves)))
    return tuple(edges)


def exhaustive_tree_tariff(values, vertex_count):
    edge_map = {(left, right): abs(value) for left, right, value in values}
    totals = []
    maxima = []
    for word in product(range(vertex_count), repeat=max(0, vertex_count - 2)):
        edges = prufer_edges(word, vertex_count)
        weights = tuple(edge_map[edge] for edge in edges)
        totals.append(sum(weights))
        maxima.append(max(weights) if weights else 0)
    return (
        vertex_count ** max(0, vertex_count - 2),
        min(totals),
        totals.count(min(totals)),
        min(maxima),
        maxima.count(min(maxima)),
    )


def scalar_hub_record(prime, degree, scalar):
    q = prime * degree
    word = hub_owner_word(prime, degree)
    delta = hub_delta(word, degree)
    require(gcd(scalar, prime) == 1 and 7 * scalar * delta < q, ("inadmissible hub", prime, degree, scalar))
    require(len(word) == prime and len(set(word)) == prime, (prime, degree, word))
    require({value % prime for value in word[1:]} == set(range(1, prime)), (prime, degree, word))
    owners = tuple(degree * value for value in word)
    physical_time = Fraction(scalar, 2 * degree * q)
    target_residues = (0,) + tuple(
        (-scalar * pow(value, -1, prime)) % prime
        for value in word[1:]
    )
    require(set(target_residues) == set(range(prime)), (prime, degree, scalar, target_residues))

    blocks = []
    centres = []
    for owner, target in zip(owners, target_residues):
        mask = direct_danger_mask(q, owner, physical_time)
        expected_mask = residue_mask(q, prime, target)
        require(mask == expected_mask, ("direct residue partition", prime, degree, scalar, owner, target, mask, expected_mask))
        block = tuple(sheet for sheet in range(q) if mask & (1 << sheet))
        blocks.append(block)
        centres.append(containing_single_atom_centre(q, owner, block, physical_time))
    require(sum(len(block) for block in blocks) == q, (prime, degree, scalar, blocks))
    expected_centres = (Fraction(0),) + tuple(
        Fraction(scalar, prime * degree * value)
        for value in word[1:]
    )
    require(tuple(centres) == expected_centres, ("hub centres", prime, degree, scalar, centres, expected_centres))

    lifts = tuple(2 * q * owner * centre for owner, centre in zip(owners, centres))
    require(all(lift.denominator == 1 for lift in lifts), (prime, degree, scalar, lifts))
    lifts = tuple(lift.numerator for lift in lifts)
    expected_lifts = (0,) + (2 * scalar * degree,) * (prime - 1)
    require(lifts == expected_lifts, ("hub lifts", prime, degree, scalar, lifts, expected_lifts))
    physical_words = tuple(2 * q * owner * physical_time for owner in owners)
    require(all(value.denominator == 1 for value in physical_words), physical_words)
    errors = tuple(lift - value.numerator for lift, value in zip(lifts, physical_words))

    projective = projective_record(owners, lifts)
    values = projective[0]
    edge_map = edge_dictionary(values)
    base_tariff = 2 * scalar * degree * degree
    require(
        tuple(edge_map[(0, index)] for index in range(1, prime))
        == (-base_tariff,) * (prime - 1),
        ("hub root star", prime, degree, scalar, values),
    )
    for left, right in combinations(range(1, prime), 2):
        require(
            edge_map[(left, right)] == base_tariff * (word[right] - word[left]),
            ("hub horizontal determinant", prime, degree, scalar, left, right, values),
        )
    exhaustive = exhaustive_tree_tariff(values, prime)
    expected_tariffs = (base_tariff * (prime - 1), base_tariff)
    require(exhaustive[1] == expected_tariffs[0] and exhaustive[3] == expected_tariffs[1], ("hub tariffs", exhaustive, expected_tariffs))
    require(projective[5] == expected_tariffs, ("hub Kruskal tariffs", projective, expected_tariffs))
    return (
        prime,
        q,
        degree,
        scalar,
        delta,
        word,
        target_residues,
        lifts,
        errors,
        tuple(value for _, _, value in values),
        exhaustive,
        expected_tariffs,
    )


def scalar_hub_audit():
    records = []
    for prime in (2, 3, 5):
        for q in range(8, 41):
            if q % prime != 0:
                continue
            degree = q // prime
            if prime == 2 and degree < 4:
                continue
            if prime == 3 and not (degree == 3 or (degree >= 5 and degree % 2 == 1)):
                continue
            if prime == 5 and not (degree >= 5 and degree % 2 == 1 and degree % 3 != 0):
                continue
            word = hub_owner_word(prime, degree)
            delta = hub_delta(word, degree)
            for scalar in range(1, q + 1):
                if gcd(scalar, prime) == 1 and 7 * scalar * delta < q:
                    records.append(scalar_hub_record(prime, degree, scalar))
    require(len(records) == 44, len(records))
    counts = tuple((prime, sum(record[0] == prime for record in records)) for prime in (2, 3, 5))
    require(counts == ((2, 30), (3, 11), (5, 3)), counts)

    # A unit failure well inside the metric radius: all q=39 ternary owners
    # collapse onto one residue class when a is divisible by three.
    q = 39
    degree = 13
    word = hub_owner_word(3, degree)
    bad_unit_time = Fraction(3, 2 * degree * q)
    bad_unit_masks = tuple(direct_danger_mask(q, degree * value, bad_unit_time) for value in word)
    bad_unit_union = 0
    for mask in bad_unit_masks:
        bad_unit_union |= mask
    require(bad_unit_union.bit_count() == 13 < q, ("unit hostile", bad_unit_masks, bad_unit_union.bit_count()))

    # A strict-radius failure with a still prime to five.
    q = 35
    degree = 7
    word = hub_owner_word(5, degree)
    bad_radius_time = Fraction(3, 2 * degree * q)
    bad_radius_masks = tuple(direct_danger_mask(q, degree * value, bad_radius_time) for value in word)
    bad_radius_union = 0
    for mask in bad_radius_masks:
        bad_radius_union |= mask
    require(bad_radius_union.bit_count() == 21 < q, ("radius hostile", bad_radius_masks, bad_radius_union.bit_count()))

    return (
        tuple(records),
        counts,
        (39, 3, tuple(mask.bit_count() for mask in bad_unit_masks), bad_unit_union.bit_count()),
        (35, 3, tuple(mask.bit_count() for mask in bad_radius_masks), bad_radius_union.bit_count()),
    )


def q15_regressions():
    ternary = scalar_hub_record(3, 5, 1)
    require(ternary[9] == (-50, -50, 100), ternary)
    require(ternary[11] == (100, 50), ternary)

    speeds = (1, 7, 8, 9, 11, 13)
    values_flat = (2, 2, 3, 3, 4, -2, 3, -1, 2, 6, 2, 6, -6, -3, 5)
    values = tuple(
        (left, right, value)
        for (left, right), value in zip(pair_order(len(speeds)), values_flat)
    )
    lifts = (0,) + tuple(-values_flat[index] for index in range(5))
    require(lifts == (0, -2, -2, -3, -3, -4), lifts)
    require(wedge_values(speeds, lifts) == values, ("q15 exceptional wedge", speeds, lifts, values))
    projective = projective_record(speeds, lifts)
    exhaustive = exhaustive_tree_tariff(values, len(speeds))
    require(exhaustive == (1296, 10, 15, 3, 104), exhaustive)
    require(projective[5] == (10, 3) and projective[6] == (8, 3), projective)

    shear_records = []
    for shift in (-2, -1, 0, 1, 2):
        sheared = tuple(lift + 30 * shift * speed for speed, lift in zip(speeds, lifts))
        sheared_projective = projective_record(speeds, sheared)
        require(sheared_projective[0] == values, ("integer source-lift gauge", shift, sheared_projective[0]))
        require(sheared_projective[1] == projective[1], ("content gauge", shift, sheared_projective[1]))
        shear_records.append((shift, sheared_projective[1], sheared_projective[5]))

    dilation_records = []
    for scale in (1, 2, 3, 7):
        dilated_speeds = tuple(scale * speed for speed in speeds)
        dilated_lifts = tuple(scale * lift for lift in lifts)
        dilated = projective_record(dilated_speeds, dilated_lifts)
        require(
            dilated[5] == (scale * scale * 10, scale * scale * 3),
            ("dilation tariff", scale, dilated[5]),
        )
        dilation_records.append(
            (
                scale,
                dilated[5],
                Fraction(dilated[5][0], (15 * scale) ** 2),
                Fraction(dilated[5][1], (15 * scale) ** 2),
            )
        )
    require(
        {record[2:] for record in dilation_records} == {(Fraction(2, 45), Fraction(1, 75))},
        dilation_records,
    )
    return (ternary, lifts, projective, exhaustive, tuple(shear_records), tuple(dilation_records))


def main():
    for name, path, expected in PINS:
        require(lf_hash(path) == expected, ("dependency hash", name, path, lf_hash(path), expected))

    synthetic = synthetic_projective_audit()
    hubs = scalar_hub_audit()
    q15 = q15_regressions()
    hub_records = hubs[0]
    hub_digest = digest_repr(hub_records)
    semantic_payload = (
        synthetic,
        hubs,
        q15,
        "P=A_wedge_u",
        "ray_edge=c_i*c_j*primitive_determinant",
        "tau_hub=(2*a*(p-1)*d^2,2*a*d^2)",
        "q2_dilation_scale_q3_scalar_fibre_diameter",
    )
    semantic = digest_repr(semantic_payload)
    if EXPECTED_SEMANTIC_DIGEST != "__FILL_AFTER_FIRST_EXACT_RUN__":
        require(semantic == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", semantic, EXPECTED_SEMANTIC_DIGEST))

    first_last = tuple(
        (
            prime,
            next(record for record in hub_records if record[0] == prime),
            next(record for record in reversed(hub_records) if record[0] == prime),
        )
        for prime in (2, 3, 5)
    )
    print("THM-3410 PROJECTIVE COCHAIN WEDGE AND RESIDUE SCALAR HUBS")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINS)}")
    print("status=PROVED analytic integral-wedge/ray-tree theorem and p=2,3,5 scalar-hub family;VERIFIED-EXACT companion;INDEPENDENT_AUDIT_REQUESTED")
    print("wedge=A_i=2q*u_i*x_i in Z;P_ij=A_i*u_j-u_i*A_j;integer_source_lift_is_unimodular_shear")
    print("ray_quotient=c_i=gcd(u_i,A_i);r_i=(u_i,A_i)/c_i;edge_weight=c_C*c_D*abs(det(r_C,r_D));contract_equal_rays")
    print("tree_tariff=tau1_exact_MST_on_ray_quotient;tauinf_exact_connectivity_threshold;content_bounds=c_min*sum_other_c_and_c_min*c_max")
    print(f"synthetic_projective_audit=(cases,one_ray,multiray,L1_bound_equal,Linf_bound_equal)={synthetic[:5]}")
    print(f"equal_ray_contraction_hostile={synthetic[6]}")
    print("scalar_hub=p_in_2,3,5;q=p*d;gcd(a,p)=1;Delta=max(1,max_abs(2d-v));7aDelta<q;c=a/(2dq);A=(0,2ad,...,2ad)")
    print("scalar_hub_cochain=P_0j=-2ad^2;P_ij=2ad^2(v_j-v_i);tariffs=(2a(p-1)d^2,2ad^2)")
    print(f"scalar_hub_q8_q40=(count,counts_by_p,sha256)=({len(hub_records)},{hubs[1]},{hub_digest})")
    print(f"scalar_hub_first_last={first_last}")
    print(f"hostiles=(unit_failure_q39_a3,radius_failure_q35_a3)=({hubs[2]},{hubs[3]})")
    print(f"q15_halfgrid_regression=(P,tau)={(q15[0][9], q15[0][11])}")
    print(f"q15_exceptional=(A,contents,rays,tau,content_bounds,tree_enum)={(q15[1], q15[2][1], q15[2][2], q15[2][5], q15[2][6], q15[3])}")
    print(f"q15_source_lift_shears={q15[4]}")
    print(f"q15_pure_dilations=(scale,tau,tau1_over_q2,tauinf_over_q2)={q15[5]}")
    print("scale_boundary=fixed_projective_scalar_pure_dilation_is_q^2;admissible_scalar_fibres_have_a=Theta(q),so_uniform_family_diameter_is_Theta(q^3);q^-2_not_compact,q^-3_order_sharp")
    print("scope=realized_THM3398_packets_and_displayed_halfgrid_families;ray_graph_does_not_imply_cover;no_capped_owner_or_LRC14_ledger_decrement")
    print(f"semantic_sha256={semantic}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()

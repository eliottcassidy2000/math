#!/usr/bin/env python3
"""Exact all-D q=D/7 screen for the k=3/four-drift no-small residual.

The bounded q-profile scout uses a literal lcm-state multiset DP.  That is
an excellent independent control, but repeating it at all 217 resolving
denominators is unnecessary.  Here the exact-lcm coefficient is extracted
by divisor-lattice Mobius inversion.  For each downward alphabet E|D, its
denominators d>4 are grouped by the complete q-profile feature

    (ambient capacity, baseline load, increment, high-fibre count).

The arity-four multiset product is then propagated only in grouped feature
space.  Mobius inversion recovers the exact lcm-D distribution.  This never
materializes a denominator four-tuple.

For the canonical local modulus seven, q=D/7, there is a further exact
collapse.  If g=gcd(d,q), then d/g divides 7 and is coprime to q/g, hence
d/g is 1 or 7.  Consequently every denominator d>4 is exactly one of

* ``d|q``: a seven-spike of load 7 on
  R=(q/d)ceil(d/7) fibres and load 0 elsewhere;
* ``d not|q``: a uniform-one mask, of load 1 on every q-fibre.

Thus a four-mask profile is only ``(b,R)``, where b is the number of
uniform-one masks and R is the sum of the spike cardinalities.  Since a
q-fibre has at most seven points, the entire fractional-cover profile test
reduces to the one exact inequality

    #{target fibres of load > b} <= min(q,R).

The canonical generating function is therefore propagated in this smaller
feature space.  This inequality already implies scalar capacity: outside
the heavy set target load is at most b, while on it target load is at most
seven.  The implementation nevertheless writes the stage as an explicit
intersection with scalar capacity and audits this implication feature by
feature.  Counts for any extra ``--local-moduli`` are reported separately:
minima across screens are only single-screen upper bounds and are not
asserted to be their intersection.
"""

from argparse import ArgumentParser
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb, gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = Path(__file__).with_name(
    "lrc14_k3_four_drift_divisor_status_gf_thm2928.py"
)
EXPECTED_BASE_SHA256 = (
    "2fcd1fa7f122517feff3d3e0b3a21a6664fefaa12588e4db008572078989d6eb"
)
if sha256(BASE_PATH.read_bytes()).hexdigest() != EXPECTED_BASE_SHA256:
    raise RuntimeError("frozen four-drift base dependency changed")
spec = spec_from_file_location("lrc14_k3_status_gf", BASE_PATH)
base = module_from_spec(spec)
spec.loader.exec_module(base)
support = base.support
combined = base.combined


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def q_type(D, d, q):
    """Exact baseline/increment/high-count profile of a d-mask mod q."""
    ell = (d + 6) // 7
    common = gcd(d, q)
    quotient, remainder = divmod(ell, common)
    height = D // lcm(d, q)
    low = height * quotient
    increment = height if remainder else 0
    high_count = (q // common) * remainder
    require(0 <= high_count <= q, "high marginal count changed")
    return low, increment, high_count


def no_small_shape_count(D):
    total = 0
    for E in base.divisors_of(D):
        alphabet = sum(d > 4 for d in base.divisors_of(E))
        total += base.mobius(D // E) * comb(alphabet + 3, 4)
    require(total >= 0, ("negative no-small shape count", D))
    return total


@lru_cache(maxsize=None)
def literal_profile_distribution(D, q):
    """Independent divisor-symbol/current-lcm q-profile DP."""
    require(D % q == 0, "q does not divide D")
    states = {(0, 1, 0, 0, ()): 1}
    for divisor in (d for d in base.divisors_of(D) if d > 4):
        ambient = (D // divisor) * ((divisor + 6) // 7)
        low, increment, high_count = q_type(D, divisor, q)
        additions = Counter()
        for state, multiplicity in tuple(states.items()):
            used, current_lcm, old_ambient, baseline, high_profile = state
            for copies in range(1, 5 - used):
                new_profile = list(high_profile)
                if increment:
                    new_profile.extend([(increment, high_count)] * copies)
                    new_profile.sort()
                additions[
                    (
                        used + copies,
                        lcm(current_lcm, divisor),
                        old_ambient + copies * ambient,
                        baseline + copies * low,
                        tuple(new_profile),
                    )
                ] += multiplicity
        for state, multiplicity in additions.items():
            states[state] = states.get(state, 0) + multiplicity
    result = Counter()
    for state, multiplicity in states.items():
        used, current_lcm, ambient, baseline, high_profile = state
        if used == 4 and current_lcm == D:
            result[(ambient, baseline, high_profile)] += multiplicity
    require(sum(result.values()) == no_small_shape_count(D), "literal lcm total changed")
    return result


def target_thresholds(histogram):
    result = []
    running = 0
    for load, count in reversed(histogram):
        running += count
        if load:
            result.append((load, running))
    return tuple(result)


@lru_cache(maxsize=None)
def generic_profile_passes(q, baseline, high_profile, thresholds):
    increments = tuple(item[0] for item in high_profile)
    marginals = tuple(Q(item[1], q) for item in high_profile)
    for target_load, target_count in thresholds:
        need = target_load - baseline
        if need <= 0:
            continue
        maximum_count = q * base.maximum_upward_mass(
            increments, marginals, need
        )
        if Q(target_count) > maximum_count:
            return False
    return True


def exact_needle(D, d, phase, step):
    width = (d + 6) // 7
    classes = {(phase + index * step) % d for index in range(width)}
    return sum(1 << x for x in range(D) if x % d in classes)


def bitset_arcs(mask, D):
    arcs = []
    left = None
    for point in range(D + 1):
        occupied = point < D and bool((mask >> point) & 1)
        if occupied and left is None:
            left = point
        elif not occupied and left is not None:
            arcs.append((left, point))
            left = None
    return tuple(arcs)


def positive_profile_controls():
    """A literal four-needle union passes every generic q relaxation."""
    D = 84
    denominators = (6, 7, 12, 14)
    phases = (1, 2, 3, 4)
    steps = (1, 2, 5, 3)
    masks = tuple(
        exact_needle(D, d, phase, step)
        for d, phase, step in zip(denominators, phases, steps)
    )
    union = 0
    for mask in masks:
        union |= mask
    arcs = bitset_arcs(union, D)
    checks = 0
    for q in base.divisors_of(D):
        baseline = 0
        high_profile = []
        for divisor in denominators:
            low, increment, high_count = q_type(D, divisor, q)
            baseline += low
            if increment:
                high_profile.append((increment, high_count))
        high_profile.sort()
        require(
            generic_profile_passes(
                q,
                baseline,
                tuple(high_profile),
                target_thresholds(combined.residue_load_histogram(arcs, q)),
            ),
            ("positive q-profile control rejected", q),
        )
        checks += 1
    return checks, union.bit_count()


@lru_cache(maxsize=None)
def mobius_profile_distribution(D, q):
    """Exact-lcm profile distribution by grouped Mobius inversion."""
    require(D % q == 0, ("q does not divide D", D, q))
    result = Counter()
    for E in base.divisors_of(D):
        sign = base.mobius(D // E)
        if not sign:
            continue
        groups = Counter()
        for divisor in (d for d in base.divisors_of(E) if d > 4):
            ambient = (D // divisor) * ((divisor + 6) // 7)
            low, increment, high_count = q_type(D, divisor, q)
            groups[(ambient, low, increment, high_count)] += 1

        # State: used, ambient, baseline, sorted nonzero (increment,R) list.
        states = {(0, 0, 0, ()): 1}
        for feature, alphabet_size in groups.items():
            ambient, low, increment, high_count = feature
            additions = Counter()
            for state, multiplicity in tuple(states.items()):
                used, old_ambient, baseline, high_profile = state
                for copies in range(1, 5 - used):
                    new_profile = list(high_profile)
                    if increment:
                        new_profile.extend(
                            [(increment, high_count)] * copies
                        )
                        new_profile.sort()
                    additions[
                        (
                            used + copies,
                            old_ambient + copies * ambient,
                            baseline + copies * low,
                            tuple(new_profile),
                        )
                    ] += multiplicity * comb(
                        alphabet_size + copies - 1, copies
                    )
            for state, multiplicity in additions.items():
                states[state] = states.get(state, 0) + multiplicity

        for state, multiplicity in states.items():
            used, ambient, baseline, high_profile = state
            if used == 4:
                result[(ambient, baseline, high_profile)] += (
                    sign * multiplicity
                )

    require(
        all(multiplicity >= 0 for multiplicity in result.values()),
        ("negative exact-lcm q-profile coefficient", D, q),
    )
    result += Counter()
    require(
        sum(result.values()) == no_small_shape_count(D),
        ("q-profile shape coefficient changed", D, q),
    )
    return result


@lru_cache(maxsize=None)
def canonical_q7_distribution(D):
    """Exact-lcm distribution of the sufficient q=D/7 feature (b,R)."""
    require(D % 7 == 0, ("canonical quotient is not integral", D))
    q = D // 7
    result = Counter()
    for E in base.divisors_of(D):
        sign = base.mobius(D // E)
        if not sign:
            continue
        groups = Counter()
        for divisor in (d for d in base.divisors_of(E) if d > 4):
            ell = (divisor + 6) // 7
            if q % divisor == 0:
                feature = (0, (q // divisor) * ell)
                expected = (0, 7, feature[1])
            else:
                feature = (1, 0)
                expected = (1, 0, 0)
            require(
                q_type(D, divisor, q) == expected,
                ("uniform-one/seven-spike identity failed", D, divisor),
            )
            ambient = (D // divisor) * ell
            require(
                ambient == feature[0] * q + 7 * feature[1],
                ("canonical ambient identity failed", D, divisor),
            )
            groups[feature] += 1

        # State: used, uniform baseline b, total seven-spike cardinality R.
        states = {(0, 0, 0): 1}
        for feature, alphabet_size in groups.items():
            unit_b, unit_R = feature
            additions = Counter()
            for state, multiplicity in tuple(states.items()):
                used, baseline, spike_count = state
                for copies in range(1, 5 - used):
                    additions[
                        (
                            used + copies,
                            baseline + copies * unit_b,
                            spike_count + copies * unit_R,
                        )
                    ] += multiplicity * comb(
                        alphabet_size + copies - 1, copies
                    )
            for state, multiplicity in additions.items():
                states[state] = states.get(state, 0) + multiplicity
        for (used, baseline, spike_count), multiplicity in states.items():
            if used == 4:
                result[(baseline, spike_count)] += sign * multiplicity

    require(
        all(multiplicity >= 0 for multiplicity in result.values()),
        ("negative exact-lcm canonical q7 coefficient", D),
    )
    result += Counter()
    require(
        sum(result.values()) == no_small_shape_count(D),
        ("canonical q7 shape coefficient changed", D),
    )
    return result


def canonical_q7_passes(q, baseline, spike_count, histogram):
    """Closed form of the full q-profile relaxation at q=D/7."""
    require(0 <= baseline <= 4, "uniform baseline out of range")
    require(spike_count >= 0, "negative spike cardinality")
    require(
        all(0 <= load <= 7 for load, _count in histogram),
        "a q=D/7 fibre does not have size seven",
    )
    heavy_count = sum(
        count for load, count in histogram if load > baseline
    )
    return heavy_count <= min(q, spike_count)


def build_rows():
    """Reconstruct the complete k=3 projected-support row universe."""
    by_divisor = defaultdict(list)
    body_count = 0
    body_divisor_rows = 0
    for body in combinations(range(1, 15), 6):
        body_count += 1
        L, ranges = support.safe_cell_ranges(body)
        for D in base.divisors_of(L):
            body_divisor_rows += 1
            support_count = support.support_size_bitset(D, ranges)
            if Q(support_count, D) > base.SUPPORT_CUTOFF:
                continue
            arcs = combined.projected_support_arcs(D, ranges)
            require(
                sum(right - left for left, right in arcs) == support_count,
                ("support projection changed", body, D),
            )
            by_divisor[D].append((support_count, body, L, arcs))

    require(body_count == 3003, "body universe changed")
    require(body_divisor_rows == 251536, "body/divisor universe changed")
    require(sum(map(len, by_divisor.values())) == 26970, "row universe changed")
    require(len(by_divisor) == 217, "resolving-D universe changed")
    require(all(D % 7 == 0 for D in by_divisor), "nonseptimal D appeared")
    for records in by_divisor.values():
        records.sort(key=lambda row: (row[0], row[1], row[2]))
    return by_divisor, body_count, body_divisor_rows


def stage_digest_update(digest, D, record, count):
    support_count, body, L, _arcs = record
    digest.update(
        f"{body}|{L}|{D}|{support_count}|{count}\n".encode()
    )


def audit_direct_controls():
    """Compare grouped Mobius with the independent literal lcm-state DP."""
    cases = ((56, 8), (84, 12), (27720, 3960))
    for D, q in cases:
        literal = literal_profile_distribution(D, q)
        require(
            mobius_profile_distribution(D, q)
            == literal,
            ("Mobius/lcm-state q-profile disagreement", D, q),
        )
        collapsed = Counter()
        for (ambient, baseline, high_profile), multiplicity in literal.items():
            require(
                all(increment == 7 for increment, _R in high_profile),
                ("canonical increment is not seven", D),
            )
            spike_count = sum(R for _increment, R in high_profile)
            require(
                ambient == baseline * q + 7 * spike_count,
                ("canonical profile ambient identity failed", D),
            )
            collapsed[(baseline, spike_count)] += multiplicity
        require(
            canonical_q7_distribution(D) == collapsed,
            ("closed q7 profile/lcm-state disagreement", D),
        )
    return len(cases)


def main(local_moduli):
    require(7 in local_moduli, "canonical local modulus seven is mandatory")
    base.activity_controls()
    positive_checks, positive_union_size = positive_profile_controls()
    by_divisor, body_count, body_divisor_rows = build_rows()

    stages = ("raw_no_small", "scalar_no_small")
    occurrences = Counter()
    shapes = Counter()
    rows = {stage: set() for stage in stages}
    divisors = {stage: set() for stage in stages}
    semantic = {stage: sha256() for stage in stages}
    q_occurrences = Counter()
    q_shapes = Counter()
    q_rows = {modulus: set() for modulus in local_moduli}
    q_divisors = {modulus: set() for modulus in local_moduli}
    q_semantic = {modulus: sha256() for modulus in local_moduli}
    q_features = Counter()
    per_D = []
    feature_semantic = sha256()

    for D in sorted(by_divisor):
        records = by_divisor[D]
        applicable = tuple(
            modulus for modulus in local_moduli if D % modulus == 0
        )
        distributions = {
            modulus: (
                canonical_q7_distribution(D)
                if modulus == 7
                else mobius_profile_distribution(D, D // modulus)
            )
            for modulus in applicable
        }
        canonical = distributions[7]
        raw_shape_count = sum(canonical.values())
        raw_occurrence_count = raw_shape_count * len(records)
        occurrences["raw_no_small"] += raw_occurrence_count
        shapes["raw_no_small"] += raw_shape_count
        if raw_shape_count:
            divisors["raw_no_small"].add(D)
        for record in records:
            if raw_shape_count:
                rows["raw_no_small"].add((record[1], D))
                stage_digest_update(
                    semantic["raw_no_small"], D, record, raw_shape_count
                )

        scalar_occurrence_count = 0
        scalar_surviving_features = set()
        for record in records:
            support_count = record[0]
            row_count = sum(
                multiplicity
                for (baseline, spike_count), multiplicity
                in canonical.items()
                if baseline * (D // 7) + 7 * spike_count >= support_count
            )
            if not row_count:
                continue
            scalar_occurrence_count += row_count
            rows["scalar_no_small"].add((record[1], D))
            divisors["scalar_no_small"].add(D)
            stage_digest_update(
                semantic["scalar_no_small"], D, record, row_count
            )
            for feature, multiplicity in canonical.items():
                baseline, spike_count = feature
                if baseline * (D // 7) + 7 * spike_count >= support_count:
                    scalar_surviving_features.add(feature)

        scalar_shape_count = sum(
            canonical[feature] for feature in scalar_surviving_features
        )
        occurrences["scalar_no_small"] += scalar_occurrence_count
        shapes["scalar_no_small"] += scalar_shape_count

        q_D_summaries = []
        for modulus in applicable:
            q = D // modulus
            distribution = distributions[modulus]
            q_features[modulus] += len(distribution)
            passing_features = set()
            occurrence_count = 0
            passing_row_count = 0
            for record in records:
                support_count, body, _L, arcs = record
                histogram = combined.residue_load_histogram(arcs, q)
                thresholds = target_thresholds(histogram)
                row_count = 0
                for feature, multiplicity in distribution.items():
                    if modulus == 7:
                        baseline, spike_count = feature
                        ambient = baseline * q + 7 * spike_count
                        passes = canonical_q7_passes(
                            q, baseline, spike_count, histogram
                        )
                    else:
                        ambient, baseline, high_profile = feature
                        passes = generic_profile_passes(
                            q, baseline, high_profile, thresholds
                        )
                    if modulus == 7 and passes:
                        require(
                            ambient >= support_count,
                            (
                                "canonical q7 condition failed to imply "
                                "scalar capacity",
                                D,
                                body,
                                feature,
                            ),
                        )
                    # Keep the intersection semantics explicit: this is the
                    # scalar necessary condition AND the chosen q screen.
                    if ambient < support_count:
                        continue
                    if not passes:
                        continue
                    row_count += multiplicity
                    passing_features.add(feature)
                if not row_count:
                    continue
                occurrence_count += row_count
                passing_row_count += 1
                q_rows[modulus].add((body, D))
                q_divisors[modulus].add(D)
                stage_digest_update(
                    q_semantic[modulus], D, record, row_count
                )
            shape_count = sum(
                distribution[feature] for feature in passing_features
            )
            q_occurrences[modulus] += occurrence_count
            q_shapes[modulus] += shape_count
            q_D_summaries.append(
                (
                    modulus,
                    len(distribution),
                    occurrence_count,
                    shape_count,
                    passing_row_count,
                )
            )
            for feature, multiplicity in sorted(distribution.items()):
                feature_semantic.update(
                    f"{D}|{modulus}|{feature}|{multiplicity}\n".encode()
                )

        per_D.append(
            (
                D,
                len(records),
                raw_shape_count,
                raw_occurrence_count,
                scalar_shape_count,
                scalar_occurrence_count,
                tuple(q_D_summaries),
            )
        )

    direct_control_count = audit_direct_controls()
    require(
        q_occurrences[7] <= occurrences["scalar_no_small"],
        "q=7 occurrence screen enlarged scalar stage",
    )
    require(
        q_shapes[7] <= shapes["scalar_no_small"],
        "q=7 shape screen enlarged scalar stage",
    )
    # Previously audited bounded terminal cases.
    D_summary = {row[0]: row for row in per_D}
    require(D_summary[56][5] == 110, "D=56 scalar control changed")
    require(D_summary[56][6][0][2] == 0, "D=56 q=D/7 control changed")
    require(D_summary[84][5] == 20, "D=84 scalar control changed")
    require(D_summary[84][6][0][2] == 0, "D=84 q=D/7 control changed")
    require(
        D_summary[27720][5] == 108929658,
        "D=27720 scalar control changed",
    )
    require(
        D_summary[27720][6][0][2] == 11523257,
        "D=27720 q=D/7 control changed",
    )

    print("LRC14 k=3/four-drift all-D exact q-profile GF scout")
    print(f"body_count={body_count}")
    print(f"body_divisor_rows={body_divisor_rows}")
    print(f"support_rows={sum(map(len, by_divisor.values()))}")
    print(f"support_divisors={len(by_divisor)}")
    print(f"local_moduli={local_moduli}")
    print(f"positive_profile_checks={positive_checks}")
    print(f"positive_union_size={positive_union_size}")
    print(f"direct_lcm_state_control_count={direct_control_count}")
    print(f"feature_semantic_sha256={feature_semantic.hexdigest()}")
    for stage in stages:
        print(
            f"stage={stage},occurrences={occurrences[stage]},"
            f"shapes={shapes[stage]},rows={len(rows[stage])},"
            f"divisors={len(divisors[stage])},"
            f"semantic_sha256={semantic[stage].hexdigest()}"
        )
    for modulus in local_moduli:
        print(
            f"stage=q_local_{modulus}_AND_scalar,"
            f"occurrences={q_occurrences[modulus]},"
            f"shapes={q_shapes[modulus]},rows={len(q_rows[modulus])},"
            f"divisors={len(q_divisors[modulus])},"
            f"aggregated_features={q_features[modulus]},"
            f"semantic_sha256={q_semantic[modulus].hexdigest()}"
        )
    print(f"per_D={per_D}")
    print("screen_intersection_warning=separate q screens are not intersected")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--local-moduli", default="7")
    arguments = parser.parse_args()
    main(
        tuple(
            sorted(
                {
                    int(value)
                    for value in arguments.local_moduli.split(",")
                    if value
                }
            )
        )
    )

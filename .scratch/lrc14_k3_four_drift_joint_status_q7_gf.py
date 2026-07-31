#!/usr/bin/env python3
"""Exact joint common-u activity x q=D/7 screen for k3/four drifts.

All 217 projected-support resolving denominators are divisible by seven.
Put q=D/7.  A large denominator d>4 has the exact q-fibre type proved in
``lrc14_k3_four_drift_q7_all_D_gf.py``:

* d not dividing q contributes uniform baseline one;
* d dividing q contributes a 0/7 spike on
  R=(q/d)ceil(d/7) fibres.

The small denominators d=2,3,4 all divide q.  Conditional on common aligned
phase u, such a mask is either empty or a 0/7 spike on q/d fibres, and its
activity bit has exact marginal d/7.

Let c be the number of uniform large masks, R the sum of the large-spike
cardinalities, and

    N_c(S_D) = #{b mod q: lambda_q(S_D,b) > c}.

For an active small-label set E, the arbitrary-location q relaxation can
cover the target only if

    N_c(S_D) <= R + sum_{i in E} q/d_i.                 (*)

This is an upward threshold event in the small activity bits.  The exact
fractional-cover theorem bounds its possible common-u mass from the
one-marginals d_i/7.  A counterexample would force this event to contain the
compact aligned-safe carrier of mass at least 55/91.  The event is open, so
the equality endpoint is also impossible.  Therefore a feature survives
only when the exact maximum upward-event mass is strictly greater than
55/91.

The exact-lcm multiset GF is compressed to

    (m2,m3,m4,c,R).

Before applying (*), the program intersects with the already-audited
row-specific common-u support-status screen.  Thus the reported joint stage
dominates that prior screen and specializes exactly to the no-small q7
screen when (m2,m3,m4)=(0,0,0).  No denominator tuple is materialized.
"""

from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
Q7_PATH = Path(__file__).with_name(
    "lrc14_k3_four_drift_q7_all_D_gf.py"
)
spec = spec_from_file_location("lrc14_k3_q7", Q7_PATH)
q7 = module_from_spec(spec)
spec.loader.exec_module(q7)
base = q7.base
support = q7.support
combined = q7.combined


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def joint_distribution(D):
    """Exact-lcm distribution on (m2,m3,m4,c,R) by Mobius inversion."""
    require(D % 7 == 0, ("nonseptimal resolving denominator", D))
    q = D // 7
    result = Counter()
    for E in base.divisors_of(D):
        sign = base.mobius(D // E)
        if not sign:
            continue
        groups = Counter()
        for divisor in (d for d in base.divisors_of(E) if d > 1):
            if divisor in (2, 3, 4):
                feature = (
                    int(divisor == 2),
                    int(divisor == 3),
                    int(divisor == 4),
                    0,
                    0,
                )
                require(q % divisor == 0, ("small d does not divide q", D, divisor))
                require(
                    q7.q_type(D, divisor, q)
                    == (0, 7, q // divisor),
                    ("small seven-spike identity failed", D, divisor),
                )
            else:
                ell = (divisor + 6) // 7
                if q % divisor == 0:
                    feature = (0, 0, 0, 0, (q // divisor) * ell)
                else:
                    feature = (0, 0, 0, 1, 0)
            groups[feature] += 1

        # State: used,m2,m3,m4,c,R.
        states = {(0, 0, 0, 0, 0, 0): 1}
        for feature, alphabet_size in groups.items():
            unit2, unit3, unit4, unit_c, unit_R = feature
            additions = Counter()
            for state, multiplicity in tuple(states.items()):
                used, m2, m3, m4, c, R = state
                for copies in range(1, 5 - used):
                    additions[
                        (
                            used + copies,
                            m2 + copies * unit2,
                            m3 + copies * unit3,
                            m4 + copies * unit4,
                            c + copies * unit_c,
                            R + copies * unit_R,
                        )
                    ] += multiplicity * comb(
                        alphabet_size + copies - 1, copies
                    )
            for state, multiplicity in additions.items():
                states[state] = states.get(state, 0) + multiplicity
        for state, multiplicity in states.items():
            used, m2, m3, m4, c, R = state
            if used == 4:
                result[((m2, m3, m4), c, R)] += sign * multiplicity

    require(
        all(multiplicity >= 0 for multiplicity in result.values()),
        ("negative joint exact-lcm coefficient", D),
    )
    result += Counter()

    # All-D weighted-feature control: forgetting (c,R) recovers the
    # independently audited primary distribution coefficientwise.
    collapsed = Counter()
    for (pattern, c, R), multiplicity in result.items():
        collapsed[(pattern, c * q + 7 * R)] += multiplicity
    require(
        collapsed == base.shape_distribution(D),
        ("joint GF does not collapse to weighted primary GF", D),
    )
    return result


@lru_cache(maxsize=None)
def joint_activity_limit(pattern, D):
    """Largest q-spike demand allowed above R by common-u transport."""
    q = D // 7
    weights = []
    marginals = []
    for divisor, copies in zip((2, 3, 4), pattern):
        if not copies:
            continue
        require(q % divisor == 0, ("small d does not divide q", D, divisor))
        weights.extend([q // divisor] * copies)
        marginals.extend([Q(divisor, 7)] * copies)
    return base.status_need_limit(tuple(weights), tuple(marginals))


def main():
    base.activity_controls()
    by_divisor, body_count, body_divisor_rows = q7.build_rows()
    stages = ("raw", "support_status", "joint_status_q7")
    occurrences = Counter()
    shapes = Counter()
    rows = {stage: set() for stage in stages}
    bodies = {stage: set() for stage in stages}
    divisors = {stage: set() for stage in stages}
    pattern_occurrences = {stage: Counter() for stage in stages}
    pattern_shapes = {stage: Counter() for stage in stages}
    semantic = {stage: sha256() for stage in stages}
    feature_semantic = sha256()
    feature_count = 0

    for D in sorted(by_divisor):
        q = D // 7
        distribution = joint_distribution(D)
        feature_count += len(distribution)
        for feature, multiplicity in sorted(distribution.items()):
            feature_semantic.update(f"{D}|{feature}|{multiplicity}\n".encode())

        by_pattern_c = defaultdict(Counter)
        for (pattern, c, R), multiplicity in distribution.items():
            by_pattern_c[(pattern, c)][R] += multiplicity
        suffix_data = {
            key: base.suffix_counter(counter)
            for key, counter in by_pattern_c.items()
        }
        raw_shape_count = sum(distribution.values())
        shapes["raw"] += raw_shape_count
        if raw_shape_count:
            divisors["raw"].add(D)
        raw_pattern_shapes = Counter()
        for (pattern, _c), counter in by_pattern_c.items():
            raw_pattern_shapes[pattern] += sum(counter.values())
        pattern_shapes["raw"].update(raw_pattern_shapes)
        pattern_occurrences["raw"].update(
            {
                pattern: count * len(by_divisor[D])
                for pattern, count in raw_pattern_shapes.items()
            }
        )

        minimum_threshold = {
            "support_status": {},
            "joint_status_q7": {},
        }
        for record in by_divisor[D]:
            support_count, body, L, arcs = record
            small_loads = {}
            for divisor in (2, 3, 4):
                if D % divisor == 0:
                    histogram_d = combined.residue_load_histogram(arcs, divisor)
                    small_loads[divisor] = combined.top_class_load(
                        histogram_d, 1
                    )
            histogram_q = combined.residue_load_histogram(arcs, q)
            heavy_by_c = {
                c: sum(count for load, count in histogram_q if load > c)
                for c in range(5)
            }
            raw_row_count = raw_shape_count
            support_row_count = 0
            joint_row_count = 0
            support_pattern_count = Counter()
            joint_pattern_count = Counter()
            support_limits = {}
            activity_limits = {}
            for pattern, _c in by_pattern_c:
                if pattern not in support_limits:
                    local_weights, marginals = base.small_vectors(
                        pattern, D, small_loads
                    )
                    support_limits[pattern] = base.status_need_limit(
                        local_weights, marginals
                    )
                    activity_limits[pattern] = joint_activity_limit(
                        pattern, D
                    )

            for key, (keys, suffix) in suffix_data.items():
                pattern, c = key
                support_limit = support_limits[pattern]
                # c*q+7R+support_limit >= support_count.
                support_R = (
                    support_count - support_limit - c * q + 6
                ) // 7
                support_count_here = base.count_at_least(
                    keys, suffix, support_R
                )
                if not support_count_here:
                    continue
                support_row_count += support_count_here
                support_pattern_count[pattern] += support_count_here
                old = minimum_threshold["support_status"].get(key)
                if old is None or support_R < old:
                    minimum_threshold["support_status"][key] = support_R

                # N_c <= R + active-small spike capacity.
                joint_R = max(
                    support_R,
                    heavy_by_c[c] - activity_limits[pattern],
                )
                joint_count_here = base.count_at_least(keys, suffix, joint_R)
                if not joint_count_here:
                    continue
                joint_row_count += joint_count_here
                joint_pattern_count[pattern] += joint_count_here
                old = minimum_threshold["joint_status_q7"].get(key)
                if old is None or joint_R < old:
                    minimum_threshold["joint_status_q7"][key] = joint_R

            row_key = (body, D)
            for stage, count in (
                ("raw", raw_row_count),
                ("support_status", support_row_count),
                ("joint_status_q7", joint_row_count),
            ):
                if not count:
                    continue
                occurrences[stage] += count
                rows[stage].add(row_key)
                bodies[stage].add(body)
                divisors[stage].add(D)
                semantic[stage].update(
                    f"{body}|{L}|{D}|{support_count}|{count}\n".encode()
                )
            pattern_occurrences["support_status"].update(support_pattern_count)
            pattern_occurrences["joint_status_q7"].update(joint_pattern_count)

        for stage in ("support_status", "joint_status_q7"):
            local = Counter()
            for key, threshold in minimum_threshold[stage].items():
                pattern, _c = key
                keys, suffix = suffix_data[key]
                count = base.count_at_least(keys, suffix, threshold)
                local[pattern] += count
                shapes[stage] += count
            pattern_shapes[stage].update(local)

    require(occurrences["raw"] == 21357714101, "raw occurrence total changed")
    require(shapes["raw"] == 694921995, "raw shape total changed")
    require(
        occurrences["support_status"] == 13280722299,
        "support-status occurrence total changed",
    )
    require(
        shapes["support_status"] == 694254050,
        "support-status shape total changed",
    )
    require(
        len(rows["support_status"]) == 18599,
        "support-status row total changed",
    )
    require(
        len(divisors["support_status"]) == 186,
        "support-status D total changed",
    )
    require(
        occurrences["joint_status_q7"]
        <= occurrences["support_status"],
        "joint stage enlarged support-status occurrences",
    )
    require(
        shapes["joint_status_q7"] <= shapes["support_status"],
        "joint stage enlarged support-status shapes",
    )
    # Exact specialization control: with no small labels the joint screen
    # must reproduce the frozen all-D no-small q7 census.
    no_small_pattern = (0, 0, 0)
    require(
        pattern_occurrences["support_status"][no_small_pattern]
        == 12852450428,
        "no-small support-status control changed",
    )
    require(
        pattern_occurrences["joint_status_q7"][no_small_pattern]
        == 1437089787,
        "no-small joint q7 occurrence control changed",
    )
    require(
        pattern_shapes["joint_status_q7"][no_small_pattern]
        == 171159117,
        "no-small joint q7 shape control changed",
    )

    print("LRC14 k=3/four-drift joint common-u x q=D/7 GF scout")
    print(f"body_count={body_count}")
    print(f"body_divisor_rows={body_divisor_rows}")
    print(f"support_rows={sum(map(len, by_divisor.values()))}")
    print(f"support_divisors={len(by_divisor)}")
    print(f"aggregated_features={feature_count}")
    print(f"feature_semantic_sha256={feature_semantic.hexdigest()}")
    for stage in stages:
        print(
            f"stage={stage},occurrences={occurrences[stage]},"
            f"shapes={shapes[stage]},rows={len(rows[stage])},"
            f"bodies={len(bodies[stage])},divisors={len(divisors[stage])},"
            f"semantic_sha256={semantic[stage].hexdigest()}"
        )
        print(
            f"stage={stage},pattern_occurrences="
            f"{Counter(dict(sorted(pattern_occurrences[stage].items())))}"
        )
        print(
            f"stage={stage},pattern_shapes="
            f"{Counter(dict(sorted(pattern_shapes[stage].items())))}"
        )
    print(f"joint_activity_limit_cache={joint_activity_limit.cache_info()}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()

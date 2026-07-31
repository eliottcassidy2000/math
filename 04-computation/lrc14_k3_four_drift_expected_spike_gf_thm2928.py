#!/usr/bin/env python3
"""All-D expected-spike and exact one-spike screens for k=3/four drifts.

Every projected-support resolving denominator ``D`` is divisible by seven.
Put ``q=D/7``.  For a denominator ``d|D`` one has

    d/gcd(d,q) in {1,7}.

We grant a mask with ``d not| q`` the uniform baseline one on every
q-fibre.  This is an upper relaxation of its literal strict address mask.
If ``c`` of the four masks have this type, the remaining ``m=4-c`` masks
have ``d|q`` and hit whole q-fibres with load seven.  Fubini gives the
exact expected number ``q/7`` of filled fibres for every such mask.

For a projected-safe target let

    N_c = #{b mod q : lambda_q(S_D,b) > c}.

A cover forces the total spike support ``T`` to satisfy ``T>=N_c`` on the
compact aligned-safe carrier, whose measure is at least 55/91.  The event
``T>=N_c`` is open.  If it is the whole circle its mass is 1>55/91;
otherwise compact/open separation makes its mass strictly larger than the
safe mass.  Markov therefore gives the necessary condition

    N_c=0  or  55*N_c < 13*(4-c)*q.                 (mean screen)

The equality endpoint is killed.  This program extracts the exact-lcm
arity-four GF on ``(m2,m3,m4,c,large_capacity)`` and intersects that mean
screen with the independently established row-specific support-status
screen.  No denominator four-tuple is materialized.

The next stage sharpens exactly the ``c=3`` slice.  There is then one spike
denominator ``d|q``.  Write ``d=7*a+r`` and ``w=q/d``.  On open phase
cells, and hence almost everywhere, its number of filled q-fibres is

    a*w + w*X_d(u),       Haar(X_d=1)=r/7.

At a strict-mask equality phase the actual count can be smaller.  The
displayed two-level law is therefore also a valid pointwise upper envelope,
which is the direction used by every necessary screen below.

Thus a threshold event can have mass strictly above 55/91 exactly through

    N_3 <= a*w + 1_{r in {5,6}}*w.                 (one-spike screen)

The c=3 exact-lcm distribution is derived independently by fixing the one
spike denominator and choosing a three-multiset from the uniform alphabet,
then applying divisor-lattice Mobius inversion.  It is checked against the
c=3 slice of the all-c GF coefficientwise.

Both new screens are necessary upper relaxations.  They are intersected
with support status only at the level of separate necessary conditions;
the shared common-phase bits are not optimized jointly.
"""

from bisect import bisect_left
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations_with_replacement
from math import comb, gcd, lcm
from pathlib import Path


Q7_PATH = Path(__file__).with_name(
    "lrc14_k3_four_drift_q7_all_D_gf_thm2928.py"
)
EXPECTED_Q7_SHA256 = (
    "3cc07195d580c5c5c01457ea95b58837a25c2176d326d12feaccc8e0bfa28dcc"
)
if sha256(Q7_PATH.read_bytes()).hexdigest() != EXPECTED_Q7_SHA256:
    raise RuntimeError("frozen q7 four-drift dependency changed")
spec = spec_from_file_location("lrc14_k3_q7", Q7_PATH)
q7 = module_from_spec(spec)
spec.loader.exec_module(q7)
base = q7.base
combined = q7.combined


STAGES = (
    "raw",
    "support_status",
    "expected_spike_AND_support",
    "c3_exact_one_spike_AND_expected_AND_support",
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def denominator_feature(D, divisor):
    """Return (small pattern,c,large ceiling capacity) for one symbol."""
    q = D // 7
    pattern = (
        int(divisor == 2),
        int(divisor == 3),
        int(divisor == 4),
    )
    if q % divisor:
        require(divisor % 7 == 0, ("bad uniform denominator", D, divisor))
        require(not any(pattern), ("small denominator became uniform", D, divisor))
        return pattern, 1, q
    large_capacity = 0
    if not any(pattern):
        large_capacity = (D // divisor) * ((divisor + 6) // 7)
    return pattern, 0, large_capacity


@lru_cache(maxsize=None)
def expected_distribution(D):
    """Exact-lcm GF on (small pattern,c,large ceiling capacity)."""
    require(D % 7 == 0, ("nonseptimal D", D))
    result = Counter()
    for E in base.divisors_of(D):
        sign = base.mobius(D // E)
        if not sign:
            continue
        groups = Counter(
            denominator_feature(D, divisor)
            for divisor in base.divisors_of(E)
            if divisor > 1
        )
        # State: used,m2,m3,m4,c,large_capacity.
        states = {(0, 0, 0, 0, 0, 0): 1}
        for feature, alphabet_size in groups.items():
            pattern, unit_c, unit_capacity = feature
            additions = Counter()
            for state, multiplicity in tuple(states.items()):
                used, m2, m3, m4, c, capacity = state
                for copies in range(1, 5 - used):
                    additions[
                        (
                            used + copies,
                            m2 + copies * pattern[0],
                            m3 + copies * pattern[1],
                            m4 + copies * pattern[2],
                            c + copies * unit_c,
                            capacity + copies * unit_capacity,
                        )
                    ] += multiplicity * comb(alphabet_size + copies - 1, copies)
            for state, multiplicity in additions.items():
                states[state] = states.get(state, 0) + multiplicity
        for state, multiplicity in states.items():
            used, m2, m3, m4, c, capacity = state
            if used == 4:
                result[((m2, m3, m4), c, capacity)] += sign * multiplicity

    require(all(value >= 0 for value in result.values()), ("negative GF", D))
    result += Counter()
    collapsed = Counter()
    for (pattern, _c, capacity), multiplicity in result.items():
        collapsed[(pattern, capacity)] += multiplicity
    require(
        collapsed == base.shape_distribution(D),
        ("expected GF does not collapse to primary", D),
    )
    return result


def one_spike_allowance(q, divisor):
    """Largest positive-mass-safe fibre demand for one d=7a+r spike."""
    require(q % divisor == 0, ("one-spike denominator does not divide q", q, divisor))
    a, r = divmod(divisor, 7)
    width = q // divisor
    return a * width + (width if r >= 5 else 0)


@lru_cache(maxsize=None)
def c3_denominator_distribution(D):
    """Independent exact-lcm c=3 GF retaining the unique spike d."""
    q = D // 7
    result = Counter()
    for E in base.divisors_of(D):
        sign = base.mobius(D // E)
        if not sign:
            continue
        alphabet = tuple(d for d in base.divisors_of(E) if d > 1)
        uniform_count = sum(q % d != 0 for d in alphabet)
        uniform_multisets = comb(uniform_count + 2, 3)
        if not uniform_multisets:
            continue
        for divisor in alphabet:
            if q % divisor:
                continue
            pattern, unit_c, spike_capacity = denominator_feature(D, divisor)
            require(unit_c == 0, ("spike marked uniform", D, divisor))
            large_capacity = 3 * q + spike_capacity
            allowance = one_spike_allowance(q, divisor)
            result[(pattern, divisor, large_capacity, allowance)] += (
                sign * uniform_multisets
            )

    require(all(value >= 0 for value in result.values()), ("negative c3 GF", D))
    result += Counter()
    collapsed = Counter()
    for (pattern, _d, capacity, _allowance), multiplicity in result.items():
        collapsed[(pattern, capacity)] += multiplicity
    expected_c3 = Counter()
    for (pattern, c, capacity), multiplicity in expected_distribution(D).items():
        if c == 3:
            expected_c3[(pattern, capacity)] += multiplicity
    require(collapsed == expected_c3, ("c3 denominator GF mismatch", D))
    return result


def mean_passes(N, c, q):
    return N == 0 or 55 * N < 13 * (4 - c) * q


def suffix(counter):
    keys = sorted(counter)
    totals = [0] * (len(keys) + 1)
    for index in range(len(keys) - 1, -1, -1):
        totals[index] = totals[index + 1] + counter[keys[index]]
    return keys, totals


def at_least(data, threshold):
    keys, totals = data
    return totals[bisect_left(keys, threshold)]


def distribution_controls():
    """Literal multiset controls for small D, plus the residue endpoint."""
    cases = 0
    shapes = 0
    for D in (14, 28, 42, 56, 70, 84):
        q = D // 7
        alphabet = tuple(d for d in base.divisors_of(D) if d > 1)
        brute = Counter()
        brute_c3 = Counter()
        for values in combinations_with_replacement(alphabet, 4):
            if lcm(*values) != D:
                continue
            pattern = tuple(values.count(d) for d in (2, 3, 4))
            c = 0
            capacity = 0
            spike = []
            for divisor in values:
                _pattern, unit_c, unit_capacity = denominator_feature(D, divisor)
                c += unit_c
                capacity += unit_capacity
                if not unit_c:
                    spike.append(divisor)
            brute[(pattern, c, capacity)] += 1
            if c == 3:
                require(len(spike) == 1, ("c3 does not have one spike", D, values))
                divisor = spike[0]
                brute_c3[
                    (pattern, divisor, capacity, one_spike_allowance(q, divisor))
                ] += 1
            shapes += 1
        require(brute == expected_distribution(D), ("literal expected GF failed", D))
        require(brute_c3 == c3_denominator_distribution(D), ("literal c3 GF failed", D))
        cases += 1

    for r in range(7):
        limit = base.status_need_limit((1,), (Q(r, 7),))
        require(limit == int(r >= 5), ("one-spike residue endpoint failed", r, limit))
    return cases, shapes, 7


def update_minimum(old, candidate):
    return candidate if old is None or candidate < old else old


def main():
    base.activity_controls()
    literal_D_cases, literal_shapes, residue_controls = distribution_controls()
    by_divisor, body_count, body_divisor_rows = q7.build_rows()

    occurrences = Counter()
    shapes = Counter()
    rows = {stage: set() for stage in STAGES}
    bodies = {stage: set() for stage in STAGES}
    divisors = {stage: set() for stage in STAGES}
    pattern_occurrences = {stage: Counter() for stage in STAGES}
    pattern_shapes = {stage: Counter() for stage in STAGES}
    c_occurrences = {stage: Counter() for stage in STAGES}
    c_shapes = {stage: Counter() for stage in STAGES}
    semantic = {stage: sha256() for stage in STAGES}
    feature_semantic = sha256()
    c3_feature_semantic = sha256()
    feature_count = 0
    c3_feature_count = 0
    equality_pairs = 0
    equality_by_c = Counter()
    equality_support_occurrences = 0
    equality_support_by_c = Counter()
    minimum_equality = None
    c3_endpoint_occurrences = 0
    c3_endpoint_features = set()
    minimum_c3_endpoint = None
    final_c3_occurrences_by_residue = Counter()
    final_c3_shapes_by_residue = Counter()
    minimum_final = None
    minimum_final_c3 = None

    for D in sorted(by_divisor):
        q = D // 7
        distribution = expected_distribution(D)
        c3_distribution = c3_denominator_distribution(D)
        feature_count += len(distribution)
        c3_feature_count += len(c3_distribution)
        for feature, multiplicity in sorted(distribution.items()):
            feature_semantic.update(f"{D}|{feature}|{multiplicity}\n".encode())
        for feature, multiplicity in sorted(c3_distribution.items()):
            c3_feature_semantic.update(f"{D}|{feature}|{multiplicity}\n".encode())

        by_pattern_c = defaultdict(Counter)
        for (pattern, c, capacity), multiplicity in distribution.items():
            by_pattern_c[(pattern, c)][capacity] += multiplicity
        suffix_data = {key: suffix(counter) for key, counter in by_pattern_c.items()}
        raw_shapes = sum(distribution.values())
        raw_pattern = Counter()
        raw_c = Counter()
        for (pattern, c), counter in by_pattern_c.items():
            count = sum(counter.values())
            raw_pattern[pattern] += count
            raw_c[c] += count
        shapes["raw"] += raw_shapes
        pattern_shapes["raw"].update(raw_pattern)
        c_shapes["raw"].update(raw_c)
        if raw_shapes:
            divisors["raw"].add(D)

        survivor_features = {
            "support_status": set(),
            "expected_spike_AND_support": set(),
        }
        final_non_c3_features = set()
        final_c3_features = set()

        for support_count, body, L, arcs in by_divisor[D]:
            histogram_q = combined.residue_load_histogram(arcs, q)
            N_by_c = {
                c: sum(count for load, count in histogram_q if load > c)
                for c in range(5)
            }
            small_loads = {}
            for divisor in (2, 3, 4):
                if D % divisor == 0:
                    histogram_d = combined.residue_load_histogram(arcs, divisor)
                    small_loads[divisor] = combined.top_class_load(histogram_d, 1)

            support_limits = {}
            for pattern, _c in by_pattern_c:
                if pattern in support_limits:
                    continue
                weights, marginals = base.small_vectors(pattern, D, small_loads)
                support_limits[pattern] = base.status_need_limit(weights, marginals)

            row_stage_count = Counter(raw=raw_shapes)
            row_stage_pattern = {stage: Counter() for stage in STAGES}
            row_stage_c = {stage: Counter() for stage in STAGES}
            row_stage_pattern["raw"] = raw_pattern.copy()
            row_stage_c["raw"] = raw_c.copy()
            support_c_here = Counter()

            for (pattern, c), data in suffix_data.items():
                threshold = support_count - support_limits[pattern]
                count = at_least(data, threshold)
                if not count:
                    continue
                support_c_here[c] += count
                row_stage_count["support_status"] += count
                row_stage_pattern["support_status"][pattern] += count
                row_stage_c["support_status"][c] += count
                for capacity in data[0][bisect_left(data[0], threshold):]:
                    survivor_features["support_status"].add((pattern, c, capacity))

                if not mean_passes(N_by_c[c], c, q):
                    continue
                row_stage_count["expected_spike_AND_support"] += count
                row_stage_pattern["expected_spike_AND_support"][pattern] += count
                row_stage_c["expected_spike_AND_support"][c] += count
                for capacity in data[0][bisect_left(data[0], threshold):]:
                    survivor_features["expected_spike_AND_support"].add(
                        (pattern, c, capacity)
                    )
                if c != 3:
                    row_stage_count[
                        "c3_exact_one_spike_AND_expected_AND_support"
                    ] += count
                    row_stage_pattern[
                        "c3_exact_one_spike_AND_expected_AND_support"
                    ][pattern] += count
                    row_stage_c[
                        "c3_exact_one_spike_AND_expected_AND_support"
                    ][c] += count
                    for capacity in data[0][bisect_left(data[0], threshold):]:
                        final_non_c3_features.add((pattern, c, capacity))
                        minimum_final = update_minimum(
                            minimum_final,
                            (D, body, pattern, c, capacity, 0),
                        )

            # Replace, rather than further aggregate, the c=3 slice by its
            # denominator-resolved exact distribution.
            if mean_passes(N_by_c[3], 3, q):
                for feature, multiplicity in c3_distribution.items():
                    pattern, spike_d, capacity, allowance = feature
                    if capacity + support_limits[pattern] < support_count:
                        continue
                    if allowance < N_by_c[3]:
                        continue
                    require(
                        mean_passes(N_by_c[3], 3, q),
                        ("exact c3 survivor failed mean screen", D, body, feature),
                    )
                    stage = "c3_exact_one_spike_AND_expected_AND_support"
                    row_stage_count[stage] += multiplicity
                    row_stage_pattern[stage][pattern] += multiplicity
                    row_stage_c[stage][3] += multiplicity
                    final_c3_features.add(feature)
                    final_c3_occurrences_by_residue[spike_d % 7] += multiplicity
                    candidate = (D, body, pattern, 3, capacity, spike_d)
                    minimum_final = update_minimum(minimum_final, candidate)
                    minimum_final_c3 = update_minimum(
                        minimum_final_c3,
                        (
                            D,
                            body,
                            L,
                            support_count,
                            pattern,
                            spike_d,
                            capacity,
                            allowance,
                            N_by_c[3],
                        ),
                    )
                    if N_by_c[3] > 0 and allowance == N_by_c[3]:
                        c3_endpoint_occurrences += multiplicity
                        c3_endpoint_features.add((D, feature))
                        minimum_c3_endpoint = update_minimum(
                            minimum_c3_endpoint,
                            (
                                D,
                                body,
                                L,
                                support_count,
                                pattern,
                                spike_d,
                                allowance,
                            ),
                        )

            for c in range(5):
                N = N_by_c[c]
                if N > 0 and 55 * N == 13 * (4 - c) * q:
                    equality_pairs += 1
                    equality_by_c[c] += 1
                    killed = support_c_here[c]
                    equality_support_occurrences += killed
                    equality_support_by_c[c] += killed
                    minimum_equality = update_minimum(
                        minimum_equality,
                        (D, body, L, support_count, c, N, q, killed),
                    )

            row_key = (body, D)
            for stage in STAGES:
                count = row_stage_count[stage]
                if not count:
                    continue
                occurrences[stage] += count
                rows[stage].add(row_key)
                bodies[stage].add(body)
                divisors[stage].add(D)
                pattern_occurrences[stage].update(row_stage_pattern[stage])
                c_occurrences[stage].update(row_stage_c[stage])
                semantic[stage].update(
                    f"{body}|{L}|{D}|{support_count}|{count}\n".encode()
                )

        for stage in ("support_status", "expected_spike_AND_support"):
            for feature in survivor_features[stage]:
                pattern, c, _capacity = feature
                multiplicity = distribution[feature]
                shapes[stage] += multiplicity
                pattern_shapes[stage][pattern] += multiplicity
                c_shapes[stage][c] += multiplicity
        final_stage = "c3_exact_one_spike_AND_expected_AND_support"
        for feature in final_non_c3_features:
            pattern, c, _capacity = feature
            multiplicity = distribution[feature]
            shapes[final_stage] += multiplicity
            pattern_shapes[final_stage][pattern] += multiplicity
            c_shapes[final_stage][c] += multiplicity
        for feature in final_c3_features:
            pattern, spike_d, _capacity, _allowance = feature
            multiplicity = c3_distribution[feature]
            shapes[final_stage] += multiplicity
            pattern_shapes[final_stage][pattern] += multiplicity
            c_shapes[final_stage][3] += multiplicity
            final_c3_shapes_by_residue[spike_d % 7] += multiplicity

    require(occurrences["raw"] == 21357714101, "raw occurrence changed")
    require(shapes["raw"] == 694921995, "raw shapes changed")
    require(
        occurrences["support_status"] == 13280722299,
        "support-status occurrence changed",
    )
    require(shapes["support_status"] == 694254050, "support-status shape changed")
    expected_stage_census = {
        "raw": (21357714101, 694921995, 26970, 3003, 217),
        "support_status": (13280722299, 694254050, 18599, 3003, 186),
        "expected_spike_AND_support": (2934202044, 400005870, 2120, 2037, 115),
        "c3_exact_one_spike_AND_expected_AND_support": (
            2548901482,
            398241574,
            1904,
            1823,
            107,
        ),
    }
    for stage, expected in expected_stage_census.items():
        actual = (
            occurrences[stage],
            shapes[stage],
            len(rows[stage]),
            len(bodies[stage]),
            len(divisors[stage]),
        )
        require(actual == expected, ("frozen stage census changed", stage, actual))
    require(
        feature_semantic.hexdigest()
        == "027510f7b7e0aada9e43860384f89ad0d745ed4aee9f60aae47c5d631a9ef71f",
        "all-c feature semantic changed",
    )
    require(
        c3_feature_semantic.hexdigest()
        == "dcc4288895b95600a923ce9ef2149257aace6f2b9ce4dc234e7831af4f0a8e8f",
        "c3 feature semantic changed",
    )
    expected_stage_semantic = {
        "raw": "964cebef8f435f0964807831f5a7cfdb262e11c886e06d87713573c188fd6fc7",
        "support_status": "7576fc60722ef889ef541a3c0ede7011593fda0c71899960fb67a8a4156a24ac",
        "expected_spike_AND_support": "2dde966a4dfe379c2192689160a7ad117d93f130b4b3e1cfc74fa780a59b720a",
        "c3_exact_one_spike_AND_expected_AND_support": (
            "06b3e5a7a05d5c9a2f1633d74061f1605c2aab856b315026ad697245bac84964"
        ),
    }
    for stage, expected in expected_stage_semantic.items():
        require(semantic[stage].hexdigest() == expected, ("stage semantic changed", stage))
    for left, right in zip(STAGES, STAGES[1:]):
        require(occurrences[right] <= occurrences[left], ("occurrence monotonicity", left, right))
        require(shapes[right] <= shapes[left], ("shape monotonicity", left, right))
        require(rows[right] <= rows[left], ("row monotonicity", left, right))
    # Equality is killed by the theorem, and the finite row universe happens
    # not to attain its arithmetic endpoint at all.
    require(
        (equality_pairs, equality_support_occurrences, minimum_equality)
        == (0, 0, None),
        "mean equality census changed",
    )
    require(c3_endpoint_occurrences > 0, "no positive c3 allowance endpoint")
    require(minimum_final is not None, "final stage unexpectedly empty")
    require(minimum_final_c3 is not None, "final c3 slice unexpectedly empty")
    c3_endpoint_shape_count = sum(
        c3_denominator_distribution(D)[feature]
        for D, feature in c3_endpoint_features
    )

    print("LRC14 k=3/four-drift expected-spike x c3 exact GF scout")
    print(f"body_count={body_count}")
    print(f"body_divisor_rows={body_divisor_rows}")
    print(f"support_rows={sum(map(len, by_divisor.values()))}")
    print(f"support_divisors={len(by_divisor)}")
    print(f"literal_distribution_D_controls={literal_D_cases}")
    print(f"literal_distribution_shapes={literal_shapes}")
    print(f"one_spike_residue_controls={residue_controls}")
    print(f"all_c_features={feature_count}")
    print(f"c3_denominator_features={c3_feature_count}")
    print(f"all_c_feature_semantic_sha256={feature_semantic.hexdigest()}")
    print(f"c3_feature_semantic_sha256={c3_feature_semantic.hexdigest()}")
    print(f"mean_equality_row_c_pairs={equality_pairs}")
    print(f"mean_equality_by_c={equality_by_c}")
    print(f"mean_equality_support_occurrences_killed={equality_support_occurrences}")
    print(f"mean_equality_support_by_c={equality_support_by_c}")
    print(f"minimum_mean_equality={minimum_equality}")
    print(f"c3_positive_allowance_endpoint_occurrences={c3_endpoint_occurrences}")
    print(
        "c3_positive_allowance_endpoint_feature_types="
        f"{len(c3_endpoint_features)}"
    )
    print(f"c3_positive_allowance_endpoint_shapes={c3_endpoint_shape_count}")
    print(f"minimum_c3_positive_allowance_endpoint={minimum_c3_endpoint}")
    print(f"final_c3_occurrences_by_d_mod_7={final_c3_occurrences_by_residue}")
    print(f"final_c3_shapes_by_d_mod_7={final_c3_shapes_by_residue}")
    print(f"minimum_final_survivor={minimum_final}")
    print(f"minimum_final_c3_survivor={minimum_final_c3}")
    for stage in STAGES:
        print(
            f"stage={stage},occurrences={occurrences[stage]},"
            f"shapes={shapes[stage]},rows={len(rows[stage])},"
            f"bodies={len(bodies[stage])},divisors={len(divisors[stage])},"
            f"semantic_sha256={semantic[stage].hexdigest()}"
        )
        print(f"stage={stage},c_occurrences={c_occurrences[stage]}")
        print(f"stage={stage},c_shapes={c_shapes[stage]}")
        print(
            f"stage={stage},pattern_occurrences="
            f"{Counter(dict(sorted(pattern_occurrences[stage].items())))}"
        )
        print(
            f"stage={stage},pattern_shapes="
            f"{Counter(dict(sorted(pattern_shapes[stage].items())))}"
        )
    print("mean_predicate=N_c==0 OR 55*N_c<13*(4-c)*q")
    print("c3_predicate=N_3<=floor(d/7)*(q/d)+I[d_mod_7>=5]*(q/d)")
    print("intersection_scope=separate necessary screens; shared-bit intersection not optimized")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()

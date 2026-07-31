#!/usr/bin/env python3
"""UNFINISHED all-c floor/exception refinement of the k3 expected-spike screen.

STATUS: exploratory scratch sidecar; there is no completed all-D census.
The exact feature GF still has a high-D state explosion, so this file is not
a proved/computed dependency and must not be promoted in its present form.

This extends the exact c=3 one-spike law to every c.  For every spike
denominator d|q=D/7 write d=7a+r and w=q/d.  Its filled-fibre support is

    a*w + w*X_d(u),                  Haar(X_d=1)=r/7.

For a feature, ``R0`` is the sum of the guaranteed a*w terms and the
optional list consists of the labelled ``(w,r)`` pairs.  The exact general
upward-event theorem computes the largest threshold whose arbitrary joint
law, with marginals r/7, can occupy mass strictly above 55/91.  Hence the
necessary all-c condition is

    N_c <= R0 + status_need_limit((w_i),(r_i/7)).

The exact-lcm GF retains small pattern, uniform count c, *large-only*
ceiling capacity, R0, and the sorted optional list.  Keeping large capacity
separate is essential: d=2,3,4 capacity is paid only through the independent
row-specific support-status screen, never twice.
"""

from collections import Counter
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb, gcd
from pathlib import Path


PRIMARY_PATH = Path(__file__).with_name(
    "lrc14_k3_four_drift_full_activity_q7_gf.py"
)
spec = spec_from_file_location("lrc14_k3_expected_spike", PRIMARY_PATH)
primary = module_from_spec(spec)
spec.loader.exec_module(primary)
base = primary.base
q7 = primary.q7
combined = primary.combined
Q = primary.Q


STAGES = (
    "raw",
    "support_status",
    "expected_spike_AND_support",
    "c3_exact_one_spike_AND_expected_AND_support",
    "all_c_floor_exception_AND_support",
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def full_distribution(D):
    """Exact-lcm GF on (pattern,c,large_cap,R0,optional)."""
    q = D // 7
    result = Counter()
    for E in base.divisors_of(D):
        sign = base.mobius(D // E)
        if not sign:
            continue
        groups = Counter()
        for divisor in (d for d in base.divisors_of(E) if d > 1):
            pattern = (
                int(divisor == 2),
                int(divisor == 3),
                int(divisor == 4),
            )
            if q % divisor:
                require(divisor % 7 == 0, ("bad uniform divisor", D, divisor))
                require(not any(pattern), ("small divisor became uniform", D, divisor))
                feature = (pattern, 1, q, 0, ())
            else:
                a, r = divmod(divisor, 7)
                width = q // divisor
                large_capacity = 0
                if not any(pattern):
                    large_capacity = 7 * (a + bool(r)) * width
                optional = ((width, r),) if r else ()
                feature = (pattern, 0, large_capacity, a * width, optional)
                require(
                    Q(a * width) + Q(r, 7) * width == Q(q, 7),
                    ("spike expectation changed", D, divisor),
                )
            groups[feature] += 1

        # State: used,m2,m3,m4,c,large_cap,R0,sorted optional list.
        states = {(0, 0, 0, 0, 0, 0, 0, ()): 1}
        for feature, alphabet_size in groups.items():
            pattern, unit_c, unit_cap, unit_R0, unit_optional = feature
            additions = Counter()
            for state, multiplicity in tuple(states.items()):
                used, m2, m3, m4, c, capacity, R0, optional = state
                for copies in range(1, 5 - used):
                    new_optional = list(optional)
                    new_optional.extend(unit_optional * copies)
                    new_optional.sort()
                    additions[
                        (
                            used + copies,
                            m2 + copies * pattern[0],
                            m3 + copies * pattern[1],
                            m4 + copies * pattern[2],
                            c + copies * unit_c,
                            capacity + copies * unit_cap,
                            R0 + copies * unit_R0,
                            tuple(new_optional),
                        )
                    ] += multiplicity * comb(alphabet_size + copies - 1, copies)
            for state, multiplicity in additions.items():
                states[state] = states.get(state, 0) + multiplicity
        for state, multiplicity in states.items():
            used, m2, m3, m4, c, capacity, R0, optional = state
            if used == 4:
                result[((m2, m3, m4), c, capacity, R0, optional)] += (
                    sign * multiplicity
                )

    require(all(value >= 0 for value in result.values()), ("negative full GF", D))
    result += Counter()
    collapsed = Counter()
    for (pattern, c, capacity, _R0, _optional), multiplicity in result.items():
        collapsed[(pattern, c, capacity)] += multiplicity
    require(
        collapsed == primary.expected_distribution(D),
        ("full GF does not collapse to expected GF", D),
    )
    return result


@lru_cache(maxsize=None)
def clutter_mass(minimal_masks, marginals):
    """Exact upward-event maximum, cached by Boolean type not raw weights."""
    count = len(marginals)
    rows = tuple(
        tuple((mask >> index) & 1 for index in range(count))
        for mask in minimal_masks
    )
    constraints = [(row, Q(1)) for row in rows]
    constraints.extend(
        (
            tuple(int(index == coordinate) for index in range(count)),
            Q(0),
        )
        for coordinate in range(count)
    )
    optimum = None
    for chosen in combinations(range(len(constraints)), count):
        matrix = tuple(constraints[index][0] for index in chosen)
        rhs = tuple(constraints[index][1] for index in chosen)
        point = base.solve_square(matrix, rhs)
        if point is None or any(value < 0 for value in point):
            continue
        if any(
            sum(value * coefficient for value, coefficient in zip(point, row)) < 1
            for row in rows
        ):
            continue
        objective = sum(
            marginal * value for marginal, value in zip(marginals, point)
        )
        if optimum is None or objective < optimum:
            optimum = objective
    require(optimum is not None, ("clutter cover has no vertex", minimal_masks, marginals))
    return min(Q(1), optimum)


@lru_cache(maxsize=None)
def normalized_optional_limit(optional):
    """Largest allowed threshold, caching LPs by the induced Boolean clutter."""
    weights = tuple(width for width, _r in optional)
    marginals = tuple(Q(r, 7) for _width, r in optional)
    count = len(weights)
    require(1 <= count <= 3, ("unexpected optional arity", optional))
    thresholds = sorted(
        {
            sum(weights[index] for index in range(count) if (mask >> index) & 1)
            for mask in range(1, 1 << count)
        }
    )
    answer = 0
    for threshold in thresholds:
        minimal = []
        for mask in range(1, 1 << count):
            weight = sum(
                weights[index] for index in range(count) if (mask >> index) & 1
            )
            if weight < threshold:
                continue
            if any(
                weight - weights[index] >= threshold
                for index in range(count)
                if (mask >> index) & 1
            ):
                continue
            minimal.append(mask)
        mass = clutter_mass(tuple(minimal), marginals)
        if mass > base.THREE_SAFE_FLOOR:
            answer = threshold
        else:
            break
    return answer


def optional_limit(optional):
    """Use homogeneity to share exact LPs across scaled denominators."""
    if not optional:
        return 0
    scale = 0
    for width, _r in optional:
        scale = gcd(scale, width)
    normalized = tuple((width // scale, r) for width, r in optional)
    return scale * normalized_optional_limit(normalized)


def optional_limit_controls():
    """Compare the Boolean-type engine with the inherited exact LP."""
    cases = 0
    bank = (
        ((1, 1),),
        ((2, 4),),
        ((3, 6),),
        ((1, 2), (1, 4)),
        ((2, 1), (3, 5)),
        ((1, 2), (2, 3), (4, 6)),
        ((2, 6), (3, 1), (5, 5)),
    )
    for optional in bank:
        weights = tuple(width for width, _r in optional)
        marginals = tuple(Q(r, 7) for _width, r in optional)
        require(
            optional_limit(optional) == base.status_need_limit(weights, marginals),
            ("Boolean-type optional limit mismatch", optional),
        )
        cases += 1
    return cases


def proof_distribution(D, full=None):
    """Collapse to (pattern,c,large_cap,exact spike allowance)."""
    result = Counter()
    if full is None:
        full = full_distribution(D)
    for (pattern, c, capacity, R0, optional), multiplicity in full.items():
        allowance = R0 + optional_limit(optional)
        result[(pattern, c, capacity, allowance)] += multiplicity

    # The c=3 slice must be exactly the independently extracted one-spike GF.
    c3_direct = Counter()
    for (pattern, _d, capacity, allowance), multiplicity in (
        primary.c3_denominator_distribution(D).items()
    ):
        c3_direct[(pattern, 3, capacity, allowance)] += multiplicity
    c3_here = Counter(
        {feature: multiplicity for feature, multiplicity in result.items() if feature[1] == 3}
    )
    require(c3_here == c3_direct, ("c3 floor law mismatch", D))
    return result


def main():
    base.activity_controls()
    optional_control_cases = optional_limit_controls()
    by_divisor, body_count, body_divisor_rows = q7.build_rows()
    occurrences = Counter()
    shapes = Counter()
    rows = {stage: set() for stage in STAGES}
    bodies = {stage: set() for stage in STAGES}
    divisors = {stage: set() for stage in STAGES}
    c_occurrences = {stage: Counter() for stage in STAGES}
    c_shapes = {stage: Counter() for stage in STAGES}
    semantic = {stage: sha256() for stage in STAGES}
    full_semantic = sha256()
    proof_semantic = sha256()
    full_features = 0
    proof_features = 0
    survivor_features = {stage: set() for stage in STAGES[1:]}
    minimum_all_c = None

    for D in sorted(by_divisor):
        q = D // 7
        full = full_distribution(D)
        proof = proof_distribution(D, full)
        full_features += len(full)
        proof_features += len(proof)
        for feature, multiplicity in sorted(full.items()):
            full_semantic.update(f"{D}|{feature}|{multiplicity}\n".encode())
        for feature, multiplicity in sorted(proof.items()):
            proof_semantic.update(f"{D}|{feature}|{multiplicity}\n".encode())

        raw_shape_count = sum(proof.values())
        shapes["raw"] += raw_shape_count
        raw_c = Counter()
        for (_pattern, c, _capacity, _allowance), multiplicity in proof.items():
            raw_c[c] += multiplicity
        c_shapes["raw"].update(raw_c)

        local_survivors = {stage: set() for stage in STAGES[1:]}
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
            for pattern, _c, _capacity, _allowance in proof:
                if pattern not in support_limits:
                    weights, marginals = base.small_vectors(pattern, D, small_loads)
                    support_limits[pattern] = base.status_need_limit(weights, marginals)

            row_counts = Counter(raw=raw_shape_count)
            row_c = {stage: Counter() for stage in STAGES}
            row_c["raw"] = raw_c.copy()
            for feature, multiplicity in proof.items():
                pattern, c, capacity, allowance = feature
                if capacity + support_limits[pattern] < support_count:
                    continue
                row_counts["support_status"] += multiplicity
                row_c["support_status"][c] += multiplicity
                local_survivors["support_status"].add(feature)

                if not primary.mean_passes(N_by_c[c], c, q):
                    continue
                row_counts["expected_spike_AND_support"] += multiplicity
                row_c["expected_spike_AND_support"][c] += multiplicity
                local_survivors["expected_spike_AND_support"].add(feature)

                if c != 3 or allowance >= N_by_c[c]:
                    row_counts[
                        "c3_exact_one_spike_AND_expected_AND_support"
                    ] += multiplicity
                    row_c[
                        "c3_exact_one_spike_AND_expected_AND_support"
                    ][c] += multiplicity
                    local_survivors[
                        "c3_exact_one_spike_AND_expected_AND_support"
                    ].add(feature)

                if allowance < N_by_c[c]:
                    continue
                require(
                    primary.mean_passes(N_by_c[c], c, q),
                    ("floor survivor failed mean theorem", D, body, feature),
                )
                row_counts["all_c_floor_exception_AND_support"] += multiplicity
                row_c["all_c_floor_exception_AND_support"][c] += multiplicity
                local_survivors["all_c_floor_exception_AND_support"].add(feature)
                candidate = (
                    D,
                    body,
                    L,
                    support_count,
                    pattern,
                    c,
                    capacity,
                    allowance,
                    N_by_c[c],
                )
                if minimum_all_c is None or candidate < minimum_all_c:
                    minimum_all_c = candidate

            row_key = (body, D)
            for stage in STAGES:
                count = row_counts[stage]
                if not count:
                    continue
                occurrences[stage] += count
                c_occurrences[stage].update(row_c[stage])
                rows[stage].add(row_key)
                bodies[stage].add(body)
                divisors[stage].add(D)
                semantic[stage].update(
                    f"{body}|{L}|{D}|{support_count}|{count}\n".encode()
                )

        for stage in STAGES[1:]:
            for feature in local_survivors[stage]:
                multiplicity = proof[feature]
                shapes[stage] += multiplicity
                c_shapes[stage][feature[1]] += multiplicity
                survivor_features[stage].add((D, feature))

    require(occurrences["raw"] == 21357714101, "raw occurrence changed")
    require(shapes["raw"] == 694921995, "raw shapes changed")
    require(occurrences["support_status"] == 13280722299, "support occurrence changed")
    require(shapes["support_status"] == 694254050, "support shapes changed")
    require(occurrences["expected_spike_AND_support"] == 2934202044, "mean occurrence changed")
    require(shapes["expected_spike_AND_support"] == 400005870, "mean shapes changed")
    require(
        occurrences["c3_exact_one_spike_AND_expected_AND_support"] == 2548901482,
        "isolated c3 occurrence changed",
    )
    require(
        shapes["c3_exact_one_spike_AND_expected_AND_support"] == 398241574,
        "isolated c3 shapes changed",
    )
    for left, right in zip(STAGES, STAGES[1:]):
        require(occurrences[right] <= occurrences[left], ("occurrence monotonicity", left, right))
        require(shapes[right] <= shapes[left], ("shape monotonicity", left, right))
        require(rows[right] <= rows[left], ("row monotonicity", left, right))

    print("LRC14 k3/four-drift all-c floor/exception GF scout")
    print(f"body_count={body_count}")
    print(f"body_divisor_rows={body_divisor_rows}")
    print(f"support_rows={sum(map(len, by_divisor.values()))}")
    print(f"support_divisors={len(by_divisor)}")
    print(f"optional_limit_control_cases={optional_control_cases}")
    print(f"full_features={full_features}")
    print(f"proof_features={proof_features}")
    print(f"full_feature_semantic_sha256={full_semantic.hexdigest()}")
    print(f"proof_feature_semantic_sha256={proof_semantic.hexdigest()}")
    print(f"minimum_all_c_survivor={minimum_all_c}")
    for stage in STAGES:
        print(
            f"stage={stage},occurrences={occurrences[stage]},"
            f"shapes={shapes[stage]},rows={len(rows[stage])},"
            f"bodies={len(bodies[stage])},divisors={len(divisors[stage])},"
            f"semantic_sha256={semantic[stage].hexdigest()}"
        )
        print(f"stage={stage},c_occurrences={c_occurrences[stage]}")
        print(f"stage={stage},c_shapes={c_shapes[stage]}")
    print(f"normalized_optional_limit_cache={normalized_optional_limit.cache_info()}")
    print(f"clutter_mass_cache={clutter_mass.cache_info()}")
    print("all_c_predicate=N_c<=R0+exact_fractional_cover_optional_limit")
    print("intersection_scope=separate necessary screens; shared-bit intersection not optimized")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()

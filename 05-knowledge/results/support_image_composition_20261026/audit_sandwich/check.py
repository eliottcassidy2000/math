#!/usr/bin/env python3
"""Independent literal five-factor sandwich audit; no producer imports."""
from collections import Counter, defaultdict
from itertools import permutations, product
import json

CHECKS = Counter()


def require(ok, tag):
    CHECKS[tag] += 1
    if not ok:
        raise RuntimeError(tag)


def tournament(n, mask):
    rows = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            a, b = (i, j) if mask >> k & 1 else (j, i)
            rows[a] |= 1 << b
            k += 1
    return tuple(rows)


def induced(rows, vertices):
    return tuple(sum(1 << j for j, v in enumerate(vertices) if rows[u] >> v & 1)
                 for u in vertices)


def arc(rows, u, v):
    return rows[u] >> v & 1


def hinge(left, right):
    """Literal binary definition, rather than marker-degree formula."""
    lrows, ls, lt = left
    rrows, rs, rt = right
    n, m = len(lrows), len(rrows)
    out = list(lrows) + [(row << n) | ((1 << n) - 1) for row in rrows]
    out[lt] |= 1 << (n + rs)
    out[n + rs] &= ~(1 << lt)
    return tuple(out), ls, n + rt


E = ((0,), 0, 0)
A2 = hinge(E, E)
A2OP = ((0, 1), 0, 1)


def encode(state):
    result = E
    for factor in (E, state, E, E):
        result = hinge(result, factor)
    return result[0]


def rooted_key(state):
    rows, s, t = state
    middle = tuple(v for v in range(len(rows)) if v not in (s, t))
    return min(induced(rows, (s,) + order + (t,)) for order in permutations(middle))


def degree_orders(rows):
    groups = defaultdict(list)
    for v, row in enumerate(rows):
        groups[row.bit_count()].append(v)
    for blocks in product(*(permutations(groups[d]) for d in sorted(groups))):
        yield tuple(v for block in blocks for v in block)


def unrooted_key(rows):
    return min(induced(rows, order) for order in degree_orders(rows))


def degree_preserving_maps(rows):
    groups = defaultdict(list)
    for v, row in enumerate(rows):
        groups[row.bit_count()].append(v)
    fixed_groups = [tuple(groups[d]) for d in sorted(groups)]
    for images in product(*(permutations(group) for group in fixed_groups)):
        result = list(range(len(rows)))
        for source, target in zip(fixed_groups, images):
            for v, w in zip(source, target):
                result[v] = w
        yield tuple(result)


def is_auto(rows, p):
    return all(arc(rows, u, v) == arc(rows, p[u], p[v])
               for u in range(len(rows)) for v in range(len(rows)))


def autos(rows):
    return {p for p in degree_preserving_maps(rows) if is_auto(rows, p)}


def compose(a, b):
    return tuple(a[b[v]] for v in range(len(a)))


def rotation(n, vertices):
    result = list(range(n))
    for v, w in zip(vertices, vertices[1:] + vertices[:1]):
        result[v] = w
    return tuple(result)


def decode_candidates(rows, reencode=True):
    """Read chains literally, then demand all prescribed incidences."""
    n = len(rows)
    if n < 6:
        return []
    L = [v for v, row in enumerate(rows) if row.bit_count() == 1]
    R = [v for v, row in enumerate(rows) if row.bit_count() == n - 2]
    if len(L) not in (2, 3) or len(R) not in (2, 3):
        return []
    out = []
    for a0 in L:
        a1 = rows[a0].bit_length() - 1
        if a1 not in L:
            continue
        s = rows[a1].bit_length() - 1
        for b1 in R:
            b0 = next(v for v in range(n) if v != b1 and arc(rows, v, b1))
            if b0 not in R:
                continue
            t = next(v for v in range(n) if v != b0 and arc(rows, v, b0))
            markers = {a0, a1, b0, b1}
            if len(markers) != 4 or s in markers or t in markers or s == t:
                continue
            core = tuple(v for v in range(n) if v not in markers)
            state = induced(rows, core), core.index(s), core.index(t)
            order = (a0, a1) + core + (b0, b1)
            if not reencode or encode(state) == induced(rows, order):
                out.append(state)
    return out


def intrinsic_quotient(rows):
    """Separate choice-free-chain / homogeneous-triangle quotient route."""
    N = len(rows)
    low = {v for v, row in enumerate(rows) if row.bit_count() == 1}
    high = {v for v, row in enumerate(rows) if row.bit_count() == N - 2}
    if len(low) == 3:
        s = min(low)
        left_deleted = low - {s}
        require(all(arc(rows, w, v) for v in low for w in set(range(N)) - low),
                "left_triangle_homogeneity")
        require(all(sum(arc(rows, u, v) for v in low) == 1 for u in low),
                "left_triangle_cycle")
    else:
        require(len(low) == 2, "left_chain_size")
        a1 = next(v for v in low if not any(arc(rows, v, w) for w in low))
        s = rows[a1].bit_length() - 1
        left_deleted = low
    if len(high) == 3:
        t = min(high)
        right_deleted = high - {t}
        require(all(arc(rows, v, w) for v in high for w in set(range(N)) - high),
                "right_triangle_homogeneity")
        require(all(sum(arc(rows, u, v) for v in high) == 1 for u in high),
                "right_triangle_cycle")
    else:
        require(len(high) == 2, "right_chain_size")
        b0 = next(v for v in high if any(arc(rows, v, w) for w in high))
        t = next(v for v in range(N) if v != b0 and arc(rows, v, b0))
        right_deleted = high
    core = tuple(v for v in range(N) if v not in left_deleted | right_deleted)
    return induced(rows, core), core.index(s), core.index(t)


def verify_state(state):
    core, s, t = state
    n = len(core)
    output = encode(state)
    N = n + 4
    entrance_sink = (core[s] == 0)
    exit_source = (core[t].bit_count() == n - 1)
    exponent = entrance_sink + exit_source
    degrees = [row.bit_count() for row in output]
    require(degrees[0] == degrees[1] == 1, "left_marker_degrees")
    require(degrees[-2] == degrees[-1] == n + 2, "corrected_right_marker_degrees")
    for v in range(n):
        require(degrees[v + 2] == core[v].bit_count() + 2 - (v == s) + (v == t),
                "core_degree_formula")
    candidates = decode_candidates(output)
    require(len(candidates) == 3 ** exponent, "all_marker_choices")
    require(all(rooted_key(candidate) == rooted_key(state) for candidate in candidates),
            "decoded_joint_type")
    require(rooted_key(intrinsic_quotient(output)) == rooted_key(state),
            "independent_intrinsic_quotient")
    coregroup = {p for p in autos(core) if p[s] == s and p[t] == t}
    extended = []
    for p in coregroup:
        extended.append((0, 1) + tuple(v + 2 for v in p) + (N - 2, N - 1))
    identity = tuple(range(N))
    left = rotation(N, (0, 1, s + 2)) if entrance_sink else identity
    right = rotation(N, (t + 2, N - 2, N - 1)) if exit_source else identity
    leftpowers = {identity, left, compose(left, left)}
    rightpowers = {identity, right, compose(right, right)}
    expected = {compose(p, compose(l, r)) for p in extended for l in leftpowers for r in rightpowers}
    actual = autos(output)
    require(actual == expected, "full_automorphism_group_equality")
    require(len(actual) == len(coregroup) * 3 ** exponent, "automorphism_order")
    for p in extended:
        require(compose(p, left) == compose(left, p) and
                compose(p, right) == compose(right, p), "direct_product_commutation")
    require(compose(left, right) == compose(right, left), "disjoint_rotation_commutation")
    return output


def main():
    states = {}
    counts = []
    cases = Counter()
    for n in range(2, 6):
        bykey = {}
        for mask in range(1 << (n * (n - 1) // 2)):
            state = tournament(n, mask), 0, n - 1
            bykey.setdefault(rooted_key(state), state)
        states[n] = list(bykey.values())
        counts.append(len(bykey))
        images = {}
        for key, state in bykey.items():
            rows, s, t = state
            cases[int(rows[s] == 0), int(rows[t].bit_count() == n - 1)] += 1
            encoded = verify_state(state)
            imagekey = unrooted_key(encoded)
            require(imagekey not in images or images[imagekey] == key, "all_core_injectivity")
            images[imagekey] = key
    require(counts == [2, 8, 32, 176], "complete_rooted_counts")
    truth = set()
    for state in states[2]:
        output = encode(state)
        truth.update(induced(output, p) for p in permutations(range(6)))
    accepted = 0
    degree_candidates = 0
    step3_candidates = 0
    for mask in range(1 << 15):
        rows = tournament(6, mask)
        degree_only = (sum(row.bit_count() == 1 for row in rows) in (2, 3) and
                       sum(row.bit_count() == 4 for row in rows) in (2, 3))
        degree_candidates += degree_only
        candidates = decode_candidates(rows)
        accepted += bool(candidates)
        step3_candidates += bool(decode_candidates(rows, False))
        require(bool(candidates) == (rows in truth), "exact_order6_image_membership")
    require(accepted == len(truth) == 800, "order6_image_cardinality")
    hostile = (32, 1, 3, 7, 15, 30)
    require(sum(row.bit_count() == 1 for row in hostile) == 2 and
            sum(row.bit_count() == 4 for row in hostile) == 2 and
            not decode_candidates(hostile), "degree_cardinality_hostile")
    C, X, Y = ((2, 4, 1), 0, 2), ((4, 5, 0), 0, 2), ((6, 0, 2), 0, 2)
    require(len({rooted_key(C), rooted_key(X), rooted_key(Y)}) == 3,
            "distinct_minimality_source_types")
    contexts = [(E, E), (A2, E), (A2OP, E), (E, A2), (E, A2OP)]
    fibres = [(C, X, Y), (C, X), (C, X, Y), (C, Y), (C, X, Y)]
    collision_witnesses = []
    for (left, right), fibre in zip(contexts, fibres):
        literal_outputs = [hinge(hinge(left, core), right)[0] for core in fibre]
        outputs = [unrooted_key(rows) for rows in literal_outputs]
        require(len(set(outputs)) == 1, "each_smaller_context_collision")
        source = literal_outputs[0]
        maps = []
        for target in literal_outputs[1:]:
            witness = next(p for p in permutations(range(len(source)))
                           if all(arc(source, u, v) == arc(target, p[u], p[v])
                                  for u in range(len(source)) for v in range(len(source))))
            maps.append(witness)
        collision_witnesses.append(maps)
        order2 = {unrooted_key(hinge(hinge(left, core), right)[0]) for core in states[2]}
        require(len(order2) == 2, "no_smaller_core_collision")
    # Literal hinge associativity and actual encode-decode transport, all small source pairs.
    small = states[2] + states[3]
    for h in small:
        for k in small:
            joined = hinge(h, k)
            recovered = decode_candidates(encode(joined))
            require(bool(recovered) and all(rooted_key(x) == rooted_key(joined) for x in recovered),
                    "decode_operate_encode_composition")
            require(hinge(hinge(E, h), k) == hinge(E, hinge(h, k)), "literal_hinge_associativity")
    print(json.dumps({"status": "PASS", "checks": sum(CHECKS.values()),
                      "checks_by_kind": dict(sorted(CHECKS.items())),
                      "rooted_classes_n2_through5": counts,
                      "boundary_cases": {str(k): v for k, v in sorted(cases.items())},
                      "order6_orientations": 32768, "order6_degree_cardinality_candidates": degree_candidates,
                      "order6_distinct_valid_chain_candidates": step3_candidates,
                      "order6_exact_image": accepted,
                      "smaller_context_isomorphisms_from_C": collision_witnesses,
                      "scope": "independent literal hinge/decoder/group census; no producer import"},
                     indent=2, sort_keys=True))


if __name__ == '__main__':
    main()

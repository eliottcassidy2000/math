#!/usr/bin/env python3
"""Exact controls for THM-2608's alternative-rail transition audit.

The computation rebuilds THM-2592's common-x tensor, but chooses the
root-six theta-zero rail whenever it is available.  It then distinguishes:

* the future-clock marginal, which is a dense but noncomposable q->r cospan;
* every C7-equivariant clock matching ell5=ell4+c;
* the thirteen coefficient-level root-reference choices exposed by THM-2601;
* the full formal 91-state relation obtained by allowing arbitrary output
  clocks, which is not asserted to be one physical seven-edge carrier.

All overlap arithmetic is inherited exact integer arithmetic from THM-2592.
All relation products below are Boolean products of exact positive supports.
"""

from collections import Counter
from math import gcd

import lrc14_fallback_rail_digit_diagonal_pullback_thm2592 as old


P = old.P
Q7 = old.Q7
T = old.T
R = old.R

T_VALUES = (9, 3, 10, 1, 8, 0, 7, 2, 12, 11, 4, 6, 5)
SUCCESSOR = (7, 8, 12, 10, 6, 9, 5, 2, 0, 3, 1, 4, 11)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def build_alternative_joint():
    """Return the exact root-six-priority THM-2592 tensor and rail labels."""
    module = old.load_present_module()
    require(module.W == old.base.W, "typed row changed")

    word = module.build_word_Ta()
    word_prefix = [[None] * P for _ in range(Q7)]
    for ell in range(Q7):
        q_ell = module.subtract_comb(
            word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        for h in range(P):
            digit = old.sat.intersect_interval(
                q_ell, h * T // P, (h + 1) * T // P
            )
            word_prefix[ell][h] = module.make_prefix(digit)
    require(all(word_prefix[ell][0][2][-1] == 0 for ell in range(Q7)),
            "future digit-zero hole changed")
    require(all(word_prefix[ell][6][2][-1] > 0 for ell in range(Q7)),
            "future digit-six support changed")

    e4 = old.base.build_set(old.base.PAT_E, old.base.ZELL)
    qb = old.base.build_set(old.host.PAT_QB, old.base.ZELL)
    ust, uv, vst, vv = old.rail.packet_profiles(e4, qb)
    _, _, k_tensor = old.rail.exact_tensor(e4, qb)
    owner = old.base.clock_cells(P, Q7, T, P * P)
    deep = old.rail.deep_cells()
    arrival = [[(v * (T // P), (v + 1) * (T // P))] for v in range(P)]

    rail_bank = []
    for s4 in range(1, P):
        rvst, rvv = old.base.rotate_profile(vst, vv, s4 * (T // P), T)
        ps, pv, _, _ = old.base.product_cum(ust, uv, rvst, rvv, T)
        for ell4 in range(Q7):
            choices = [
                (v, t) for v, t in ((6, 12), (0, 0))
                if k_tensor[s4][v][t][ell4] > 0
            ]
            require(choices, "alternative theta-zero rail disappeared")
            v, t = choices[0]
            cell = old.intersect_sorted(
                old.intersect_sorted(owner[ell4], deep[t]), arrival[v]
            )
            pieces = old.profile_on_intervals(cell, ps, pv)
            numerator = P * sum(w * (b - a) for a, b, w in pieces)
            require(numerator == k_tensor[s4][v][t][ell4],
                    "alternative rail mass mismatch")
            rail_bank.append((s4, ell4, v, t, pieces))
    require(len(rail_bank) == 84, "alternative rail census changed")

    f_cache = {}
    f_starts = {}
    for s5 in range(P):
        for ell5 in range(Q7):
            present = module.build_F(ell5, s5)
            f_cache[ell5, s5] = present
            f_starts[ell5, s5] = [a for a, _ in present]

    joint = [[[[0] * P for _ in range(Q7)] for _ in range(P)]
             for _ in rail_bank]
    content = 0
    for j, (_, _, v, _, rail_pieces) in enumerate(rail_bank):
        for q in range(P):
            s5 = (-q) % P
            for ell5 in range(Q7):
                present = old.intersect_weighted_union(
                    rail_pieces, f_cache[ell5, s5], f_starts[ell5, s5]
                )
                if not present:
                    continue
                prefix = word_prefix[ell5][v]
                for r in range(1, P):
                    probed = old.intersect_weighted_comb(
                        present, module.C3, 182,
                        14 * r - 13, 14 * r + 13,
                    )
                    value = old.delayed_weighted_numerator(probed, prefix)
                    joint[j][q][ell5][r] = value
                    content = gcd(content, value)
    require(content > 0, "alternative joint bank vanished")
    return joint, rail_bank, content


def relation_product(left, right):
    """Boolean product; matrices are tuples of integer row bitmasks."""
    out = []
    for row in left:
        image = 0
        active = row
        while active:
            bit = active & -active
            image |= right[bit.bit_length() - 1]
            active -= bit
        out.append(image)
    return tuple(out)


def relation_power(matrix, exponent):
    identity = tuple(1 << i for i in range(len(matrix)))
    out = identity
    base = matrix
    n = exponent
    while n:
        if n & 1:
            out = relation_product(out, base)
        base = relation_product(base, base)
        n //= 2
    return out


def support_size(matrix):
    return sum(row.bit_count() for row in matrix)


def rotate13(mask, step):
    step %= P
    all_bits = (1 << P) - 1
    if step == 0:
        return mask
    return ((mask << step) | (mask >> (P - step))) & all_bits


def root_matrix(joint, rail_index, ell5=None, reference_offset=0):
    """q-row to q'-column, with q'=reference_offset+r."""
    rows = []
    for q in range(P):
        mask = 0
        clocks = range(Q7) if ell5 is None else (ell5,)
        for future_clock in clocks:
            for r in range(P):
                if joint[rail_index][q][future_clock][r] > 0:
                    mask |= 1 << r
        rows.append(rotate13(mask, reference_offset))
    return tuple(rows)


def difference_returns(matrix):
    return tuple(
        sum((matrix[q] >> ((q + delta) % P)) & 1 for q in range(P))
        for delta in range(P)
    )


def scalar_reference_atlas():
    require(set(T_VALUES) == set(range(P)), "THM-2601 scalar lost bijectivity")
    q_by_t = [None] * P
    for q, value in enumerate(T_VALUES):
        q_by_t[value] = q
    require(tuple(SUCCESSOR[T_VALUES[q]] for q in range(P))
            == tuple(T_VALUES[(q + 1) % P] for q in range(P)),
            "THM-2601 successor covariance changed")

    atlas = []
    for c in range(P):
        q0 = q_by_t[c]
        phi = tuple(T_VALUES[(q0 + r) % P] for r in range(P))
        require(set(phi) == set(range(P)), "reference map is not bijective")
        require(all(phi[(r + 1) % P] == SUCCESSOR[phi[r]]
                    for r in range(P)),
                "reference map is not successor-equivariant")
        atlas.append((c, q0, phi))
    require(len({phi for _, _, phi in atlas}) == P,
            "reference atlas did not have thirteen choices")
    return tuple(atlas)


def combined_clock_relation(joint, lookup, s4, reference_offset):
    """Formal 91-state relation (ell4,q)->(ell5,q0+r)."""
    rows = []
    for ell4 in range(Q7):
        j = lookup[s4, ell4]
        for q in range(P):
            mask = 0
            for ell5 in range(Q7):
                for r in range(P):
                    if joint[j][q][ell5][r] > 0:
                        qp = (reference_offset + r) % P
                        mask |= 1 << (P * ell5 + qp)
            rows.append(mask)
    return tuple(rows)


def combined_twisted_returns(matrix7):
    """Return counts at the same clock with q correction -7a."""
    out = []
    for a in range(1, P):
        count = 0
        for ell in range(Q7):
            for q in range(P):
                target_q = (q - Q7 * a) % P
                target = P * ell + target_q
                count += (matrix7[P * ell + q] >> target) & 1
        out.append(count)
    return tuple(out)


def main():
    joint, rail_bank, content = build_alternative_joint()
    lookup = {(s4, ell4): j
              for j, (s4, ell4, _, _, _) in enumerate(rail_bank)}
    atlas = scalar_reference_atlas()

    positive_pairs = sum(
        any(joint[j][q][ell5][r] > 0
            for ell5 in range(Q7) for r in range(P))
        for j in range(84) for q in range(P)
    )
    positive_rails = sum(
        any(joint[j][q][ell5][r] > 0
            for q in range(P) for ell5 in range(Q7) for r in range(P))
        for j in range(84)
    )
    require(content == 4244240, "alternative global content changed")
    require(positive_pairs == 1012, "alternative positive-pair census changed")
    require(positive_rails == 79, "alternative positive-rail census changed")

    marginal_summary = {}
    formal_reference_histogram = Counter()
    for s4 in range(1, P):
        base_matrices = [root_matrix(joint, lookup[s4, ell4])
                         for ell4 in range(Q7)]
        base_product = tuple(1 << q for q in range(P))
        for matrix in base_matrices:
            base_product = relation_product(base_product, matrix)
        marginal_summary[s4] = (
            tuple(support_size(matrix) for matrix in base_matrices),
            support_size(base_product),
            difference_returns(base_product),
        )
        for _, q0, _ in atlas:
            product = tuple(1 << q for q in range(P))
            for ell4 in range(Q7):
                matrix = root_matrix(
                    joint, lookup[s4, ell4], reference_offset=q0
                )
                product = relation_product(product, matrix)
            formal_reference_histogram[
                (support_size(product), difference_returns(product))
            ] += 1

    good_s = tuple(s for s in range(1, P) if s != 6)
    good_marginal = (
        (156, 156, 144, 156, 156, 156, 156),
        156,
        (12,) * P,
    )
    bad_marginal = ((132, 132, 0, 0, 0, 0, 0), 0, (0,) * P)
    require(all(marginal_summary[s] == good_marginal for s in good_s),
            "dense future-clock marginal changed")
    require(marginal_summary[6] == bad_marginal,
            "exceptional displacement-six marginal changed")
    require(formal_reference_histogram == Counter({
        (156, (12,) * P): 143,
        (0, (0,) * P): 13,
    }), "thirteen-reference marginal prism changed")

    clock_signature_histogram = Counter()
    lawful_zero_products = 0
    for s4 in range(1, P):
        for clock_step in range(1, Q7):
            order = tuple((k * clock_step) % Q7 for k in range(Q7))
            matrices = [
                root_matrix(
                    joint, lookup[s4, ell4],
                    ell5=(ell4 + clock_step) % Q7,
                )
                for ell4 in order
            ]
            product = tuple(1 << q for q in range(P))
            for matrix in matrices:
                product = relation_product(product, matrix)
            signature = tuple(support_size(matrix) for matrix in matrices)
            clock_signature_histogram[signature] += 1
            lawful_zero_products += int(support_size(product) == 0)

    expected_clock_histogram = Counter({
        (156, 132, 132, 0, 0, 0, 144): 11,
        (0, 0, 0, 0, 144, 0, 132): 11,
        (132, 0, 144, 0, 144, 0, 132): 11,
        (0, 132, 0, 0, 0, 132, 132): 11,
        (0, 156, 0, 0, 0, 144, 0): 11,
        (0, 0, 0, 144, 132, 0, 0): 11,
        (0, 0, 0, 0, 0, 0, 0): 4,
        (0, 0, 0, 0, 132, 0, 0): 1,
        (132, 0, 0, 0, 0, 0, 0): 1,
    })
    require(lawful_zero_products == 72,
            "a C7-equivariant clock product unexpectedly survived")
    require(clock_signature_histogram == expected_clock_histogram,
            "clock-matched edge signature census changed")

    combined_case_histogram = Counter()
    combined_trace_histogram = Counter()
    torsor_trace_by_s = {s4: [0] * 12 for s4 in range(1, P)}
    for s4 in range(1, P):
        for _, q0, _ in atlas:
            matrix = combined_clock_relation(joint, lookup, s4, q0)
            seventh = relation_power(matrix, Q7)
            edge_support = support_size(matrix)
            seventh_support = support_size(seventh)
            traces = combined_twisted_returns(seventh)
            combined_case_histogram[edge_support, seventh_support] += 1
            combined_trace_histogram.update(traces)
            torsor_trace_by_s[s4] = [x + y for x, y in
                                      zip(torsor_trace_by_s[s4], traces)]
            if s4 == 6:
                require((edge_support, seventh_support, traces)
                        == (264, 0, (0,) * 12),
                        "exceptional combined relation changed")
            else:
                require((edge_support, seventh_support) == (2760, 4320),
                        "formal combined relation support changed")
                require(set(traces) <= {47, 48} and min(traces) == 47,
                        "formal combined twisted return changed")
    require(combined_case_histogram == Counter({
        (2760, 4320): 143,
        (264, 0): 13,
    }), "combined 91-state case census changed")
    require(combined_trace_histogram == Counter({47: 1584, 48: 132, 0: 156}),
            "combined 91-state return census changed")
    require(all(tuple(torsor_trace_by_s[s]) == (612,) * 12 for s in good_s),
            "full reference-torsor trace lost gauge invariance")
    require(tuple(torsor_trace_by_s[6]) == (0,) * 12,
            "exceptional full reference-torsor trace changed")

    print("THM-2608 exact alternative-rail clock-collapse controls")
    print("global_content=", content, sep="")
    print("alternative_positive_pairs=", positive_pairs, "/1092", sep="")
    print("alternative_positive_rails=", positive_rails, "/84", sep="")
    print("marginal_good_s=", good_s,
          " edges=", good_marginal[0],
          " product=", good_marginal[1],
          " returns=", good_marginal[2], sep="")
    print("marginal_exception_s6=", bad_marginal, sep="")
    print("formal_13_reference_histogram=", dict(formal_reference_histogram), sep="")
    print("lawful_clock_products_zero=", lawful_zero_products, "/72", sep="")
    print("lawful_clock_edge_signature_histogram=",
          dict(clock_signature_histogram), sep="")
    print("combined_91_state_case_histogram=", dict(combined_case_histogram), sep="")
    print("combined_91_state_twisted_trace_histogram=",
          dict(combined_trace_histogram), sep="")
    print("combined_1183_state_good_twisted_trace=", (612,) * 12, sep="")
    print("combined_1183_state_exception_s6=", (0,) * 12, sep="")
    print("reference_c_to_q0=", tuple((c, q0) for c, q0, _ in atlas), sep="")
    print("verdict=PASS: future-clock marginal is formally dense; every "
          "C7-equivariant clock-matched product collapses")


if __name__ == "__main__":
    main()

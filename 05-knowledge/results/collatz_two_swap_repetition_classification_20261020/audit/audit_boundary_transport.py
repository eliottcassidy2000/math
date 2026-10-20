#!/usr/bin/env python3
"""Independent exact boundary audit, using integer masks and a rightward scan.

No imports from either existing verifier. All tests remain active under -O.
"""
import json
from collections import Counter
from math import gcd


checks = 0
counts = Counter()


def need(condition, label, context=None):
    global checks
    checks += 1
    if not condition:
        raise ArithmeticError((label, context))


def carry(mask, length):
    """Sum right to left, multiplying the suffix weight only at a one."""
    suffix = 1
    answer = 0
    for position in range(length - 1, -1, -1):
        if (mask >> position) & 1:
            answer += (1 << position) * suffix
            suffix *= 3
    return answer


def rotate(mask, length, places=1):
    return (mask >> places) | ((mask & ((1 << places) - 1)) << (length - places))


def swap(mask, length, edge):
    end = (edge + 1) % length
    if ((mask >> edge) ^ (mask >> end)) & 1:
        return mask ^ (1 << edge) ^ (1 << end)
    return None


def all_swaps(mask, length):
    for edge in range(length):
        output = swap(mask, length, edge)
        if output is not None:
            yield edge, output


def word(mask, length):
    return "".join(str((mask >> i) & 1) for i in range(length))


def necklace(mask, length):
    return min(rotate(mask, length, i) for i in range(length))


def least_period(mask, length):
    return next(k for k in range(1, length + 1)
                if length % k == 0 and rotate(mask, length, k) == mask)


def audit_family(base, base_length, repetitions):
    weight = base.bit_count()
    q = base_length * repetitions
    source = sum(base << (j * base_length) for j in range(repetitions))
    a, c = 1 << base_length, 3 ** weight
    h = (a ** repetitions - c ** repetitions) // (a - c)
    clock = (1 << q) - 3 ** source.bit_count()
    need(carry(source, q) == h * carry(base, base_length), "repetition carry")
    counts["sources"] += 1
    nonminimal = least_period(base, base_length) < base_length
    seen_outputs = set()
    for edge1, intermediate in all_swaps(source, q):
        for edge2, output in all_swaps(intermediate, q):
            counts["paths"] += 1
            counts["nonminimal_base_paths"] += int(nonminimal)
            counts["length_four_paths"] += int(q == 4)
            counts["original_cut_paths"] += int(q - 1 in (edge1, edge2))
            expected = edge1 == edge2 or (
                repetitions == 2 and (edge2 - edge1) % q == base_length)
            content = gcd(h, carry(output, q))
            need((content == h) == expected, "classification",
                 (word(base, base_length), repetitions, edge1, edge2))
            need((source == output) == (edge1 == edge2), "undo iff same edge")
            if expected:
                counts["retained"] += 1
                if edge1 != edge2:
                    counts["mirrored"] += 1
                    edited_base = output & (a - 1)
                    need(output == edited_base | (edited_base << base_length),
                         "literal square")
                    need(edited_base in {v for _, v in all_swaps(base, base_length)},
                         "base actual edit")
            if output not in seen_outputs:
                seen_outputs.add(output)
                # The cocycle verifies exact arithmetic, not just gcd transport.
                first = output & 1
                rotated = rotate(output, q)
                need(2 * carry(rotated, q)
                     == 3 ** first * carry(output, q) + first * clock,
                     "rotation cocycle")
                need(gcd(h, carry(rotated, q)) == content, "gcd cut transport")
                counts["distinct_output_transport_checks"] += 1
            if edge1 == edge2:
                continue
            displacement = (edge2 - edge1) % q
            if displacement in (1, q - 1):
                block = edge1 if displacement == 1 else edge2
                before = rotate(source, q, block)
                after = rotate(output, q, block)
                local = carry(after & 7, 3) - carry(before & 7, 3)
                suffix_weight = (before >> 3).bit_count()
                need(abs(local) in (3, 5), "actual overlap word")
                need(carry(after, q) - carry(before, q) == local * 3 ** suffix_weight,
                     "overlap exact outside multiplier")
                need(content == gcd(h, abs(local)) < h, "overlap depth")
                counts["overlap_factor_" + str(abs(local))] += 1
            else:
                start = edge1 if displacement <= q // 2 else edge2
                other = edge2 if start == edge1 else edge1
                length = (other - start) % q
                before = rotate(source, q, start)
                after = rotate(output, q, start)
                m = ((before >> 2) & ((1 << length) - 1)).bit_count()
                suffix = (before >> (length + 2)).bit_count()
                sign1 = (before & 1) - ((before >> 1) & 1)
                sign2 = ((before >> length) & 1) - ((before >> (length + 1)) & 1)
                binomial = sign1 * 3 ** m + sign2 * (1 << length)
                need(carry(after, q) - carry(before, q) == 3 ** suffix * binomial,
                     "actual disjoint perturbation")
                need(content == gcd(h, binomial), "disjoint exact depth")
                if repetitions >= 3:
                    need(length < (repetitions - 1) * base_length
                         and m <= (repetitions - 1) * weight,
                         "periodic bound without balance")
                    need(0 < abs(binomial) < h, "strict multirepetition bound")
                else:
                    need(0 < abs(binomial) <= h, "double bound")
                    need((abs(binomial) == h) == (
                        length == base_length and m == weight and sign1 == sign2),
                        "all equality conditions")


# A different universe from either pre-existing verifier: add d=6, d=7,
# and base length 9 for squares, while stressing every actual ordered path.
for base_length in range(2, 8):
    for base in range(1, (1 << base_length) - 1):
        counts["bases_in_full_repetition_range"] += 1
        for repetitions in range(2, 8):
            audit_family(base, base_length, repetitions)
for base in range(1, (1 << 9) - 1):
    counts["additional_nine_bit_square_bases"] += 1
    audit_family(base, 9, 2)


# Prove-facing finite checks of the separate elementary primitive-base lemma.
unit_shapes = []
for q in range(2, 201):
    for p in range(1, q):
        clock = (1 << q) - 3 ** p
        if abs(clock) == 1:
            unit_shapes.append([p, q])
        if gcd(p, q) != 1:
            continue
        v0 = sum((((i + 1) * p // q) - (i * p // q)) << i for i in range(q))
        v1 = sum(((((i + 1) * p + 1) // q) - ((i * p + 1) // q)) << i
                 for i in range(q))
        rotation = pow(p, -1, q)
        need(rotate(v0, q, rotation) == v1, "intercept one is a rotation")
        difference = v0 ^ v1
        indices = [i for i in range(q) if (difference >> i) & 1]
        need(len(indices) == 2 and indices[1] == indices[0] + 1,
             "intercept adjacent pair")
        edge = indices[0]
        need(((v0 >> edge) & 3) == 2 and ((v1 >> edge) & 3) == 1,
             "intercept direction 01 to 10")
        delta = carry(v0, q) - carry(v1, q)
        need(delta == (1 << edge) * 3 ** (v0 >> (edge + 2)).bit_count(),
             "primitive smooth difference")
        need(gcd(abs(clock), carry(v0, q)) == 1, "primitive balanced freeness")
        counts["primitive_mechanical_shapes"] += 1
need(unit_shapes == [[1, 2], [2, 3]], "finite unit-clock boundary")


hostiles = []
for base_text, edges in [("0011", (3, 2)), ("11101111000", (3, 14))]:
    base = sum(int(t) << i for i, t in enumerate(base_text))
    q0 = len(base_text)
    source = base | (base << q0)
    q = 2 * q0
    output = source
    for edge in edges:
        output = swap(output, q, edge)
        need(output is not None, "hostile actual edit")
    b = carry(output, q)
    h = (1 << q0) + 3 ** base.bit_count()
    d = (1 << q) - 3 ** output.bit_count()
    target = necklace(output, q)
    distance = 0 if target == necklace(source, q) else (
        1 if target in {necklace(w, q) for _, w in all_swaps(source, q)} else 2)
    if base_text == "0011":
        need((h, b, d, gcd(h, b)) == (25, 680, 175, 5), "five-adic hostile")
    else:
        need(b == -17 * d and distance == 2, "genuine integer distance two hostile")
        n = b // d
        for i in range(q):
            bit = (output >> i) & 1
            need(n % 2 == bit, "negative orbit parity")
            n = (3 * n + 1) // 2 if bit else n // 2
        need(n == -17, "negative orbit closes")
    hostiles.append({"base": base_text, "edges": list(edges), "H": h,
                     "B": b, "D": d, "H_gcd": gcd(h, b), "distance": distance})

print(json.dumps({"status": "PASS_EXACT", "checks": checks,
                  "universe": "q0=2..7, d=2..7; q0=9, d=2; all actual ordered paths; primitive shapes q<=200",
                  "counts": dict(sorted(counts.items())),
                  "unit_shapes_p_q": unit_shapes, "hostiles": hostiles},
                 sort_keys=True, indent=2))

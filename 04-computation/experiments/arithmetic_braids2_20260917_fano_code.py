"""Exact Fano incidence / Hamming code / E8 bridge; no external dependencies.

python 04-computation/experiments/arithmetic_braids2_20260917_fano_code.py
Truth-bearing checks use require(), so python -O performs the same audit.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path
import argparse
import json


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def dot(a, b):
    return (a & b).bit_count() & 1


def xor_sum(values):
    result = 0
    for value in values:
        result ^= value
    return result


def incidence_image(mask, lines):
    return sum((sum((mask >> (x - 1)) & 1 for x in line) % 2) << i
               for i, line in enumerate(lines))


def syndrome(mask):
    return xor_sum(a for a in range(1, 8) if (mask >> (a - 1)) & 1)


def cd_basis(i, j, dimension=8):
    # Baez convention (a,b)(c,d)=(ac-db*, a*d+cb).
    if dimension == 1:
        return 1, 0
    h = dimension // 2
    if i < h and j < h:
        return cd_basis(i, j, h)
    if i < h <= j:
        sign, k = cd_basis(i, j - h, h)
        return sign * (1 if i == 0 else -1), k + h
    if j < h <= i:
        sign, k = cd_basis(j, i - h, h)
        return sign, k + h
    sign, k = cd_basis(j - h, i - h, h)
    return -sign * (1 if i == h else -1), k


def tournament(mask, lines):
    adj = [0] * 7
    for i, (a, b, c) in enumerate(lines):
        cycle = (a, b, c) if (mask >> i) & 1 else (a, c, b)
        for x, y in zip(cycle, cycle[1:] + cycle[:1]):
            adj[x - 1] |= 1 << (y - 1)
    return adj


def hamiltonian_dp(adj):
    dp = [[0] * 7 for _ in range(128)]
    for v in range(7):
        dp[1 << v][v] = 1
    for mask in range(1, 128):
        for v in range(7):
            if not dp[mask][v]:
                continue
            choices = adj[v] & ~mask
            while choices:
                bit = choices & -choices
                choices ^= bit
                dp[mask | bit][bit.bit_length() - 1] += dp[mask][v]
    return sum(dp[-1])


def hamiltonian_brute(adj):
    return sum(all(adj[a] & (1 << b) for a, b in zip(order, order[1:]))
               for order in permutations(range(7)))


def determinant(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    result = Fraction(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            result *= -1
        value = a[i][i]
        result *= value
        for j in range(i + 1, len(a)):
            scale = a[j][i] / value
            a[j] = [x - scale * y for x, y in zip(a[j], a[i])]
    return result


def rref_binary(masks, width):
    rows = list(masks)
    pivots = []
    rank = 0
    for col in range(width):
        index = next((i for i in range(rank, len(rows)) if rows[i] & (1 << col)), None)
        if index is None:
            continue
        rows[rank], rows[index] = rows[index], rows[rank]
        for i in range(len(rows)):
            if i != rank and rows[i] & (1 << col):
                rows[i] ^= rows[rank]
        pivots.append(col)
        rank += 1
    return rows[:rank], pivots


def parity_mask(vector):
    return sum((value & 1) << i for i, value in enumerate(vector))


def isometry_doubled(vector):
    # Input z represents z/sqrt(2); output represents 2*Q(z/sqrt(2)).
    result = []
    for i in range(0, 8, 2):
        result.extend((vector[i] + vector[i + 1], vector[i] - vector[i + 1]))
    return tuple(result)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or Path(__file__).with_suffix(".json")
    # The line at position a-1 has normal a in the standard F2 dot product.
    lines = [tuple(x for x in range(1, 8) if not dot(a, x)) for a in range(1, 8)]
    require(all(len(line) == 3 and xor_sum(line) == 0 for line in lines), "Fano lines")
    image_list = [incidence_image(g, lines) for g in range(128)]
    code = set(image_list)
    kernel = {g for g in range(128) if image_list[g] == 0}
    characters = {sum(dot(u, x) << (x - 1) for x in range(1, 8)) for u in range(8)}
    affine_words = {sum((s ^ dot(v, a)) << (a - 1) for a in range(1, 8))
                    for s in range(2) for v in range(8)}
    require(code == affine_words == {c for c in range(128) if syndrome(c) == 0}, "three code definitions")
    require(kernel == characters, "kernel is character/simplex code")
    require(Counter(image_list) == {c: 8 for c in code}, "every gauge fiber has size eight")
    require(Counter(c.bit_count() for c in code) == {0: 1, 3: 7, 4: 7, 7: 1}, "Hamming enumerator")
    require(Counter(c.bit_count() for c in kernel) == {0: 1, 4: 7}, "simplex enumerator")
    for g in range(128):
        s = g.bit_count() % 2
        v = xor_sum(x for x in range(1, 8) if g & (1 << (x - 1)))
        formula = sum((s ^ dot(v, a)) << (a - 1) for a in range(1, 8))
        require(formula == image_list[g], "incidence formula")
    base = sum((cd_basis(line[0], line[1])[0] == 1) << i for i, line in enumerate(lines))
    gauges = {base ^ c for c in code}
    counts = {}
    cosets = {s: [] for s in range(8)}
    for m in range(128):
        adj = tournament(m, lines)
        require(all(row.bit_count() == 3 for row in adj), "all line orientations regular")
        counts[m] = hamiltonian_dp(adj)
        require(counts[m] == hamiltonian_brute(adj), ("independent Hamiltonian check", m))
        require(counts[m] == (189 if m in gauges else 171), "Hamiltonian locus")
        s = syndrome(m ^ base)
        cosets[s].append(m)
        corrected = m if s == 0 else m ^ (1 << (s - 1))
        nearest = [g for g in gauges if (m ^ g).bit_count() <= 1]
        require(nearest == [corrected], ("perfect one-line correction", m))
    # Check actual gauge multiplication signs, not just incidence cardinality.
    for g in range(128):
        adj = tournament(base ^ image_list[g], lines)
        for a, b in permutations(range(1, 8), 2):
            sign, index = cd_basis(a, b)
            require(index == a ^ b, "XOR product sidecar")
            sign *= (-1) ** (((g >> (a - 1)) ^ (g >> (b - 1)) ^ (g >> (index - 1))) & 1)
            require(bool(adj[a - 1] & (1 << (b - 1))) == (sign == 1), "actual sign gauge")
    extended = {(c << 1) | (c.bit_count() % 2) for c in code}
    affine_full = {sum((s ^ dot(v, a)) << a for a in range(8))
                   for s in range(2) for v in range(8)}
    require(extended == affine_full, "parity extension is full RM(1,3)")
    require(Counter(c.bit_count() for c in extended) == {0: 1, 4: 14, 8: 1}, "extended enumerator")
    require({x for x in range(256) if all(not dot(x, c) for c in extended)} == extended, "self-duality exhaustive")
    require(all((a ^ b) in extended for a in extended for b in extended), "extended linearity")
    basis, pivots = rref_binary(sorted(extended), 8)
    integer_basis = [[(c >> i) & 1 for i in range(8)] for c in basis]
    integer_basis += [[2 * (i == j) for i in range(8)] for j in range(8) if j not in pivots]
    gram = [[sum(a * b for a, b in zip(u, v)) // 2 for v in integer_basis] for u in integer_basis]
    require(abs(determinant(integer_basis)) == 16, "integer preimage index")
    require(all(sum(a*b for a,b in zip(u,v)) % 2 == 0 for u in integer_basis for v in integer_basis), "integral Gram")
    require(all(gram[i][i] % 2 == 0 for i in range(8)), "even lattice")
    require(determinant(gram) == 1, "unimodular lattice")
    roots = set()
    for i in range(8):
        for sign in (-1, 1):
            roots.add(tuple(2 * sign if j == i else 0 for j in range(8)))
    for c in extended:
        if c.bit_count() == 4:
            support = [i for i in range(8) if c & (1 << i)]
            for signs in product((-1, 1), repeat=4):
                vector = [0] * 8
                for i, sign in zip(support, signs):
                    vector[i] = sign
                roots.add(tuple(vector))
    # Independent full box covers every squared-length-four integer vector.
    brute_roots = {z for z in product(range(-2, 3), repeat=8)
                   if sum(t * t for t in z) == 4 and parity_mask(z) in extended}
    require(roots == brute_roots and len(roots) == 240, "complete root census")
    standard_doubled_roots = set()
    for i, j in combinations(range(8), 2):
        for s, t in product((-2, 2), repeat=2):
            z = [0] * 8
            z[i], z[j] = s, t
            standard_doubled_roots.add(tuple(z))
    standard_doubled_roots |= {z for z in product((-1, 1), repeat=8) if sum(z) % 4 == 0}
    require({isometry_doubled(z) for z in roots} == standard_doubled_roots, "explicit isometry onto standard E8 roots")
    for u in roots:
        for v in roots:
            inner2 = sum(a * b for a, b in zip(u, v))
            require(inner2 % 2 == 0, "root inner products integral")
            reflection = tuple(b - inner2 // 2 * a for a, b in zip(u, v))
            require(reflection in roots, "root reflections close")
    # Hostile controls: two line flips miscorrect; affine locus need not contain zero.
    two_errors = (1 << 0) ^ (1 << 1)
    error_syndrome = syndrome(two_errors)
    decoded_error = two_errors ^ (1 << (error_syndrome - 1))
    require(decoded_error != 0 and decoded_error in code and decoded_error.bit_count() == 3, "two errors do not recover original")
    require(any(c.bit_count() == 3 for c in code), "zero extension is not doubly even")
    result = {
        "status": "FINITE-EXACT, with universal algebraic proof in companion note",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "all_checks_passed": True,
        "coordinate_convention": "point x=1..7; line at bit a-1 consists of nonzero x with binary dot(a,x)=0; extended bit0 is a=0",
        "lines_in_normal_order": lines,
        "incidence_matrix": [[int(x in line) for x in range(1, 8)] for line in lines],
        "code_dimension": 4, "kernel_dimension": 3,
        "hamming_code_masks": sorted(code), "simplex_kernel_masks": sorted(kernel),
        "code_weight_enumerator": dict(Counter(c.bit_count() for c in code)),
        "kernel_weight_enumerator": dict(Counter(c.bit_count() for c in kernel)),
        "base_octonion_orientation_mask": base,
        "base_orientation_syndrome": syndrome(base),
        "octonion_gauge_orientation_masks": sorted(gauges),
        "relative_syndrome_cosets": {s: {"orientation_masks": masks, "H_histogram": dict(Counter(counts[m] for m in masks))} for s, masks in cosets.items()},
        "one_line_corrections_checked": 128,
        "hamiltonian_DP_brute_checks": 128,
        "extended_code_masks": sorted(extended),
        "extended_weight_enumerator": dict(Counter(c.bit_count() for c in extended)),
        "integer_lattice_basis_before_sqrt2_scaling": integer_basis,
        "scaled_Gram_matrix": gram, "Gram_determinant": 1,
        "integer_preimage_index": 16, "scaled_lattice_covolume": 1,
        "root_count": 240, "coordinate_roots": 16, "weight_four_roots": 224,
        "standard_E8_integer_roots": 112, "standard_E8_half_integer_roots": 128,
        "root_reflection_checks": 240 * 240,
        "root_universe": "all integer z with sum(z_i^2)=4, independently contained in [-2,2]^8; actual roots z/sqrt(2)",
        "hostile_two_error": {"line_normals": [1, 2], "syndrome": error_syndrome, "wrong_codeword_mask_after_decoding": decoded_error},
        "dimension_generalization": [{"r": r, "length": 2**r - 1, "image_dimension": r + 1,
                                       "extended_length": 2**r, "self_dual_dimension_possible": 2 * (r + 1) == 2**r}
                                      for r in range(2, 8)],
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"output": str(output), "all_checks_passed": True, "base_mask": base,
                      "code_weights": result["code_weight_enumerator"], "root_count": len(roots)}, indent=2))


if __name__ == "__main__":
    main()

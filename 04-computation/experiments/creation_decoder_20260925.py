"""Exact parity/Fano carries and all-depth proof controls.

See 05-knowledge/results/creation_decoder_20260925.md.
No floating-point arithmetic or convergence assumption is used.
"""

from itertools import combinations, product
import json
from pathlib import Path


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def step(n, sigma):
    return (3 * n + sigma) // 2 if n & 1 else n // 2


def parity(n, depth, sigma=1):
    value = 0
    for j in range(depth):
        value |= (n & 1) << j
        n = step(n, sigma)
    return value


def inverse_word(word, depth, sigma):
    # Independent affine integrality recovery from the output parity word.
    a, b, denominator = 1, 0, 1
    for j in range(depth):
        if (word >> j) & 1:
            a, b = 3 * a, 3 * b + sigma * denominator
        denominator *= 2
    return (-b * pow(a, -1, denominator)) % denominator


def anf(values, depth):
    result = list(values)
    for i in range(depth):
        for mask in range(1 << depth):
            if mask & (1 << i):
                result[mask] ^= result[mask ^ (1 << i)]
    return [mask for mask, coefficient in enumerate(result) if coefficient]


def permutation_sign(p):
    seen, cycles = set(), 0
    for start in range(len(p)):
        if start in seen:
            continue
        cycles += 1
        current = start
        while current not in seen:
            seen.add(current)
            current = p[current]
    return (-1) ** (len(p) - cycles)


def bits3(x):
    return tuple((x >> j) & 1 for j in range(3))


def polar(x, y):
    x0, x1, _ = bits3(x)
    y0, y1, _ = bits3(y)
    return 4 * ((x0 * y1) ^ (x1 * y0))


def star(x, y):
    return x ^ y ^ polar(x, y)


def main():
    degree_rows = []
    inverse_checks = 0
    for depth in range(1, 13):
        size = 1 << depth
        plus = [parity(n, depth, 1) for n in range(size)]
        minus = [parity(n, depth, -1) for n in range(size)]
        degrees, terms, signs = [], [], []
        for sigma, values in ((1, plus), (-1, minus)):
            check(sorted(values) == list(range(size)), "not a permutation")
            for n, value in enumerate(values):
                check(inverse_word(value, depth, sigma) == n, "affine inverse")
                inverse_checks += 1
            top = [(value >> (depth - 1)) & 1 for value in values]
            support = anf(top, depth)
            degree = max(mask.bit_count() for mask in support)
            if depth >= 3:
                full_lower = (1 << (depth - 1)) - 1
                check((full_lower in support) == (sigma == 1), "top monomial")
                check(degree == depth - 1 if sigma == 1 else degree <= depth - 2,
                      "degree theorem")
            degrees.append(degree)
            terms.append(len(support))
            signs.append(permutation_sign(values))
        for n in range(size):
            check(minus[n] == plus[-n % size], "signed conjugacy")
        if depth >= 3:
            check(signs == [-1, 1], "permutation signs")
        degree_rows.append({"depth": depth, "plus_minus_degrees": degrees,
                            "plus_minus_terms": terms, "signs": signs})

    q = [parity(x, 3) for x in range(8)]
    qm = [parity(x, 3, -1) for x in range(8)]
    check(q == [0, 5, 2, 3, 4, 1, 6, 7], "plus three-bit map")
    check(qm == [0, 7, 6, 1, 4, 3, 2, 5], "minus three-bit map")
    for x in range(8):
        x0, x1, x2 = bits3(x)
        check(q[x] == x0 + 2 * x1 + 4 * (x2 ^ x0 ^ (x0 * x1)), "plus ANF")
        check(qm[x] == x0 + 2 * (x1 ^ x0) + 4 * (x2 ^ x0 ^ x1), "minus ANF")
        check(q[q[x]] == x, "three-bit involution")
    for x, y in product(range(8), repeat=2):
        check(q[x] ^ q[y] ^ q[x ^ y] == polar(x, y), "polar carry")
        check(q[star(x, y)] == q[x] ^ q[y], "transported addition")
        check(qm[x ^ y] == qm[x] ^ qm[y], "minus linearity")
        check(star(x, y) == star(y, x), "commutativity")
        check(star(x, x) == 0 and star(x, 0) == x, "group laws")
    for x, y, z in product(range(8), repeat=3):
        check(star(star(x, y), z) == star(x, star(y, z)), "associativity")
    check(q[1] ^ q[2] != q[1 ^ 2], "linearity hostile")
    lines = {tuple(sorted((x, y, x ^ y))) for x, y in combinations(range(1, 8), 2)}
    transported = {tuple(sorted(q[x] for x in line)) for line in lines}
    common = sorted(lines & transported)
    check(common == [(1, 4, 5), (2, 4, 6), (3, 4, 7)], "common Fano pencil")
    check(len(lines - transported) == len(transported - lines) == 4, "four-line trade")

    section_values = []
    for k in range(1, 13):
        n = 3 ** k - 1
        section_values.append(parity(n, 20))
        for m in range(64):
            source = (1 << k) * m + (1 << k) - 1
            expected = (1 << k) - 1 + (parity(3 ** k * m + 3 ** k - 1, 16) << k)
            check(parity(source, k + 16) == expected, "all-ones residual section")
    check(len(set(section_values)) == 12, "distinct finite residual controls")

    affine_sign_controls = []
    for m in range(2, 13):
        size, half = 1 << m, 1 << (m - 1)
        flips = sum((3 * r // half) & 1 for r in range(half))
        affine = [(3 * r + 2) % size for r in range(size)]
        check(flips % 2 == 1 and permutation_sign(affine) == -1, "affine sign proof")
        affine_sign_controls.append({"bits": m, "odd_fiber_flips": flips})

    output = {
        "status": "FINITE-EXACT controls; all-depth proofs are in the companion note",
        "universe": "q plus/minus for all residues through depth12; section controls k1..12,m0..63",
        "affine_inverse_checks": inverse_checks,
        "degree_rows": degree_rows,
        "q_plus_3": q, "q_minus_3": qm,
        "fano_common": common,
        "fano_removed": sorted(lines - transported),
        "fano_added": sorted(transported - lines),
        "distinct_sections_at_tail_zero_mod_2pow20": section_values,
        "affine_sign_controls": affine_sign_controls,
    }
    target = Path(__file__).resolve().parents[2] / "05-knowledge/results/creation_decoder_20260925.out"
    target.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: {inverse_checks} independent affine inverses; parity/Fano/section controls.")


if __name__ == "__main__":
    main()

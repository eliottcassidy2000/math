"""Squarefree Collatz prefixes and a labelled Fano/octonion tournament audit.

python 04-computation/experiments/arithmetic_braids2_20260917_squarefree_symmetry.py
No nonstandard dependencies; all truth-bearing checks survive -O.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from math import isqrt, pi, prod
from pathlib import Path
import argparse
import json


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def prime_sieve(n):
    flags = bytearray(b"\1") * (n + 1)
    flags[:2] = b"\0\0"
    for p in range(2, isqrt(n) + 1):
        if flags[p]:
            flags[p * p::p] = b"\0" * ((n - p * p) // p + 1)
    return [p for p in range(2, n + 1) if flags[p]]


def squarefree_sieve(n):
    flags = bytearray(b"\1") * (n + 1)
    flags[0] = 0
    for p in prime_sieve(isqrt(n)):
        square = p * p
        flags[square::square] = b"\0" * (n // square)
    return flags


def trial_squarefree(n, primes):
    for p in primes:
        if p * p > n:
            break
        if n % (p * p) == 0:
            return False
    return True


def odd_step(n):
    raw = 3 * n + 1
    k = 0
    while raw % 2 == 0:
        raw //= 2
        k += 1
    return raw, k


def roots_growth(length, p):
    modulus = p * p
    if p == 2:
        return set()
    if p == 3:
        return {pow(2**(length + 1), -1, modulus)}
    root = pow(pow(2, length + 1, modulus), -1, modulus)
    ratio = 2 * pow(3, -1, modulus) % modulus
    roots = set()
    for _ in range(length + 1):
        roots.add(root)
        root = root * ratio % modulus
    return roots


def multiplicative_order(a, modulus):
    x, order = a % modulus, 1
    while x != 1:
        x = x * a % modulus
        order += 1
    return order


def growth_census(length, bound):
    coefficients = [3**j * 2**(length + 1 - j) for j in range(length + 1)]
    primes = prime_sieve(isqrt(max(coefficients) * bound - 1))
    valid = bytearray(b"\1") * (bound + 1)
    valid[0] = 0
    for p in primes:
        modulus = p * p
        for root in roots_growth(length, p):
            require(root != 0, ("nonzero local root", length, p))
            for q in range(root, bound + 1, modulus):
                valid[q] = 0
    counts, first = Counter(), {}
    for q in range(1, bound + 1):
        values = [a * q - 1 for a in coefficients]
        if q <= 300:
            direct = all(trial_squarefree(n, primes) for n in values)
            require(direct == bool(valid[q]), ("independent squarefree prefix", length, q))
        if not valid[q]:
            continue
        row = values[0] % 6
        counts[row] += 1
        if row not in first:
            for a, b in zip(values, values[1:]):
                require(odd_step(a) == (b, 1), "exponent-one prefix")
            first[row] = {"q": q, "nodes": values}
    require(set(first) == {1, 3, 5}, "all rows occur with squarefree growing prefixes")
    return {"q_bound": bound, "counts_by_source_row": dict(counts), "first_by_source_row": first,
            "total": sum(counts.values()), "all_intermediate_nodes_squarefree": True}


def cd_basis(i, j, dimension=8):
    """Baez convention (a,b)(c,d)=(ac-d*b*,a* d+c b)."""
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


def hamiltonian_count(adj):
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v]:
                for w in range(n):
                    if not (mask >> w & 1) and adj[v] >> w & 1:
                        dp[mask | 1 << w][w] += dp[mask][v]
    return sum(dp[-1])


def hamiltonian_brute(adj):
    return sum(all(adj[a] >> b & 1 for a, b in zip(path, path[1:]))
               for path in permutations(range(len(adj))))


def relabel(adj, mapping):
    out = [0] * len(adj)
    for a in range(len(adj)):
        for b in range(len(adj)):
            if adj[a] >> b & 1:
                out[mapping[a]] |= 1 << mapping[b]
    return tuple(out)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or Path(__file__).with_suffix(".json")
    bound = 600_000
    sf = squarefree_sieve(bound)
    row_counts = {r: {"total": 0, "squarefree": 0} for r in (1, 3, 5)}
    words = ((1,), (2,), (1, 1), (1, 1, 1), (2, 3), (3, 2), (4, 1, 2))
    word_counts = {str(word): {r: {"total": 0, "squarefree": 0} for r in (1, 3, 5)} for word in words}
    for n in range(1, bound + 1, 2):
        row = n % 6
        row_counts[row]["total"] += 1
        row_counts[row]["squarefree"] += sf[n]
        actual, m = [], n
        for _ in range(3):
            m, k = odd_step(m)
            actual.append(k)
        for word in words:
            if tuple(actual[:len(word)]) == word:
                word_counts[str(word)][row]["total"] += 1
                word_counts[str(word)][row]["squarefree"] += sf[n]
    small_primes = prime_sieve(100)
    for n in range(1, 10_001):
        require(bool(sf[n]) == trial_squarefree(n, small_primes), ("squarefree direct", n))

    growth = {str(length): growth_census(length, 10000) for length in (1, 2, 3, 6, 10)}
    z = 101
    truncated = {}
    primes = prime_sieve(z)
    for length in (1, 2, 3, 6, 10, 20, 30):
        density = Fraction(1)
        records = []
        for p in primes:
            roots = roots_growth(length, p)
            if p >= 5:
                order = multiplicative_order(3 * pow(2, -1, p * p), p * p)
                require(len(roots) == min(length + 1, order), "order formula")
                if p <= 19:
                    direct = {q for q in range(p * p)
                              if any((3**j * 2**(length + 1 - j) * q - 1) % (p * p) == 0
                                     for j in range(length + 1))}
                    require(roots == direct, "direct root enumeration")
            density *= 1 - Fraction(len(roots), p * p)
            records.append({"p": p, "excluded_roots": len(roots)})
        lower = density * (1 - Fraction(length + 1, z))
        require(lower > 0, "strict positive lower certificate")
        truncated[str(length)] = {"cutoff": z, "local_counts": records,
                                  "upper_bound": str(density), "lower_bound": str(lower),
                                  "upper_decimal_display": float(density), "lower_decimal_display": float(lower)}

    lines = sorted({tuple(sorted((a, b, a ^ b))) for a, b in combinations(range(1, 8), 2)})
    require(len(lines) == 7, "seven Fano lines")
    pairs = [pair for line in lines for pair in combinations(line, 2)]
    require(len(pairs) == len(set(pairs)) == 21, "Fano pair partition")
    table = {(a, b): cd_basis(a, b) for a in range(8) for b in range(8)}
    for a in range(1, 8):
        require(table[a, a] == (-1, 0), "imaginary squares")
        for b in range(1, 8):
            require(table[a, b][1] == a ^ b, "xor product support")
            if a != b:
                require(table[a, b][0] == -table[b, a][0], "anticommutation")
    gauges, gauge_hist = set(), Counter()
    for mask in range(128):
        adj = [0] * 7
        for a, b in combinations(range(1, 8), 2):
            sign, c = table[a, b]
            exponent = ((mask >> (a - 1)) & 1) + ((mask >> (b - 1)) & 1) + ((mask >> (c - 1)) & 1)
            sign *= (-1)**exponent
            if sign == 1:
                adj[a - 1] |= 1 << (b - 1)
            else:
                adj[b - 1] |= 1 << (a - 1)
        require(all(row.bit_count() == 3 for row in adj), "regular Fano tournament")
        gauges.add(tuple(adj))
        gauge_hist[hamiltonian_count(adj)] += 1
    all_hist, high = Counter(), set()
    for mask in range(128):
        adj = [0] * 7
        for bit, (a, b, c) in enumerate(lines):
            triple = (a, b, c) if mask >> bit & 1 else (a, c, b)
            for u, v in zip(triple, triple[1:] + triple[:1]):
                adj[u - 1] |= 1 << (v - 1)
        count = hamiltonian_count(adj)
        require(count == hamiltonian_brute(adj), "independent Hamiltonian count")
        all_hist[count] += 1
        if count == 189:
            high.add(tuple(adj))
    require(gauge_hist == {189: 128}, "octonion sign gauges")
    require(all_hist == {171: 112, 189: 16}, "all line orientations")
    require(high == gauges and len(gauges) == 16, "exact gauge locus")
    reference = min(gauges)
    linear_maps = []
    for a, b, c in permutations(range(1, 8), 3):
        if c in (a, b, a ^ b):
            continue
        linear_maps.append(tuple(((a if v & 1 else 0) ^ (b if v & 2 else 0)
                                 ^ (c if v & 4 else 0)) - 1 for v in range(1, 8)))
    first_orbit = {relabel(reference, mapping) for mapping in linear_maps}
    converse = tuple(127 ^ (1 << a) ^ reference[a] for a in range(7))
    second_orbit = {relabel(converse, mapping) for mapping in linear_maps}
    full_orbit = {relabel(reference, mapping) for mapping in permutations(range(7))}
    require(len(linear_maps) == 168, "GL3F2 order")
    require(len(first_orbit) == len(second_orbit) == 8, "two labelled Fano chiral orbits")
    require(not first_orbit & second_orbit and first_orbit | second_orbit == gauges, "chirality partition")
    require(gauges <= full_orbit and len(full_orbit) == 240, "unlabelled chirality collapse")
    odd_cycles, odd_cycle_counts = [], {}
    for length in (3, 5, 7):
        cycles = [cycle for cycle in permutations(range(7), length)
                  if cycle[0] == min(cycle)
                  and all(reference[a] >> b & 1 for a, b in zip(cycle, cycle[1:] + cycle[:1]))]
        odd_cycle_counts[length] = len(cycles)
        odd_cycles.extend(cycles)
    cycle_masks = [sum(1 << v for v in cycle) for cycle in odd_cycles]
    disjoint_pairs = sum(not (a & b) for a, b in combinations(cycle_masks, 2))
    require(odd_cycle_counts == {3: 14, 5: 42, 7: 24} and disjoint_pairs == 7, "full odd-cycle sidecar")
    require(1 + 2 * len(odd_cycles) + 4 * disjoint_pairs == 189, "independent OCF count")
    associator_witness = None
    for a, b, c in permutations(range(1, 8), 3):
        s1, i = table[a, b]
        s2, left = table[i, c]
        s3, j = table[b, c]
        s4, right = table[a, j]
        if (s1 * s2, left) != (s3 * s4, right):
            associator_witness = {"a_b_c": [a, b, c], "left": [s1 * s2, left], "right": [s3 * s4, right]}
            break
    require(associator_witness is not None, "nonassociative hostile")

    result = {
        "status": "FINITE-EXACT controls; universal density proofs in companion markdown",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "squarefree_census_bound": bound, "squarefree_all": sum(sf),
        "row_counts": row_counts,
        "asymptotic_relative_density_row_1_5": "9/pi^2",
        "asymptotic_relative_density_row_3": "6/pi^2",
        "decimal_display_only": {"9/pi^2": 9 / pi**2, "6/pi^2": 6 / pi**2},
        "selected_word_counts": word_counts,
        "simultaneously_squarefree_growth": growth,
        "certified_density_intervals": truncated,
        "Fano_lines_xor_labels": lines, "Fano_point_count": 7, "Fano_pair_count": 21,
        "octonion_gauge_masks": 128, "octonion_distinct_gauge_tournaments": len(gauges),
        "octonion_gauge_H_histogram": dict(gauge_hist), "all_Fano_line_orientation_H_histogram": dict(all_hist),
        "octonion_gauge_locus_equals_H189_locus": True,
        "octonion_GL3F2_orbit_sizes": [len(first_orbit), len(second_orbit)],
        "octonion_GL3F2_orbits_exchanged_by_converse": True,
        "octonion_S7_orbit_size": len(full_orbit),
        "octonion_reference_adjacency_bitmasks": reference,
        "octonion_reference_odd_cycle_counts": odd_cycle_counts,
        "octonion_reference_disjoint_odd_cycle_pairs": disjoint_pairs,
        "octonion_multiplication_table_sign_index": {f"{a},{b}": value for (a, b), value in table.items()},
        "octonion_associator_hostile": associator_witness,
        "independent_controls": {"squarefree_integers": [1, 10000], "growth_q_each_length": [1, 300],
                                 "local_root_enumeration_primes_through": 19, "Hamiltonian_bruteforce_tournaments": 128},
        "all_checks_passed": True,
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"output": str(output), "all_checks_passed": True, "growth": growth,
                      "octonion_gauges": dict(gauge_hist), "all_orientations": dict(all_hist)}, indent=2))


if __name__ == "__main__":
    main()

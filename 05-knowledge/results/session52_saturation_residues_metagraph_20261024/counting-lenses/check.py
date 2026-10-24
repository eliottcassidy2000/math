#!/usr/bin/env python3
"""Independent exact comparison of class counts, signed traces, and native wiring.
Only Python's standard library is required; all checks survive python -O.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, gcd, lcm
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0

def require(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)

def parts(n, lo=1):
    if n == 0:
        yield ()
    for j in range(lo, n + 1):
        for tail in parts(n - j, j):
            yield (j,) + tail

def z_type(lam):
    ans = 1
    for a, m in Counter(lam).items():
        ans *= a ** m * factorial(m)
    return ans

def representative(lam):
    perm = []
    start = 0
    for a in lam:
        perm.extend(list(range(start + 1, start + a)) + [start])
        start += a
    return tuple(perm)

def signed_pair_cycles(perm):
    """Literal action on ordered representatives i<j; no cycle-type formula."""
    pairs = list(combinations(range(len(perm)), 2))
    unseen = set(pairs)
    result = []
    while unseen:
        start = min(unseen)
        cur = start
        length = 0
        sign = 1
        while cur in unseen:
            unseen.remove(cur)
            a, b = perm[cur[0]], perm[cur[1]]
            sign *= 1 if a < b else -1
            cur = (min(a, b), max(a, b))
            length += 1
        require(cur == start, 'pair action failed to close')
        result.append((length, sign))
    return sorted(result)

def closed_pair_cycles(lam):
    """Independent closed signed-cycle inventory from vertex cycle lengths."""
    result = []
    for a in lam:
        if a % 2:
            result += [(a, 1)] * ((a - 1) // 2)
        else:
            result += [(a, 1)] * (a // 2 - 1)
            result.append((a // 2, -1))
    for a, b in combinations(lam, 2):
        result += [(lcm(a, b), 1)] * gcd(a, b)
    return sorted(result)

def trace_poly(cycles):
    out = [1]
    for length, sign in cycles:
        nxt = out + [0] * length
        for i, x in enumerate(out):
            nxt[i + length] += sign * x
        out = nxt
    return out

def fixed_from_cycles(cycles, anti=False):
    for length, sign in cycles:
        if sign * ((-1) ** length if anti else 1) != 1:
            return 0
    return 1 << len(cycles)

def sc_allowed(lam):
    return lam.count(1) <= 1 and all(a == 1 or a % 4 == 2 for a in lam)

def partition_fixed(lam, anti=False):
    """Direct partition exponent, independent of signed pair-cycle code."""
    if anti:
        if not sc_allowed(lam):
            return 0
        evens = [a for a in lam if a != 1]
        e = sum(a // 2 for a in evens)
        e += sum(gcd(a, b) for a, b in combinations(evens, 2))
        e += lam.count(1) * len(evens)
    else:
        if any(a % 2 == 0 for a in lam):
            return 0
        e = sum((a - 1) // 2 for a in lam)
        e += sum(gcd(a, b) for a, b in combinations(lam, 2))
    return 1 << e

def partition_counts(n):
    ordinary = Fraction(0)
    converse = Fraction(0)
    for lam in parts(n):
        ordinary += Fraction(partition_fixed(lam), z_type(lam))
        converse += Fraction(partition_fixed(lam, True), z_type(lam))
    require(ordinary.denominator == converse.denominator == 1, 'nonintegral count')
    return int(ordinary), int(converse)

def signed_polynomial(n):
    d = n * (n - 1) // 2
    weighted = [0] * (d + 1)
    nf = factorial(n)
    for lam in parts(n):
        cycles = signed_pair_cycles(representative(lam))
        require(cycles == closed_pair_cycles(lam), ('closed cycles', lam))
        require(fixed_from_cycles(cycles) == partition_fixed(lam), ('ordinary filter', lam))
        require(fixed_from_cycles(cycles, True) == partition_fixed(lam, True), ('anti filter', lam))
        weight = nf // z_type(lam)
        for k, c in enumerate(trace_poly(cycles)):
            weighted[k] += weight * c
    require(all(c % nf == 0 for c in weighted), ('nonintegral coefficients', n))
    coeffs = [c // nf for c in weighted]
    require(all(c >= 0 for c in coeffs), ('negative multiplicity', n))
    return coeffs

def literal_permutation_average(n):
    """All n! permutations, not cycle-type representatives."""
    d = n * (n - 1) // 2
    total = [0] * (d + 1)
    positive = negative = 0
    for perm in permutations(range(n)):
        cycles = signed_pair_cycles(perm)
        for k, c in enumerate(trace_poly(cycles)):
            total[k] += c
        positive += fixed_from_cycles(cycles)
        negative += fixed_from_cycles(cycles, True)
    nf = factorial(n)
    require(all(v % nf == 0 for v in total), 'literal average divisibility')
    return [v // nf for v in total], (positive // nf, negative // nf)

def relabel(mask, perm, pairs, pair_index):
    out = 0
    for k, (a, b) in enumerate(pairs):
        u, v = perm[a], perm[b]
        bit = ((mask >> k) & 1) ^ (u > v)
        out |= bit << pair_index[(min(u, v), max(u, v))]
    return out

def literal_mask_counts(n):
    """No orbit-sign logic: enumerate all masks, relabel, compare bits."""
    pairs = list(combinations(range(n), 2))
    pi = {e: i for i, e in enumerate(pairs)}
    full = (1 << len(pairs)) - 1
    ordinary = converse = 0
    for p in permutations(range(n)):
        fixed = anti = 0
        for mask in range(full + 1):
            image = relabel(mask, p, pairs, pi)
            fixed += image == mask
            anti += image == (mask ^ full)
        ordinary += fixed
        converse += anti
    nf = factorial(n)
    require(ordinary % nf == converse % nf == 0, 'literal masks divisibility')
    return ordinary // nf, converse // nf

def rows_from_tiling(n, mask):
    rows = [0] * n
    free = list((i, j) for i in range(n) for j in range(i + 2, n))
    bits = {e: (mask >> k) & 1 for k, e in enumerate(free)}
    for i, j in combinations(range(n), 2):
        if j == i + 1 or bits[(i, j)]:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return tuple(rows)

def directed(rows, u, v):
    return bool(rows[u] & (1 << v))

def canonical_rows(rows, merged=False):
    n = len(rows)
    vals = []
    for p in permutations(range(n)):
        vals.append(tuple(int(directed(rows, p[i], p[j])) for i, j in combinations(range(n), 2)))
    if merged:
        vals += [tuple(1 - b for b in x) for x in vals]
    return min(vals)

def operator_witnesses():
    # Native all-free-bit complement at n=3: both tilings are reflection-fixed.
    reps = [rows_from_tiling(3, 1), rows_from_tiling(3, 0)]  # transitive, cyclic
    keys = [canonical_rows(r) for r in reps]
    native = [[0, 0], [0, 0]]
    for mask in (0, 1):
        a = keys.index(canonical_rows(rows_from_tiling(3, mask)))
        b = keys.index(canonical_rows(rows_from_tiling(3, mask ^ 1)))
        native[a][b] += 1
    arc = [[0, 0], [0, 0]]
    for i, rows in enumerate(reps):
        for u, v in combinations(range(3), 2):
            new = list(rows)
            new[u] ^= 1 << v
            new[v] ^= 1 << u
            arc[i][keys.index(canonical_rows(new))] += 1
    require(native == [[0, 1], [1, 0]], 'native witness')
    require(arc == [[2, 1], [3, 0]], 'arc witness')
    # Characteristic polynomial of 2x2 is x^2 - tr*x + determinant.
    def charpoly(a):
        return [1, -a[0][0] - a[1][1], a[0][0] * a[1][1] - a[0][1] * a[1][0]]
    require(charpoly(native) == [1, 0, -1], 'native eigenvalues')
    require(charpoly(arc) == [1, -2, -3], 'arc eigenvalues')
    # Minimal n=4 non-descent witness, and exhaustive absence at smaller n.
    minimal = None
    for n in range(1, 5):
        m = (n - 1) * (n - 2) // 2
        by_class = {}
        for t in range(1 << m):
            source = canonical_rows(rows_from_tiling(n, t), True)
            target = canonical_rows(rows_from_tiling(n, t ^ ((1 << m) - 1)), True)
            by_class.setdefault(source, set()).add(target)
        if any(len(x) > 1 for x in by_class.values()):
            minimal = n
            break
    require(minimal == 4, 'minimal complement non-descent')
    r0, r1 = rows_from_tiling(4, 0), rows_from_tiling(4, 1)
    require(r0 == (2, 4, 9, 3) and r1 == (6, 4, 8, 3), 'literal n4 rows')
    p = (1, 2, 3, 0)
    require(all(directed(r0, u, v) == directed(r1, p[u], p[v]) for u, v in permutations(range(4), 2)), 'n4 isomorphism')
    targets = [rows_from_tiling(4, t) for t in (7, 6)]
    scores = [sorted(row.bit_count() for row in rows) for rows in targets]
    require(scores == [[0, 1, 2, 3], [0, 2, 2, 2]], 'n4 distinct merged targets')
    require(canonical_rows(targets[0], True) != canonical_rows(targets[1], True), 'n4 merged distinction')
    return {'n3_order': ['transitive', 'cyclic'], 'native': native, 'arc_flip': arc,
            'native_characteristic_polynomial_descending': charpoly(native),
            'arc_characteristic_polynomial_descending': charpoly(arc),
            'n4_tilings': [0, 1], 'n4_rows': [r0, r1], 'n4_complement_scores': scores,
            'minimal_non_descent_order': minimal}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--receipt', default=str(HERE / 'receipt.json'))
    args = parser.parse_args()
    baseline = json.loads((HERE / 'inputs' / 'counts-0-20.json').read_text())
    table = []
    for n in range(21):
        coeffs = signed_polynomial(n)
        a, sc = partition_counts(n)
        require(sum(coeffs) == a, ('P(1)', n))
        require(sum((-1) ** k * c for k, c in enumerate(coeffs)) == sc, ('P(-1)', n))
        require((a, sc) == (baseline['ordinary'][n], baseline['self_converse'][n]), ('frozen input', n))
        if n <= 8:
            literal_coeffs, literal_counts = literal_permutation_average(n)
            require(literal_coeffs == coeffs and literal_counts == (a, sc), ('all permutations', n))
        if n <= 5:
            require(literal_mask_counts(n) == (a, sc), ('all labelled masks', n))
        table.append({'n': n, 'ordinary': a, 'self_converse': sc, 'merged': (a + sc) // 2,
                      'unpaired_class_count': sc, 'paired_merged_nodes': (a - sc) // 2,
                      'signed_polynomial_coefficients': coeffs})
    require(table[3]['signed_polynomial_coefficients'] == [1, 0, 1, 0], 'P3')
    hostiles = []
    for lam in ((1, 1), (4,), (1, 1, 2), (3,), (2,), (2, 2), (6,), (1, 6)):
        cycles = signed_pair_cycles(representative(lam))
        fixed = fixed_from_cycles(cycles, True)
        require(fixed == partition_fixed(lam, True), ('hostile', lam))
        hostiles.append({'cycle_type': lam, 'signed_pair_cycles': cycles, 'anti_fixed_labelled': fixed})
    witness = operator_witnesses()
    out = {'status': 'FINITE-EXACT independent replay of inherited formulas; no new counting algorithm',
           'scope': {'partition_coefficients_max_n': 20, 'all_permutations_max_n': 8,
                     'all_labelled_masks_and_permutations_max_n': 5},
           'checks': CHECKS, 'table': table, 'cycle_type_controls': hostiles,
           'operator_witnesses': witness,
           'inputs_sha256': hashlib.sha256((HERE / 'inputs' / 'counts-0-20.json').read_bytes()).hexdigest()}
    Path(args.receipt).write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'checks': CHECKS, 'status': 'PASS', 'A20': table[20]['ordinary'],
                      'SC20': table[20]['self_converse'], 'receipt': str(args.receipt)}, sort_keys=True))

if __name__ == '__main__':
    main()

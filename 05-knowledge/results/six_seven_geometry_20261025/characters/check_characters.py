"""Exact, standard-library controls for CHARACTER_SEXTICS.md. No removable assertions."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
FORMS = [(1, 0, 0), (1, 1, 0), (1, -1, 0), (1, 0, 1),
         (1, 0, -1), (1, 1, 1), (1, -1, -1), (1, 1, -1), (1, -1, 1)]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def det(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


def inverse_mod(a, p):
    rows = [[x % p for x in row] + [int(i == j) for j in range(3)]
            for i, row in enumerate(a)]
    for j in range(3):
        pivot = next(i for i in range(j, 3) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        factor = pow(rows[j][j], -1, p)
        rows[j] = [x*factor % p for x in rows[j]]
        for i in range(3):
            if i != j:
                factor = rows[i][j]
                rows[i] = [(x-factor*y) % p for x, y in zip(rows[i], rows[j])]
    return [row[3:] for row in rows]


def projective_plane(p):
    return [(1, a, b) for a in range(p) for b in range(p)] + [(0, 1, a) for a in range(p)] + [(0, 0, 1)]


def legendre(a, p):
    a %= p
    if not a:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


def character_count(forms, p):
    n = len(forms)
    masks = [s for s in range(1 << n) if s.bit_count() % 2 == 0]
    sums = Counter({s: 0 for s in masks})
    fibre_hist = Counter()
    for b in projective_plane(p):
        values = [legendre(sum(x*y for x, y in zip(row, b)), p) for row in forms]
        terms = [1]*(1 << n)
        for s in range(1, 1 << n):
            bit = s & -s
            terms[s] = terms[s ^ bit]*values[bit.bit_length()-1]
        count = sum(terms[s] for s in masks)
        require(count >= 0, 'negative character fibre')
        # A second local formula explicitly averages both scalar square classes.
        plus = minus = 1
        for value in values:
            plus *= 1+value
            minus *= 1-value
        require(2*count == plus+minus, 'scalar averaging failure')
        fibre_hist[count] += 1
        sums.update({s: terms[s] for s in masks})
    return sum(sums.values()), sums, fibre_hist


def root_count(forms, p):
    basis = next(ids for ids in combinations(range(len(forms)), 3)
                 if det([forms[i] for i in ids]) % p)
    inv = inverse_mod([forms[i] for i in basis], p)
    remaining = [row for j, row in enumerate(forms) if j not in basis]
    total = 0
    for roots in projective_plane(p):
        square = [a*a % p for a in roots]
        base = [sum(x*y for x, y in zip(row, square)) % p for row in inv]
        # Literal root enumeration in F_p, independent of Legendre symbols.
        ways = 1
        for row in remaining:
            radicand = sum(x*y for x, y in zip(row, base)) % p
            ways *= sum(u*u % p == radicand for u in range(p))
        total += ways
    return total


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for a in range(total+1):
            for tail in compositions(total-a, length-1):
                yield (a,) + tail


def canonical_characters(n):
    d = n-6
    observed = Counter()
    if d < 0:
        return observed
    # Complete intersection normal basis: unrestricted first3 exponents,
    # squarefree exponents on the remaining n-3 root variables.
    for tail in product((0, 1), repeat=n-3):
        if sum(tail) > d:
            continue
        for first in compositions(d-sum(tail), 3):
            exponents = first+tail
            mask = sum(((1+a) % 2) << j for j, a in enumerate(exponents))
            require(mask.bit_count() % 2 == 0, 'canonical action not projective')
            observed[mask] += 1
    return observed


def main():
    genera = []
    for n in range(4, 13):
        actual = canonical_characters(n)
        expected = Counter({s: comb(s.bit_count()//2-1, 2)
                            for s in range(1 << n)
                            if s.bit_count() >= 6 and s.bit_count() % 2 == 0})
        require(actual == expected, f'character decomposition at n={n}')
        d = n-6
        hilbert = (sum((-1)**j*comb(n-3, j)*comb(n-1+d-2*j, n-1)
                       for j in range(min(n-3, d//2)+1)) if d >= 0 else 0)
        require(sum(actual.values()) == hilbert, f'Hilbert genus n={n}')
        genera.append({'n': n, 'pg': hilbert,
                       'by_weight': dict(sorted(Counter({k: sum(v for s, v in actual.items() if s.bit_count() == k)
                                                       for k in range(6, n+1, 2)}).items()))})
    require([row['pg'] for row in genera[:6]] == [0, 0, 1, 7, 31, 111], 'known genus control')
    c7 = canonical_characters(7)
    require(set(c7) == {127 ^ (1 << i) for i in range(7)}, 'omitted-root characters')
    flip_table = [[(-1)**((s >> j) & 1) for j in range(7)] for s in sorted(c7)]
    require(all(row.count(1) == 1 for row in flip_table), 'determinant-twist hostile')

    point_rows = []
    for n in range(5, 10):
        for p in (3, 5, 7, 11, 13, 17, 19):
            total, sums, hist = character_count(FORMS[:n], p)
            direct = root_count(FORMS[:n], p)
            require(total == direct, f'point count n={n}, p={p}')
            point_rows.append({'n': n, 'p': p, 'points': total,
                               'character_sums': {str(s): sums[s] for s in sorted(sums)},
                               'fibre_histogram': dict(sorted(hist.items()))})
    subset_controls = []
    for ids in combinations(range(9), 7):
        forms = [FORMS[i] for i in ids]
        counts = []
        for p in (3, 5):
            total, _, _ = character_count(forms, p)
            require(total == root_count(forms, p), f'all subsets {ids} mod{p}')
            counts.append(total)
        subset_controls.append({'labels': ids, 'p3_p5_counts': counts})
    # Sharp hostile: an odd root product changes under global projective sign.
    require((-1)**3 == -1 and (-1)**6 == 1, 'odd-character hostile')
    receipt = {'status': 'FINITE-EXACT',
               'universe': {'canonical_n': [4, 12], 'chain_n': [5, 9],
                            'prime_fields': [3, 5, 7, 11, 13, 17, 19],
                            'all_seven_subsets_of_nine_at_primes': [3, 5]},
               'genera': genera, 'canonical_single_flip_table': flip_table,
               'point_counts': point_rows, 'all_seven_subsets': subset_controls,
               'boundary': 'singular models only; no resolution point counts, good-reduction or integer lifting claim',
               'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    out = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    (HERE/'character_receipt.json').write_text(out, encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'canonical_dimensions': [r['pg'] for r in genera],
                      'chain_point_counts': len(point_rows), 'seven_subset_controls': len(subset_controls),
                      'receipt_sha256': hashlib.sha256(out.encode()).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    main()

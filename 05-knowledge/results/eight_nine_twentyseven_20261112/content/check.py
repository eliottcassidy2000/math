"""Small exact controls for the content proofs; no large tournament census."""
from collections import Counter
from fractions import Fraction
from functools import reduce
import hashlib
from itertools import permutations
import json
from math import comb, factorial, gcd
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
BASE = 'c022b73c3cb5'
BANK = '05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/extended_census.json'
CORE = (56,117,129,22,100,76,141,59)
EXPECTED = (513,6102,32508,99432,187200,216000,141120,40320)
FPOLY = (0,0,246,0,232,0,34,0,1)


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def vp(value, prime):
    if value == 0:
        return 10**9
    result = 0
    while value % prime == 0:
        result += 1
        value //= prime
    return result


def prime_factors(n):
    out, divisor = [], 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            out.append(divisor)
            while n % divisor == 0:
                n //= divisor
        divisor += 1
    if n > 1:
        out.append(n)
    return out


def content(covers):
    h = gcd(*covers)
    g = gcd(*(factorial(k) * x for k, x in enumerate(covers)))
    invariants = []
    E = 1
    for p in prime_factors(h):
        a = vp(h, p)
        kappa = next(k for k, x in enumerate(covers, 1) if vp(x, p) == a)
        b = vp(factorial(kappa - 1), p)
        invariants.append(dict(p=p, a=a, kappa=kappa, b=b))
        E *= p**b
    return dict(raw=h, contextual=g, E=E, prime_data=invariants)


def covers_dp(rows):
    n, size = len(rows), 1 << len(rows)
    h = [0] * size
    h[0] = 1
    ends = [[0] * n for _ in range(size)]
    for v in range(n):
        ends[1 << v][v] = 1
    for mask in range(1, size):
        for v, count in enumerate(ends[mask]):
            if not count:
                continue
            todo = rows[v] & ~mask
            while todo:
                bit = todo & -todo
                todo -= bit
                ends[mask | bit][bit.bit_length() - 1] += count
        h[mask] = sum(ends[mask])
    parts = [[0] * (n + 1) for _ in range(size)]
    parts[0][0] = 1
    for mask in range(1, size):
        root, block = mask & -mask, mask
        while block:
            if block & root:
                remainder = mask ^ block
                for k in range(1, n + 1):
                    parts[mask][k] += h[block] * parts[remainder][k - 1]
            block = (block - 1) & mask
    return tuple(parts[-1][k] * factorial(k) for k in range(1, n + 1))


def literal_covers(rows):
    """Independent permutation/cut enumeration via the backward-adjacency count."""
    n = len(rows)
    descents = Counter()
    for path in permutations(range(n)):
        backwards = sum(not (rows[a] >> b & 1) for a, b in zip(path, path[1:]))
        descents[backwards] += 1
    covers = [sum(count * comb(n - 1 - b, k - 1 - b)
                  for b, count in descents.items() if b <= k - 1)
              for k in range(1, n + 1)]
    return tuple(covers), [descents[b] for b in range(n)]


def from_covers(covers):
    result = [Fraction(0)] * (len(covers) + 1)
    falling = [1]
    for k, count in enumerate(covers, 1):
        new = [0] * (len(falling) + 1)
        for i, coefficient in enumerate(falling):
            new[i] -= (k - 1) * coefficient
            new[i + 1] += coefficient
        falling = new
        for i, coefficient in enumerate(falling):
            result[i] += Fraction(count * coefficient, factorial(k))
    require(all(x.denominator == 1 for x in result), 'integer monomial coefficients')
    return tuple(int(x) for x in result)


def evaluate(poly, x):
    return sum(c * x**k for k, c in enumerate(poly))


def newton(values):
    out, values = [], list(values)
    while values:
        out.append(values[0])
        values = [b - a for a, b in zip(values, values[1:])]
    return tuple(out[1:])


def descent_vector(poly):
    n = len(poly) - 1
    out = []
    for x in range(1, n + 1):
        out.append(evaluate(poly, x) -
                   sum(count * comb(x + n - 1 - b, n) for b, count in enumerate(out)))
    return out


def is_strong(rows):
    n = len(rows)
    reach = [row | (1 << i) for i, row in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return all(x == (1 << n) - 1 for x in reach)


def main():
    covers = covers_dp(CORE)
    literal, descents = literal_covers(CORE)
    require(covers == literal == EXPECTED, 'two independent full core profiles')
    require(from_covers(covers) == FPOLY, 'core slot polynomial')
    require(is_strong(CORE), 'core is strong')
    automorphisms = sum(all((CORE[i] >> j & 1) == (CORE[p[i]] >> p[j] & 1)
                            for i in range(8) for j in range(8))
                        for p in permutations(range(8)))
    require(automorphisms == 1, 'asymmetric core')
    source = dict(rows=CORE, covers=covers, descents=descents,
                  polynomial_coefficients=FPOLY, strong=True,
                  automorphisms=automorphisms, **content(covers))
    require((source['raw'], source['contextual'], source['E']) == (9,27,3), 'sharp source')

    power_rows = []
    for r in range(1, 13):
        c = newton(evaluate(FPOLY, x)**r for x in range(8 * r + 1))
        data = content(c)
        require(data['raw'] == 9**r, 'raw ordinal power identity')
        require(data['contextual'] == 3 * 9**r, 'context ordinal power identity')
        require(data['prime_data'][0]['kappa'] == 4, 'first minimum persists')
        power_rows.append(dict(r=r, order=8*r, H=513**r, **data))

    formal = (0,18,0,4,0,4,0,1)
    formal_covers = newton(evaluate(formal, x) for x in range(8))
    formal_descents = descent_vector(formal)
    require(formal_covers == (27,270,2430,9360,17280,15120,5040), 'formal seven profile')
    require(formal_descents == [27,108,1485,1800,1485,108,27], 'formal descent vector')
    require(sum(formal_descents) == factorial(7), 'formal total permutation mass')
    no_parity = (0,3,11,-6,1)
    no_parity_covers = newton(evaluate(no_parity, x) for x in range(5))
    require(no_parity_covers == (9,0,0,24), 'nonparity hostile')
    require((content(no_parity_covers)['raw'], content(no_parity_covers)['contextual']) == (3,9),
            'parity is necessary for strict square bound')

    bank_bytes = subprocess.check_output(['git', 'show', BASE + ':' + BANK])
    records = json.loads(bank_bytes)['records']
    order_counts, ordinary_counts, pairs = Counter(), Counter(), Counter()
    fibres = {}
    gap_hits = []
    for record in records:
        if not 3 <= record['n'] <= 7:
            continue
        record_covers = covers_dp(record['adjacency_outmasks'])
        data = content(record_covers)
        H = record_covers[0]
        witness = dict(order=record['n'], rows=record['adjacency_outmasks'],
                       ordered_covers=record_covers, context_content=data['contextual'])
        if H not in fibres or data['contextual'] < fibres[H]['minimum_context_content']:
            fibres[H] = dict(H=H, parity_lower_bound=gcd(H,3),
                             minimum_context_content=data['contextual'], witness=witness)
        order_counts[record['n']] += 1
        ordinary_counts[record['n']] += record['ordinary_multiplicity']
        pairs[data['raw'], data['contextual']] += record['ordinary_multiplicity']
        if data['raw'] != data['contextual']:
            gap_hits.append(record)
        require(data['contextual'] % data['raw'] == 0, 'h divides g')
        require((data['raw'] * data['E']) % data['contextual'] == 0, 'g divides hE')
        for entry in data['prime_data']:
            require(entry['b'] <= entry['a'] - 1, 'parity first-index bound')
    require(not gap_hits, 'no amplification through order seven')
    require(sum(order_counts.values()) == 321 and sum(ordinary_counts.values()) == 530,
            'inherited exact small universe')
    for rows in ((0,), (2,0)):
        require(content(covers_dp(rows))['contextual'] == 1, 'orders one and two')

    receipt = dict(status='PASS; all-order quantifiers come from PROOF.md',
                   universe=dict(core_permutations=factorial(8), ordinal_power_exponents=list(range(1,13)),
                                 small_bank_revision=BASE, small_bank_path=BANK,
                                 small_bank_sha256=hashlib.sha256(bank_bytes).hexdigest(),
                                 converse_representatives_by_order=dict(order_counts),
                                 ordinary_classes_by_order=dict(ordinary_counts),
                                 added_transitive_orders=[1,2], filters='none within named orders'),
                   source=source, ordinal_powers=power_rows,
                   formal_seven_hostile=dict(polynomial_coefficients=formal, covers=formal_covers,
                                             descents=formal_descents, **content(formal_covers)),
                   nonparity_hostile=dict(polynomial_coefficients=no_parity, covers=no_parity_covers,
                                          **content(no_parity_covers)),
                   small_bank_content_pairs=[dict(h=h,g=g,ordinary_classes=count)
                                             for (h,g),count in sorted(pairs.items())],
                   bounded_H_fibres=[fibres[H] for H in sorted(fibres)],
                   bounded_normalization_misses=[fibres[H] for H in sorted(fibres)
                                                if fibres[H]['minimum_context_content']>gcd(H,3)],
                   source_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                                  for name in ('PROOF.md','REPAIR_LADDER.md','check.py')})
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(status='PASS', core_content=[9,27], power_controls=12,
                          small_converse_representatives=321, ordinary_classes=532,
                          bounded_H_fibres=len(fibres),
                          bounded_normalization_misses=[(H,fibres[H]['minimum_context_content'])
                              for H in sorted(fibres)
                              if fibres[H]['minimum_context_content']>gcd(H,3)])))


if __name__ == '__main__':
    main()

"""Exact controls; no assertions, so python -O performs identical checks."""
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent

def require(ok, detail):
    if not ok:
        raise RuntimeError(detail)

def oddpart(n):
    while n % 2 == 0:
        n //= 2
    return n

def sign(n):
    return 1 if n % 3 == 1 else -1

def g(n):
    if n < 2:
        return 1
    q, r = divmod(n, 2)
    return g(q) if r else g(q) + q + q % 2

def paths(rows):
    n = len(rows)
    for u in range(n):
        require(not ((rows[u] >> u) & 1), ('loop', u))
        for v in range(u):
            require(((rows[u] >> v) & 1) + ((rows[v] >> u) & 1) == 1,
                    ('not tournament', u, v))
    ends = [[0] * n for _ in range(1 << n)]
    counts = [0] * (1 << n)
    counts[0] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not ((mask >> v) & 1):
                continue
            prev = mask ^ (1 << v)
            ends[mask][v] = (sum(ends[prev][u] for u in range(n)
                                    if (rows[u] >> v) & 1) if prev else 1)
        counts[mask] = sum(ends[mask])
    return ends, counts

def state(rows, root):
    ends, counts = paths(rows)
    n = len(rows)
    full = (1 << n) - 1
    reverse = [sum(1 << u for u in range(n) if (rows[u] >> v) & 1)
               for v in range(n)]
    revends, _ = paths(reverse)
    f = sum(ends[mask][root] * counts[full ^ mask]
            for mask in range(1, full + 1) if (mask >> root) & 1)
    return [counts[full], revends[full][root], counts[full ^ (1 << root)], f]

def birth(rows, root):
    n = len(rows)
    return [row | ((1 << n) if u == root else 0)
            for u, row in enumerate(rows)] + [((1 << n) - 1) ^ (1 << root)]

inverse_checks = 0
for m in range(1, 202, 2):
    if m % 3 == 0:
        continue
    for a in range(1, 25):
        amp = (1 << a) * m
        e = sign(amp)
        n = (amp - e) // 3
        require(n > 0 and n % 2 == 1, ('inverse domain', m, a))
        require(oddpart(3 * n + e) == m, ('target', m, a))
        require((n % 3 != 0) == (amp % 9 != e % 9), ('guard', m, a))
        for j in range(7):
            nn, ee = n, e
            for _ in range(j):
                nn, ee = 2 * nn + ee, -ee
            require(nn == (1 << j) * n + ((1 << j) - (-1) ** j) // 3 * e,
                    ('iterate', m, a, j))
            require(ee == (-1) ** j * e and oddpart(3 * nn + ee) == m,
                    ('iterate target', m, a, j))
        inverse_checks += 1

first_return = []
for c, e in [(1, 1), (2, 1), (2, -1), (1, -1)]:
    cc, ee, gap = c, e, 0
    while True:
        cc, ee = (2 * cc + ee) % 3, -ee
        gap += 1
        if cc:
            break
    first_return.append([c, e, cc, ee, gap])
require(first_return == [[1,1,2,1,2],[2,1,2,-1,1],
                         [2,-1,1,-1,2],[1,-1,1,1,1]], 'return table')
require([oddpart(3 * z - 1) for z in [11, 9, 7]] == [1, 13, 5], 'eight hostile')

stencil_checks = 0
for width in range(6):
    for suffix in range(1 << width):
        for j in range(25):
            vals = [g((1 << width) * (4 * j + i) + suffix) for i in range(4)]
            require(vals[0] - vals[1] - vals[2] + vals[3] == -2,
                    ('stencil', width, suffix, j, vals))
            stencil_checks += 1

receipt_path = RESULTS / 'changing_core_returns_20261113/returns/birth_budget_receipt.json'
budget = json.loads(receipt_path.read_text())
inherited = budget['neighboring_ray']['literal_levels']
physical = []
for i, item in enumerate(inherited):
    rows, root = item['rows'], item['root']
    actual = state(rows, root)
    require(actual == item['state'], ('actual state', i, actual))
    h, t, q, f = actual
    require(t == q == 71 and h == f, ('balanced', i))
    amp = h - t
    e = sign(amp)
    n = (amp - e) // 3
    require(oddpart(3*n+e) == 101, ('physical odd target', i))
    if i:
        require(rows == birth(inherited[i-1]['rows'], root), ('physical R', i))
    physical.append({'h':h, 't':t, 'amplitude':amp, 'n':n, 'sign':e, 'mod3':n%3})
require([x['mod3'] for x in physical] == [2,2,0,1], 'physical pruning hostile')

cover_values=[]
for core in budget['strong_core_presentations']:
    counts=paths(core['rows'])[1]
    full=len(counts)-1
    p2=sum(counts[s]*counts[full^s] for s in range(full+1))
    require(counts[full]==71 and p2==core['P2'], 'literal71 core cover count')
    require(405 % oddpart(p2) != 0, 'module readout obstruction')
    cover_values.append(p2)
require(len(cover_values)==11 and set(cover_values)=={1012,1072,1976}, 'complete input head')

def ordinal(left, right):
    a, b = len(left), len(right)
    return [row | (((1 << b)-1) << a) for row in left] + [row << a for row in right]

small = [[], [0], [2,0], [2,4,1], [6,4,0], [14,4,8,2]]
product_checks = 0
for left in small:
    for right in small:
        ih = set(paths(left)[1]); jh = set(paths(right)[1])
        both = set(paths(ordinal(left, right))[1])
        require(both == {a*b for a in ih for b in jh}, ('hereditary product', left, right))
        product_checks += 1

# C3[TT3,1,1]: the first block beats vertex3, which beats4, which beats the block.
inflated = [14, 12, 8, 16, 7]
first = ordinal(inflated, [0])
second = ordinal([2,4,1], [2,4,1])
first_counts, second_counts = paths(first)[1], paths(second)[1]
require(first_counts[-1] == second_counts[-1] == 9, 'same-H hostile')
require(set(first_counts) == {1,3,5,9} and set(second_counts) == {1,3,9},
        'hereditary sidecar hostile')

source_names = [
 'thirtysix_digits_20260925.md',
 'collatz_mod6_20260917_reverse_tree_pieces.md',
 'inverse_predecessor_sections_20261005.md',
 'central_modular_chanpak_20261110/compiler/PROOF.md',
 'h13_paired_descent_20261111/square_graceful/PROOF.md',
 'multiplicand_saturation_chanpak_20261111/chanpak/OFFSET_TRANSPORT.md',
 'changing_core_returns_20261113/returns/FIXED_BALANCED.md',
 'changing_core_returns_20261113/returns/BIRTH_BUDGET.md',
 'changing_core_returns_20261113/returns/birth_budget_receipt.json',
]
pins = [{'path':'05-knowledge/results/'+name,
         'sha256':hashlib.sha256((RESULTS/name).read_bytes()).hexdigest()}
        for name in source_names]
output = {'status':'FINITE-EXACT controls for scoped proved statements',
          'universe':{'odd_targets':'1..201, excluding multiples of3',
                      'exponents':'1..24', 'iterate_depth':'0..6',
                      'binary_width':'0..5, all suffixes', 'stencil_j':'0..24',
                      'physical_graphs':'four inherited404-ray graphs, orders10..13',
                      'module_cores':'eleven literal strong71 cores from inherited complete head',
                      'hereditary_products':'36 ordered pairs of six literal graphs'},
          'inverse_cases':inverse_checks, 'first_return_table':first_return,
          'stencil_cases':stencil_checks, 'physical_404_ray':physical,
          'strong71_P2_values':cover_values,
          'hereditary_product_cases':product_checks,
          'same_H_hostile':{'first_rows':first,'second_rows':second,
                            'H':9,'first_I':sorted(set(first_counts)),
                            'second_I':sorted(set(second_counts))}}
suffix = '.optimized' if sys.flags.optimize else ''
(HERE/('check'+suffix+'.json')).write_text(json.dumps(output,indent=2)+'\n')
(HERE/'source_pins.json').write_text(json.dumps(pins,indent=2)+'\n')
print(json.dumps({'inverse_cases':inverse_checks,'stencil_cases':stencil_checks,
                  'literal_graphs':len(physical),'hereditary_products':product_checks}))

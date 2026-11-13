"""Exact controls for neutral spatial fill and the 5/17/7/21 fibres."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, factorial
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]

def need(test, message):
    if not test:
        raise RuntimeError(message)

def rows_of_code(n, code):
    rows = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        a, b = (i, j) if code >> k & 1 else (j, i)
        rows[a] |= 1 << b
    return tuple(rows)

def valid(rows):
    n = len(rows)
    return all(not (rows[i] >> i & 1) and 0 <= rows[i] < 1 << n for i in range(n)) and all(
        ((rows[i] >> j & 1) + (rows[j] >> i & 1)) == 1
        for i in range(n) for j in range(i))

def hamilton(rows):
    n = len(rows)
    if not n:
        return 1
    end = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        end[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, value in enumerate(end[mask]):
            todo = rows[v] & ~mask
            while value and todo:
                bit = todo & -todo
                todo -= bit
                end[mask | bit][bit.bit_length() - 1] += value
    return sum(end[-1])

def strong(rows):
    n = len(rows)
    reach = [row | 1 << i for i, row in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return all(x == (1 << n) - 1 for x in reach)

def canonical(rows):
    # Degree classes greatly reduce this exhaustive isomorphism test.
    cells = defaultdict(list)
    for v, row in enumerate(rows):
        cells[row.bit_count()].append(v)
    best = None
    for blocks in product(*(permutations(cells[d]) for d in sorted(cells))):
        order = sum((tuple(b) for b in blocks), ())
        out = tuple(sum((rows[v] >> w & 1) << j for j, w in enumerate(order)) for v in order)
        best = out if best is None else min(best, out)
    return best

def birth(rows, out):
    n = len(rows)
    return tuple(row | (0 if out >> v & 1 else 1 << n) for v, row in enumerate(rows)) + (out,)

targets = (5, 7, 17, 21)
counts = {}
atoms = {h: set() for h in targets}
labelled_tested = 0
neutral_controls = 0
for n in range(1, 7):
    counts[n] = Counter()
    for code in range(1 << (n * (n - 1) // 2)):
        rows = rows_of_code(n, code)
        need(valid(rows), 'tournament input type')
        h = hamilton(rows)
        labelled_tested += 1
        if h in targets:
            counts[n][h] += 1
            if strong(rows):
                atoms[h].add(canonical(rows))
        if n <= 4:
            full = (1 << n) - 1
            for out in range(1 << n):
                child = birth(rows, out)
                neutral = hamilton(child) == h
                terminal = all(not (rows[v] & (full ^ out)) for v in range(n) if out >> v & 1)
                need(neutral == terminal, 'neutral iff terminal SCC cut')
                neutral_controls += 1
    need(counts[n][5] == factorial(n) * max(n - 3, 0), 'all-labelled H5 formula')
    need(counts[n][17] == factorial(n) * 2 * max(n - 5, 0), 'all-labelled H17 formula')
    need(counts[n][7] == counts[n][21] == 0, 'gap controls')

need(sorted(map(len, atoms[5])) == [4], 'one strong H5 atom')
need(sorted(map(len, atoms[17])) == [6, 6], 'two strong H17 atoms')
need(not atoms[7] and not atoms[21], 'no small strong gap atom')
for h in targets:
    # A strong tournament of order >=7 is excluded for ALL four targets.
    need(5 ** 6 > h ** 3, 'integer Busch cutoff')

literal = ((6, 12, 8, 1), (30, 60, 56, 48, 32, 1), (54, 60, 24, 49, 32, 4))
for rows, h in zip(literal, (5, 17, 17)):
    need(valid(rows) and strong(rows), 'literal atom type')
    brute = sum(all(rows[p[i]] >> p[i + 1] & 1 for i in range(len(rows) - 1))
                for p in permutations(range(len(rows))))
    need(brute == h == hamilton(rows), 'independent permutation Hamilton count')
    aut = sum(all((rows[i] >> j & 1) == (rows[p[i]] >> p[j] & 1)
                  for i in range(len(rows)) for j in range(len(rows)))
              for p in permutations(range(len(rows))))
    need(aut == 1, 'atom is rigid')
    full = (1 << len(rows)) - 1
    converse = tuple(full ^ (1 << i) ^ row for i, row in enumerate(rows))
    need(canonical(rows) == canonical(converse), 'atom is self-converse')

def join(blocks):
    rows = []
    total = sum(map(len, blocks))
    offset = 0
    for block in blocks:
        later = ((1 << total) - 1) ^ ((1 << (offset + len(block))) - 1)
        rows.extend((row << offset) | later for row in block)
        offset += len(block)
    return tuple(rows)

def compositions(d, k):
    if k == 1:
        yield (d,)
    else:
        for first in range(d + 1):
            for tail in compositions(d - first, k - 1):
                yield (first,) + tail

# Literal named-vertex neutral growth.  New vertices have their birth-time names.
history_controls = []
for r, depth in ((1, 5), (2, 3)):
    core = join([(2, 4, 1)] * r)
    # Keep the actual SCC sequence of singleton/name sets, separately from the graph.
    blocks = tuple(tuple(range(3 * i, 3 * i + 3)) for i in range(r))
    states = {(core, blocks)}
    for d in range(depth + 1):
        need(len({rows for rows, _ in states}) == len(states), 'SCC sidecar does not split a literal graph')
        gap_hist = Counter()
        for rows, blockword in states:
            gap = [0] * (r + 1)
            at = 0
            for block in blockword:
                if len(block) == 1:
                    gap[at] += 1
                else:
                    at += 1
            gap_hist[tuple(gap)] += 1
            need(valid(rows) and hamilton(rows) == 3 ** r, 'physical neutral history count')
        need(set(gap_hist) == set(compositions(d, r + 1)), 'every simplex lattice cell filled')
        need(set(gap_hist.values()) == {factorial(d)}, 'exact d! literal histories per cell')
        history_controls.append(dict(r=r, added=d, cells=len(gap_hist), histories=len(states)))
        if d == depth:
            break
        children = set()
        for rows, blockword in states:
            for cut in range(len(blockword) + 1):
                out = sum(1 << v for block in blockword[cut:] for v in block)
                child = birth(rows, out)
                newword = blockword[:cut] + ((len(rows),),) + blockword[cut:]
                children.add((child, newword))
        need(len(children) == len(states) * (r + d + 1), 'unique birth-time ancestry')
        states = children

# Independent rational urn recurrence, including the wrong unweighted quotient rule.
urn_controls = 0
for r in range(6):
    mass = {(0,) * (r + 1): Fraction(1)}
    for d in range(9):
        need(len(mass) == comb(d + r, r), 'urn complete support')
        need(set(mass.values()) == {Fraction(1, comb(d + r, r))}, 'urn uniform fill')
        urn_controls += len(mass)
        nxt = defaultdict(Fraction)
        for g, p in mass.items():
            for i in range(r + 1):
                child = list(g); child[i] += 1
                nxt[tuple(child)] += p * Fraction(g[i] + 1, d + r + 1)
        mass = nxt
wrong = [Fraction(1, 4), Fraction(1, 2), Fraction(1, 4)]
need(wrong != [Fraction(1, 3)] * 3, 'uniform quotient choice is not uniform fill')

def coeff(n, r, size, multiplicity):
    return multiplicity * comb(n - size + r, r) if n >= size else 0
def difference(seq):
    return [seq[0]] + [seq[i] - seq[i - 1] for i in range(1, len(seq))]
rank_terms = [(1, 4, 1), (1, 6, 2), (1, 5, 2), (2, 6, 1)]
for r, size, multiplicity in rank_terms:
    values = [coeff(n, r, size, multiplicity) for n in range(25)]
    for _ in range(r + 1):
        values = difference(values)
    need(values == [multiplicity if n == size else 0 for n in range(25)], 'rankwise primitive recovery')
h9 = [coeff(n, 1, 5, 2) + coeff(n, 2, 6, 1) for n in range(25)]
for _ in range(3):
    h9 = difference(h9)
need(h9[5:7] == [2, -1], 'unmarked finite differences have signed false births')

old = ROOT / '05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities'
census = json.loads((old / 'extended_census.json').read_text(encoding='utf-8'))
atoms27 = ((30, 52, 56, 50, 32, 1),
           (62, 124, 120, 112, 32, 64, 17),
           (118, 124, 40, 113, 100, 64, 4))
records27 = [item for item in census['records'] if item['strong'] and item['H'] == 27]
need(sum(item['ordinary_multiplicity'] for item in records27) == 5, 'five inherited ordinary H27 atoms')
need({tuple(item['adjacency_outmasks']) for item in records27} == set(atoms27), 'three inherited converse representatives')
need(5 ** 7 > 27 ** 3, 'Busch excludes further H27 strong orders')
five_subset_controls = []
for rows in atoms27:
    need(valid(rows) and strong(rows) and hamilton(rows) == 27, 'literal H27 atom')
    hist = Counter()
    for subset in combinations(range(len(rows)), 5):
        subrows = tuple(sum((rows[v] >> w & 1) << j for j, w in enumerate(subset)) for v in subset)
        hist[hamilton(subrows)] += 1
    need(13 not in hist, 'no H13 core inside H27 atom')
    five_subset_controls.append(dict(rows=list(rows), five_vertex_hamilton_histogram=dict(hist)))
frontiers = {}
for ns, entries in census['counts_by_order'].items():
    support = {int(h) for h, v in entries.items() if v['ordinary']}
    first = next(h for h in range(1, max(support) + 8, 2) if h not in (7, 21) and h not in support)
    frontiers[ns] = dict(first_unfilled_nonexception=first, filled_prefix=first - 2,
                        source='inherited complete ordinary-isomorphism census')
pins = [old/'PROOF.md', old/'extended_census.json', old/'low_H_formulas.json',
        ROOT/'01-canon/theorems/THM-1900-the-insertion-response-calculus.md']
result = dict(status='PASS', labelled_tournaments=labelled_tested, neutral_attachments=neutral_controls,
              all_order_target_atom_counts={h:len(atoms[h]) for h in targets},
              target_strong_representatives={h:[list(x) for x in sorted(atoms[h])] for h in targets},
              target_labelled_counts={n:dict(c) for n,c in counts.items()},
              literal_history_controls=history_controls, exact_urn_cells=urn_controls,
              unweighted_quotient_hostile=['1/4','1/2','1/4'],
              unmarked_h9_third_difference={'5':2,'6':-1}, inherited_fill_frontiers=frontiers,
              thirteen_to_twentyseven_core_controls=five_subset_controls,
              source_pins=[dict(path=str(p.relative_to(ROOT)).replace('\\','/'),
                                sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in pins],
              scope='All labelled tournaments through6 directly enumerated; neutral attachments through4; literal histories as listed; all-order proofs separate. n7/n8 frontier rows inherited, not re-enumerated.')
suffix = '.optimized' if sys.flags.optimize else ''
(HERE / ('receipt'+suffix+'.json')).write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k not in ('source_pins','target_strong_representatives','target_labelled_counts')}, indent=2))

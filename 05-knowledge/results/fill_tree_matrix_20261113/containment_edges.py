"""All-order count-only induced containment for targets through73.

The inherited complete class bank is not regenerated. Every subset count
is independently evaluated in forward and backward endpoint coordinates.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys
from tree_blueprint import contracted_parent, HOLES

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BANK = HERE.parent / 'actual_hamilton_multiplicities_20261021/session24-multiplicities'

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def counts(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    totals = [0] * (1 << n)
    totals[0] = 1
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            val = dp[mask][v]
            if not val:
                continue
            nxt = rows[v] & ~mask
            while nxt:
                bit = nxt & -nxt
                w = bit.bit_length() - 1
                dp[mask | bit][w] += val
                nxt ^= bit
        totals[mask] = sum(dp[mask])
    return totals

def main():
    need(5**8 > 73**3, 'all-order strong cutoff')
    records = json.loads((BANK / 'census.json').read_text(encoding='utf-8'))['records']
    cores = [(r['H'], tuple(r['adjacency_outmasks']), 'census:' + r['node'])
             for r in records if r['strong'] and r['H'] <= 73]
    with (BANK / 'n8_classes.csv').open(newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            if int(r['strong']) and int(r['H']) <= 73:
                code = int(r['code'])
                cores.append((int(r['H']), tuple((code >> (8*i)) & 255 for i in range(8)), 'n8:' + r['code']))
    inherited_atoms = {}
    subset_tests = 0
    union = {h: set() for h in range(1, 74, 2)}
    witness = {h: {} for h in union}
    union[1] = {1}
    witness[1][1] = dict(kind='empty')
    for h, rows, source in cores:
        n = len(rows)
        full = (1 << n) - 1
        for i in range(n):
            need(not rows[i] >> i & 1, 'loop')
            for j in range(i):
                need(((rows[i] >> j) ^ (rows[j] >> i)) & 1, 'tournament')
        a = counts(rows)
        b = counts([full ^ (1 << i) ^ rows[i] for i in range(n)])
        need(a == b and a[-1] == h, ('subset and bank counts', source))
        subset_tests += len(a)
        inherited_atoms[h] = inherited_atoms.get(h, 0) + 1
        for mask, value in enumerate(a):
            union[h].add(value)
            witness[h].setdefault(value, dict(kind='strong', rows=rows, subset_mask=mask, bank=source))
    for h in range(3, 74, 2):
        for a in range(3, h, 2):
            if h % a:
                continue
            b = h // a
            for x in union[a]:
                for y in union[b]:
                    union[h].add(x*y)
                    witness[h].setdefault(x*y, dict(kind='join', factors=[a,b], subset_counts=[x,y]))
    missing_edges = []
    edge_witnesses = {}
    for h in range(3, 74, 2):
        if h in HOLES:
            need(not union[h], 'forbidden target')
            continue
        need(h in union[h], ('ordinary target present', h))
        p = contracted_parent(h)
        if p not in union[h]:
            missing_edges.append([p, h])
        else:
            edge_witnesses[h] = dict(parent=p, witness=witness[h][p])
    report = dict(status='PASS', scope='All orders for target H<=73, using inherited complete strong census through8 and Busch cutoff',
                  inherited_core_representatives=len(cores), subset_tests=subset_tests,
                  missing_arithmetic_tree_edges=missing_edges,
                  all_order_induced_count_unions={h: sorted(s) for h,s in union.items()},
                  arithmetic_tree_edge_witnesses=edge_witnesses,
                  source_pins={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in (BANK/'census.json', BANK/'n8_classes.csv')})
    suffix = '.optimized' if sys.flags.optimize else ''
    (HERE / ('containment_edges' + suffix + '.json')).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('all_order_induced_count_unions','arithmetic_tree_edge_witnesses','source_pins')}))

if __name__ == '__main__':
    main()

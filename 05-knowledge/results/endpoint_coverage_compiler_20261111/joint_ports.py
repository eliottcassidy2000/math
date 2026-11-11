"""Exact two-port quotient words; no interpolation from a finite size grid."""
from collections import defaultdict, Counter
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
import json
import composition_probe as cp
import inflation_bank as ib

HERE = Path(__file__).resolve().parent


@lru_cache(maxsize=None)
def ends(rows, marked):
    n = len(rows)
    old = [v for v in range(n) if v not in marked]
    order = old + list(marked)
    m = len(old)
    # A two-vertex tournament is acyclic. Each visit to either marked vertex
    # after its first on a path needs another old vertex: the cap n is safe
    # also for a pair of paths. Increasing masks/visit tuples is a DAG.
    dp = [dict() for _ in range(1 << m)]
    for mask in range(1 << m):
        for counts in product(range(n+1), repeat=2):
            value = [0]*n
            for v in range(m):
                if mask >> v & 1:
                    rest = mask ^ (1 << v)
                    if rest == 0 and counts == (0,0):
                        value[v] = 1
                    else:
                        prior = dp[rest].get(counts, [0]*n)
                        value[v] = sum(prior[u] for u in range(n) if rows[order[u]] >> order[v] & 1)
            for b in range(2):
                if counts[b]:
                    before = tuple(counts[i] - (i == b) for i in range(2))
                    v = m+b
                    value[v] = int(mask == 0 and before == (0,0))
                    prior = dp[mask].get(before, [0]*n)
                    value[v] += sum(prior[u] for u in range(n) if rows[order[u]] >> order[v] & 1)
            if any(value):
                dp[mask][counts] = value
    return order, dp


def defect_coefficients(rows, marked, selected):
    order, end = ends(tuple(rows), tuple(marked))
    reverse = tuple(sum(1 << b for b in range(len(rows)) if rows[b] >> a & 1) for a in range(len(rows)))
    _, start = ends(reverse, tuple(marked))
    full = len(end)-1
    selected_positions = [i for i,v in enumerate(order) if v in selected]
    d = defaultdict(int)
    for counts, values in end[full].items():
        d[counts] += sum(values)
    for mask in range(full+1):
        for a, left in end[mask].items():
            x = sum(left[v] for v in selected_positions)
            if not x:
                continue
            for b, right in start[full^mask].items():
                y = sum(right[v] for v in selected_positions)
                if y:
                    counts = (a[0]+b[0], a[1]+b[1])
                    d[counts] -= x*y
    # Both substituted modules are nonempty. P_0(size)=0 there.
    return {k:v for k,v in d.items() if all(k) and v}


def main():
    data = json.loads((HERE / 'inflation_bank_7.json').read_text())['records']
    groups = defaultdict(list)
    for item in data:
        groups[(item['n'], item['code'], item['selected'])].append(item)
    checked, failures, size_controls = 0, [], 0
    levels = Counter()
    for key, items in groups.items():
        for left, right in combinations(items, 2):
            rows = left['rows']
            selected = [v for v in range(len(rows)) if key[2] >> v & 1]
            marked = [left['inflated_vertex'], right['inflated_vertex']]
            coeff = defect_coefficients(rows, marked, selected)
            checked += 1
            levels[key[0]] += 1
            if coeff != {(1,1):1}:
                failures.append(dict(n=key[0], code=key[1], rows=rows, selected=key[2], marked=marked,
                                     defect_terms=[[a,b,d] for (a,b),d in sorted(coeff.items())]))
            # Separately implemented literal graph at a mixed size.
            if key[0] <= 5 or checked % 17 == 0 or coeff != {(1,1):1}:
                expanded, es = cp.expand(rows, selected, dict(zip(marked, (2,3))))
                actual, _ = cp.defect(expanded, es)
                predicted = sum(d*ib.onto(2,a)*ib.onto(3,b) for (a,b),d in coeff.items())
                ib.tc.require(actual == predicted, 'two-port mixed-size control')
                size_controls += 1
        # Avoid retaining DP tables for every base across the whole census.
        ends.cache_clear()
    report = dict(separately_stable_pairs=checked, jointly_unstable= len(failures),
                  pairs_by_order=[levels[n] for n in range(1,8)],
                  mixed_size_literal_controls=size_controls, failures=failures)
    (HERE / 'joint_ports.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()

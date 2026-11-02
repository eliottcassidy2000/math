"""Exact endpoints and an independent native audit of discarded-degree cancellation."""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import json

GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(why)


def graph(n, weights):
    edges = list(combinations(range(n), 2))
    remaining = {e: w for e, w in zip(edges, weights) if w}
    degree, discarded_degree = [0]*n, [0]*n
    for (u, v), w in remaining.items():
        degree[u] += w; degree[v] += w
    twice_wedge = sum(x*x for x in degree)-2*sum(w*w for w in remaining.values())
    total_discarded = 0
    while remaining:
        chosen = max(remaining, key=lambda e: (remaining[e], e))
        chosen_weight = remaining[chosen]
        group = [e for e in remaining if set(e) & set(chosen)]
        for e in group:
            w = remaining.pop(e)
            if e != chosen:
                need(chosen_weight >= w, 'greedy square invoice is paid')
                total_discarded += w
                for v in e:
                    discarded_degree[v] += w
    need(twice_wedge >= sum(x*x for x in discarded_degree), 'native discarded-degree cancellation')
    need(n*twice_wedge >= 4*total_discarded**2, 'native stronger global bound')
    return twice_wedge//2, total_discarded


def main():
    c, t, d, r = F(999, 4000), F(4, 119), F(1793, 10000), 119
    lam = t*t/d
    q = r+F(33, 2)+F(8, r)
    need(t < c/2, 'small variance exclusion')
    need(c-4*t > 0, 'heavy-edge lower bound increasing')
    need((c-4*t)*r-4-F(2, r) == F(4618839, 476000) > 9, 'ten heavy edges endpoint')
    slack = c-d-q/(r-10)*F(3, 4)**10
    need(slack == F(6985054879, 17001349120000) > 0, 'tight geometric tail')
    need(q*F(3, 4)**(r-1) < c-d, 'other geometric endpoint')
    need(F(3, 4)*(1+F(1, r)) < 1, 'tail endpoint monotonicity')
    need(d-4*t > 0, 'positive light deficit')
    need(lam == F(160000, 25390673), 'exact dual parameter')
    need(3*t/4-4*lam == F(101, 25390673) > 0, 'nonnegative dual remainder')
    need(t-4*lam == F(213468, 25390673) > 0, 'positive total coefficient')
    need(t*119 == 4, 'final contradiction including exact endpoint')
    # Independent coefficient expansion in the variables S,w.
    square = {(1, 0): F(1, 2), (0, 0): -t, (0, 1): F(1, 2)}
    rhs = {}
    for e, x in square.items():
        for f, y in square.items():
            key = tuple(a+b for a, b in zip(e, f))
            rhs[key] = rhs.get(key, 0)+x*y
    rhs[(0, 1)] += t-4*lam
    rhs[(0, 2)] -= F(1, 4)
    lhs = {(0, 0): lam*d, (1, 0): -t, (2, 0): F(1, 4),
           (1, 1): F(1, 2), (0, 1): -4*lam}
    need({e: x for e, x in rhs.items() if x} == lhs, 'dual polynomial identity')
    for n in range(238, 1500):
        rr = n//2
        qq = rr+F(33, 2)+F(8, rr)
        need((c-4*t)*rr-4-F(2, rr) > 9, 'finite heavy-edge floor replay')
        need(max(qq*F(3, 4)**(rr-1), qq/(rr-10)*F(3, 4)**10) < c-d,
             'finite tail replay')
        need(F(4, rr) <= t, 'finite closing coefficient replay')
    graph_count = 0
    for n in range(2, 6):
        for weights in product((0, 1), repeat=n*(n-1)//2):
            graph(n, weights)
            graph_count += 1
    for weights in product((0, 1, 2), repeat=6):
        graph(4, weights)
        graph_count += 1
    sharp = []
    for rr in range(2, 21):
        weights = [int(u < rr <= v) for u, v in combinations(range(2*rr), 2)]
        wedge, mass = graph(2*rr, weights)
        ratio = F(2*rr*wedge, mass*mass)
        need(ratio == F(2*rr, rr-1), 'asymptotic sharpness family of general graph lemma')
        sharp.append({'r': rr, 'nW_over_D_squared': str(ratio)})
    report = {'status': 'PASS', 'gates': GATES, 'threshold': 238,
              'remaining_open': [10, 237], 'closed_since_original_1500': 1262,
              'graph_universe_count': graph_count, 'bipartite_controls': sharp,
              'constants': {'c': str(c), 't': str(t), 'd': str(d), 'lambda': str(lam)},
              'tight_tail_slack': str(slack)}
    Path(__file__).with_name('stability_238_checks.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'gates': GATES, 'threshold': 238, 'graphs': graph_count}))


if __name__ == '__main__':
    main()

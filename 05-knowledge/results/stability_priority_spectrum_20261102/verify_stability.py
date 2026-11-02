"""Exact endpoint audit and native controls; the infinite proof is in STABILITY_466.md."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, factorial
from pathlib import Path
import json
import random

GATES = 0


def need(test, message):
    global GATES
    GATES += 1
    if not test:
        raise RuntimeError(message)


def kappa(n):
    r = n // 2
    return F(4, r) * (1 - F(3, 4)**r - (F(1, 4)**r if n % 2 == 0 else 0))


def det(a):
    a = [[F(x) for x in row] for row in a]
    answer = F(1)
    for i in range(len(a)):
        j = next((j for j in range(i, len(a)) if a[j][i]), None)
        if j is None:
            return F(0)
        if i != j:
            a[i], a[j] = a[j], a[i]
            answer = -answer
        pivot = a[i][i]
        answer *= pivot
        for j in range(i+1, len(a)):
            ratio = a[j][i] / pivot
            for k in range(i+1, len(a)):
                a[j][k] -= ratio * a[i][k]
    return answer


def native(a):
    n = len(a)
    dp = {(1 << v, v): F(1) for v in range(n)}
    for mask in range(1, 1 << n):
        for v in range(n):
            value = dp.get((mask, v), 0)
            if not value:
                continue
            for u in range(n):
                if (mask >> u) & 1 or not a[v][u]:
                    continue
                key = (mask | (1 << u), u)
                dp[key] = dp.get(key, 0) + value * a[v][u]
    h = sum(dp.get(((1 << n)-1, v), 0) for v in range(n))
    u = det([[a[i][j] + (i == j) for j in range(n)] for i in range(n)])
    d = h - u - abs(det(a))
    qs = {(i, j): a[i][j]*a[j][i] for i, j in combinations(range(n), 2)}
    wedges = {frozenset((e, f)): qs[e]*qs[f]
              for e, f in combinations(qs, 2) if set(e) & set(f)}
    return h, d, sum(qs.values()), qs, wedges


def fair_pairs(n):
    a = [[F(i < j) for j in range(n)] for i in range(n)]
    for i in range(0, n-1, 2):
        a[i][i+1] = a[i+1][i] = F(1, 2)
    return a


def greedy_audit(qs, wedges, n):
    remaining = {e for e in qs if qs[e]}
    assigned_wedges, discarded = set(), F(0)
    while remaining:
        edge = max(remaining, key=lambda e: (qs[e], e))
        group = {f for f in remaining if set(f) & set(edge)}
        lost = group - {edge}
        mass = sum(qs[f] for f in lost)
        used = {frozenset((e, f)) for e, f in combinations(sorted(group), 2)
                if set(e) & set(f)}
        need(not (used & assigned_wedges), 'group contributions disjoint')
        assigned_wedges |= used
        budget = sum(wedges[pair] for pair in used)
        need(budget >= qs[edge]*mass, 'group selected/discarded wedges')
        need(budget >= mass*mass/4, 'group square mass bound')
        discarded += mass
        remaining -= group
    need(sum(wedges.values()) >= discarded**2/F(2*n), 'global discarded bound')


def endpoint_checks():
    c, eta, delta = F(999, 4000), F(1, 58), F(11, 100)
    astar = c - 4*eta - delta
    need(sum(F(8**j, factorial(j)) for j in range(8)) == F(47259, 35), 'exp8 partial sum')
    need(F(47259, 35) > 1000, 'exp8 lower bound')
    need(F(4, 233) < c/2, 'small variance exclusion')
    g233 = (c-4*eta)*233 - 16 - F(8, 233)
    need(g233-26 == F(2390619, 27028000) > 0, 'retained matching endpoint')
    need(c-4*eta-F(1, 6) == F(4913, 348000) > 0, 'retained matching slope')
    need(25850*3**27 < 11*4**27, 'epsilon endpoint')
    need(F(3, 4)*(F(517, 2)+6) < F(517, 2), 'epsilon step; slope 3/4 < 1')
    need(astar == F(8211, 116000), 'astar exact')
    need(astar-4*eta == F(211, 116000) > 0, 'positive left coefficient')
    need(astar-eta == F(6211, 116000) > 0, 'small right remainder')
    need(8000*3**233 < 4**233, 'rho endpoint')
    need(F(232, 233)+F(1, 1000) < 1, 'final contradiction')
    need(eta > F(4, 233), 'no-light strict boundary')
    for n in range(466, 1500):
        r, k = n//2, kappa(n)
        j = (r-233)//6
        eps = F(3, 4)**(27+j)
        need(0 < k <= F(4, r), 'candidate bounds')
        need(eps*(F(n, 2)+20) < delta, 'finite tail replay')
        need(k/eta+eta*(4-k*r)/astar < 1, 'finite contradiction replay')
    # Independently expand the deletion identity in a sparse polynomial ring.
    zero = (0,)*8
    def const(x):
        return {zero: F(x)}
    def var(i):
        ex = list(zero); ex[i] = 1
        return {tuple(ex): F(1)}
    def add(*ps):
        out = {}
        for p in ps:
            for e, x in p.items():
                out[e] = out.get(e, 0)+x
        return {e: x for e, x in out.items() if x}
    def mul(p, q):
        out = {}
        for e, x in p.items():
            for f, y in q.items():
                z = tuple(a+b for a, b in zip(e, f))
                out[z] = out.get(z, 0)+x*y
        return {e: x for e, x in out.items() if x}
    def neg(p):
        return {e: -x for e, x in p.items()}
    k, v, b, r, ell, eps, d, eta_p = map(var, range(8))
    four, one = const(4), const(1)
    om = add(one, neg(eps))
    raw = add(mul(add(mul(k, v), neg(b)), add(r, neg(ell))),
              neg(mul(mul(four, om), add(v, neg(d), neg(mul(eta_p, ell))))))
    rhs = add(mul(add(mul(k, r), neg(four), mul(four, eps)), v),
              mul(mul(four, om), d), neg(mul(b, r)))
    bracket = add(mul(k, v), neg(b), neg(mul(mul(four, eta_p), om)))
    need(raw == add(rhs, neg(mul(ell, bracket))), 'deletion identity as polynomial')
    return {'threshold': 466, 'formerly_open_orders_closed': 1034,
            'remaining_open': [10, 465], 'constants': list(map(str, (c, eta, delta, astar)))}


def native_checks():
    pair_data = {}
    for n in range(4, 13):
        h, d, v, qs, wedges = native(fair_pairs(n))
        need(h == 1 and d == kappa(n)*v, 'native fair-pair formula')
        greedy_audit(qs, wedges, n)
        pair_data[n] = {'V': str(v), 'Delta': str(d), 'ratio': str(d/v),
                        'global_optimality_proved_here': n <= 9}
    rng = random.Random(466)
    for n in range(4, 7):
        edges = list(combinations(range(n), 2))
        for trial in range(8):
            a = [[F(0) for _ in range(n)] for _ in range(n)]
            uncertain = rng.sample(edges, 4)
            for e in edges:
                i, j = e
                p = rng.choice([F(1, 4), F(1, 2), F(3, 4)]) if e in uncertain else F(rng.randrange(2))
                a[i][j], a[j][i] = p, 1-p
            h, d, v, qs, wedges = native(a)
            baseline, rounded_h = F(0), F(0)
            for bits in product(range(2), repeat=len(uncertain)):
                t, prob = [row[:] for row in a], F(1)
                for (i, j), bit in zip(uncertain, bits):
                    prob *= a[i][j] if bit else a[j][i]
                    t[i][j], t[j][i] = F(bit), F(1-bit)
                bh, bd, _, _, _ = native(t)
                baseline += prob*bd
                rounded_h += prob*bh
            need(rounded_h == h, 'native independent rounding H identity')
            need(baseline >= sum(wedges.values()), 'native binary baseline dominates wedges')
            need(sum(wedges.values()) >= F(2, n)*v*v-v/4, 'native global wedge lower bound')
            need(d >= kappa(n)*v, 'proved small-order native stability')
            greedy_audit(qs, wedges, n)
    # Exact hostile to copying the five-vertex inequality into the full order.
    local_data = {}
    for n in (6, 11):
        a = fair_pairs(n)
        _, full_gap, v, _, _ = native(a)
        local_sum = sum(native([[a[i][j] for j in sub] for i in sub])[1]
                        for sub in combinations(range(n), 5))
        r = n//2
        formula = comb(n-2, 3)*v-F((n-4)*comb(r, 2), 16)
        need(local_sum == formula, 'five-vertex incidence exact formula')
        normalized = local_sum/comb(n-2, 3)
        need(full_gap < normalized, 'false direct local-to-global transfer rejected')
        local_data[n] = {'global_gap': str(full_gap), 'edge_normalized_five_gap': str(normalized),
                         'defect': str(normalized-full_gap)}
    need(local_data[6]['defect'] == '3/32', 'minimal named transfer hostile')
    return {'fair_pair_data': pair_data, 'local_transfer_hostiles': local_data,
            'random_native_controls': 24, 'random_seed': 466}


def main():
    report = {'scope': 'Exact finite certificate; infinite quantifiers proved in STABILITY_466.md',
              'endpoints': endpoint_checks(), 'native': native_checks()}
    report['gates'] = GATES
    Path(__file__).with_name('stability_checks.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'gates': GATES, 'threshold': 466}))


if __name__ == '__main__':
    main()

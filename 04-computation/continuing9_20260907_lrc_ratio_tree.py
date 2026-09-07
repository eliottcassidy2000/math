#!/usr/bin/env python3
"""Prescribed-clock ratio-tree gluing and a composite-path7200 consumer.

All arithmetic is exact. No mathematical producer is imported or executed.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations, combinations_with_replacement, product
from math import gcd, lcm
from hashlib import sha256
import json
import sys
sys.stdout.reconfigure(encoding='utf-8', newline='\n')
HERE = Path(__file__).resolve()
ROOT = HERE.parent.parent if HERE.parent.name == '04-computation' else Path('C:/w/s0905')
DEST = ROOT / '05-knowledge/results' if HERE.parent.name == '04-computation' else HERE.parent
GATES = 0


def need(x, why):
    global GATES
    GATES += 1
    if not x:
        raise ArithmeticError(why)


def canonical(x):
    return json.dumps(x, sort_keys=True, separators=(',', ':')).encode()


def factors(n):
    out = {}; p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0)+1; n //= p
        p += 1
    if n > 1:
        out[n] = 1
    return out


def vp(n, p):
    v = 0
    while n % p == 0:
        n //= p; v += 1
    return v


def projective(n, edges):
    adj = [[] for _ in range(n)]
    for i, j, ratio in edges:
        adj[i].append((j, ratio)); adj[j].append((i, 1/ratio))
    values = [None]*n; values[0] = F(1); todo = [0]
    for i in todo:
        for j, ratio in adj[i]:
            if values[j] is None:
                values[j] = values[i]*ratio; todo.append(j)
            else:
                need(values[j] == values[i]*ratio, 'consistent directed ratio propagation')
    need(len(todo) == n, 'connected tree or consistent connected ratio graph')
    D = lcm(*(x.denominator for x in values))
    integers = [int(D*x) for x in values]; g = gcd(*integers)
    return values, tuple(x//g for x in integers)


def valuation_iff(t, d, values):
    if len(set(values)) < len(values) or any(t % a for a in d):
        return False
    for p, k in factors(t).items():
        h = [vp(x.numerator, p)-vp(x.denominator, p) for x in values]
        m = min(h)
        if [min(k, a-m) for a in h] != [vp(a, p) for a in d]:
            return False
    return True


def atlas_sum(s):
    return all(p % 3 == 2 and e <= 2 for p, e in factors(s).items())


def geometry(p, q):
    L = 14*p*q
    walls = {0, L}
    for v, z in [(p, q), (q, p)]:
        for k in range(v+1):
            for sign in [-1, 1]:
                x = (14*k+sign)*z
                if 0 < x < L:
                    walls.add(x)
    walls = sorted(walls)
    def danger(x, v):
        r = v*x % (2*L)
        return 14*min(r, 2*L-r) < 2*L
    intervals = [(a, b) for a, b in zip(walls, walls[1:]) if danger(a+b, p) and danger(a+b, q)]
    need(intervals[0][0] == 0 and intervals[-1][1] == L, 'origin split is genuine')
    return L, [(intervals[-1][0]-L, intervals[0][1])] + intervals[1:-1]


def grid_count(n, L, intervals, phase2):
    den = 2*L
    return sum(-((-(2*n*b-phase2))//den)-(2*n*a-phase2)//den-1 for a, b in intervals)


def profile(n, p, q):
    L, intervals = geometry(p, q)
    walls = sorted({n*x % L for ab in intervals for x in ab})
    phases = [2*x for x in walls] + [walls[i] + (walls[i+1] if i+1 < len(walls) else walls[0]+L) for i in range(len(walls))]
    rows = [[r, grid_count(n, L, intervals, r)] for r in phases]
    return L, intervals, walls, rows


def literal_pair(n, p, q, phase):
    A, B = phase.numerator, phase.denominator
    den = B*n
    return sum(all(14*min((v*(B*j+A)) % den, den-(v*(B*j+A)) % den) < den for v in (p, q)) for j in range(n))


def mst_tree(d, weights):
    parents = list(range(7)); edges = []
    def find(i):
        while parents[i] != i:
            i = parents[i]
        return i
    for cost, i, j in sorted((weights[tuple(sorted((d[i], d[j])))][0], i, j) for i, j in combinations(range(7), 2)):
        a, b = find(i), find(j)
        if a != b:
            parents[a] = b
            e, p, q = weights[tuple(sorted((d[i], d[j])))][1:4]
            edges.append((i, j, e, p, q))
    return edges


def orientations(d, edges):
    choices = []
    for i, j, e, p, q in edges:
        n = 7200//e; local = []
        if (e*gcd(n, p), e*gcd(n, q)) == (d[i], d[j]):
            local.append((i, j, F(q, p)))
        if (e*gcd(n, q), e*gcd(n, p)) == (d[i], d[j]):
            local.append((i, j, F(p, q)))
        need(bool(local), 'selected ratio has a local margin orientation')
        choices.append(local)
    return product(*choices)


def main():
    raw = (ROOT/'05-knowledge/results/continuing8_20260906_lrc_minimum_tree_certificate.json').read_bytes()
    need(sha256(raw).hexdigest() == '580a7c930103aab3bea867ad463a90b0e0208323a90ee95a685ff811a761582d', 'frozen complete7200 source certificate')
    inherited = json.loads(raw); R = next(r for r in inherited['clocks'] if r['t'] == 7200)
    weights = {tuple(ab): v for ab, v in R['weights']}
    need(len(R['survivors']) == 15, 'complete15 projected survivor list')
    rawP = (ROOT/'04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_bytes()
    need(sha256(rawP).hexdigest() == '935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f', 'complete hereditary profiles pin')
    J = json.loads(rawP)
    P = {int(k): {(c, tuple(w)) for c, w in row['profiles']} for k, row in J['levels'].items()}
    sets = [(I, tuple(j for j in range(7) if j not in I)) for size in range(6, 0, -1) for I in combinations(range(7), size)]
    need(len(sets) == 126, 'complete positional profile universe')
    def valid(d):
        if gcd(*d) != 1:
            return False
        for I, K in sets:
            c = gcd(*(d[i] for i in I))
            if (c, tuple(sorted(gcd(c, d[j]) for j in K))) not in P[len(K)]:
                return False
        return True
    def excess(d):
        return sum(a*((7200+7*a-1)//(7*a)) for a in d)-7200

    # General iff controls, including rejection of all false proposed margins.
    generic_rows = 0
    for row in combinations(range(1, 13), 3):
        if gcd(*row) != 1:
            continue
        for T in [[(0, 1), (1, 2)], [(0, 1), (0, 2)]]:
            edges = [(i, j, F(row[j], row[i])) for i, j in T]
            v, u = projective(3, edges)
            need(u == row, 'rational LCM normalization recovers primitive row')
            generic_rows += 1
            for t in range(1, 25):
                d = tuple(gcd(t, x) for x in row)
                need(valuation_iff(t, d, v), 'prime-potential formula on actual rows')
    v, u = projective(3, [(0, 1, F(4)), (1, 2, F(3, 2))])
    need(u == (1, 4, 6) and not valuation_iff(2, (1, 2, 1), v), 'saturated middle loses absolute depth')
    need(atlas_sum(5), 'both minimal hostile edges belong to strict sum5 atlas')
    need((gcd(2, 1), gcd(2, 4)) == (1, 2) and (gcd(2, 2), gcd(2, 3)) == (2, 1), 'both hostile edges locally compatible')
    for t in range(1, 25):
        divisors = [a for a in range(1, t+1) if t % a == 0]
        for d in product(divisors, repeat=3):
            need(valuation_iff(t, d, v) == (d == tuple(gcd(t, x) for x in u)), 'iff including every false divisor labelling')
    vv, uu = projective(3, [(0, 1, F(2)), (1, 2, F(1, 2))])
    need(uu == (1, 2, 1) and not valuation_iff(2, (1, 2, 1), vv), 'distinctness is a separate required condition')

    # This is a complete fixed-tree/first-attainer orientation bank, not every
    # low-credit ratio choice or every possible connecting tree.
    fixed = []; realized = 0; orientation_count = 0; duplicates = 0; valuation_failures = 0
    for d, E, M in R['survivors']:
        need(valid(d) and excess(d) == E, 'survivor retains every inherited profile')
        T = mst_tree(d, weights); rows = []
        need(sum(weights[tuple(sorted((d[i], d[j])))][0] for i, j, e, p, q in T) == M, 'deterministic minimum tree cost')
        for edges in orientations(d, T):
            v, u = projective(7, edges); orientation_count += 1
            margins = tuple(gcd(7200, x) for x in u)
            good = margins == tuple(d) and len(set(u)) == 7
            need(valuation_iff(7200, d, v) == good, 'valuation iff versus independent rational normalization')
            if good:
                realized += 1
            elif len(set(u)) < 7:
                duplicates += 1
            else:
                valuation_failures += 1
            rows.append(dict(ratios=[[i, j, r.numerator, r.denominator] for i, j, r in edges], row=list(u), margins=list(margins), valid=good))
        fixed.append(dict(word=d, E=E, M=M, tree=T, orientations=rows))
    need((orientation_count, realized, duplicates, valuation_failures) == (29, 20, 8, 1), 'complete fixed-bank boundaries')
    conflict = next(row for record in fixed for row in record['orientations'] if not row['valid'] and len(set(row['row'])) == 7)
    need(conflict['margins'] == [45, 8, 9, 90, 32, 36, 48] and vp(conflict['row'][0], 3)-vp(conflict['row'][1], 3) == 2,
         'fixed-tree conflict occurs between two prescribed3-free anchors')

    U = (89012, 1154241, 27596853, 171248, 123354, 20328, 19179776)
    d = tuple(gcd(7200, x) for x in U)
    need(gcd(*U) == 1 and len(set(U)) == 7, 'actual primitive distinct seven-tail row')
    need(d == (4, 9, 9, 16, 18, 24, 32) and valid(d), 'actual row retains worst projected word')
    T = mst_tree(d, weights); actual = []
    for i, j in combinations(range(7), 2):
        D = gcd(U[i], U[j]); p, q = sorted((U[i]//D, U[j]//D)); e = gcd(7200, D)
        if p+q <= 356 and atlas_sum(p+q):
            actual.append((i, j, e, p, q, D//e))
    need([x[:5] for x in actual] == T, 'entire actual strict graph is precisely the selected tree')
    masks = []
    for j in range(7200):
        mask = 0
        for i, speed in enumerate(U):
            r = speed*(2*j+1) % 14400
            if 14*min(r, 14400-r) < 14400:
                mask |= 1 << i
        masks.append(mask)
    paircounts = [[i, j, sum(bool(m & (1 << i)) and bool(m & (1 << j)) for m in masks)] for i, j in combinations(range(7), 2)]
    pc = {(i, j): c for i, j, c in paircounts}
    need(all(pc[i, j] == 0 for i, j, e, p, q in T), 'all six atlas minima coexist at the same alpha1/2')
    marginals = [sum(bool(m & (1 << i)) for m in masks) for i in range(7)]
    histogram = [sum(m.bit_count() == k for m in masks) for k in range(8)]
    need(marginals == [1032, 1026, 1026, 1024, 1044, 1008, 1056], 'literal native danger marginals')
    need(histogram == [2020, 3394, 1548, 226, 12, 0, 0, 0], 'exact7200-cell danger histogram')
    need(sum(marginals)-sum(bool(m) for m in masks) == 2036, 'actual multiplicity excess is in non-atlas overlaps')
    need(pc[0, 4] == 232, 'two-edge composite pair sees the missing overlap at alpha1/2')

    # The shortest composite path has ratio578/801 and e=2. Every real phase,
    # including every wall, is certified by the following complete arrangement.
    need(U[0] == 154*578 and U[5] == 154*132 and U[4] == 154*801, 'actual primitive two-arm pattern')
    need(gcd(578, 801) == 1 and gcd(7200, 154) == 2, 'composite primitive ratio and sheet multiplicity')
    need(F(578, 132) == F(289, 66) and F(801, 132) == F(267, 44), 'composite path preserves both directed arms')
    L, intervals, walls, rows = profile(3600, 578, 801)
    lo = min(c for r, c in rows); hi = max(c for r, c in rows)
    need(len(intervals) == 197 and len(walls) == 344 and len(rows) == 688, 'complete unpruned composite wall/chamber arrangement')
    need(lo == 57, 'uniform composite pair credit114')
    for r, c in rows:
        need(0 <= c <= 3600, 'native count bounds')
    low = next(F(r, 2*L) for r, c in rows if c == lo)
    for phase in [F(0), F(1, 2), low, F(1, 7), F(3, 623)]:
        count = grid_count(3600, L, intervals, 2*L*phase)
        need(count.denominator == 1 if isinstance(count, F) else True, 'rational phase count integral')
        need(count == literal_pair(3600, 578, 801, phase), 'independent literal pair-grid control')
    need(2*lo > 74 and 2*grid_count(3600, L, intervals, L) == 232, 'uniform and displayed-phase credits distinguished')

    # Complete conditional word universe, with one occurrence of4,18,24 forced.
    alphabet = sorted(a for a in J['levels']['6']['gcds'] if 7200 % a == 0)
    conditional = []; attempted = 0
    for rest in combinations_with_replacement(alphabet, 4):
        attempted += 1; word = tuple(sorted((4, 18, 24)+rest))
        if valid(word):
            conditional.append([list(word), excess(word)])
    need(attempted == 23751 and len(conditional) == 1033, 'complete four-free-position profile census')
    conditional.sort()
    maximum = max(E for word, E in conditional)
    need(maximum == 108 and 2*lo-maximum == 6, 'uniform two-arm full-profile closure margin')
    need([[1, 4, 9, 18, 24, 32, 60], 108] in conditional, 'attained conditional maximum')
    body = tuple(7200*i for i in range(1, 7))
    tails = (1156, 264, 1602, 9081, 9117, 16304, 32672)
    need(gcd(*body) == 7200 and gcd(*(body+tails)) == 1 and len(set(body+tails)) == 13, 'nonunit primitive full13 positive scope control')
    need(valid(tuple(sorted(gcd(7200, x) for x in tails))), 'full13 control projected profiles verified')
    tail_edges = []
    for i, j in combinations(range(7), 2):
        g = gcd(tails[i], tails[j]); p, q = sorted((tails[i]//g, tails[j]//g))
        if p+q <= 356 and atlas_sum(p+q):
            tail_edges.append((i, j))
    need(tail_edges == [(0, 1), (1, 2)], 'full13 positive control has disconnected seven-complement graph')
    safe_lifts = [j for j in range(7200) if all(14*min((v*(7*j+1)) % 50400, 50400-(v*(7*j+1)) % 50400) >= 50400 for v in body+tails)]
    need(len(safe_lifts) >= 114-excess(tuple(gcd(7200, x) for x in tails)), 'literal thirteen-speed lift preserves six-body phase and composite credit')

    certificate = dict(status='FINITE-EXACT complete stated universes; analytic promotion belongs to independent audit',
        inherited_sha256=sha256(raw).hexdigest(), profile_sha256=sha256(rawP).hexdigest(), fixed_tree_bank=fixed,
        actual_row=list(U), actual_margins=list(d), actual_graph=actual, alpha=[1, 2], paircounts=paircounts,
        danger_marginals=marginals, danger_histogram=histogram,
        composite=dict(n=3600, p=578, q=801, L=L, intervals=intervals, walls=walls, phase2_counts=rows,
                       minimum=lo, maximum=hi, minimizer=[low.numerator, low.denominator]),
        conditional_words=conditional, conditional_attempted=attempted, conditional_maximum=maximum,
        full13_control=dict(body=body, tails=tails, complement_edges=tail_edges, alpha=[1, 7], safe_lifts=safe_lifts),
        conditional_words_sha256=sha256(canonical(conditional)).hexdigest())
    data = canonical(certificate)+b'\n'
    (DEST/(HERE.stem+'_certificate.json')).write_bytes(data)
    print('GENERAL_IFF primitive projective normalization plus clipped prime potentials and distinctness; generic row/tree controls', generic_rows)
    print('FIXED_BANK 15 words;29 locally oriented first-attainer trees;20 distinct realizations on10 words;8 duplicate failures;1 valuation failure')
    print('ACTUAL_U', U, 'MARGINS', d, 'STRICT_GRAPH_EDGES', len(actual))
    print('COMMON_ALPHA 1/2; ATLAS_PAIR_COUNTS',(0,)*6, 'DANGER_MARGINALS', tuple(marginals), 'SAFE', histogram[0])
    print('ALL21_PAIR_COUNTS', paircounts)
    print('COMPOSITE ratio578/801; e2; intervals', len(intervals), 'walls', len(walls), 'chambers', len(walls), 'min', lo, 'max', hi, 'minimizer', str(low))
    print('CONDITIONAL_PROFILE 23751 complete candidates;1033 admissible; maxE108; uniform overlap114; strict integer margin6')
    print('FULL13_CONTROL nonunit primitive row; disconnected complement with2 edges; alpha1/7 safe lifts', len(safe_lifts))
    print('SCOPE named actual two-arm family closes; no all7200 closure, no alleged unsafe row')
    print('CERTIFICATE_SHA256', sha256(data).hexdigest())
    print('PASS', GATES, 'always-active exact gates; raw LF')


if __name__ == '__main__':
    main()

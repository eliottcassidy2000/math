"""Exact cycle ledger and native square-guard complement quotient.

No floating point, assertions, external graph atlas, or cached research outputs.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import importlib.util
import json
import random

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def parts(n, lo=2):
    yield ()
    for a in range(lo, n+1):
        for rest in parts(n-a, a):
            yield (a,)+rest


def add(*items):
    ans = defaultdict(F)
    for c, p in items:
        for key, value in p.items():
            ans[key] += c*value
    return {key: value for key, value in ans.items() if value}


def u(n):
    return {p: F((-1)**sum(a % 2 == 0 for a in p)) for p in parts(n)}


def v(n):
    return {p: F((-1)**len(p)) for p in parts(n) if sum(p) == n}


def selected(a, poly):
    # Sum over selected cycles. Repeated equal-sized cycles have that many
    # choices of distinguished occurrence; the M basis is unordered.
    return {tuple(sorted((a,)+p)): c*(p.count(a)+1) for p, c in poly.items()}


def ledger():
    n = 10
    H = {p: F(2**len(p)) for p in parts(n) if all(a % 2 for a in p)}
    P = {(2,): 1, (2, 2): -1, (2, 2, 2): 1, (2, 2, 2, 2): -1}
    tri = {(3,): 1, (3, 3): 3, (2, 2, 3): -1, (2, 2, 3, 3): -2}
    common = add((1, P), (1, tri), (1, {(2, 3): 1, (5,): 1,
        (2, 2, 5): -1, (2, 5): 1, (3, 5): 3,
        (7,): 1, (8,): 1, (2, 8): -2, (9,): 1}))
    receipts = {}
    for sign in (-1, 1):
        G = add((1, H), (-1, u(n)), (sign, v(n)))
        E4 = selected(4, add((1, u(6)), (-sign, v(6)),
                            (1, {(4,): F(1, 2)}),
                            (int(sign == -1), {(2, 4): -1})))
        E6 = selected(6, {(): 1, (2,): -1, (3,): 1,
                          (2, 2): 2 if sign == -1 else 0})
        rest = add((1, G), (-1, common), (-1, E4), (-1, E6))
        need(all(c >= 0 for c in rest.values()), (sign, rest))
        receipts[str(sign)] = {str(p): str(c) for p, c in sorted(rest.items())}
    A, B = F(9, 56), F(37, 64)
    need(B*B/(4*A) == F(9583, 18432), "seven-vertex completed square")
    need(1-B*B/(4*A) == F(8849, 18432), "triangle reserve")
    # Check the matching equality boundary independently of the graph census.
    from math import comb
    values = [25*r-16*comb(r, 2)+4*comb(r, 3)-comb(r, 4) for r in range(6)]
    need(values == [0, 25, 34, 31, 19, 0], "matching endpoints")
    m1, m2, m3, m4, m5, p21, wedge = 6, 14, 16, 9, 2, 4, 1
    hostile11 = 499*m1-320*m2+80*m3-20*m4+5*m5+80*p21
    need(hostile11 == -56 and F(-hostile11, 320*wedge) == F(7, 40),
         'order eleven bare graph obstruction')
    return {"all_cycle_collections": len(list(parts(n))),
            "nonnegative_unused_terms": receipts, "matching_values": values,
            "triangle_reserve": "8849/18432",
            "order11_bare_graph_witness": {"graph": "P3 disjoint union 4K2",
                "value": hostile11, "necessary_wedge_coefficient": "7/40"}}


def load_paths():
    source = HERE.parent / 'joint_incidence_clocked_frames_20261102/path_kernel.py'
    spec = importlib.util.spec_from_file_location('inherited_path_kernel', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def complement_carrier(n):
    source = load_paths()
    edges, paths = source.graph_paths(n, True)
    pats = source.patterns(edges, paths)
    full = (1 << len(edges))-1
    active = set()
    for pos, neg in pats:
        free = full ^ (pos | neg)
        sub = free
        while True:
            active.add(pos | sub)
            if not sub:
                break
            sub = (sub-1) & free
    states = sorted(active)
    index = {x: i for i, x in enumerate(states)}
    def mask(x):
        return sum(1 << j for j, (p, m) in enumerate(pats)
                   if x & p == p and not x & m)
    guards = [mask(x) for x in states]
    transitions = []
    for x, guard in zip(states, guards):
        row = []
        for j, (pos, neg) in enumerate(pats):
            if not guard >> j & 1:
                row.append(-1)
                continue
            y = x ^ (full ^ (pos | neg))
            need(y in index and mask(y) >> j & 1, "native complement legality")
            # Independent bit-by-bit native relation: path arcs kept, every
            # other guard arc reversed. Outside-guard arcs are irrelevant.
            need(all(((x ^ y) >> e & 1) == (not ((pos | neg) >> e & 1))
                     for e in range(len(edges))), "literal target relation")
            row.append(index[y])
        transitions.append(tuple(row))
    def partition(keys):
        labels = {x: i for i, x in enumerate(sorted(set(keys)))}
        return [labels[x] for x in keys]
    classes = partition(guards)
    rounds = [len(set(classes))]
    while True:
        refined = partition([(guards[i], tuple(-1 if j < 0 else classes[j]
                             for j in row)) for i, row in enumerate(transitions)])
        rounds.append(len(set(refined)))
        if len(set(refined)) == len(set(classes)):
            classes = refined
            break
        classes = refined
    group = defaultdict(list)
    for i, c in enumerate(classes):
        group[c].append(i)
    # Replay all transitions of every representative in every quotient fibre.
    for ids in group.values():
        signatures = {(guards[i], tuple(-1 if j < 0 else classes[j]
                      for j in transitions[i])) for i in ids}
        need(len(signatures) == 1, "closed quotient fibre")
    histogram = Counter(g.bit_count() for g in guards)
    histogram[0] = (1 << len(edges))-len(states)
    comp = [(i, j) for i, j in combinations(range(len(pats)), 2)
            if not (pats[i][0] & pats[j][1] or pats[i][1] & pats[j][0])]
    need(all(sum(i in e for e in comp) <= 1 for i in range(len(pats))),
         "matching compatibility hypothesis")
    sizes = Counter(map(len, group.values()))
    if n == 25:
        need(len(edges) == 32 and len(paths) == 20 and len(states) == 4864,
             "inherited square25 controls")
        need(rounds == [24, 32, 32] and sizes == {64: 12, 128: 8, 256: 12},
             "minimal actual-complement carrier")
        need(histogram == {0: 4294962432, 1: 4608, 2: 256}, "path spectrum")
    return {"n": n, "guard_edges": len(edges), "directed_paths": len(paths),
            "compatible_pairs": comp, "active_orientations": len(states),
            "refinement_class_counts": rounds,
            "active_class_size_histogram": dict(sorted(sizes.items())),
            "minimal_states_including_dead": len(group)+int(histogram[0] > 0),
            "path_count_histogram": dict(sorted(histogram.items())),
            "full_tournament_free_bits": n*(n-1)//2-len(edges)}


def fibonacci():
    # Path matchings, inherited separator warning: carry endpoint occupancy.
    polys = [[1], [1]]
    for n in range(2, 21):
        p = polys[-1][:]+[0]
        for k, c in enumerate(polys[-2]):
            p[k+1] += c
        while p[-1] == 0:
            p.pop()
        polys.append(p)
    values = []
    for n, p in enumerate(polys):
        val = sum(c*F(-1, 4)**k for k, c in enumerate(p))
        need(val == F(n+1, 2**n), "double-root transfer")
        values.append(sum(p))
    return {"matching_count_P0_through_P20": values,
            "weighted_identity": "Z(P_n,-1/4)=(n+1)/2^n"}


def determinant(matrix):
    a = [row[:] for row in matrix]
    n, old, sign = len(a), 1, 1
    for k in range(n-1):
        if not a[k][k]:
            row = next((i for i in range(k+1, n) if a[i][k]), None)
            if row is None:
                return 0
            a[k], a[row] = a[row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                top = a[i][j]*pivot-a[i][k]*a[k][j]
                need(top % old == 0, 'Bareiss exact division')
                a[i][j] = top//old
            a[i][k] = 0
        old = pivot
    return sign*a[-1][-1]


def hamilton(a):
    n = len(a)
    dp = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for s in range(1, 1 << n):
        for i in range(n):
            if not s >> i & 1:
                continue
            t = s ^ (1 << i)
            if t:
                dp[s][i] = sum(dp[t][j]*a[j][i] for j in range(n) if t >> j & 1)
    return sum(dp[-1])


def cycles(a):
    n = len(a)
    counts = [0]*(n+1)
    for root in range(n):
        m = n-root
        dp = [[0]*m for _ in range(1 << m)]
        dp[1][0] = 1
        for s in range(1, 1 << m, 2):
            for i in range(1, m):
                if not s >> i & 1:
                    continue
                t = s ^ (1 << i)
                dp[s][i] = sum(dp[t][j]*a[root+j][root+i]
                               for j in range(m) if t >> j & 1)
                counts[s.bit_count()] += dp[s][i]*a[root+i][root]
    return counts


def numerical_controls():
    rng = random.Random(10258127)
    n = 10
    matrices = []
    for _ in range(50):
        a = [[0]*n for _ in range(n)]
        for i, j in combinations(range(n), 2):
            a[i][j] = rng.randrange(5)
            a[j][i] = 4-a[i][j]
        matrices.append(a)
    for pairs in range(6):
        a = [[4*int(i < j) for j in range(n)] for i in range(n)]
        for i in range(pairs):
            a[2*i][2*i+1] = a[2*i+1][2*i] = 2
        matrices.append(a)
    margins = []
    for a in matrices:
        H, det = hamilton(a), determinant(a)
        u = determinant([[a[i][j]+4*int(i == j) for j in range(n)] for i in range(n)])
        delta = F(4*H-u-abs(det), 4**n)
        variance = F(sum(a[i][j]*a[j][i] for i, j in combinations(range(n), 2)), 16)
        tau = F(sum(a[i][j]*a[j][k]*a[k][i]+a[i][k]*a[k][j]*a[j][i]
                    for i, j, k in combinations(range(n), 3)), 64)
        margin = delta-F(39, 64)*variance-F(8849, 18432)*tau
        cyc = cycles(a)
        need(F(cyc[3], 64) == tau, 'cycle DP versus triple enumeration')
        margin -= sum(F(cyc[j], 4**j)*weight for j, weight in
                      ((4, F(1, 32)), (5, 1), (6, F(1, 8)),
                       (7, 1), (8, F(1, 2)), (9, 1)))
        need(margin >= 0, 'direct determinant/Hamilton control')
        margins.append(margin)
    need(margins[-1] == 0 and margins[-6] == 0 and all(x > 0 for x in margins[-5:-1]),
         'fair-pair equality controls')
    # Scope hostile: same restricted carrier, different full Hamilton values.
    source = load_paths()
    edges, paths = source.graph_paths(15, True)
    p = paths[0]
    rank = {x: i for i, x in enumerate(p)}
    used = {tuple(sorted(e)) for e in zip(p, p[1:])}
    e = next(e for e in edges if e not in used)
    a = [[int(rank[i] < rank[j]) for j in range(15)] for i in range(15)]
    b = [row[:] for row in a]
    i, j = e
    b[i][j], b[j][i] = b[j][i], b[i][j]
    guards_a = [all(a[u][v] for u, v in zip(path, path[1:])) for path in paths]
    guards_b = [all(b[u][v] for u, v in zip(path, path[1:])) for path in paths]
    need(guards_a == guards_b and sum(guards_a) == 1, 'same restricted guard state')
    old, new = hamilton(a), hamilton(b)
    distance = abs(rank[i]-rank[j])
    need(old == 1 and new == 1+2**(distance-1) and new > old,
         'unrestricted Hamilton observer hostile')
    return {'exact_order10_samples': len(matrices),
            'square15_scope_hostile': {'extra_guard_edge': [i+1, j+1],
                'path_distance': distance, 'full_H': [old, new],
                'same_guarded_carrier': True}}


if __name__ == '__main__':
    receipt = {"ledger": ledger(), "complement_carriers": [
        complement_carrier(n) for n in (15, 16, 17, 23, 25)],
        "fibonacci": fibonacci(), "independent_controls": numerical_controls()}
    (HERE/'audit.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf8')
    print(json.dumps(receipt))

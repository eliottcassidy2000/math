#!/usr/bin/env python3
"""S38 audit (independent referee), exploratory: the saturated relation lattice of {zeta55^k, zeta55^k v_1,
zeta55^k v_3, zeta55^k v_4} contains a relation of L1 length 20 (audit_short_coupling_relations.py), i.e. two
10-step unit walks with a common endpoint built only from rotated FIRST-POWER conjugates of v (compare
clock_duality_oct15's C2: two 32-step walks in the powers v^i xi^j).  This script builds walk banks from that
20-step circuit, in the manner of C2's circuit assemblies, and tests them:

  bank(K): the partial sums of K cyclic orderings of the 20 signed steps (seeded RNG), all unit pairs found exactly
  (three-prime modular filter on every pair + exact norm), non-generator edges counted, triangles counted,
  three-colourability by SAT (colouring re-verified on every edge).
  bank(K) + translates: the bank translated by every point of the bank's first ordering (Minkowski-type closure).

Exploratory, FINITE-EXACT for the graphs as built; it decides nothing about larger assemblies.
"""
import itertools, json, math, os, random, sys, time
import numpy as np
from pysat.solvers import Solver
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import audit_cubic_probes_independent as A

HERE = os.path.dirname(os.path.abspath(__file__))
OUTJ = os.path.join(HERE, 'audit_short_circuit_walk_banks.json')


def log(*a):
    print(*a, flush=True)


rel = json.load(open(os.path.join(HERE, 'audit_short_coupling_relations.json')))['triples'][0]
steps_spec = rel['saturated_relation_steps']
vq, N, Nb = A.build_v()
ms = [m for m in range(1, 55) if math.gcd(m, 55) == 1 and m % 5 == 1]
vconj = {m % 11: A.sigma_q(vq, m) for m in ms}
DEN = 3082
one = A.rou_q(0)
steps = []
for w, k, cf in steps_spec:
    wel = one if w == '1' else vconj[int(w[1:])]
    vec = A.q_to_intvec(A.qmul(wel, A.rou_q(k)), DEN)
    for _ in range(abs(cf)):
        steps.append(vec * (1 if cf > 0 else -1))
assert len(steps) == 20 and not np.any(sum(steps))
gens = {tuple(s.tolist()) for s in steps} | {tuple((-s).tolist()) for s in steps}
primes = A.find_primes(3, avoid=[DEN])


def analyse(points, label):
    P = np.array(points)
    n = len(P)
    cand = A.modular_unit_pairs(P, DEN, primes)
    E = [(a, b) for a, b in cand if A.exact_unit(P[a] - P[b], DEN)]
    nongen = sum(1 for a, b in E if tuple((P[b] - P[a]).tolist()) not in gens)
    adj = [set() for _ in range(n)]
    for a, b in E:
        adj[a].add(b); adj[b].add(a)
    tri = sum(len(adj[a] & adj[b]) for a, b in E) // 3
    var = lambda v, c: 3 * v + c + 1
    with Solver(name='cadical195') as S:
        for v in range(n):
            S.add_clause([var(v, c) for c in range(3)])
        for a, b in E:
            for c in range(3):
                S.add_clause([-var(a, c), -var(b, c)])
        ok = S.solve()
        col = None
        if ok:
            m = set(l for l in S.get_model() if l > 0)
            col = [next(c for c in range(3) if var(v, c) in m) for v in range(n)]
            assert all(col[a] != col[b] for a, b in E)
    rec = dict(assembly=label, vertices=n, unit_edges=len(E), non_generator_edges=nongen, triangles=tri,
               three_colourable=bool(ok), max_degree=max((len(s) for s in adj), default=0))
    log(json.dumps(rec))
    return rec


rng = random.Random(20261016)
out = dict(relation=steps_spec, banks=[])
for K in (1, 4, 16):
    pts, seen = [], set()
    orders = []
    for _ in range(K):
        order = list(range(20))
        rng.shuffle(order)
        orders.append(order)
        cur = np.zeros(120, dtype=np.int64)
        for i in order:
            t = tuple(cur.tolist())
            if t not in seen:
                seen.add(t); pts.append(cur.copy())
            cur = cur + steps[i]
        assert not np.any(cur)
    out['banks'].append(analyse(pts, f'bank of {K} cyclic orderings'))
    if K == 4:
        base = [p for p in pts]
        first = []
        cur = np.zeros(120, dtype=np.int64)
        for i in orders[0]:
            first.append(cur.copy()); cur = cur + steps[i]
        pts2, seen2 = [], set()
        for t0 in first:
            for p in base:
                q = p + t0
                t = tuple(q.tolist())
                if t not in seen2:
                    seen2.add(t); pts2.append(q)
        out['banks'].append(analyse(pts2, 'bank of 4 orderings translated along the first ordering'))
json.dump(out, open(OUTJ, 'w'), indent=1, default=str)
log('wrote', OUTJ)

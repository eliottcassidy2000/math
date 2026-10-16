#!/usr/bin/env python3
"""S38 audit (independent referee): re-derive, from the House of Graphs adjacency lists stored in
ei19_fish_placement_q_zeta55_20261012.py, the two graph facts that THM-4661 statement 4 inherits:

(a) EI19: in every faithful (injective, unit-edge) plane realization, every unit 4-cycle is a parallelogram, and the
    resulting direction classes satisfy 2 u(2->6) + 3 u(6->13) + 3 u(18->2) = 0 and u5 = 3(u4 + u9) (support three),
    forcing the ratios (-1 +- 2 sqrt(-2))/3 and (1 +- sqrt(-35))/6.  Checked: the displayed relation lies in the
    Q-span of the cycle-space relations among the direction classes (my own union-find and elimination).
(b) Fish: the five 4-cycles (0,10,8,20), (2,8,10,19), (1,11,9,21), (3,9,11,18), (18,19,21,20) exist, 22 is adjacent
    to 0..3, 0 ~ 20, 1 ~ 21, so v0 + v2 = v19 + v20 = v18 + v21 = v1 + v3 and the derivation of equations (1) of the
    fish note applies (algebra checked in audit_field_and_counting_checks.py).
"""
import itertools, json, os, re
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
OUTJ = os.path.join(HERE, 'audit_ei19_fish_structure.json')
src = open(os.path.join(ROOT, '04-computation/experiments/ei19_fish_placement_q_zeta55_20261012.py')).read()
block = src[src.index('GRAPHS = {'):src.index('}', src.index('GRAPHS = {')) + 1]
GRAPHS = eval(block[len('GRAPHS = '):])
out = {}


def log(*a):
    print(*a, flush=True)


def edges_of(ADJ):
    E = sorted({(min(i, j), max(i, j)) for i, nb in enumerate(ADJ) for j in nb})
    for i, nb in enumerate(ADJ):
        for j in nb:
            assert i in ADJ[j], 'adjacency not symmetric'
    return E


# ---------------------------------------------------------------- (a) EI19
ADJ = GRAPHS['EI19']
E = edges_of(ADJ)
n = len(ADJ)
adj = {v: set(a) for v, a in enumerate(ADJ)}
assert n == 19 and len(E) == 35
# all 4-cycles (a,b,c,d), as cyclic sequences up to rotation/reflection
cyc = set()
for a in range(n):
    for b in adj[a]:
        for c in adj[b] - {a}:
            for d in (adj[c] & adj[a]) - {b}:
                if len({a, b, c, d}) == 4:
                    seq = (a, b, c, d)
                    rots = [seq[i:] + seq[:i] for i in range(4)]
                    rots += [tuple(reversed(r)) for r in rots]
                    cyc.add(min(rots))
cyc = sorted(cyc)
# union-find over oriented edges: variable for edge (i<j) is v_j - v_i; parallelogram (a,b,c,d): v_b - v_a = v_c - v_d,
# v_c - v_b = v_d - v_a
eidx = {e: k for k, e in enumerate(E)}
parent = list(range(len(E)))
sign = [1] * len(E)         # u_e = sign[e] * class_vector(root)


def find(x):
    if parent[x] == x:
        return x, 1
    r, s = find(parent[x])
    parent[x] = r
    sign[x] = sign[x] * s
    return r, sign[x]


def oriented(i, j):          # returns (edge index, +1 if v_j - v_i equals the edge variable)
    return (eidx[(i, j)], 1) if i < j else (eidx[(j, i)], -1)


def union(i1, j1, i2, j2):    # v_j1 - v_i1 = v_j2 - v_i2
    e1, s1 = oriented(i1, j1)
    e2, s2 = oriented(i2, j2)
    r1, t1 = find(e1)
    r2, t2 = find(e2)
    # u_e1 = t1 * U_r1, u_e2 = t2 * U_r2 ; s1 * u_e1 = s2 * u_e2  ->  U_r1 = (s2 t2)/(s1 t1) U_r2
    if r1 == r2:
        assert s1 * t1 == s2 * t2, 'inconsistent parallelogram orientation'
        return
    parent[r1] = r2
    sign[r1] = (s2 * t2) * (s1 * t1)


for (a, b, c, d) in cyc:
    union(a, b, d, c)
    union(b, c, a, d)
roots = sorted({find(k)[0] for k in range(len(E))})
cls = {r: i for i, r in enumerate(roots)}
m = len(roots)


def evec(i, j):
    e, s = oriented(i, j)
    r, t = find(e)
    v = [0] * m
    v[cls[r]] = s * t
    return v


# cycle space: fundamental cycles of a BFS tree
import collections
par = {0: None}
dq = collections.deque([0])
while dq:
    x = dq.popleft()
    for y in sorted(adj[x]):
        if y not in par:
            par[y] = x
            dq.append(y)


def path_to_root(x):
    p = [x]
    while par[p[-1]] is not None:
        p.append(par[p[-1]])
    return p


rels = []
tree = {(min(x, par[x]), max(x, par[x])) for x in par if par[x] is not None}
for (i, j) in E:
    if (i, j) in tree:
        continue
    # cycle: root -> i, edge i->j, j -> root
    pi, pj = path_to_root(i), path_to_root(j)
    vec = [0] * m
    # walk from j up to root and from root down to i, then i->j
    for a, b in zip(pi[:-1], pi[1:]):      # i -> ... -> root : add (b - a) steps reversed later
        w = evec(a, b)
        vec = [x + y for x, y in zip(vec, w)]  # sum of (v_b - v_a) along i -> root
    for a, b in zip(pj[:-1], pj[1:]):
        w = evec(a, b)
        vec = [x - y for x, y in zip(vec, w)]  # minus sum along j -> root
    w = evec(i, j)
    vec = [x + y for x, y in zip(vec, w)]      # (v_root - v_i) - (v_root - v_j) + (v_j - v_i) = ... consistent cycle
    # the above equals 2(v_j - v_i)?  recompute cleanly: closed walk i -> j -> ... -> root -> ... -> i
    vec = [0] * m
    walk = [i, j] + pj[1:] + list(reversed(pi[:-1]))
    for a, b in zip(walk[:-1], walk[1:]):
        w = evec(a, b)
        vec = [x + y for x, y in zip(vec, w)]
    rels.append(vec)
Rm = sp.Matrix(rels)
rank = Rm.rank()
target1 = [2 * x + 3 * y + 3 * z for x, y, z in zip(evec(2, 6), evec(6, 13), evec(18, 2))]
in_span1 = sp.Matrix(rels + [target1]).rank() == rank
log(f'EI19: {len(cyc)} unit 4-cycles, {m} direction classes, cycle-space relation rank {rank}; '
    f'2 u(2->6) + 3 u(6->13) + 3 u(18->2) = 0 in span: {in_span1}')
# second relation: classes of (2,3) =: u4 and (3,16) =: u9 in the S24 numbering; test u = 3(u_a + u_b) for the third
# class among all class pairs: search support-three relations n_a u_a + n_b u_b + n_c u_c = 0 in the span
sup3 = []
for a_, b_, c_ in itertools.combinations(range(m), 3):
    # does the span contain a vector supported on {a,b,c} with all three nonzero?
    sub = sp.Matrix(rels)
    # nullspace approach: vectors x with x^T rels-span ... simpler: project span onto coordinates outside {a,b,c}
    others = [k for k in range(m) if k not in (a_, b_, c_)]
    # find combinations of relations vanishing on 'others'
    A_ = sp.Matrix([[r[k] for k in others] for r in rels]).T
    ns = A_.nullspace()
    for v in ns:
        comb = (sp.Matrix(rels).T * v)
        vals = [comb[k] for k in (a_, b_, c_)]
        if all(x != 0 for x in vals):
            g = sp.gcd_list([int(sp.nsimplify(x * sp.lcm_list([sp.fraction(sp.nsimplify(y))[1] for y in vals]))) for x in vals])
            l = sp.lcm_list([sp.fraction(sp.nsimplify(y))[1] for y in vals])
            iv = [int(sp.nsimplify(x * l / g)) for x in vals]
            sup3.append(((a_, b_, c_), iv))
            break
log('EI19 support-three relations (class indices, coefficients):', sup3)
out['EI19'] = dict(four_cycles=len(cyc), classes=m, relation_rank=int(rank), relation_2_3_3_in_span=in_span1,
                   support_three=[[list(k), v] for k, v in sup3])
assert in_span1
# the forced ratios: |n_b u_b + n_c u_c| = |n_a| with unit u's: ratio u_c/u_b has real part (n_a^2 - n_b^2 - n_c^2)/(2 n_b n_c)
for (k, (na, nb, nc)) in sup3:
    for (x, y, z) in ((na, nb, nc), (nb, na, nc), (nc, na, nb)):
        re_ = sp.Rational(x * x - y * y - z * z, 2 * y * z)
        disc = sp.simplify(re_ ** 2 - 1)
        log(f'  ratio of the classes with coefficients {y},{z}: real part {re_}, imaginary part sqrt({disc}) '
            f'-> field Q(sqrt({sp.factorint(sp.Rational(disc).p * sp.Rational(disc).q)}))')
        break

# ---------------------------------------------------------------- (b) fish
ADJ = GRAPHS['FISH']
Ef = set(edges_of(ADJ))
has = lambda i, j: (min(i, j), max(i, j)) in Ef
cycles = [(0, 10, 8, 20), (2, 8, 10, 19), (1, 11, 9, 21), (3, 9, 11, 18), (18, 19, 21, 20)]
ok_cycles = all(has(c[i], c[(i + 1) % 4]) for c in cycles for i in range(4))
ok_hub = all(has(22, i) for i in range(4))
ok_other = has(0, 20) and has(1, 21) and has(10, 0) and has(10, 19)
log(f'fish: 23 vertices {len(ADJ) == 23}, 42 edges {len(Ef) == 42}; the five 4-cycles present {ok_cycles}; '
    f'hub 22 ~ 0..3 {ok_hub}; 0~20, 1~21, 10~0, 10~19 {ok_other}')
out['fish'] = dict(cycles_present=ok_cycles, hub=ok_hub, other_edges=ok_other, vertices=len(ADJ), edges=len(Ef))
assert ok_cycles and ok_hub and ok_other
json.dump(out, open(OUTJ, 'w'), indent=1, default=str)
log('wrote', OUTJ)

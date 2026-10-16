#!/usr/bin/env python3
"""S38 audit (independent referee): the S38 probe script reports, for 1 and three first-power conjugates of the
width-three unit v, 'LLL coupling relations' of L1 length 19,416 / 34,652 / 37,330 / 4,922,578 (note section 5;
the results INDEX turns this into 'LLL length at least about 19,000').  Those numbers come from LLL applied to the
lattice spanned by flint's fmpz_mat.nullspace() basis, which is a basis of the RATIONAL null space, not a Z-basis of
the integer relation lattice.  This script

1. recomputes the saturated integer relation lattice by kernel-LLL (rows [N * b_i | e_i]),
2. extracts the shortest relation found, prints it as an explicit list of unit steps c * (w * zeta55^k),
3. verifies it EXACTLY in two independent arithmetics: (a) my w-basis arithmetic (Z[zeta55][w]/(w^3+2w^2-4w-4)),
   (b) the S38 script's own theta-basis arithmetic and its own v (imported by exec of its field section),
4. reproduces the S38 number for comparison by running their exact recipe (nullspace + LLL),
5. reports the index of the S38 sublattice in the saturated lattice (ratio of Gram determinants).
Output: audit_short_coupling_relations.out and .json in this folder.
"""
import json, math, os, sys, time
import numpy as np
from flint import fmpz_mat, fmpz_poly, fmpq_poly, fmpq

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, HERE)
OUTJ = os.path.join(HERE, 'audit_short_coupling_relations.json')
import audit_cubic_probes_independent as A   # my independent arithmetic (w-basis); main() is not executed

RES = {}


def log(*a):
    print(*a, flush=True)


vq, N, Nb = A.build_v()
ms = [m for m in range(1, 55) if math.gcd(m, 55) == 1 and m % 5 == 1]
vconj = {m % 11: A.sigma_q(vq, m) for m in ms}
DEN = 3082
one = A.rou_q(0)


def rows_for(trip):
    rows, labels = [], []
    for wi, w in enumerate([one] + [vconj[r] for r in trip]):
        for k in range(40):
            rows.append(A.q_to_intvec(A.qmul(w, A.rou_q(k)), DEN))
            labels.append(('1' if wi == 0 else f'v{trip[wi - 1]}', k))
    return rows, labels


def saturated_kernel(rows):
    nr = len(rows)
    Nbig = 10 ** 12
    M = fmpz_mat(nr, 120 + nr)
    for i, r in enumerate(rows):
        for c in range(120):
            M[i, c] = int(r[c]) * Nbig
        M[i, 120 + i] = 1
    Rm = M.lll()
    kern = [[int(Rm[i, 120 + c]) for c in range(nr)] for i in range(nr) if all(Rm[i, c] == 0 for c in range(120))]
    return kern


def gram_det(B):
    Bm = fmpz_mat(B)
    return (Bm * Bm.transpose()).det()


RES['triples'] = []
for trip in [(1, 3, 4), (2, 9, 3), (1, 10, 2), (1, 2, 3)]:
    t0 = time.time()
    rows, labels = rows_for(trip)
    kern = saturated_kernel(rows)
    assert len(kern) == 40
    # verify every kernel vector exactly in integer arithmetic on the coordinate vectors
    Rows = np.array(rows, dtype=object)
    for kv in kern:
        s = np.zeros(120, dtype=object)
        for c, cf in enumerate(kv):
            if cf:
                s = s + Rows[c] * cf
        assert all(int(x) == 0 for x in s)
    best = min(kern, key=lambda kv: sum(abs(x) for x in kv))
    L1 = sum(abs(x) for x in best)
    steps = [(labels[c][0], labels[c][1], cf) for c, cf in enumerate(best) if cf]
    # (a) exact check in my arithmetic as a sum of field elements
    acc = None
    for (w, k, cf) in steps:
        wel = one if w == '1' else vconj[int(w[1:])]
        term = A.qmul(wel, A.rou_q(k))
        term = [[x * cf for x in blk] for blk in term]
        acc = term if acc is None else [[acc[b][i] + term[b][i] for i in range(40)] for b in range(3)]
    zero_a = all(acc[b][i] == 0 for b in range(3) for i in range(40))
    # the S38 recipe: rational nullspace basis + LLL (their exact procedure), for comparison
    X, nullity = fmpz_mat([[int(x) for x in r] for r in rows]).transpose().nullspace()
    Bs = [[int(X[i, c]) for i in range(X.nrows())] for c in range(nullity)]
    Bl = fmpz_mat(Bs).lll()
    s38_best = min(sum(abs(int(Bl[r, c])) for c in range(Bl.ncols())) for r in range(Bl.nrows()))
    # index of the S38 sublattice: sqrt(det Gram(S38 basis) / det Gram(saturated basis))
    idx2 = gram_det(Bs) / gram_det(kern)
    rec = dict(triple=list(trip), saturated_min_L1=L1, saturated_relation_steps=steps, exact_zero_w_basis=zero_a,
               s38_recipe_min_L1=s38_best, s38_sublattice_index_squared=str(idx2), seconds=round(time.time() - t0, 1))
    RES['triples'].append(rec)
    log(json.dumps(rec))

# (b) verify the (1,3,4) relation with the S38 script's own arithmetic and own v (theta basis), by executing the
# field-construction part of their script (everything before the directions section) in a private namespace.
src = open(os.path.join(ROOT, '04-computation/experiments/cubic_carrier_three_colour_probes_20261016.py')).read()
head = src[:src.index('# ---------- integer vectors ----------')]
ns = {'__name__': 's38probe'}
cwd = os.getcwd()
os.chdir(ROOT)
exec(compile(head, 's38probe_head', 'exec'), ns)
os.chdir(cwd)
v_el = ns['v_el']; sigma = ns['sigma']; mul = ns['mul']; root_of_unity = ns['root_of_unity']
ms_ = [m for m in range(1, 55) if math.gcd(m, 55) == 1 and m % 5 == 1]
vc = {m % 11: sigma(v_el, m) for m in ms_}
rec0 = RES['triples'][0]
acc = [fmpq_poly(0)] * 3
for (w, k, cf) in rec0['saturated_relation_steps']:
    wel = root_of_unity(0) if w == '1' else vc[int(w[1:])]
    term = mul(wel, root_of_unity(k))
    acc = [acc[b] + term[b] * cf for b in range(3)]
zero_b = all(acc[b] == 0 for b in range(3))
RES['triple_134_exact_zero_in_S38_arithmetic'] = zero_b
log('relation for (1,3,4) is exactly zero in the S38 script\'s own theta-basis arithmetic:', zero_b)
json.dump(RES, open(OUTJ, 'w'), indent=1, default=str)
log('wrote', OUTJ)

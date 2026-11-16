#!/usr/bin/env python3
"""Wave 39: exact verifier for the syllable-adaptive cone certificate of literal span-3 ribbon words.

Claim verified: every bottleneck-free literal span-3 word (a first vertex followed by letters, a letter being the set
of gaps 1, 2, 3 at which the new vertex beats an older one) with n vertices has
        H >= (p/q)^(n - 1 - Linit),
where p/q, Linit, the automaton and the functionals are read from the certificate file.

Ingredients (all rebuilt here from the fragment rule of ribbon.py; nothing is trusted from the search):
  * the 119-state nonnegative profile transfer A_l (THM-4880), the start vector v0 after the first vertex, readout h;
  * the reachable span (basis Rb) and the observable rows O = (h A_w) with O Rb of full row rank 61: the minimal
    representation after the first vertex.  "Bottleneck-free" means: every factor has rank >= 2 there, i.e.
    rank(O A_w Rb) >= 2.  All entries of O A_w Rb are nonnegative integers, so rank <= 1 iff all nonzero columns share
    one primitive direction (exact test);
  * the automaton: state = longest suffix of the word in PI, PI = {strings of length <= L0} U {u_1...u_j p: u_i floor
    units from the listed set, p a proper unit prefix, total length <= LMAX}.  PI is closed under deleting the last
    letter, so the state after appending l is the longest PI-suffix of (state + l).  The verifier recomputes every
    transition and every rank decision (an edge is absent exactly when state + l has rank <= 1, a bottleneck factor);
  * cones K_s = {O A_s X : X >= 0 supported on the raw support R}, and integer functionals Phi_s (denominator D).

Checks, all in exact integer arithmetic:
  (growth)  q * Phi_t . (O A_l A_s e_j) >= p * Phi_s . (O A_s e_j)  for every live edge s -l-> t and j in R,
  (readout) Phi_s . (O A_s e_j) <= D * h . (A_s e_j)                  for every live s and j in R,
  (initial) Phi_s . (O A_w v0) >= D                                  for every bottleneck-free word w of length Linit,
  (closure) every state reachable from the initial states is live.
Then for a bottleneck-free word with n vertices, Phi grows by p/q per letter after Linit letters and H >= Phi/D.
Usage: python3 w39_verify_literal_cert.py CERT.pkl
"""
import itertools
import os
import pickle
import sys
import time
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RIB = os.path.join(HERE, '..', '..', '..', '..', '04-computation', 'experiments', 'ribbon_spectrum_20261104')
sys.path.insert(0, RIB)
from ribbon import transfer_system  # noqa: E402

UNITS = ['0.2', '1.2', '0.123', '2.13', '0.1.3', '0.3.1', '0.13.2', '0.1.123', '0.3.12', '0.12.3', '1.2.3',
         '1.12.3', '0.1.23']


def tob(s):
    return sum(1 << (int(c) - 1) for c in s if c in '123')


UL = [tuple(tob(t) for t in u.split('.')) for u in UNITS]

t_start = time.time()
cert = pickle.load(open(sys.argv[1], 'rb'))
p, q, D, Linit = cert['p'], cert['q'], cert['D'], cert['Linit']
L0, LMAX, unit_ids = cert['L0'], cert['LMAX'], cert['unit_ids']
order, edges, Phi = cert['order'], cert['edges'], cert['Phi']      # Phi: dict state index -> list of 61 ints

# ---- the transfer ----
order_raw, states_raw, mats, final = transfer_system(3)
N = len(order_raw)
A = []
for b in range(8):
    M = np.zeros((N, N), dtype=np.int64)
    for (j, i), m in mats[frozenset(g for g in (1, 2, 3) if (b >> (g - 1)) & 1)].items():
        M[j, i] += m
    A.append(M)
e0 = np.zeros(N, dtype=np.int64)
e0[0] = 1
v0 = A[0] @ e0
h = np.zeros(N, dtype=np.int64)
for j in final:
    h[j] = 1


def rank_mod(rows, pmod=2147483629):
    M = [[int(x) % pmod for x in r] for r in rows]
    rk = 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c]), None)
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][c], pmod - 2, pmod)
        M[rk] = [x * inv % pmod for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c]:
                f = M[i][c]
                M[i] = [(a - f * b) % pmod for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk


# reachable span basis and observable rows (greedy, independence mod a prime; the final rank claims are exact below)
B = [v0.copy()]
fr = [v0.copy()]
while fr:
    nf = []
    for v in fr:
        for b in range(8):
            w = A[b] @ v
            if w.any() and rank_mod([list(x) for x in B] + [list(w)]) > len(B):
                B.append(w)
                nf.append(w)
    fr = nf
Rb = np.array(B, dtype=np.int64).T
Orows = [h.copy()]
fr = [h.copy()]
while fr:
    nf = []
    for u in fr:
        for b in range(8):
            w = u @ A[b]
            if rank_mod([list(o @ Rb) for o in Orows] + [list(w @ Rb)]) > len(Orows):
                Orows.append(w)
                nf.append(w)
    fr = nf
O = np.array(Orows, dtype=np.int64)
OR = O @ Rb
assert OR.min() >= 0
# exact full row rank of O Rb: Fraction elimination on an explicit set of columns
r = O.shape[0]
print('transfer: %d raw states, reachable span %d, observable rows %d' % (N, Rb.shape[1], r))
assert r == 61, r
assert all(len(Phi[a]) == r for a in Phi)
Mx = [[Fraction(int(x)) for x in row] for row in OR.tolist()]
rk = 0
for c in range(len(Mx[0])):
    piv = next((i for i in range(rk, r) if Mx[i][c] != 0), None)
    if piv is None:
        continue
    Mx[rk], Mx[piv] = Mx[piv], Mx[rk]
    for i in range(r):
        if i != rk and Mx[i][c] != 0:
            f = Mx[i][c] / Mx[rk][c]
            Mx[i] = [a - f * b for a, b in zip(Mx[i], Mx[rk])]
    rk += 1
assert rk == 61, rk
print('exact: O Rb has full row rank 61, so rank(O A_w Rb) is the rank of the minimal representation')

# ---- PI and the automaton ----
units = [UL[i] for i in unit_ids]
PI = set()
for L in range(L0 + 1):
    for t in itertools.product(range(8), repeat=L):
        PI.add(t)
prefixes = {u[:k] for u in units for k in range(len(u))}
fr = [()]
seen = {()}
while fr:
    nf = []
    for base in fr:
        for pr in prefixes:
            if len(base) + len(pr) <= LMAX:
                PI.add(base + pr)
        for u in units:
            nb = base + u
            if len(nb) <= LMAX and nb not in seen:
                seen.add(nb)
                nf.append(nb)
                PI.add(nb)
    fr = nf
for s in PI:
    assert s[:-1] in PI or not s     # closed under deleting the last letter
sidx = {s: i for i, s in enumerate(order)}
outd = {}
for (a, l, b) in edges:
    assert (a, l) not in outd
    outd[(a, l)] = b

lc = {}


def OAw(s):
    if s in lc:
        return lc[s]
    M = O.copy() if not s else OAw(s[1:]) @ A[s[0]]
    lc[s] = M
    return M


def rank_le1(s):
    M = OAw(s) @ Rb
    assert M.max() < 2 ** 61
    cols = M[:, M.any(axis=0)]
    if cols.shape[1] <= 1:
        return True
    g = np.gcd.reduce(cols, axis=0)
    P = cols // g
    return bool((P == P[:, :1]).all())


# initial states and live closure
R = set(np.nonzero(v0)[0])
fr = list(R)
while fr:
    nf = []
    for i in fr:
        for b in range(8):
            for j in np.nonzero(A[b][:, i])[0]:
                if j not in R:
                    R.add(j)
                    nf.append(j)
    fr = nf
R = np.array(sorted(R))
init_ok = 0
live = set()
for w in itertools.product(range(8), repeat=Linit):
    s = ()
    ok = True
    for k, l in enumerate(w):
        if rank_le1(s + (l,)):
            ok = False
            break
        s = next((s + (l,))[m:] for m in range(len(s) + 2) if (s + (l,))[m:] in PI)
    if not ok:
        continue
    a = sidx[s]
    live.add(a)
    X = v0.copy()
    for l in w:
        X = A[l] @ X
    y = [int(x) for x in O @ X]
    val = sum(int(c) * v for c, v in zip(Phi[a], y))
    assert val >= D, ('initial', w, val, D)
    init_ok += 1
fr = list(live)
while fr:
    nf = []
    for a in fr:
        for l in range(8):
            if (a, l) in outd and outd[(a, l)] not in live:
                live.add(outd[(a, l)])
                nf.append(outd[(a, l)])
    fr = nf
assert all(a in Phi for a in live), 'a live state has no functional'
print('initial words checked: %d (Linit = %d); live states %d' % (init_ok, Linit, len(live)))

# ---- edges: rank decisions, targets, growth and readout inequalities ----
n_edge = n_gen = n_removed = 0
for a in sorted(live):
    s = order[a]
    Ps = np.identity(N, dtype=np.int64)
    for l in s:
        Ps = A[l] @ Ps
    G = Ps[:, R]
    Ys = O @ G
    assert np.abs(Ys).max() < 2 ** 62
    phis = np.array([int(c) for c in Phi[a]], dtype=object)
    vs = phis.dot(Ys.astype(object))                          # Phi_s(O A_s e_j)
    Hs = (h @ G).astype(object)
    assert all(vs[j] <= D * Hs[j] for j in range(G.shape[1])), ('readout', s)
    for l in range(8):
        w = s + (l,)
        if rank_le1(w):
            assert (a, l) not in outd, ('edge present on a bottleneck factor', w)
            n_removed += 1
            continue
        t = next(w[m:] for m in range(len(w) + 1) if w[m:] in PI)
        assert (a, l) in outd and order[outd[(a, l)]] == t, ('transition', w, t)
        b = outd[(a, l)]
        Y2 = O @ (A[l] @ G)
        assert np.abs(Y2).max() < 2 ** 62
        phit = np.array([int(c) for c in Phi[b]], dtype=object)
        vt = phit.dot(Y2.astype(object))
        for j in range(G.shape[1]):
            if q * vt[j] < p * vs[j]:
                raise AssertionError(('growth', s, l, j, vt[j], vs[j]))
        n_edge += 1
        n_gen += G.shape[1]
print('edges verified %d (%d growth inequalities), bottleneck transitions removed %d' % (n_edge, n_gen, n_removed))
print('VERIFIED: every bottleneck-free literal span-3 word with n vertices has H >= (%d/%d)^(n-1-%d) = %.6f^(n-%d)'
      % (p, q, Linit, p / q, 1 + Linit))
print('time %.0fs' % (time.time() - t_start))

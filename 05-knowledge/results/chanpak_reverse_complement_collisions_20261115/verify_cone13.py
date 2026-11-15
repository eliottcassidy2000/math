#!/usr/bin/env python3
"""Exact verifier for the 13-unit bottleneck-free growth certificate (span-3 floor units).

Usage: python3 verify_cone13.py <ribbon_dir> <hankel.pkl> <cert.pkl>

Proves: every bottleneck-free unit word W (letter 0 followed by span-3 floor units) has
    H(W) >= C * rho0^|W|,     rho0 = cert['rho0'] (rational),
where |W| is the number of vertices and C is printed.  Exact integer/rational arithmetic throughout.
"""
import os, sys, pickle, itertools, random, math
import numpy as np
from fractions import Fraction as Fr
RIB, HPKL, CPKL = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, RIB)
from ribbon import transfer_system, tournament_from_word, ham_count_bruteforce
D = pickle.load(open(HPKL, 'rb')); CERT = pickle.load(open(CPKL, 'rb'))
U = D['units']; n = len(U); Lw = [len(u.split('.')) for u in U]
RHO = Fr(CERT['rho0']); DEPTH = CERT['depth']; Q = CERT['Q']; Phi = CERT['Phi']

# 1. nonnegative profile realization (ribbon transfer, THM-4880)
order, states, mats, final = transfer_system(3)
NP = len(order)
def lets(s): return [frozenset(int(c) for c in t if c in '123') for t in s.split('.')]
def letmat(a):
    A = np.zeros((NP, NP), dtype=np.int64)
    for (j, i), m in mats[a].items(): A[j, i] += m
    return A
LM = {a: letmat(a) for a in mats}
AU = []
for u in U:
    A = np.identity(NP, dtype=np.int64)
    for a in lets(u): A = LM[a] @ A
    AU.append(A)
assert all((A >= 0).all() for A in AU)
v0 = LM[frozenset()] @ np.eye(NP, dtype=np.int64)[0]
fin = np.array([1 if j in set(final) else 0 for j in range(NP)], dtype=np.int64)
REACH = set(int(j) for j in np.nonzero(v0)[0]); fr = list(REACH)
while fr:
    nf = []
    for j in fr:
        for A in AU:
            for k in np.nonzero(A[:, j])[0]:
                if int(k) not in REACH: REACH.add(int(k)); nf.append(int(k))
    fr = nf
REACH = sorted(REACH)
SUF = [None, 0, 1, 2, 3, 5, 6, 7]
O = np.array([fin if j is None else fin @ AU[j] for j in SUF], dtype=np.int64)
print('1. realization: %d profile states, %d reachable at unit boundaries, nonnegative unit matrices' % (NP, len(REACH)))

# 2. Hankel representation: correct on the reachable span V, reachable, observable (hence minimal)
M = [[[Fr(int(e.p), int(e.q)) for e in D['red'][u].row(i)] for i in range(8)] for u in range(n)]
class Ech:
    def __init__(s): s.rows = []
    def reduce(s, v):
        r = [Fr(int(a)) for a in v]
        for p, b in s.rows:
            if r[p]:
                f = r[p] / b[p]; r = [x - f * y for x, y in zip(r, b)]
        return r
    def add(s, v):
        r = s.reduce(v); p = next((i for i, x in enumerate(r) if x), None)
        if p is None: return False
        s.rows.append((p, r)); return True
E = Ech(); E.add(v0); basis = [v0]; fr = [v0]
while fr:
    nf = []
    for x in fr:
        for A in AU:
            y = A @ x
            if E.add(y): basis.append(y); nf.append(y)
    fr = nf
for x in basis:
    Ox = [Fr(int(a)) for a in O @ x]
    for u in range(n):
        lhs = [Fr(int(a)) for a in O @ (AU[u] @ x)]
        rhs = [sum(M[u][i][k] * Ox[k] for k in range(8)) for i in range(8)]
        assert lhs == rhs
E8 = Ech(); r8 = sum(1 for x in basis if E8.add(O @ x)); assert r8 == 8
for idx, j in enumerate(SUF):           # observability: fin . M_{s_j} = e_j
    row = [Fr(int(i == 0)) for i in range(8)] if j is None else M[j][0]
    assert row == [Fr(int(i == idx)) for i in range(8)]
# every unit maps every reachable profile state into V
for A in AU:
    for j in REACH:
        assert not any(E.reduce(A[:, j])), 'A_u e_j outside the reachable span'
print('2. reachable span dim %d; O A_u = M_u O on it; O(V) = Q^8; fin M_(s_j) = e_j; A_u e_j in V for all reachable j'
      % len(basis))

# 3. exact bottleneck-free triples and 4-words (rank >= 2 in the minimal representation)
den = [math.lcm(*[e.denominator for row in M[u] for e in row]) for u in range(n)]
NMa = [np.array([[int(e * den[u]) for e in row] for row in M[u]], dtype=np.int64) for u in range(n)]
def rk2(seq):
    P = np.identity(8, dtype=np.int64)
    for u in seq: P = NMa[u] @ P
    P = [list(map(int, r)) for r in P if any(r)]
    if not P: return False
    r0 = P[0]
    return any(r0[i] * r[j] != r0[j] * r[i] for r in P[1:] for i in range(8) for j in range(i + 1, 8))
BF2 = set(w for w in itertools.product(range(n), repeat=2) if rk2(w))
BF3 = set(w for w in itertools.product(range(n), repeat=3) if w[:2] in BF2 and w[1:] in BF2 and rk2(w))
BF4 = set(w for w in itertools.product(range(n), repeat=4) if w[:3] in BF3 and w[1:] in BF3 and rk2(w))
st = CERT['states']; sidx = {s: i for i, s in enumerate(st)}
assert set(st) == BF3
print('3. exact ranks: %d / %d / %d bottleneck-free unit pairs / triples / 4-words; all others have rank <= 1'
      % (len(BF2), len(BF3), len(BF4)))

# 4. cone generators (primitive Hankel directions) and the edge inequalities
def prim(v):
    g = 0
    for a in v: g = math.gcd(g, abs(int(a)))
    return tuple(int(a) // g for a in v)
GEN = {}
for s in st:
    u1, u2, u3 = s
    B = O @ AU[u3] @ AU[u2] @ AU[u1]
    mats_ = [B] if DEPTH == 3 else [B @ AU[u0] for u0 in range(n) if (u0,) + s in BF4]
    cols = set()
    for G in mats_:
        for j in REACH:
            if G[:, j].any(): cols.add(prim(G[:, j]))
    GEN[s] = sorted(cols)
def phi(s, v): return sum(Phi[sidx[s]][k] * v[k] for k in range(8))
fails = checks = 0
for w in sorted(BF4):
    s, t, u = w[:3], w[1:], w[3]
    num, dnm = RHO.numerator ** Lw[u], RHO.denominator ** Lw[u]
    for g in GEN[s]:
        mg = [sum(M[u][i][k] * g[k] for k in range(8)) for i in range(8)]
        lhs = phi(t, mg) * dnm; rhs = phi(s, g) * num
        checks += 1
        if lhs < rhs: fails += 1
print('4. depth-%d cones: %d (edge, generator) inequalities phi_t(M_u g) >= rho0^|u| phi_s(g): %d failures'
      % (DEPTH, checks, fails))
assert fails == 0

# 5. readout, initial values, constant
cmin = None
for s in st:
    for g in GEN[s]:
        assert g[0] >= 0
        p = phi(s, g)
        if p > 0:
            c = Fr(g[0] * Q, p); cmin = c if cmin is None or c < cmin else cmin
init_words = BF3 if DEPTH == 3 else BF4
imin = None
for w in init_words:
    x = v0
    for u in w: x = AU[u] @ x
    p = Fr(phi(w[-3:], [int(a) for a in O @ x]), Q)
    assert p > 0
    val = p / RHO ** (1 + sum(Lw[u] for u in w)); imin = val if imin is None or val < imin else imin
C = cmin * imin
# short words (fewer units than the chain start) checked directly
short = [()] + [(u,) for u in range(n)] + sorted(BF2) + (sorted(BF3) if DEPTH == 4 else [])
for w in short:
    x = v0
    for u in w: x = AU[u] @ x
    assert int(fin @ x) >= C * RHO ** (1 + sum(Lw[u] for u in w))
print('5. readout c = %.6g, initial minimum %.6g, short words checked directly' % (float(cmin), float(imin)))
print('THEOREM: every bottleneck-free word over the 13 span-3 floor units has H(W) >= %.5g * (%s)^|W|'
      % (float(C), RHO))

# 6. control: brute-force Hamiltonian path counts on random bottleneck-free words
rng = random.Random(11); ok = 0
B4 = sorted(BF4)
for trial in range(40):
    w = list(rng.choice(B4))
    while len(w) < 5:
        nxt = [u for u in range(n) if tuple(w[-3:]) + (u,) in BF4]
        if not nxt: break
        w.append(rng.choice(nxt))
    word = [frozenset()]
    for u in w: word += lets(U[u])
    word = [frozenset(x for x in a if x <= q) for q, a in enumerate(word)]
    h = ham_count_bruteforce(tournament_from_word(3, word))
    x = v0
    for u in w: x = AU[u] @ x
    assert h == int(fin @ x) and h >= C * RHO ** len(word)
    ok += 1
print('6. brute force agrees and obeys the bound on %d random bottleneck-free words' % ok)

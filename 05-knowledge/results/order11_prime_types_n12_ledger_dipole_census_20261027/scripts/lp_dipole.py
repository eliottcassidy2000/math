#!/usr/bin/env python3
"""lp_dipole.py n [idx ...] -- the fractional dipole distance to the transitive tournament.
A native word of length l from T to T' satisfies (K' - (-1)^l K)/2 = sum_i (-1)^(l-i) A_{Q_i} with
A_{Q_i} the signed matrix of the i-th directed path.  Relaxing the orientation constraint and
integrality gives, for each parity, LP(D) = min sum c_Q  s.t.  sum_Q c_Q A_Q = D, c >= 0, over all
n! directed Hamiltonian paths Q of K_n.  Then d(T,T') >= min(ceil_even LP(D_even), ceil_odd LP(D_odd)),
minimised over the labellings of T (the dual certificates are exactly the linear negation-invariant
potentials).  Compares with the exact class distance to TT from analyze_mg/movegraph."""
import sys, itertools, collections
import numpy as np
from scipy.optimize import linprog
n = int(sys.argv[1]); sel = [int(x) for x in sys.argv[2:]]
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]; P = len(pairs); pidx = {p: k for k, p in enumerate(pairs)}
# path matrices as columns (vector over pairs: +1 if arc i->j with i<j, -1 if j->i)
cols = []
for perm in itertools.permutations(range(n)):
    v = np.zeros(P)
    for a, b in zip(perm, perm[1:]):
        if a < b: v[pidx[(a, b)]] = 1
        else: v[pidx[(b, a)]] = -1
    cols.append(v)
Aeq = np.array(cols).T  # P x n!
cost = np.ones(Aeq.shape[1])
def lp(D):
    r = linprog(cost, A_eq=Aeq, b_eq=D, bounds=(0, None), method='highs')
    assert r.status == 0, r.message
    return r.fun
def kvec(A):
    return np.array([1 if A[i, j] else -1 for (i, j) in pairs], dtype=float)
V = {}; E = collections.defaultdict(list)
for line in open(f'mg{n}.txt'):
    t = line.split()
    if t[0] == 'V': V[int(t[1])] = dict(key=int(t[2], 16), s=list(map(int, t[3].split(','))), sc=int(t[5]), aut=int(t[6]))
    else: a, b, blue, c = map(int, t[1:]); E[a].append(b)
def adj_from_key(key):
    A = np.zeros((n, n), dtype=int); bp = 0
    for i in range(n):
        for j in range(i+1, n):
            if key >> bp & 1: A[i, j] = 1
            else: A[j, i] = 1
            bp += 1
    return A
TT = [i for i, v in V.items() if sorted(v['s']) == list(range(n))][0]
dist = {TT: 0}; q = collections.deque([TT])
while q:
    a = q.popleft()
    for b in E[a]:
        if b not in dist: dist[b] = dist[a] + 1; q.append(b)
Att = np.triu(np.ones((n, n), dtype=int), 1); Ktt = kvec(Att)
def ceil_par(x, par):
    l = int(np.ceil(x - 1e-7)); 
    if l % 2 != par: l += 1
    return l
def L_pot(s):
    s = sorted(s); return sum(abs(s[i] - i) for i in range(n)) // 2
classes = sel if sel else sorted(V)
summary = collections.Counter()
for c in classes:
    A = adj_from_key(V[c]['key']); best = 10**9; best_lp = None
    seen = set()
    for perm in itertools.permutations(range(n)):
        B = A[np.ix_(perm, perm)]; kb = kvec(B); key = kb.tobytes()
        if key in seen: continue
        seen.add(key)
        De = (Ktt - kb) / 2; Do = (Ktt + kb) / 2
        le = lp(De) if np.any(De) else 0.0; lo = lp(Do) if np.any(Do) else 0.0
        bound = min(ceil_par(le, 0), ceil_par(lo, 1))
        if bound < best: best = bound; best_lp = (le, lo)
    d = dist[c]; summary[(d, best, d == best)] += 1
    print(f'class {c} scores={V[c]["s"]} sc={V[c]["sc"]} aut={V[c]["aut"]} dist={d} L={L_pot(V[c]["s"])} LPbound={best} (lp_even={best_lp[0]:.3f}, lp_odd={best_lp[1]:.3f}) {"EXACT" if d == best else "GAP"}', flush=True)
print('summary (dist, LPbound, exact): count')
for k, v in sorted(summary.items()): print('  ', k, v)

"""Audit S50 - a0: self-tests of the independent helpers (alib.py) before they are used.
  * ham_dp == ham_ie (inclusion-exclusion) == brute force on random tournaments;
  * ear_values (quadratic-form extension formula) == direct DP on every extension of random parents;
  * canon: random relabelings give the same code; code equality agrees with networkx (VF2) isomorphism on random
    pairs and on 'one arc flipped' pairs; |Aut| from the leaf count agrees with networkx automorphism counts."""
import random, time, itertools
import numpy as np
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher
import alib as L

T0 = time.time()
FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


rng = random.Random(12345)


def rand_t(n):
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
    return adj


def relabel(adj, perm):
    n = len(adj)
    new = [0] * n
    for u in range(n):
        for w in range(n):
            if (adj[u] >> w) & 1:
                new[perm[u]] |= 1 << perm[w]
    return new


def to_nx(adj):
    G = nx.DiGraph()
    n = len(adj)
    G.add_nodes_from(range(n))
    G.add_edges_from((u, w) for u in range(n) for w in range(n) if (adj[u] >> w) & 1)
    return G


# 1. path counters
ok = True
for _ in range(300):
    n = rng.randint(1, 8)
    t = rand_t(n)
    a = np.array(t, dtype=np.int64)
    ok &= L.ham_dp(n, a) == L.ham_ie(n, a) == L.ham_brute(t)
for _ in range(40):
    n = rng.randint(9, 12)
    t = rand_t(n)
    a = np.array(t, dtype=np.int64)
    ok &= L.ham_dp(n, a) == L.ham_ie(n, a)
check(ok, 'ham_dp = ham_ie = brute force (300 random tournaments n <= 8), ham_dp = ham_ie (40 random, n = 9..12)')

# 2. ear formula
ok = True; cnt = 0
for _ in range(60):
    n = rng.randint(3, 9)
    t = rand_t(n)
    vals = L.ear_values(n, np.array(t, dtype=np.int64), np.array(L.converse(t), dtype=np.int64))
    for x in range(1 << n):
        cnt += 1
        ok &= int(vals[x]) == L.H(L.extend(t, x))
check(ok, f'ear_values (quadratic form) = direct DP on all {cnt} one-vertex extensions of 60 random parents')

# 3. canonical form
ok = True
for _ in range(400):
    n = rng.randint(3, 12)
    t = rand_t(n)
    c = L.canon(t)
    for _ in range(3):
        p = list(range(n)); rng.shuffle(p)
        ok &= L.canon(relabel(t, p)) == c
    ok &= L.canon(L.decode(n, c)) == c
check(ok, 'canon invariant under random relabeling and idempotent on decoded codes (400 random tournaments, n <= 12)')

ok = True; tested = 0
for _ in range(1500):
    n = rng.randint(4, 9)
    t1 = rand_t(n)
    if rng.random() < 0.5:
        t2 = rand_t(n)
    else:                                   # flip one arc: a hard near-miss
        t2 = list(t1); i, j = rng.sample(range(n), 2)
        if (t2[i] >> j) & 1:
            t2[i] ^= 1 << j; t2[j] |= 1 << i
        else:
            t2[j] ^= 1 << i; t2[i] |= 1 << j
    same = L.canon(t1) == L.canon(t2)
    iso = nx.is_isomorphic(to_nx(t1), to_nx(t2))
    ok &= same == iso; tested += 1
check(ok, f'canon(T1) == canon(T2) iff networkx VF2 says isomorphic ({tested} random and arc-flip pairs, n = 4..9)')


def circ(n, S):
    return [sum(1 << ((x + d) % n) for d in S) for x in range(n)]


def lex(A, B):
    na, nb = len(A), len(B)
    T = [0] * (na * nb)
    for a in range(na):
        for b in range(nb):
            for a2 in range(na):
                for b2 in range(nb):
                    if (a != a2 and (A[a] >> a2) & 1) or (a == a2 and (B[b] >> b2) & 1):
                        T[a * nb + b] |= 1 << (a2 * nb + b2)
    return T


ok = True; rows = []
tests = [('C3', circ(3, {1})), ('Z5{1,2}', circ(5, {1, 2})), ('P7', circ(7, {1, 2, 4})), ('Z7{1,2,3}', circ(7, {1, 2, 3})),
         ('Z9{1,2,3,5}', circ(9, {1, 2, 3, 5})), ('C3[C3]', lex(circ(3, {1}), circ(3, {1}))), ('P11', circ(11, {1, 3, 4, 5, 9}))]
for _ in range(20):
    tests.append(('rand', rand_t(rng.randint(5, 9))))
for name, t in tests:
    c, a = L.canon(t, want_aut=True)
    G = to_nx(t)
    a2 = sum(1 for _ in DiGraphMatcher(G, G).isomorphisms_iter())
    ok &= a == a2
    if name != 'rand':
        rows.append(f'{name}: |Aut| = {a}')
check(ok, '|Aut| from canon leaf count = networkx automorphism count: ' + ', '.join(rows) + ' (+ 20 random)')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')

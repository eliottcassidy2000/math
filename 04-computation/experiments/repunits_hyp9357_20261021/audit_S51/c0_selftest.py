"""S51 audit c0: self-test of alib51 (the auditor's own helpers) against brute force and networkx.
Prints ALL CHECKS PASSED or the failures."""
import itertools, random, time
import numpy as np
import networkx as nx
import alib51 as L

FAIL = []
T0 = time.time()


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


def rand_tour(n, rng):
    adj = np.zeros(n, dtype=np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
    return adj


def relabel(adj, p):
    n = len(adj)
    new = np.zeros(n, dtype=np.int64)
    for u in range(n):
        for w in range(n):
            if (int(adj[u]) >> w) & 1:
                new[p[u]] |= 1 << p[w]
    return new


def aut_brute(adj):
    n = len(adj)
    return sum(1 for p in itertools.permutations(range(n)) if all(((int(adj[u]) >> w) & 1) == ((int(adj[p[u]]) >> p[w]) & 1) for u in range(n) for w in range(n)))


rng = random.Random(20261021)
# Hamiltonian paths
ok = True
for _ in range(300):
    n = rng.randint(1, 8)
    a = rand_tour(n, rng)
    ok &= L.is_tournament(a) and int(L.ham(n, a)) == L.ham_brute(a)
check(ok, 'ham (subset DP) = brute force on 300 random tournaments, n <= 8')
ok = True
for _ in range(40):
    n = rng.randint(3, 8)
    a = rand_tour(n, rng)
    Hs = L.ham_subsets(n, a)
    for S in range(1, 1 << n):
        vs = [v for v in range(n) if (S >> v) & 1]
        sub = np.array([sum(1 << vs.index(w) for w in vs if (int(a[u]) >> w) & 1) for u in vs], dtype=np.int64)
        ok &= int(Hs[S]) == int(L.ham(len(vs), sub))
check(ok, 'ham_subsets = ham of every induced subtournament (40 random tournaments)')
# strongness
ok = True
for _ in range(300):
    n = rng.randint(1, 9)
    a = rand_tour(n, rng)
    G = nx.DiGraph()
    G.add_nodes_from(range(n))
    for u in range(n):
        for w in range(n):
            if (int(a[u]) >> w) & 1:
                G.add_edge(u, w)
    S = rng.randint(1, (1 << n) - 1)
    vs = [v for v in range(n) if (S >> v) & 1]
    ok &= bool(L.is_strong_mask(n, a, S)) == nx.is_strongly_connected(G.subgraph(vs))
check(ok, 'is_strong_mask = networkx strong connectivity on 300 random induced subtournaments')
# feedback arc number
ok = True
for _ in range(150):
    n = rng.randint(1, 7)
    a = rand_tour(n, rng)
    best = min(sum(1 for i in range(n) for j in range(i + 1, n) if (int(a[p[j]]) >> p[i]) & 1) for p in itertools.permutations(range(n)))
    ok &= int(L.fas(n, a)) == best
check(ok, 'fas (subset DP) = brute force over all orderings, 150 random tournaments n <= 7')
# canonical form: invariance, |Aut|, and class counts by extension
ok = True
for _ in range(300):
    n = rng.randint(1, 12)
    a = rand_tour(n, rng)
    p = list(range(n)); rng.shuffle(p)
    c1 = L.canon(n, a); c2 = L.canon(n, relabel(a, p))
    d = L.decode(n, c1[0], c1[1])
    c3 = L.canon(n, d)
    ok &= c1 == c2 and c3 == c1 and int(L.ham(n, d)) == int(L.ham(n, a))
check(ok, 'canon is relabelling-invariant and decode(canon) is a representative with the same code and H (300 random, n <= 12)')
classes = {1: {L.canon(1, np.zeros(1, dtype=np.int64))}}
for n in range(2, 8):
    s = set()
    for (w1, w0, _a) in classes[n - 1]:
        par = L.decode(n - 1, w1, w0)
        for x in range(1 << (n - 1)):
            s.add(L.canon(n, L.extend(n - 1, par, x)))
    classes[n] = s
cnt = [len(classes[n]) for n in range(1, 8)]
check(cnt == [1, 1, 2, 4, 12, 56, 456], f'class counts by extension n = 1..7: {cnt} (OEIS A000568: 1, 1, 2, 4, 12, 56, 456)')
ok = True
for n in range(1, 7):
    for (w1, w0, a) in classes[n]:
        ok &= aut_brute(L.decode(n, w1, w0)) == a
check(ok, '|Aut| from the canonical search = brute-force automorphism count, every class n <= 6')
ok = True
for (w1, w0, a) in rng.sample(sorted(classes[7]), 60):
    ok &= aut_brute(L.decode(7, w1, w0)) == a
check(ok, '|Aut| = brute force on 60 random classes at n = 7')
# ear values
ok = True
for _ in range(60):
    n = rng.randint(2, 8)
    a = rand_tour(n, rng)
    ev = L.ear_all(n, a)
    for x in range(1 << n):
        ok &= int(ev[x]) == int(L.ham(n + 1, L.extend(n, a, x)))
check(ok, 'ear_all (quadratic form) = DP of every one-vertex extension, 60 random tournaments n <= 8')
# one-apex closed form
ok = True
for m in range(2, 11):
    for r in range(0, m):
        for A in itertools.combinations(range(1, m), r):
            ok &= L.one_apex_paths(m, A) == int(L.ham(m, L.one_apex(m, A)))
check(ok, 'one_apex_paths (path decomposition, any A) = DP for every in-set A, m <= 10')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')

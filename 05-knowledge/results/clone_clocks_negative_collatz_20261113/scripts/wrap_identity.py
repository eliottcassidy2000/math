"""Wrap identity test.  For a two-cherry skeleton S (in-set I, out-set O of u; all I->O except two cherries
s1,s2->t1, s3,s4->t2) and a vertex v in I whose only in-neighbour is the target t1, the two-module sectors satisfy
X_1 - 2 X_2 = H(B) - s_t1(B) - E_t1(B) + 2 H(B - t1),  B = S - v  (proved by the first-run bijection).
Question: is X_1 - 2 X_2 = H(M) (the wrap law) for EVERY such skeleton, clock or not?  Checked through the
two-module visit table (c_11 = 1 and c_21 = c_31 = 0 in the reduced form X_1 - 2 X_2 = sum_i (a_i1 - 2 a_i2) K_i)."""
import random, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from twomod import visit_table
rng = random.Random(13)
stats = {}
for trial in range(400):
    p = rng.randint(3, 4); q = rng.randint(4, 5); n = p + q + 1; u = p
    adj = [0]*n
    # random T[I], T[O]
    for a in range(p):
        for b in range(a+1, p):
            if rng.random() < 0.5: adj[a] |= 1 << b
            else: adj[b] |= 1 << a
    for a in range(q):
        for b in range(a+1, q):
            x, y = p+1+a, p+1+b
            if rng.random() < 0.5: adj[x] |= 1 << y
            else: adj[y] |= 1 << x
    for a in range(p):
        adj[a] |= 1 << u
        for b in range(q): adj[a] |= 1 << (p+1+b)
    for b in range(q): adj[u] |= 1 << (p+1+b)
    t1, t2 = rng.sample(range(p), 2)
    S = rng.sample(range(q), 4)
    for k, (s, t) in enumerate(zip(S, [t1, t1, t2, t2])):
        o = p+1+s; adj[t] &= ~(1 << o); adj[o] |= 1 << t
    # choose v in I, non-target, force its only in-neighbour to be t1
    cands = [x for x in range(p) if x not in (t1, t2)]
    if not cands: continue
    v = rng.choice(cands)
    for x in range(n):
        if x == v: continue
        if x == t1:
            adj[x] |= 1 << v; adj[v] &= ~(1 << x)
        else:
            adj[v] |= 1 << x; adj[x] &= ~(1 << v)
    tab = visit_table(adj, n, u, v, CU=8, CV=6)
    # reduced sector D_i = a_i1 - 2 a_i2 ; wrap law <=> D = (1,0,0,...) i.e. X1 - 2X2 = K_1 = H(M)
    D = [tab.get((i, 1), 0) - 2*tab.get((i, 2), 0) for i in range(1, 6)]
    key = 'wrap law holds' if D == [1, 0, 0, 0, 0] else 'fails'
    stats[key] = stats.get(key, 0) + 1
    if key == 'fails' and stats[key] <= 4:
        print('fails: p q', p, q, 'D =', D, ' ternary at u?', tab.get((4,1),0)==0 and tab.get((2,1),0)==3*tab.get((3,1),0))
print(stats)

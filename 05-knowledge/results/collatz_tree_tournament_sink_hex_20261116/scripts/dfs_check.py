#!/usr/bin/env python3
"""dfs_check.py: for random rooted trees on n <= 8 nodes and random linear extensions l (build orders), the tree tournament
T(P, l) (ancestor -> descendant; incomparable pairs: later-built beats earlier-built) — test the claims
  (i)  T is transitive  <=>  l is a preorder (depth-first) linear extension, i.e. every subtree is built contiguously;
  (ii) the 3-cycles of T are exactly the triples (x, y, z) with x an ancestor of y, z incomparable to both, l(x) < l(z) < l(y);
  (iii) H(T) = 1 iff T transitive (Redei), and H >= 1 + 2 * #3-cycles.
Also prints, for the Collatz BFS orders, the 3-cycle counts c3(N) and the bound 1 + 2 c3 against H."""
import random, itertools
random.seed(20261116)
def random_tree(n):
    par = [None] + [random.randrange(i) for i in range(1, n)]
    return par
def anc(par, a, b):
    while b is not None:
        if b == a: return True
        b = par[b]
    return False
def linear_extensions(par, n):
    # all orders where parents precede children
    res = []
    def rec(order, placed):
        if len(order) == n: res.append(tuple(order)); return
        for v in range(n):
            if v not in placed and (par[v] is None or par[v] in placed):
                placed.add(v); order.append(v); rec(order, placed); order.pop(); placed.discard(v)
    rec([], set()); return res
def tournament(par, n, l):
    pos = {v: i for i, v in enumerate(l)}
    def beats(a, b):
        if anc(par, a, b): return True
        if anc(par, b, a): return False
        return pos[a] > pos[b]
    return beats
def is_preorder(par, n, l):
    # every subtree contiguous: for each v, the positions of v's descendants (incl v) form an interval starting at pos[v]
    pos = {v: i for i, v in enumerate(l)}
    for v in range(n):
        desc = [u for u in range(n) if anc(par, v, u)]
        ps = sorted(pos[u] for u in desc)
        if ps != list(range(ps[0], ps[0] + len(ps))) or ps[0] != pos[v]: return False
    return True
def ham(beats, n):
    return sum(1 for p in itertools.permutations(range(n)) if all(beats(p[i], p[i+1]) for i in range(n-1)))
bad = 0; tested = 0; bad3 = 0
for trial in range(60):
    n = random.randint(4, 7); par = random_tree(n)
    exts = linear_extensions(par, n); random.shuffle(exts)
    for l in exts[:12]:
        beats = tournament(par, n, l); pos = {v: i for i, v in enumerate(l)}
        trans = all(not (beats(x, y) and beats(y, z) and beats(z, x)) for x in range(n) for y in range(n) for z in range(n) if len({x, y, z}) == 3)
        pre = is_preorder(par, n, l)
        if trans != pre: bad += 1
        cyc = set()
        for x, y, z in itertools.permutations(range(n), 3):
            if beats(x, y) and beats(y, z) and beats(z, x): cyc.add(frozenset((x, y, z)))
        pred = set()
        for x in range(n):
            for y in range(n):
                if x != y and anc(par, x, y):
                    for z in range(n):
                        if z not in (x, y) and not anc(par, x, z) and not anc(par, z, x) and not anc(par, y, z) and not anc(par, z, y) and pos[x] < pos[z] < pos[y]:
                            pred.add(frozenset((x, y, z)))
        if cyc != pred: bad3 += 1
        H = ham(beats, n)
        if (H == 1) != trans or H < 1 + 2*len(cyc): bad += 1
        tested += 1
print(f"random trees/orders tested: {tested}; failures of (i)/(iii): {bad}; failures of (ii): {bad3}")
# Collatz BFS orders: 3-cycle counts vs H
def collatz_tree(order_even_first, want):
    lab = [1]; par = [None]; depth = [0]; q = [0]
    while len(lab) < want + 4:
        v = q.pop(0); n = lab[v]
        kids = [2*n] + ([(n-1)//3] if n % 6 == 4 and (n-1)//3 > 1 else [])
        for k in kids: lab.append(k); par.append(v); depth.append(depth[v]+1); q.append(len(lab)-1)
    idx = sorted(range(len(lab)), key=lambda i: (depth[i], (lab[i] % 2) if order_even_first else (1 - lab[i] % 2), lab[i]))
    return lab, par, idx
for ef in (True, False):
    lab, par, idx = collatz_tree(ef, 20); pos = {v: i for i, v in enumerate(idx)}
    out = []
    for N in range(1, 21):
        nodes = idx[:N]
        c3 = 0
        for x in nodes:
            for y in nodes:
                if x != y and anc(par, x, y):
                    for z in nodes:
                        if z not in (x, y) and not anc(par, x, z) and not anc(par, z, x) and not anc(par, y, z) and not anc(par, z, y) and pos[x] < pos[z] < pos[y]: c3 += 1
        out.append(c3)
    print(("even-first" if ef else "odd-first"), "c3(N) for N=1..20:", out)

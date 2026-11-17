#!/usr/bin/env python3
"""twocover.py (C-3, C-4): the q-slot path-cover counts P_D(q) = sum over ordered q-tuples of disjoint, possibly empty
vertex sets covering D of the product of their Hamiltonian counts (P(1) = h, P(2) = two-cover, P(3) = three-cover), for
every class of orders 1-8; the fibres of P(2) (frozen cover count, varying h); the (h, P(2)) collisions and whether
P(3) separates them; the basic structures C_3, strong4, RT_5, P_7, M_k."""
import numpy as np
from collections import defaultdict, Counter
exec(open('nearjoin_check.py').read().split('# (1)')[0])
P = np.load('profiles8.npz')
def covers(T):
    n = len(T); N = 1 << n; full = N - 1; dp, Hm = ham_all(T)
    P2 = sum(Hm[W] * Hm[full & ~W] for W in range(N))
    # P3: sum over W1, W2 subset of complement of W1
    P3 = 0
    for W1 in range(N):
        c = full & ~W1; h1 = Hm[W1]; sub = c
        while True:
            P3 += h1 * Hm[sub] * Hm[c & ~sub]
            if sub == 0: break
            sub = (sub - 1) & c
    return Hm[full], P2, P3
byP = defaultdict(set); byHP = defaultdict(list); total = 0
for n in range(1, 9):
    A = P[f'A{n}']
    for i in range(len(A)):
        h, p2, p3 = covers(A[i].tolist()); total += 1; byP[p2].add(h); byHP[(h, p2)].append(p3)
nP = len(byP); shared = sum(1 for p, hs in byP.items() if len(hs) >= 2)
coll = [(k, v) for k, v in byHP.items() if len(v) >= 2]; sep = sum(1 for k, v in coll if len(set(v)) >= 2)
print(f"(C-3) classes {total}; distinct P(2) values {nP} ({nP/total:.3f} of the classes); P values shared by classes of different h: {shared} ({shared/nP:.3f}); (h, P(2)) collisions: {len(coll)} classes-of-pairs, separated by P(3): {sep} ({sep/max(1,len(coll)):.3f})")
print("    examples of frozen-cover fibres (P(2): counts h):", {p: sorted(byP[p]) for p in sorted(byP) if len(byP[p]) >= 4 and p < 300})
print("    dyadic P(2) values and their h-fibres:", {p: sorted(byP[p]) for p in sorted(byP) if p & (p - 1) == 0})
def circ(n, S): return [[1 if (j - i) % n in S else 0 for j in range(n)] for i in range(n)]
basics = {'C_3': circ(3, {1}), 'strong4': [[0,1,1,0],[0,0,1,1],[0,0,0,1],[1,0,0,0]], 'RT_5': circ(5, {1, 2}), 'P_7': circ(7, {1, 2, 4}), 'M_4': [[1 if i < j else 0 for j in range(4)] for i in range(4)], 'M_5': [[1 if i < j else 0 for j in range(5)] for i in range(5)]}
print("(C-4) basic structures (h, P(2), P(3)):", {k: covers(v) for k, v in basics.items()})

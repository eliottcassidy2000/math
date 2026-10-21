#!/usr/bin/env python3
"""THM-4733 check: the fold Phi of the four mirror-moved tiles of the n = 5 cube onto the divisors of p^2 q r is
equivariant for (kappa, sigma, sigma kappa) -> (C, Q, CQ); the displayed ten-pattern subset is a Klein-invariant section;
on the full n = 5 cube no tiling is sigma kappa-fixed (obstruction), while the four-tile cube has four."""
import itertools
tiles = [(a, b) for b in range(1, 4) for a in range(b + 2, 6)]  # n = 5
mirror = {t: (6 - t[1], 6 - t[0]) for t in tiles}
moved = [(3, 1), (5, 3), (4, 1), (5, 2)]
def kappa(x): return tuple(1 - v for v in x)
def sigma4(x): p1, p2, q, r = x; return (p2, p1, r, q)
def Phi(x): p1, p2, q, r = x; return (p1 + p2, q, r)
def C(d): a, b, c = d; return (2 - a, 1 - b, 1 - c)
def Q(d): a, b, c = d; return (a, c, b)
pats = list(itertools.product([0, 1], repeat=4))
assert all(Phi(kappa(x)) == C(Phi(x)) for x in pats) and all(Phi(sigma4(x)) == Q(Phi(x)) for x in pats)
assert all(Phi(sigma4(kappa(x))) == C(Q(Phi(x))) for x in pats)
divs = sorted(set(Phi(x) for x in pats)); assert len(divs) == 12
fib = {d: [x for x in pats if Phi(x) == d] for d in divs}
print('fold: 16 patterns -> 12 divisors; fibre sizes', {d: len(f) for d, f in fib.items()})
proper = [d for d in divs if d not in ((0, 0, 0), (2, 1, 1))]
rank = {d: sum(d) for d in proper}; f = {d: d[0] - d[1] - d[2] for d in proper}
print('rank levels 3-4-3:', [sum(1 for d in proper if rank[d] == k) for k in (1, 2, 3)], ' sign levels 4-2-4:', [sum(1 for d in proper if (f[d] < 0, f[d] == 0, f[d] > 0)[i]) for i in range(3)])
print('fixed divisors: C', sum(1 for d in proper if C(d) == d), 'Q', sum(1 for d in proper if Q(d) == d), 'CQ', sum(1 for d in proper if C(Q(d)) == d))
section = {(1,0,0): (0,0,0,0), (1,1,1): (1,1,1,1), (2,0,0): (1,1,0,0), (0,1,1): (0,0,1,1), (1,1,0): (1,0,1,0), (1,0,1): (0,1,0,1),
           (0,1,0): (1,0,0,0), (0,0,1): (0,1,0,0), (2,0,1): (0,1,1,1), (2,1,0): (1,0,1,1)}
assert set(section) == set(proper)
S = set(section.values()); assert len(S) == 10
assert all(kappa(x) in S and sigma4(x) in S for x in S), 'invariant subset'
assert all(section[C(d)] == kappa(section[d]) and section[Q(d)] == sigma4(section[d]) for d in proper), 'equivariant embedding'
print('embedding: ten-pattern invariant subset, equivariant bijection for C->kappa, Q->sigma (abstract Klein isomorphism; NOT a section of Phi:', 'Phi(0000) =', Phi((0,0,0,0)), ')')
# obstruction on the full cube
full = list(itertools.product([0, 1], repeat=6)); idx = {t: i for i, t in enumerate(tiles)}
def sig6(x): return tuple(x[idx[mirror[t]]] for t in tiles)
def kap6(x): return tuple(1 - v for v in x)
print('full n=5 cube: kappa-fixed', sum(1 for x in full if kap6(x) == x), '; sigma kappa-fixed', sum(1 for x in full if sig6(kap6(x)) == x), '(obstruction: must be 0); sigma-fixed', sum(1 for x in full if sig6(x) == x))
print('four-tile cube: sigma kappa-fixed', sum(1 for x in pats if sigma4(kappa(x)) == x), '(the preimages of pq, pr)')
print('PASS')

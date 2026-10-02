"""burnside.py -- number of Aut(Q_n)-orbits of a-subsets of V(Q_n), via Burnside / cycle index.
Independent of all C code.  Aut(Q_n) = {x -> pi(x xor t)}: 2^n n! elements.
Usage: python3 burnside.py n   -> prints a and #orbits for a = 0..2^n, writes data/burnside_n.txt
"""
import itertools, sys, os
from math import comb, factorial

def group_cycle_types(n):
    N = 1 << n
    for pi in itertools.permutations(range(n)):
        for t in range(N):
            perm = [0] * N
            for v in range(N):
                w = v ^ t; u = 0
                for i in range(n):
                    if (w >> i) & 1: u |= 1 << pi[i]
                perm[v] = u
            seen = [False] * N; cyc = []
            for v in range(N):
                if not seen[v]:
                    L = 0; x = v
                    while not seen[x]:
                        seen[x] = True; x = perm[x]; L += 1
                    cyc.append(L)
            yield tuple(sorted(cyc))

def burnside(n):
    N = 1 << n
    from collections import Counter
    types = Counter(group_cycle_types(n))
    G = sum(types.values()); assert G == (1 << n) * factorial(n)
    tot = [0] * (N + 1)
    for cyc, mult in types.items():
        poly = [1] + [0] * N
        for L in cyc:
            new = poly[:]
            for i in range(N + 1 - L):
                if poly[i]: new[i + L] += poly[i]
            poly = new
        for i in range(N + 1): tot[i] += mult * poly[i]
    assert all(x % G == 0 for x in tot)
    return [x // G for x in tot], G, len(types)

if __name__ == '__main__':
    n = int(sys.argv[1])
    orb, G, ntypes = burnside(n)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'burnside_%d.txt' % n)
    with open(out, 'w') as f:
        for a, x in enumerate(orb): f.write('%d %d\n' % (a, x))
    print('|Aut(Q%d)| = %d, distinct cycle types = %d' % (n, G, ntypes))
    for a in range(min(len(orb), 20)): print(a, orb[a], '  (C(N,a)/|G| = %.1f)' % (comb(1 << n, a) / G))

#!/usr/bin/env python3
"""Independent exhaustive isoclass and exact-rational up/down audit, n<=5.
No production canonicalizer or tournament banks imported. Checks survive -O.
"""
import itertools as it
import json
import math
from collections import Counter
from pathlib import Path
import sympy as s
from sympy.matrices.normalforms import smith_normal_form


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def edges(n):
    return list(it.combinations(range(n), 2))


def arc(code, n, u, v):
    if u > v:
        return 1 - arc(code, n, v, u)
    return (code >> edges(n).index((u, v))) & 1


def induced(code, n, vertices):
    return sum(arc(code, n, vertices[i], vertices[j]) << k
               for k, (i, j) in enumerate(edges(len(vertices))))


def main():
    N = 5
    canon, banks, sizes, aut = {}, {}, {}, {}
    for n in range(N + 1):
        canon[n] = {}
        permutations = list(it.permutations(range(n)))
        for code in range(1 << math.comb(n, 2)):
            if code in canon[n]:
                continue
            orbit = {induced(code, n, p) for p in permutations}
            representative = min(orbit)
            for member in orbit:
                canon[n][member] = representative
        banks[n] = sorted(set(canon[n].values()))
        sizes[n] = Counter(canon[n].values())
        aut[n] = {t: math.factorial(n) // sizes[n][t] for t in banks[n]}
    counts = [len(banks[n]) for n in range(N + 1)]
    require(counts == [1, 1, 1, 2, 4, 12], 'independent class counts')
    D, U, E = {}, {}, {}
    for n in range(1, N + 1):
        D[n] = s.zeros(len(banks[n-1]), len(banks[n]))
        U[n-1] = s.zeros(len(banks[n]), len(banks[n-1]))
        for col, code in enumerate(banks[n]):
            for deleted in range(n):
                small = induced(code, n, [v for v in range(n) if v != deleted])
                D[n][banks[n-1].index(canon[n-1][small]), col] += 1
        for col, code in enumerate(banks[n-1]):
            for mask in range(1 << (n-1)):
                extension = sum((arc(code, n-1, i, j) if j < n-1 else (mask >> i) & 1) << k
                                for k, (i, j) in enumerate(edges(n)))
                U[n-1][banks[n].index(canon[n][extension]), col] += 1
        E[n-1] = U[n-1] / 2**(n-1)
        require(U[n-1].rank() == len(banks[n-1]), 'U injective')
        require(D[n].rank() == len(banks[n-1]), 'D surjective')
        for i, t in enumerate(banks[n-1]):
            for j, u in enumerate(banks[n]):
                require(U[n-1][j,i]*aut[n][u] == D[n][i,j]*aut[n-1][t], 'marked orbit identity')
    delta, sourcefree = {}, {}
    sourcefree[0] = 1
    for n in range(1,N+1):
        delta[n] = s.zeros(len(banks[n-1]),len(banks[n]))
        sourcefree[n] = 0
        for col, t in enumerate(banks[n]):
            sources = [v for v in range(n) if all(arc(t,n,v,w) for w in range(n) if w != v)]
            require(len(sources) <= 1, 'unique universal source')
            if sources:
                deleted = induced(t,n,[w for w in range(n) if w != sources[0]])
                delta[n][banks[n-1].index(canon[n-1][deleted]),col] = 1
            else:
                sourcefree[n] += 1
        require(sourcefree[n] == len(banks[n])-len(banks[n-1]), 'sourcefree count harmonic dimension')
    for n in range(N):
        source_prior = U[n-1]*delta[n] if n else s.zeros(1)
        require(delta[n+1]*U[n]-source_prior == s.eye(len(banks[n])), 'integer source-deletion Weyl relation')
    for n in range(N):
        prior = E[n-1]*D[n] if n else s.zeros(1)
        require(D[n+1]*E[n] - prior == s.eye(len(banks[n])), 'Weyl commutator')
    result = {'class_counts': counts, 'harmonic_dimensions': [], 'smith_deletion': {}, 'sourcefree_counts': [sourcefree[n] for n in range(N+1)], 'levels': []}
    for n in range(1,N+1):
        smith = smith_normal_form(D[n], domain=s.ZZ)
        diagonal = [abs(int(smith[i,i])) for i in range(smith.rows)]
        require(all(diagonal) and all(diagonal[i+1] % diagonal[i] == 0 for i in range(len(diagonal)-1)), 'Smith full rank and divisibility')
        require(s.ones(1,len(banks[n-1]))*D[n] == n*s.ones(1,len(banks[n])), 'integral augmentation obstruction')
        require(math.prod(diagonal) % n == 0, 'cokernel order divisible by n')
        inherited_entries = [n-r for r in range(n) for _ in range(counts[r]-(counts[r-1] if r else 0))]
        inherited_smith = smith_normal_form(s.diag(*inherited_entries), domain=s.ZZ)
        require(diagonal == [abs(int(inherited_smith[i,i])) for i in range(inherited_smith.rows)], 'inherited all-degree deletion cokernel theorem')
        result['smith_deletion'][n] = diagonal
    parity_witnesses = [list(v) for v in it.product(range(2), repeat=D[5].rows) if any(v) and all(int(x)%2 == 0 for x in s.Matrix([v])*D[5])]
    require(bool(parity_witnesses), 'degree-five extra mod-two obstruction')
    result['degree5_mod2_left_kernel'] = parity_witnesses
    result['degree4_canonical_codes'] = banks[4]
    for n in range(N + 1):
        dim = len(banks[n])
        h = dim - (len(banks[n-1]) if n else 0)
        result['harmonic_dimensions'].append(h)
        P = s.eye(dim)
        for k in range(1, n+1):
            down, up = s.eye(dim), s.eye(len(banks[n-k]))
            for degree in range(n, n-k, -1):
                down = D[degree]*down
            for degree in range(n-k, n):
                up = E[degree]*up
            P += (-1)**k * up*down / math.factorial(k)
        require(P*P == P and P.rank() == h, 'harmonic projector')
        if n:
            require(D[n]*P == s.zeros(len(banks[n-1]), dim), 'projector image kernel')
            require(P*E[n-1] == s.zeros(dim,len(banks[n-1])), 'projector kills E')
        blocks = []
        for j in range(n+1):
            basis = D[j].nullspace() if j else [s.ones(1,1)]
            for v in basis:
                for degree in range(j,n):
                    v = E[degree]*v
                blocks.append(v)
        require(s.Matrix.hstack(*blocks).rank() == dim, 'direct sum harmonic lifts')
        if not n:
            continue
        K = E[n-1]*D[n]/n
        require(n*s.trace(K) == sum(counts[:n]), 'trace cumulative count identity')
        # Independent labeled star-refresh operation: no deletion or extension matrices.
        for col, code in enumerate(banks[n]):
            frequencies = Counter()
            for vertex in range(n):
                star = [k for k,(i,j) in enumerate(edges(n)) if vertex in (i,j)]
                for mask in range(1 << (n-1)):
                    changed = code
                    for bit, edge in enumerate(star):
                        changed = (changed & ~(1 << edge)) | (((mask >> bit)&1) << edge)
                    frequencies[canon[n][changed]] += 1
            for row, target in enumerate(banks[n]):
                require(K[row,col] == s.Rational(frequencies[target], n*2**(n-1)), 'direct labeled refresh agreement')
        pi = s.Matrix([s.Rational(sizes[n][t], 2**math.comb(n,2)) for t in banks[n]])
        require(K*pi == pi and sum(pi) == 1, 'stationarity')
        require(K*s.diag(*pi) == s.diag(*pi)*K.T, 'detailed balance')
        if n >= 3:
            centered_cycles = s.Matrix([sum(all(sum(arc(t,n,u,v) for v in triple if v != u) == 1 for u in triple) for triple in it.combinations(range(n),3)) - s.Rational(math.comb(n,3),4) for t in banks[n]])
            require(K.T*centered_cycles == s.Rational(n-3,n)*centered_cycles, 'cyclic-triangle eigenfunction')
        expected = {s.Rational(n-j,n):len(banks[j])-(len(banks[j-1]) if j else 0) for j in range(n+1)}
        expected = {val:mult for val,mult in expected.items() if mult}
        require(K.eigenvals() == expected, 'exact spectrum')
        result['levels'].append({'n':n,'spectrum':{str(v):m for v,m in expected.items()}, 'projector_rank':h})
    result['status'] = 'PASS: exhaustive isoclasses, ranks, commutators, projections, direct sums, direct labeled refresh, stationary law, detailed balance, exact spectra'
    print(json.dumps(result, indent=2))

if __name__ == '__main__':
    main()

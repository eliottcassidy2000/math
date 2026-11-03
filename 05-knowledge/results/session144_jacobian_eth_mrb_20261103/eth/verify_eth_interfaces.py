#!/usr/bin/env python3
"""Exact finite audits of interfaces used by the imported ETH claim.

This does not run the ETH algorithm or verify its complete proof. Python 3
standard library only. No floating-point comparisons and no asserts, so -O
does not disable checks. Outputs a JSON receipt to stdout.
"""
from fractions import Fraction as Q
from itertools import combinations, permutations
from math import comb
import json


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def rank(a):
    a = [[Q(v) for v in row] for row in a]
    r = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        c = a[r][col]
        a[r] = [v / c for v in a[r]]
        for i in range(len(a)):
            if i != r and a[i][col]:
                c = a[i][col]
                a[i] = [x - c*y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def solve(a, b):
    n = len(a)
    a = [[Q(v) for v in row] + [Q(y)] for row, y in zip(a, b)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        c = a[j][j]
        a[j] = [v/c for v in a[j]]
        for i in range(n):
            if i != j:
                c = a[i][j]
                a[i] = [x-c*y for x,y in zip(a[i],a[j])]
    return [row[-1] for row in a]


def permanent_polynomial(a, b):
    """Coefficients of per(a + T*b), via subset DP."""
    n = len(a)
    dp = {0: [Q(1)]}
    for i in range(n):
        nxt = {}
        for mask, p in dp.items():
            for j in range(n):
                if mask >> j & 1:
                    continue
                out = nxt.setdefault(mask | 1 << j, [Q(0)] * (i+2))
                for d, c in enumerate(p):
                    out[d] += c * a[i][j]
                    out[d+1] += c * b[i][j]
        dp = nxt
    return dp[(1 << n)-1]


def perm_direct(a):
    total = Q(0)
    for p in permutations(range(len(a))):
        term = Q(1)
        for i,j in enumerate(p):
            term *= a[i][j]
            if not term:
                break
        total += term
    return total


def jet(p, t, h):
    return sum((Q(comb(d,h))*p[d]*Q(t)**(d-h)
                for d in range(h,len(p))), Q(0))


def cube_audit():
    # Three commuting coordinate flips. p(x)=-x matches parity at
    # Hamming weights 1 and 2; the two extreme eigenspaces form E.
    q, n = 3, 8
    x = [[[Q(int(v == (u ^ (1 << j)))) for v in range(n)]
          for u in range(n)] for j in range(q)]
    s = [[sum(x[j][u][v] for j in range(q)) for v in range(n)] for u in range(n)]
    jmat = [[Q(int(v == (u ^ 7))) for v in range(n)] for u in range(n)]
    k = [[sum(x[j][u][v]/(j+1) for j in range(q)) for v in range(n)] for u in range(n)]
    e = [[jmat[u][v]+s[u][v] for v in range(n)] for u in range(n)]
    require(rank(e) == 2, 'correction rank')
    for t in (1,2,3):
        b = [[-s[u][v]+t*k[u][v] for v in range(n)] for u in range(n)]
        bit = 1 << (t-1)
        require(all(not b[u][v] for u in range(n) for v in range(n)
                    if (u & bit) != (v & bit)), 'structured coordinate block')
        require(all(jmat[u][v]+t*k[u][v] == b[u][v]+e[u][v]
                    for u in range(n) for v in range(n)), 'exact correction identity')
    p = permanent_polynomial(jmat,k)
    dropped = permanent_polynomial([[-v for v in row] for row in s], k)
    require(p[0] == 1 and dropped[0] == 81, 'sidecar hostile')
    for t in (0,1,2,3,4):
        actual = perm_direct([[jmat[u][v]+t*k[u][v] for v in range(n)] for u in range(n)])
        require(actual == jet(p,t,0), 'independent permanent enumeration')
    rows, data = [], []
    for t in (1,2,3):
        for h in range(3):
            rows.append([Q(comb(d,h))*Q(t)**(d-h) if d>=h else Q(0) for d in range(9)])
            data.append(jet(p,t,h))
    recovered = solve(rows,data)
    require(recovered == p, 'Hermite exact recovery')
    # With one derivative missing, the degree-eight vanishing witness
    # (T-1)^3(T-2)^3(T-3)^2 changes P(0) by 72 yet keeps every retained jet.
    v = [Q(1)]
    for t, multiplicity in ((1,3),(2,3),(3,2)):
        for _ in range(multiplicity):
            w = [Q(0)]*(len(v)+1)
            for d,c in enumerate(v):
                w[d] -= t*c
                w[d+1] += c
            v = w
    require(v[0] == Q(72), 'missing-jet witness constant')
    for t,m in ((1,3),(2,3),(3,2)):
        require(all(jet(v,t,h)==0 for h in range(m)), 'missing-jet invisibility')
    return {'cube_dimension': q, 'matrix_order': n, 'correction_rank': rank(e),
            'constant_with_correction': str(p[0]), 'constant_after_dropping_correction': str(dropped[0]),
            'permanent_coefficients': [str(x) for x in p],
            'independent_permutation_evaluations': 5, 'Hermite_samples': 9,
            'insufficient_jet_hostile_constant_shift': str(v[0])}


def encode_nat(n):
    return [1]*n+[0]


def frame(bits):
    return [v for b in bits for v in (1,b)]+[0]


def encode_list(enc, xs):
    return [v for x in xs for v in ([1]+frame(enc(x)))]+[0]


def encoding_audit():
    checked = 0
    for n in range(4):
        for k in range(4):
            clauses = [c for w in range(min(k,2*n)+1) for c in combinations(range(2*n),w)]
            # Complete universe of finite CNFs with <=2 clauses over n<=3,
            # including empty clauses, both signs, and unused variables.
            for m in range(3):
                for formula in combinations(clauses,m):
                    f = sorted(formula)
                    bits = encode_nat(n)+encode_list(lambda c: encode_list(encode_nat,c),f)
                    exact = n+2+4*m+sum(4*c+8 for clause in f for c in clause)
                    upper = n+2+(4+(8*n+8)*k)*m
                    require(len(bits)==exact, 'closed-form encoding length')
                    require(len(bits)<=upper, 'source encoding bound')
                    checked += 1
    return {'universe':'n=0..3, k=0..3, every distinct-clause CNF with m<=2',
            'formulas_checked':checked,
            'exact_formula':'L=n+2+4m+sum_{literal codes c}(4c+8)',
            'source_upper_bound':'L<=n+2+[4+(8n+8)k]m'}


def interpolation_grid_audit():
    exponents = [(a,b) for a in range(3) for b in range(3)]
    full = [[Q(x**a*y**b) for a,b in exponents]
            for x in range(3) for y in range(3)]
    small = [[Q(x**a*y**b) for a,b in exponents]
             for x in range(2) for y in range(2)]
    require(rank(full)==9, 'complete coefficient interpolation grid')
    require(rank(small)==4, 'undersized grid rank')
    require(all(x*(x-1)==0 for x in range(2)) and 2*(2-1)!=0,
            'undersized-grid hostile')
    return {'variables':2,'per_variable_degree_bound':2,
            'full_grid':'{0,1,2}^2','full_grid_rank':9,
            'undersized_grid':'{0,1}^2','undersized_grid_rank':4,
            'hostile_polynomial':'x(x-1)'}


def main():
    # The target lower bound concerns a separate problem, not the SAT task
    # directly contradicted by the asserted SAT algorithm.
    logical_hostiles = [dict(SAT_algorithm=a,ETH=h,separate_target_lower_bound=l)
                        for a in (False,True) for h in (False,True) for l in (False,True)
                        if (not a or not h) and (not h or l) and a and l]
    require(logical_hostiles == [dict(SAT_algorithm=True,ETH=False,separate_target_lower_bound=True)],
            'negating a hardness hypothesis does not negate its conditional conclusion')
    receipt = {'status':'FINITE-EXACT INTERFACE CHECKS PASSED',
               'not_established':'No replay of the complete ETH algorithm, Lean dependency build, or asymptotic runtime proof.',
               'cube_replacement':cube_audit(), 'encoding':encoding_audit(),
               'conditional_lower_bound_hostile':logical_hostiles,
               'Keller_interpolation_interface':interpolation_grid_audit()}
    print(json.dumps(receipt,indent=2,sort_keys=True))


if __name__ == '__main__':
    main()

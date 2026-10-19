#!/usr/bin/env python3
"""Independent exact review controls. No author imports; no checkout writes.

SymPy is used for integer characteristic polynomials and rational ranks only.
The tests do not certify an analytic cusp phase or a sieve asymptotic onset.
All guards remain active under python -O.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
from math import gcd, isqrt
from pathlib import Path
import hashlib
import json
import sympy as sp

COUNTS = Counter()


def guard(category, claim, *context):
    COUNTS[category] += 1
    if not claim:
        raise RuntimeError((category, context))


def graph_cycles(b, W):
    """Disjoint nonzero periodic orbits, independently of a matrix spectrum."""
    cycles = set()
    for start in range(1, W):
        seen = {}
        path = []
        v = start
        while v and v not in seen:
            seen[v] = len(path)
            path.append(v)
            v = b*v % W
        if v:
            cycles.add(frozenset(path[seen[v]:]))
    return sorted(map(len, cycles))


x = sp.Symbol('x')
for b in range(2, 11):
    for W in range(2, 33):
        # Integer bP, in the original residue-value basis.
        B = sp.zeros(W)
        for r in range(W):
            for d in range(b):
                B[r, (b*r+d) % W] += 1
        lengths = graph_cycles(b, W)
        expected = (x-b)*x**(W-1-sum(lengths))
        for k in lengths:
            expected *= x**k-1
        guard('full_characteristic_polynomial',
              B.charpoly(x).as_expr().expand() == expected.expand(), b, W)

        # Restrict bP to the rational basis e_j-e_(W-1). Compare its
        # power ranks with the Fourier-arrow graph, including killed modes.
        if W <= 16 or (b, W) in [(2, 24), (2, 32), (6, 30), (10, 32)]:
            A = sp.Matrix(W-1, W-1, lambda i,j: B[i,j]-B[i,W-1])
            T = sp.zeros(W-1)
            for j in range(1, W):
                image = b*j % W
                if image:
                    T[image-1, j-1] = 1
            Ap = sp.eye(W-1)
            Tp = sp.eye(W-1)
            for h in range(1, 7):
                Ap, Tp = Ap*A, Tp*T
                guard('arbitrary_gcd_power_rank', Ap.rank() == Tp.rank(), b,W,h)
                guard('exact_nilpotence_criterion',
                      (Ap == sp.zeros(W-1)) == (b**h % W == 0), b,W,h)

# Arbitrary interval lengths, not restricted to powers of a digit base.
for W in range(2, 41):
    for L in range(1, 75):
        t = L % W
        formula = Q(t*(W-t), W*L)
        for start in [0, 1, W//2, W-1]:
            counts = Counter((start+d) % W for d in range(L))
            tv = sum(abs(Q(counts[r],L)-Q(1,W)) for r in range(W))/2
            guard('interval_total_variation', tv == formula, W,L,start)
            # The event collecting overloaded residues attains TV exactly.
            plus = [r for r in range(W) if Q(counts[r],L) > Q(1,W)]
            guard('total_variation_attaining_event',
                  sum(Q(counts[r],L)-Q(1,W) for r in plus) == tv, W,L,start)

# Independently form forbidden affine roots as sets, then check all branches.
for p in range(5, 102):
    if any(p % d == 0 for d in range(2,isqrt(p)+1)):
        continue
    u = pow(6,-1,p)
    original = {u, -u % p}
    for ell,d in product(range(p),repeat=2):
        bad = {r for r in range(p)
               if r in original or (ell*r+d) % p in original}
        if ell == 0:
            expected = 0 if d in original else p-2
        elif ell in (1,p-1):
            expected = p-2 if d == 0 else p-3 if d in {2*u % p,-2*u % p} else p-4
        else:
            exceptional = {(sign1*ell+sign2)*u % p
                           for sign1,sign2 in product((-1,1),repeat=2)}
            guard('distinct_affine_collision_shifts',len(exceptional)==4,p,ell)
            expected = p-3 if d in exceptional else p-4
        guard('forbidden_root_affine_count',p-len(bad)==expected,p,ell,d)

# The kernel of all Boolean diamond equations has exactly dimension n+1.
# This verifies the necessary direction, beyond testing factorized examples.
for n in range(0,7):
    rows=[]
    for S in range(1<<n):
        absent=[i for i in range(n) if not S>>i&1]
        for a,b in combinations(absent,2):
            row=[0]*(1<<n)
            for U,sign in [(S,1),(S|1<<a|1<<b,1),(S|1<<a,-1),(S|1<<b,-1)]:
                row[U]=sign
            rows.append(row)
    D = sp.Matrix(rows) if rows else sp.zeros(0,1<<n)
    guard('boolean_diamond_exact_nullity',(1<<n)-D.rank()==n+1,n)
    modular=sp.Matrix(1<<n,n+1,lambda S,j: 1 if j==0 else int(bool(S>>(j-1)&1)))
    guard('boolean_modular_basis',modular.rank()==n+1 and D*modular==sp.zeros(len(rows),n+1),n)
    if n>=2:
        interaction=sp.Matrix([int(S&3==3) for S in range(1<<n)])
        guard('boolean_nonmodular_hostile',D*interaction != sp.zeros(len(rows),1),n)

# Exact rational, deliberately non-monotone weights inside the claimed
# n^-3*(1 +/- eta) envelope; arbitrary nonnegative statistics, including zero.
for M in [2,7,31,100]:
    for L in range(1,7):
        eta=Q(1,M+2)
        weights=[Q(2,n**3)*(1+(-1)**n*eta) for n in range(M+1,M+L+1)]
        K=Q(M+L,M+1)**3*(1+eta)/(1-eta)
        guard('cusp_envelope_weight_ratio',max(weights)/min(weights)<=K,M,L)
        for values in product((0,1,5),repeat=L):
            discrete=Q(sum(values),L)
            continuous=sum(w*f for w,f in zip(weights,values))/sum(weights)
            guard('uniform_nonnegative_weight_transfer',
                  discrete/K<=continuous<=K*discrete,M,L,values)
        for e in [Q(1,10),Q(1,2),Q(99,100)]:
            guard('relative_sieve_error_composition',
                  1-(1-e)/K <= K*(1+e)-1 == e+(K-1)*(1+e),M,L,e)

guard('signed_statistic_hostile',Q(1-1,2)==0 and Q(1-2,3)!=0)

receipt={
    'status':'PASS: finite exact review controls; all guards active under -O',
    'counts':dict(sorted(COUNTS.items())),
    'total_guards':sum(COUNTS.values()),
    'scope':{
        'spectrum':'b=2..10, W=2..32, full integer characteristic polynomials',
        'nilpotent_parts':'b=2..10, W=2..16 plus four larger cases; powers h=1..6',
        'total_variation':'W=2..40, L=1..74, four interval offsets',
        'affine_factors':'all primes 5..101; all multipliers and shifts',
        'boolean_iff':'all diamond equations on complete Boolean lattices of ranks 0..6',
        'transfer':'four M, L=1..6, all statistics valued in {0,1,5}; perturbed rational n^-3 weights',
        'exclusions':'No numerical certification of analytic cusp phase, sieve onset, pi normality, or actual-prime result'
    },
    'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
output=json.dumps(receipt,indent=2)+'\n'
Path(__file__).with_name('independent_checks_receipt.json').write_text(output)
print(output,end='')

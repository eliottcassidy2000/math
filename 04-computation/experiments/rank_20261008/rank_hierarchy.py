#!/usr/bin/env python3
"""Standard-form hierarchy for one-step Lamperti transience (independent configurations, prime d).
For an affine coupling pi(j) = a j + b != id, the scale-free covariance is D = 2*I_nf - A with I_nf = 1 on non-unit
positions not fixed by pi and A the adjacency (with multiplicity) of the cycle graph of pi on non-unit positions.
Claims: lambda_max(D) <= 4 (< 4 if pi has no even cycle inside the non-unit positions), tr D = 2 #(non-unit, non-fixed).
Hence Q = I balances every coupling when: rank >= 6 (any coupling group); rank >= 5 (odd-order group: no even cycles);
rank >= 4 (translations: no fixed point).  Check EXACTLY (characteristic-polynomial-free: lambda_max < tr/2 iff
(tr/2) I - D positive definite, Sylvester) for every unit-position set of the given rank and every affine coupling.
Also the contracting examples."""
import itertools, math
from fractions import Fraction as Fr
from balanced_lamperti import posdef_exact
def D_of(d, units, pi):
    nonunit = [p for p in range(d) if p not in units]; idx = {p: k for k, p in enumerate(nonunit)}; r = len(nonunit)
    D = [[Fr(0)] * r for _ in range(r)]
    for j in range(d):
        w = [Fr(0)] * r
        if pi[j] in idx: w[idx[pi[j]]] += 1
        if j in idx: w[idx[j]] -= 1
        for p in range(r):
            if w[p]:
                for q in range(r): D[p][q] += w[p] * w[q]
    return D
def balanced_standard(D):
    r = len(D); t = sum(D[i][i] for i in range(r))
    if t == 0: return True                       # zero step law (identity coupling) never matters
    S = [[(t / 2 if i == j else 0) - D[i][j] for j in range(r)] for i in range(r)]
    return posdef_exact(S)
for d in (7, 11, 13):
    units_all = list(range(d))
    groups = {}
    for g in range(1, d):
        H = {1}; x = g
        while x != 1: H.add(x); x = x * g % d
        groups[len(H)] = sorted(H) if len(H) not in groups else groups[len(H)]
    for rho in range(3, min(d, 9)):
        res = {}
        for name, G in (('translations', [1]), ('odd-order', max((H for o, H in groups.items() if o % 2 == 1), key=len)), ('all units', list(range(1, d)))):
            ok = True; cnt = 0
            for Z in itertools.combinations(range(d), d - rho):
                if 0 not in Z and d - rho > 0: continue        # translate so that 0 is a unit position (AGL symmetry)
                for a in G:
                    for b in range(d):
                        if a == 1 and b == 0: continue
                        pi = [(a * j + b) % d for j in range(d)]
                        cnt += 1
                        if not balanced_standard(D_of(d, set(Z), pi)): ok = False; break
                    if not ok: break
                if not ok: break
            res[name] = ok
        print(f"d = {d}, rank {rho}: standard form balances every coupling? translations {res['translations']}, "
              f"odd-order group {res['odd-order']}, all affine couplings {res['all units']}", flush=True)
for d, m in ((7, [1, 2, 3, 5, 11, 13, 17]), (7, [1, 1, 2, 11, 23, 29, 37])):
    print(f"example Z_{d} m = {m}: residues {[x % d for x in m]}, Lambda = {sum(math.log(x / d) for x in m) / d:+.4f}, r_i = {[(-m[i] * i) % d for i in range(d)]}")

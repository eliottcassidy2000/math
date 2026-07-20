#!/usr/bin/env python3
"""jacobian_dixmier_probe_klein_S323.py -- klein-2026-07-19-S323.

TWO PROBES on the verified JC(3) counterexample
F = ((1+xy)^3 z + y^2(1+xy)(4+3xy), y + 3x(1+xy)^2 z + 3xy^2(4+3xy),
     2x - 3x^2y - x^3 z),  det J = -2.

(P1) THE DIXMIER-SIDE QUANTIZATION TEST. With G := (J^{-1})^T = adj(J)^T/(-2)
     (POLYNOMIAL entries because det is a nonzero constant), set
        X_i := F_i (multiplication operators),  D_i := sum_k G_{ki} d/dx_k.
     Then [X_i,X_j] = 0 and [D_i,X_j] = (G^T J)_{ij} = delta_ij automatically.
     The ENTIRE obstruction to "F induces an endomorphism of the Weyl algebra
     A_3" is [D_i,D_j] = 0.  For invertible F the D_i are the pushed-forward
     coordinate vector fields and commute.  QUESTION: do they commute for the
     NON-invertible Keller map?
       - If YES: x_i -> X_i, d_i -> D_i is an honest endomorphism of A_3 --
         an explicit candidate counterexample to the Dixmier Conjecture DC(3)
         (its non-surjectivity would then need a separate argument).
       - If NO: the "curvature of the fake inverse" tensor
         C^l_{ij} = sum_k (G_{ki} dG_{lj}/dx_k - G_{kj} dG_{li}/dx_k)
         is a nonzero polynomial invariant of the counterexample.
     CONTROL: the same computation on a genuine automorphism (a triangular
     map and the Nagata automorphism) must give all-zero commutators.

(P2) MOD-p FIBER CENSUS (monodromy measurement, Chebotarev-style).  For primes
     p (p != 2, det unit), tabulate the fiber-size profile of F over F_p^3.
     Degree-3 cover predictions for the Galois closure's monodromy group:
        S_3:   fractions with (3,1,0) preimages ~ (1/6, 1/2, 1/3)
        Z/3:   fractions ~ (1/3, 0, 2/3)
     The observed profile measures the group (the fleet's per-prime discipline
     applied to affine geometry).  Also records image density ~ 2/3 vs 1/3.
"""
import sympy as sp
import itertools, sys
from collections import Counter

x, y, z = X = sp.symbols('x y z')

def commutator_probe(Fs, name):
    """D_i := sum_k G[k,i] d_k with G := J^{-1} (COLUMN i of the inverse:
    D_i is the pushforward of d/du_i, since dx_k/du_i = (J^{-1})_{ki}).
    Then [D_i, X_j] = sum_k G[k,i] dF_j/dx_k = (J G)_{ji} = delta_ij.
    (S323 fix: an earlier draft used the transpose -- row i -- and the
    automorphism CONTROLS caught it by failing; controls are load-bearing.)"""
    J = sp.Matrix([[sp.diff(f, v) for v in X] for f in Fs])
    det = sp.simplify(J.det())
    G = J.inv()
    G = G.applyfunc(sp.expand)
    print(f"-- {name}: det = {det}")
    allzero = True
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        for l in range(3):
            c = sp.expand(sum(G[k, i]*sp.diff(G[l, j], X[k])
                              - G[k, j]*sp.diff(G[l, i], X[k]) for k in range(3)))
            c = sp.simplify(c)
            if c != 0:
                allzero = False
                s = str(c)
                print(f"   [D_{i+1},D_{j+1}] coeff of d_{l+1}: NONZERO "
                      f"({s[:120]}{'...' if len(s) > 120 else ''})")
    print(f"   => commutators {'ALL ZERO -- honest Weyl endomorphism' if allzero else 'DO NOT VANISH'}")
    return allzero

w = 1 + x*y
F = [w**3*z + y**2*w*(4+3*x*y),
     y + 3*x*w**2*z + 3*x*y**2*(4+3*x*y),
     2*x - 3*x**2*y - x**3*z]

print("== (P1) quantization commutator test ==")
# controls first (must be ALL ZERO)
commutator_probe([x + (y + z**2)**3, y + z**2 + 1, z], "control: triangular automorphism")
nag = [x - 2*y*(y**2 + x*z) - z*(y**2 + x*z)**2, y + z*(y**2 + x*z), z]
commutator_probe(nag, "control: Nagata automorphism")
ok = commutator_probe(F, "THE COUNTEREXAMPLE F")

print("\n== (P2) mod-p fiber census ==")
Fl = [sp.lambdify((x, y, z), sp.expand(f), 'math') for f in F]
Fpoly = [sp.Poly(sp.expand(f), x, y, z) for f in F]
print(f"{'p':>3} | {'N3':>7} {'N2':>7} {'N1':>7} {'N0':>7} | frac3  frac1  frac0 | image/p^3")
for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
    cnt = Counter()
    # evaluate F over F_p^3 with integer arithmetic mod p
    coeffs = [[(m, int(c) % p) for m, c in fp.terms()] for fp in Fpoly]
    for a in range(p):
        for b in range(p):
            # precompute powers
            for c0 in range(p):
                pt = []
                for terms in coeffs:
                    s = 0
                    for (e1, e2, e3), co in terms:
                        s += co * pow(a, e1, p) * pow(b, e2, p) * pow(c0, e3, p)
                    pt.append(s % p)
                cnt[tuple(pt)] += 1
    sizes = Counter(cnt.values())
    n3, n2, n1 = sizes.get(3, 0), sizes.get(2, 0), sizes.get(1, 0)
    img = sum(sizes.values())
    n0 = p**3 - img
    tot = p**3
    print(f"{p:>3} | {n3:>7} {n2:>7} {n1:>7} {n0:>7} | "
          f"{n3/tot:.3f}  {n1/tot:.3f}  {n0/tot:.3f} | {img/tot:.4f}")
print("S_3 Chebotarev prediction: frac3, frac1, frac0 -> 1/6=.167, 1/2=.500, 1/3=.333 (image 2/3=.667)")
print("Z/3 prediction:            frac3, frac1, frac0 -> 1/3=.333, 0,        2/3=.667 (image 1/3=.333)")

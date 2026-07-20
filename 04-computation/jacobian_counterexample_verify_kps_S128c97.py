#!/usr/bin/env python3
"""
kind-pasteur-2026-07-19-S128c97 -- HYP-8070: independent verification of the
owner-reported counterexample to the Jacobian Conjecture, and its arithmetic
fingerprints.

The map F: C^3 -> C^3:
  F1 = (1+xy)^3 z + y^2 (1+xy)(4+3xy)
  F2 = y + 3x(1+xy)^2 z + 3x y^2 (4+3xy)
  F3 = 2x - 3x^2 y - x^3 z
Claims: det J(F) == -2 identically; F sends (0,0,-1/4), (1,-3/2,13/2),
(-1,3/2,13/2) all to (-1/4, 0, 0).  Constant nonzero Jacobian + non-injective
= counterexample to JC_3 (hence JC_n false for n >= 3 by adding identity
coordinates, hence Dixmier DC_n false for n >= 3 via not-JC_n => not-DC_n).

Battery:
 (1) det J symbolically (sympy) -- THE decisive check.
 (2) exact collision arithmetic (Fractions).
 (3) equivariance F.sigma = tau.F for sigma=(-x,-y,z), tau=(a,-b,-c)
     (the Redei-shaped 1+2 fiber structure).
 (4) generic fiber degree via Groebner at a random rational target.
 (5) mod-p fiber statistics: the (1+xy)^3 / x^3 cube structure predicts
     cubic-residue stratification (p == 1 vs 2 mod 3) -- our gate-law lens.
"""
import time
from fractions import Fraction as Fr

def Fmap(x, y, z):
    u = 1 + x*y
    return (u**3*z + y**2*u*(4 + 3*x*y),
            y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y),
            2*x - 3*x**2*y - x**3*z)

def part1():
    print("== (1) det J symbolically ==", flush=True)
    import sympy as sp
    x, y, z = sp.symbols('x y z')
    F = sp.Matrix(Fmap(x, y, z))
    J = F.jacobian([x, y, z])
    d = sp.expand(J.det())
    print(f"  det J = {d}", flush=True)
    print(f"  VERDICT: {'CONSTANT -2 -- KELLER MAP CONFIRMED' if d == -2 else 'NOT constant -2 (!!)'}", flush=True)
    return d == -2

def part2():
    print("\n== (2) exact collision arithmetic ==", flush=True)
    pts = [(Fr(0), Fr(0), Fr(-1,4)), (Fr(1), Fr(-3,2), Fr(13,2)), (Fr(-1), Fr(3,2), Fr(13,2))]
    imgs = [Fmap(*p) for p in pts]
    for p, im in zip(pts, imgs):
        print(f"  F{tuple(map(str, p))} = {tuple(map(str, im))}", flush=True)
    ok = imgs[0] == imgs[1] == imgs[2] == (Fr(-1,4), Fr(0), Fr(0))
    print(f"  VERDICT: {'THREE DISTINCT POINTS, ONE IMAGE -- NON-INJECTIVITY CONFIRMED' if ok else 'MISMATCH (!!)'}", flush=True)
    return ok

def part3():
    print("\n== (3) equivariance (the Redei-shaped fiber) ==", flush=True)
    import sympy as sp
    x, y, z = sp.symbols('x y z')
    F = Fmap(x, y, z)
    Fs = Fmap(-x, -y, z)
    eq = (sp.expand(Fs[0] - F[0]) == 0 and sp.expand(Fs[1] + F[1]) == 0 and sp.expand(Fs[2] + F[2]) == 0)
    print(f"  F(-x,-y,z) = (F1, -F2, -F3): {'CONFIRMED' if eq else 'NO'}", flush=True)
    print("  => fibers over the tau-fixed plane (b=c=0) carry the sigma-involution;", flush=True)
    print("     the collision fiber = 1 sigma-fixed point + 1 swapped pair = 3 (odd = fixed mod 2: the parity lens)", flush=True)
    return eq

def part4():
    print("\n== (4) generic fiber degree (Groebner at a random rational target) ==", flush=True)
    import sympy as sp
    x, y, z = sp.symbols('x y z')
    F = Fmap(x, y, z)
    tgt = (Fr(5,7), Fr(2,3), Fr(-1,5))
    eqs = [sp.nsimplify(F[i] - tgt[i]) for i in range(3)]
    t0 = time.time()
    G = sp.groebner(eqs, x, y, z, order='lex', domain='QQ')
    # count solutions with multiplicity = dim of quotient ring
    from sympy.polys.orderings import lex
    lead = [g.LM(order=lex) for g in G.polys]
    # quotient dimension via standard monomials under the leading ideal
    import itertools
    maxd = 12
    cnt = 0
    for a in range(maxd):
        for b in range(maxd):
            for c in range(maxd):
                m = x**a*y**b*z**c
                if not any(sp.reduced(m, [sp.Poly(l, x, y, z) for l in lead])[1] != m for l in [None]):
                    pass
    # simpler: solve numerically instead
    sols = sp.solve(eqs, [x, y, z], dict=True)
    print(f"  target {tuple(map(str, tgt))}: {len(sols)} solutions over C  [{time.time()-t0:.0f}s]", flush=True)
    print(f"  generic fiber size (geometric degree) = {len(sols)}", flush=True)
    return len(sols)

def part5():
    print("\n== (5) mod-p fiber statistics (the cubic-residue prediction) ==", flush=True)
    for p in (5, 7, 11, 13, 17, 19, 31):
        from collections import Counter
        cnt = Counter()
        for X in range(p):
            for Y in range(p):
                for Z in range(p):
                    im = tuple(v % p for v in Fmap(X, Y, Z))
                    cnt[im] += 1
        fibs = Counter(cnt.values())
        missed = p**3 - len(cnt)
        tag = "p==1 mod 3" if p % 3 == 1 else "p==2 mod 3"
        print(f"  p={p} ({tag}): fiber-size histogram {dict(sorted(fibs.items()))}, image misses {missed}", flush=True)
    print("  (prediction: cube structure => 1/3-type splitting at p==1 mod 3; cleaner at p==2 mod 3)", flush=True)

if __name__ == "__main__":
    ok1 = part1()
    ok2 = part2()
    ok3 = part3()
    if ok1 and ok2:
        print("\n*** JC_3 COUNTEREXAMPLE INDEPENDENTLY VERIFIED (det + collisions) ***", flush=True)
        print("*** => JC_n FALSE for all n >= 3;  DC_n (Dixmier) FALSE for all n >= 3 via not-JC_n => not-DC_n ***", flush=True)
        print("*** surviving open: JC_2 and DC_1 -- the bottom of the doubling tower ***", flush=True)
    try:
        part4()
    except Exception as e:
        print(f"  (part 4 skipped: {e})", flush=True)
    part5()

#!/usr/bin/env python3
"""
kind-pasteur-2026-07-19-S128c98 -- THM-1310 cross-validation: the quadratic
resolvent of the fiber cubic is Q(sqrt(-L)), L = 27a^2g^2 - 18abg + 16a + b^3g - b^2.
Pointwise prediction mod p (off the walls L=0, Q=0, Delta=0):
  fiber size 1  <=> Frobenius = transposition <=> -L NON-residue mod p
  fiber size 3  <=> Frobenius = identity      <=> -L residue
  fiber size 0  <=> Frobenius = 3-cycle       <=> -L residue
Also: w-family triviality check F_w == F(x,y,z+w).
"""
import sympy as sp
from collections import Counter

x, y, z, w = sp.symbols('x y z w')
u = 1 + x*y
F = ((u**3*z + y**2*u*(4+3*x*y)), (y + 3*x*u**2*z + 3*x*y**2*(4+3*x*y)), (2*x - 3*x**2*y - x**3*z))

print("== w-family triviality ==", flush=True)
Wp = y**2*(4+3*x*y) + w*u**2
B1 = sp.expand(u*Wp); B2 = sp.expand(y + 3*x*Wp); B3 = sp.expand(x*(2-3*x*y) - w*x**3)
Fw = (sp.expand(u**3*z + B1), sp.expand(3*u**2*x*z + B2), sp.expand(-x**3*z + B3))
Fshift = tuple(sp.expand(c.subs(z, z + w)) for c in F)
print(f"  F_w == F(x,y,z+w): {all(sp.expand(Fw[i]-Fshift[i])==0 for i in range(3))}  (the solver's 1-param family = trivial z-translation)", flush=True)

print("\n== pointwise resolvent test: fiber size vs Legendre(-L) ==", flush=True)
def Fmod(X, Y, Z, p):
    U = (1 + X*Y) % p
    return ((U**3*Z + Y*Y*U*(4+3*X*Y)) % p, (Y + 3*X*U*U*Z + 3*X*Y*Y*(4+3*X*Y)) % p, (2*X - 3*X*X*Y - X**3*Z) % p)
def legendre(t, p):
    t %= p
    if t == 0: return 0
    return 1 if pow(t, (p-1)//2, p) == 1 else -1
for p in (7, 13, 19):
    fibers = Counter()
    for X in range(p):
        for Y in range(p):
            for Z in range(p):
                fibers[Fmod(X, Y, Z, p)] += 1
    table = Counter()
    QRs = {0: Counter(), 1: Counter(), -1: Counter()}
    bad = 0
    for A in range(p):
        for Bt in range(p):
            for G in range(p):
                sz = fibers.get((A, Bt, G), 0)
                L = (27*A*A*G*G - 18*A*Bt*G + 16*A + Bt**3*G - Bt*Bt) % p
                Q = (27*A*G*G - 9*Bt*G + 8) % p
                chi = legendre(-L, p)
                if L == 0 or Q == 0:
                    bad += 1
                    continue
                QRs[chi][sz] += 1
                ok = (sz == 1 and chi == -1) or (sz in (0, 3) and chi == 1)
                table[(sz, chi, ok)] += 1
    viol = sum(v for (szc, chic, okc), v in table.items() if not okc)
    tot = sum(table.values())
    print(f"  p={p}: off-wall targets {tot}, wall targets {bad}", flush=True)
    print(f"    chi=-1 (nonres): sizes {dict(sorted(QRs[-1].items()))}   [predict all size 1]", flush=True)
    print(f"    chi=+1 (res):    sizes {dict(sorted(QRs[1].items()))}   [predict sizes 0 and 3]", flush=True)
    print(f"    VIOLATIONS: {viol}/{tot}", flush=True)
print("\nDONE.", flush=True)

#!/usr/bin/env python3
"""
death-star-2026-07-19-S59n (HYP-8080) -- the rationals DECODE: every constant in
the collision data is derived from the six-function form of the counterexample.

Six-function form (t = xy, s = x^2 z the C*-invariants):
  F1 = z*A(t) + y^2*B(t),  F2 = y*C(t) + x*z*D(t),  F3 = x*(E0(t) - s)
with A = u^3, B = u(4+3t), C = 1+12t+9t^2, D = 3u^2, E0 = 2-3t,  u = 1+t.

Decode identities (all verified exactly here):
  (I)   Six-function form reproduces F exactly.
  (II)  Phi(t) := t*C(t) + E0(t)*D(t) = 4t + 6   -- LINEAR (cubic+quadratic
        terms cancel); its unique root t* = -3/2 is the collision parameter.
  (III) Psi(t) := E0(t)*A(t) + t^2*B(t) = (1+t)(2+t); Psi(t*) = -1/4 = the
        collision value a*.
  (IV)  The units cross at t*: u(t*) = (4+3t)(t*) = -1/2.
  (V)   E0(t*) = 13/2 = s* -- the invariant coordinates of the collision orbit
        are (t*, s*) = (-3/2, 13/2); z* = s* at x = 1.
  (VI)  3AD' = 2A'D (the cube-forcing relation; c2 = 0 in the det s-grading).
"""
from fractions import Fraction as Fr

def pmul(a, b):
    r = {}
    for ka, ca in a.items():
        for kb, cb in b.items():
            k = tuple(p+q for p, q in zip(ka, kb))
            r[k] = r.get(k, 0) + ca*cb
    return {k: c for k, c in r.items() if c != 0}
def padd(*ps):
    r = {}
    for p in ps:
        for k, c in p.items():
            r[k] = r.get(k, 0) + c
    return {k: c for k, c in r.items() if c != 0}
def psc(p, s): return {k: c*s for k, c in p.items() if c*s != 0}

# --- ring Q[x,y,z]
X = {(1,0,0):1}; Y = {(0,1,0):1}; Z = {(0,0,1):1}; ONE3 = {(0,0,0):1}
U3 = padd(ONE3, pmul(X,Y))
W3 = padd(psc(ONE3,4), psc(pmul(X,Y),3))
F1 = padd(pmul(pmul(pmul(U3,U3),U3), Z), pmul(pmul(pmul(Y,Y),U3), W3))
F2 = padd(Y, psc(pmul(pmul(X,pmul(U3,U3)),Z),3), psc(pmul(pmul(X,pmul(Y,Y)),W3),3))
F3 = padd(psc(X,2), psc(pmul(pmul(X,X),Y),-3), psc(pmul(pmul(pmul(X,X),X),Z),-1))

# --- ring Q[t] (univariate as 1-tuples)
T = {(1,):1}; ONE = {(0,):1}
u  = padd(ONE, T)
A  = pmul(pmul(u,u),u)
B  = pmul(u, padd(psc(ONE,4), psc(T,3)))
C  = padd(ONE, psc(T,12), psc(pmul(T,T),9))
D  = psc(pmul(u,u),3)
E0 = padd(psc(ONE,2), psc(T,-3))

def subst_ts(p_t, extra_xyz=None):
    """map Q[t] poly -> Q[x,y,z] via t = xy, then multiply by extra."""
    out = {}
    for (k,), c in p_t.items():
        out[(k, k, 0)] = out.get((k,k,0), 0) + c
    return out if extra_xyz is None else pmul(out, extra_xyz)

# (I) six-function reconstruction
F1r = padd(subst_ts(A, Z), subst_ts(B, pmul(Y,Y)))
F2r = padd(subst_ts(C, Y), subst_ts(D, pmul(X,Z)))
F3r = padd(subst_ts(E0, X), psc(pmul(pmul(pmul(X,X),X),Z), -1))
print("(I) six-function form == F exactly:",
      F1r == F1 and F2r == F2 and F3r == F3)

# (II) Phi = tC + E0 D
Phi = padd(pmul(T, C), pmul(E0, D))
print("(II) Phi(t) = tC + E0*D =", dict(sorted(Phi.items())), " == 4t+6:",
      Phi == {(0,):6, (1,):4})
tstar = Fr(-3, 2)
def ev(p, v):
    return sum(Fr(c) * v**k for (k,), c in p.items())
print("     root t* = -3/2:", ev(Phi, tstar) == 0)

# (III) Psi = E0 A + t^2 B
Psi = padd(pmul(E0, A), pmul(pmul(T,T), B))
target = pmul(padd(ONE,T), padd(psc(ONE,2), T))
print("(III) Psi(t) = E0*A + t^2*B =", dict(sorted(Psi.items())),
      " == (1+t)(2+t):", Psi == target, "; Psi(t*) =", ev(Psi, tstar), "== -1/4:",
      ev(Psi, tstar) == Fr(-1,4))

# (IV) unit crossing
w = padd(psc(ONE,4), psc(T,3))
print("(IV) u(t*) =", ev(u, tstar), "; (4+3t)(t*) =", ev(w, tstar),
      "; crossing (both -1/2):", ev(u,tstar) == ev(w,tstar) == Fr(-1,2))

# (V) E0(t*) = 13/2
print("(V) E0(t*) =", ev(E0, tstar), "== 13/2 == z* == s*:", ev(E0,tstar) == Fr(13,2))

# (VI) cube-forcing
def pdiff(p):
    return {(k-1,): c*k for (k,), c in p.items() if k > 0}
lhs = psc(pmul(A, pdiff(D)), 3)
rhs = psc(pmul(pdiff(A), D), 2)
print("(VI) 3AD' == 2A'D:", lhs == rhs)

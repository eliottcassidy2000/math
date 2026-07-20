"""opus-2026-07-20-S417 -- THE LAST TNC QUESTION: CAN EQUAL-MODULUS SADDLES CANCEL?

From THM-1615: CT(Lambda^m) = [u^{Nm}]R^m = (1/2 pi i) contour e^{m g(u)} du/u with
      g(u) = log R(u) - N log u ,      g'(u) = 0  <=>  u R'(u) = N R(u).
At nondegenerate saddles u_j the standard expansion gives
      CT(Lambda^m)  ~  sum_j  c_j * w_j^m / sqrt(m),
      w_j = R(u_j)/u_j^N = e^{g(u_j)},     c_j = 1 / ( u_j sqrt(2 pi g''(u_j)) ).
THM-1615 proved a genuine saddle ALWAYS exists, so the sum is nonempty.

THE VANDERMONDE LEMMA (this script's point).  Suppose the DOMINANT saddle values w_j
(those of maximal modulus rho) are DISTINCT.  Vanishing for all m >= 1 means
      sum_j c_j w_j^m = 0   for m = 1, 2, ..., k .
That linear system has matrix M_{mj} = w_j^m = w_j * w_j^{m-1}, i.e. M = V diag(w_j) with V
Vandermonde, so
      det M = prod_{i<j} (w_j - w_i) * prod_j w_j  !=  0
when the w_j are distinct and nonzero.  Hence ALL c_j = 0.  But
      c_j = 1/( u_j sqrt(2 pi g''(u_j)) )  !=  0
for any nondegenerate saddle (u_j != 0, g'' finite and nonzero).  CONTRADICTION.

  ==> **IF THE DOMINANT SADDLE VALUES ARE DISTINCT, TNC HOLDS.**

So the ONLY escape is several dominant saddles sharing ONE value w, with prefactors summing
to zero.  That is a FINITE ALGEBRAIC condition -- vastly sharper than 'the ladder'.  This
script (1) proves the lemma numerically, (2) measures how often dominant values COLLIDE, and
(3) tests whether the prefactor sum can actually vanish on a collision.
"""
import sympy as sp
import numpy as np
from itertools import product

u = sp.symbols('u')

def saddle_data(Rexpr, N):
    """genuine saddles, values w = R/u^N, and prefactors c = 1/(u sqrt(2 pi g''))"""
    R = sp.expand(Rexpr)
    S = sp.expand(u*sp.diff(R, u) - N*R)
    roots = sp.nroots(sp.Poly(S, u)) if sp.Poly(S, u).degree() > 0 else []
    g2 = sp.diff(sp.log(R) - N*sp.log(u), u, 2)
    out = []
    for r in roots:
        rv = complex(r)
        if abs(complex(R.subs(u, r))) < 1e-12: continue      # spurious (root of R)
        if abs(rv) < 1e-12: continue
        w = complex(R.subs(u, r))/rv**N
        try: gg = complex(g2.subs(u, r))
        except Exception: continue
        if abs(gg) < 1e-12: continue                          # degenerate saddle
        c = 1.0/(rv*np.sqrt(2*np.pi*gg))
        out.append((rv, w, c))
    return out

def CT(Rexpr, N, m):
    e = sp.expand(Rexpr**m)
    return sp.Poly(e, u).coeff_monomial(u**(N*m)) if e != 0 else 0

print("="*78)
print("(1) THE VANDERMONDE LEMMA, numerically")
print("="*78)
ws = [1.0+0j, -1.0+0j, 1j]
V = np.array([[w**m for w in ws] for m in (1, 2, 3)])
print(f"   distinct w = {ws}; matrix [w_j^m]_{{m=1..3}} has det = {np.linalg.det(V):.6g}")
print(f"   nonsingular => sum_j c_j w_j^m = 0 for m=1..3 forces c = 0.")
print("   And c_j = 1/(u_j sqrt(2 pi g''(u_j))) is NEVER 0 at a nondegenerate saddle.")
print("   => DISTINCT dominant values  ==>  TNC HOLDS.")

print()
print("="*78)
print("(2) DO DOMINANT SADDLE VALUES COLLIDE?  sweep N=2, quartic R")
print("="*78)
coll = 0; tot = 0; examples = []
for c0, c1, c2, c3 in product([-2, -1, 1, 2], [-2, -1, 0, 1, 2], [-2, -1, 0, 1, 2], [-2, -1, 0, 1, 2]):
    R = c0 + c1*u + c2*u**2 + c3*u**3 + u**4
    d = saddle_data(R, 2)
    if not d: continue
    tot += 1
    rho = max(abs(w) for _, w, _ in d)
    dom = [(uu, w, cc) for uu, w, cc in d if abs(w) > rho*(1-1e-9)]
    vals = [w for _, w, _ in dom]
    # collision = two dominant saddles with (numerically) equal value
    hit = any(abs(vals[i]-vals[j]) < 1e-8 for i in range(len(vals)) for j in range(i+1, len(vals)))
    if hit:
        coll += 1
        if len(examples) < 4:
            s = sum(cc for _, _, cc in dom)
            examples.append((R, len(dom), abs(s)))
print(f"   quartics with saddle data: {tot}")
print(f"   with COLLIDING dominant values: {coll}  ({100.0*coll/max(tot,1):.1f}%)")
for R, nd, s in examples:
    print(f"      R={R}   #dominant={nd}   |sum of prefactors| = {s:.6g}")

print()
print("="*78)
print("(3) ON A COLLISION, CAN THE PREFACTORS SUM TO ZERO?")
print("="*78)
print("   If they can, CT(Lambda^m) loses its leading term and the next order decides.")
print("   Test: among all collision cases, find the smallest |sum of prefactors|, and")
print("   check the actual CT values there.")
best = None
for c0, c1, c2, c3 in product([-2, -1, 1, 2], [-2, -1, 0, 1, 2], [-2, -1, 0, 1, 2], [-2, -1, 0, 1, 2]):
    R = c0 + c1*u + c2*u**2 + c3*u**3 + u**4
    d = saddle_data(R, 2)
    if not d: continue
    rho = max(abs(w) for _, w, _ in d)
    dom = [(uu, w, cc) for uu, w, cc in d if abs(w) > rho*(1-1e-9)]
    if len(dom) < 2: continue
    vals = [w for _, w, _ in dom]
    if not any(abs(vals[i]-vals[j]) < 1e-8 for i in range(len(vals)) for j in range(i+1, len(vals))):
        continue
    s = abs(sum(cc for _, _, cc in dom))
    if best is None or s < best[1]: best = (R, s, dom)
if best:
    R, s, dom = best
    print(f"   smallest |prefactor sum| on a collision: {s:.6g}   at R = {R}")
    cts = [CT(R, 2, m) for m in range(1, 8)]
    print(f"   CT(Lambda^m), m=1..7 = {cts}")
    print(f"   all zero? {all(v == 0 for v in cts)}")
else:
    print("   no collision case with >= 2 dominant saddles found in this sweep.")

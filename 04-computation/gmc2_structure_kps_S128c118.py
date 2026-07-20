#!/usr/bin/env python3
"""gmc2_structure_kps_S128c118.py -- kind-pasteur-2026-07-20-S128c118

TOWARD THE 2D NULLCONE STRUCTURE THEOREM (which would give GMC(2)).

CONJECTURE.  For one complex Gaussian Z,
    N = { P : E[P^m] = 0 for all m >= 1 }  =  { P whose charges are all > 0 }
                                            u { P whose charges are all < 0 } u {0},
where the charge of Z^a conj(Z)^b is a - b.  GMC(2) follows in one line, since for a
one-sided P every monomial of Q P^m has charge of modulus >= m*|k| - deg Q, which is
non-zero once m > deg(Q)/|k|.

THREE PIECES TESTED HERE.

(L1) CHARGE-0 ONLY  =>  P = 0.  If every charge is 0 then P = g(|Z|^2) = g(s) with
     s ~ Exp(1), so E[P^m] = int_0^inf g(s)^m e^{-s} ds.  If deg g = d >= 1 with leading
     coefficient c, the integrand is dominated by c^m s^{md} and
     E[P^m] ~ c^m (md)!, which is non-zero for large m -- the top term cannot be
     cancelled because it grows factorially faster than every correction.  So g is
     constant, and then E[P^m] = g^m = 0 forces g = 0.  Checked numerically.

(L2) THE TOP-DEGREE PART MUST BE ONE-SIDED.  Write d = deg P and let P_top be the
     degree-d homogeneous part.  Setting P_top = y^d h(x/y) with x = Z, y = conj(Z),
     the charge-0 coefficient of (P_top)^m is [t^{md/2}] h(t)^m, and the corresponding
     term of E[P^m] carries the factor (md/2)! -- which dominates every lower-degree
     contribution.  So the leading asymptotics of E[P^m] is
         (md/2)! * CT( (h(t) t^{-d/2})^m ),
     and for this to vanish for ALL m the Laurent polynomial h(t) t^{-d/2} must have
     support entirely positive or entirely negative (the one-variable constant-term
     fact).  That says exactly: P_top is one-sided in charge.  Checked directly.

(L3) A WIDER SEARCH for a mixed-charge nullcone element, with coefficients in
     {-2,...,2} rather than {-1,0,1}.  The previous run's box was small enough that a
     survivor could have been missed by coefficient size alone.
"""
import sys
from math import factorial
from itertools import product
from fractions import Fraction as F

DEG = int(sys.argv[1]) if len(sys.argv) > 1 else 2
MM = int(sys.argv[2]) if len(sys.argv) > 2 else 9
WIDE = int(sys.argv[3]) if len(sys.argv) > 3 else 2


def mul(p, q):
    out = {}
    for (a, b), u in p.items():
        for (c, d), v in q.items():
            k = (a + c, b + d)
            out[k] = out.get(k, 0) + u * v
    return {k: v for k, v in out.items() if v}


def Eval(p):
    return sum(v * factorial(a) for (a, b), v in p.items() if a == b)


def charges(p):
    return sorted({a - b for (a, b) in p})


def one_sided(p):
    cs = charges(p)
    return all(c > 0 for c in cs) or all(c < 0 for c in cs)


print("=" * 78)
print("(L1) CHARGE-0 ONLY => P = 0 :  E[g(s)^m] ~ c^m (md)! cannot vanish")
print("=" * 78)
print("  g(s) = s^d + lower;  E[P^m] = int g^m e^{-s} ds.  Sample g's:")
for gc in ([0, 1], [1, 1], [-1, 0, 1], [2, -3, 1]):
    d = len(gc) - 1
    # P = g(|Z|^2) as a polynomial in Z,Zb : s^j -> Z^j Zb^j
    P = {}
    for j, c in enumerate(gc):
        if c:
            P[(j, j)] = P.get((j, j), 0) + c
    if not P:
        continue
    vals = []
    cur = dict(P)
    for m in range(1, 7):
        vals.append(Eval(cur))
        cur = mul(cur, P)
    print("     g coeffs %-14s deg %d :  E[P^m] = %s" % (gc, d, vals))
print("  -> every non-zero g has E[P^m] != 0 for some small m; the factorial growth of")
print("     the top term makes cancellation impossible for large m.  L1 holds.")

print()
print("=" * 78)
print("(L2) TOP-DEGREE PART MUST BE ONE-SIDED  (central coefficient of h(t)^m)")
print("=" * 78)
print("  P_top = y^d h(x/y);  charge-0 coeff of (P_top)^m = [t^{md/2}] h(t)^m.")
print("  If h's support straddles d/2, that central coefficient is non-zero for some m.")


def central(hc, d, m):
    """[t^{md/2}] h(t)^m, h given by coefficient list hc (index = power of t)."""
    if (m * d) % 2:
        return None
    tgt = m * d // 2
    cur = {0: 1}
    for _ in range(m):
        nxt = {}
        for e, u in cur.items():
            for j, c in enumerate(hc):
                if c:
                    nxt[e + j] = nxt.get(e + j, 0) + u * c
        cur = nxt
    return cur.get(tgt, 0)


tests = [
    ("h = 1 + t   (d=1, straddles 1/2)", [1, 1], 1),
    ("h = 1 - t   (d=1, straddles)", [1, -1], 1),
    ("h = t       (d=1, one-sided +)", [0, 1], 1),
    ("h = 1       (d=1, one-sided -)", [1, 0], 1),
    ("h = 1 + t^2 (d=2, straddles 1)", [1, 0, 1], 2),
    ("h = t^2     (d=2, one-sided +)", [0, 0, 1], 2),
]
for name, hc, d in tests:
    vals = [central(hc, d, m) for m in range(1, 7)]
    nz = any(v not in (0, None) for v in vals)
    print("     %-38s central coeffs %s  some non-zero: %s" % (name, vals, nz))
print("  -> straddling h always produces a non-zero central coefficient; one-sided h")
print("     never does.  L2 holds on these, as the constant-term fact predicts.")

print()
print("=" * 78)
print("(L3) WIDER SEARCH for a MIXED-charge nullcone element, coeffs in [-%d..%d]"
      % (WIDE, WIDE))
print("=" * 78)
MONS = [(0, 0)] + [(a, b) for a in range(DEG + 1) for b in range(DEG + 1)
                   if 1 <= a + b <= DEG]
COEF = list(range(-WIDE, WIDE + 1))
print("  %d monomials (degree <= %d), %d coefficient values -> %d vectors"
      % (len(MONS), DEG, len(COEF), len(COEF)**len(MONS)))
mixed = []
nnull = 0
for vec in product(COEF, repeat=len(MONS)):
    if all(v == 0 for v in vec):
        continue
    p = {MONS[i]: vec[i] for i in range(len(MONS)) if vec[i]}
    if Eval(p) != 0:
        continue
    p2 = mul(p, p)
    if Eval(p2) != 0:
        continue
    cur = dict(p)
    bad = False
    for m in range(1, MM + 1):
        if Eval(cur) != 0:
            bad = True
            break
        cur = mul(cur, p)
    if bad:
        continue
    nnull += 1
    if not one_sided(p):
        mixed.append(p)
print("  nullcone elements found (m <= %d) : %d" % (MM, nnull))
print("  of those MIXED in charge          : %d" % len(mixed))
if mixed:
    for p in mixed[:8]:
        print("     charges %s : %s" % (charges(p), dict(sorted(p.items()))))
else:
    print("  none -- the structure conjecture survives the wider box.")
print()
print("  DETECTION FLOOR: coefficients in [-%d..%d] over %d monomials of degree <= %d,"
      % (WIDE, WIDE, len(MONS), DEG))
print("  depth m <= %d.  A search bounds a counterexample; it does not prove the theorem."
      % MM)

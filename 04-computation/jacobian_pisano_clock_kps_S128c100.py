#!/usr/bin/env python3
"""
kind-pasteur-2026-07-20-S128c100 -- HYP-8145: THE PISANO CLOCK.
Owner: relate the JC counterexample thread to pi(10) = 60 (Fibonacci last digits).

 (1) Pisano facts + the GL_2 decode of 60:
     pi(2)=3, pi(5)=20, pi(10)=60; M = [[1,1],[1,0]] (det -1: the LINEAR
     unimodular 'Keller' map of the plane).  Mod 2: GL_2(F_2) = S_3 and M is an
     order-3 element (a 3-cycle IN THE COUNTEREXAMPLE'S MONODROMY GROUP).
     Mod 5 (the ramified prime of Q(sqrt5)): char poly = (x-3)^2, order
     4*5 = 20, and M is PROJECTIVELY an order-5 element of PGL_2(F_5), where
     PSL_2(F_5) = A_5 (|SL_2(F_5)| = 120).  60 = lcm(3-rotation dressed,
     5-rotation dressed) -- and A_5 = <3-rotation, 5-rotation> (icosahedron).
 (2) THE MECHANISM BRIDGE (exact): Lucas tripling L_{3n} = L_n^3 - 3(-1)^n L_n
     == Chebyshev T_3 composition (even n: L_{3n}/2 = T_3(L_n/2)) == THM-1335's
     trisection normal form 4T^3-3T = m.  Fibonacci side: F_{3n} = 5F_n^3 +
     3(-1)^n F_n (the U-side).  Iteration: T_3 o T_3 = T_9 <-> F o F degree 9:
     the counterexample's iteration ladder IS the index-tripling semigroup.
 (3) pi_F: THE PISANO PERIOD OF THE COUNTEREXAMPLE.  F mod n on (Z/n)^3:
     eventual core (image stabilization; F is a bijection on it), cycle lengths,
     pi_F(n) = lcm.  Table for n in {2,3,4,5,7,10,11,13,17,19,23,29,31};
     headline: pi_F(10), the map's own '60'.  CRT check pi_F(10) =
     lcm(pi_F(2), pi_F(5)).  Law hunt: pi_F(p) vs p-1, p+1 factored (Pisano's
     split/inert law analog) -- honest eyeball, no forced law.
"""
import sympy as sp
from sympy import symbols, expand, Rational, sqrt, simplify
from math import gcd
from collections import Counter
import time

print("== (1) Pisano facts and the GL_2 decode of 60 ==", flush=True)
def pisano(n):
    a, b, k = 1, 1, 0
    while True:
        a, b = b, (a + b) % n
        k += 1
        if (a, b) == (1, 1):
            return k
for n in (2, 5, 10):
    print(f"  pi({n}) = {pisano(n)}", flush=True)
M = sp.Matrix([[1, 1], [1, 0]])
print(f"  det M = {M.det()}  (unimodular: the linear plane 'Keller' map)", flush=True)
M2 = M.applyfunc(lambda t: t % 2)
P = M2
for k in range(1, 8):
    if P.applyfunc(lambda t: t % 2) == sp.eye(2):
        print(f"  order of M in GL_2(F_2) = {k};  |GL_2(F_2)| = 6, GL_2(F_2) iso S_3  => M mod 2 = a 3-CYCLE in S_3", flush=True)
        break
    P = (P * M2)
x_ = symbols('x')
diff5 = sp.Poly(sp.expand((x_-3)**2 - (x_**2 - x_ - 1)), x_)
print(f"  (x-3)^2 - (x^2-x-1) = {diff5.as_expr()}: all coeffs divisible by 5: {all(c % 5 == 0 for c in diff5.all_coeffs())}  => 5 RAMIFIED in Q(sqrt5), double eigenvalue 3", flush=True)
P = M
for k in range(1, 30):
    Pm = P.applyfunc(lambda t: t % 5)
    if Pm == sp.eye(2):
        print(f"  order of M in GL_2(F_5) = {k} = 4*5 (ord(3)=4 x ramification 5)", flush=True)
        break
    P = P * M
P = M
for k in range(1, 12):
    Pm = P.applyfunc(lambda t: t % 5)
    if Pm[0, 1] == 0 and Pm[1, 0] == 0 and Pm[0, 0] == Pm[1, 1]:
        print(f"  M^{k} = {Pm[0,0]}*I mod 5  => PROJECTIVE order of M in PGL_2(F_5) = {k} (an order-5 'pentagon' element)", flush=True)
        break
    P = P * M
print(f"  |SL_2(F_5)| = 5*(25-1) = 120, PSL_2(F_5) iso A_5 (order 60) -- klein's three sixties (T1534)", flush=True)
print(f"  60 = lcm(3, 20): a 3-cycle in S_3 (mod 2) and a dressed 5-cycle (mod 5); A_5 = <3-rot, 5-rot> (icosahedron)", flush=True)

print("\n== (2) the mechanism bridge: Lucas tripling == T_3 ==", flush=True)
nn = symbols('n', integer=True)
al = (1 + sqrt(5))/2; be = (1 - sqrt(5))/2
Ln = al**nn + be**nn
lhs = al**(3*nn) + be**(3*nn)
rhs = sp.expand(Ln**3 - 3*(al*be)**nn*Ln)
print(f"  L_(3n) - (L_n^3 - 3(-1)^n L_n) == 0: {sp.simplify(lhs - rhs) == 0}   (alpha*beta = -1)", flush=True)
w = symbols('w')
T3 = 4*w**3 - 3*w
print(f"  even n: L_(3n)/2 = T_3(L_n/2): T_3(w) = 4w^3-3w; L^3-3L = 2*T_3(L/2): {sp.expand(2*T3.subs(w, symbols('LL')/2) - (symbols('LL')**3 - 3*symbols('LL'))) == 0}", flush=True)
Fn = (al**nn - be**nn)/sqrt(5)
lhsF = (al**(3*nn) - be**(3*nn))/sqrt(5)
rhsF = sp.expand(5*Fn**3 + 3*(al*be)**nn*Fn)
print(f"  F_(3n) = 5F_n^3 + 3(-1)^n F_n: {sp.simplify(lhsF - rhsF) == 0}", flush=True)
print(f"  T_3(T_3(w)) == T_9(w) (index-tripling semigroup = the F -> F.F -> ... iteration ladder, deg 3^m): "
      f"{sp.expand(T3.subs(w, T3) - sp.chebyshevt(9, w)) == 0}", flush=True)
print("  => THM-1335's modulus law 4T^3-3T = m and the Lucas/Fibonacci tripling laws are THE SAME identity;", flush=True)
print("     F's fibers trisect the modulus angle exactly as index-tripling trisects Fibonacci/Lucas indices.", flush=True)

print("\n== (3) pi_F: the Pisano period of the counterexample ==", flush=True)
def Fmod(X, Y, Z, n):
    U = (1 + X*Y) % n
    return ((U*U*U % n)*Z % n and 0 or 0, )  # placeholder, replaced below
def Fm(X, Y, Z, n):
    U = (1 + X*Y) % n
    return ((U**3*Z + Y*Y*U*(4 + 3*X*Y)) % n,
            (Y + 3*X*U*U*Z + 3*X*Y*Y*(4 + 3*X*Y)) % n,
            (2*X - 3*X*X*Y - X**3*Z) % n)
def pisano_F(n):
    pts = [(X, Y, Z) for X in range(n) for Y in range(n) for Z in range(n)]
    idx = {p: i for i, p in enumerate(pts)}
    nxt = [idx[Fm(*p, n)] for p in pts]
    # eventual core: iterate image sets
    cur = set(range(len(pts)))
    while True:
        nxts = {nxt[i] for i in cur}
        if nxts == cur:
            break
        cur = nxts
    core = cur
    # cycle lengths on core (function restricted to core is a bijection)
    seen = set()
    lens = Counter()
    for i0 in core:
        if i0 in seen: continue
        path = []
        i = i0
        while i not in seen:
            seen.add(i); path.append(i); i = nxt[i]
        # all of path is on the same cycle (core = union of cycles)
        # find cycle length: walk from i0 until return
        j, l = nxt[i0], 1
        while j != i0:
            j = nxt[j]; l += 1
        lens[l] += 1
        # mark rest of that cycle seen already via path
    from math import lcm as _lcm
    piF = 1
    for l in lens:
        piF = _lcm(piF, l)
    return len(core), dict(sorted(lens.items())), piF

results = {}
for n in (2, 3, 4, 5, 7, 10, 11, 13, 17, 19, 23, 29, 31):
    t0 = time.time()
    c, lens, piF = pisano_F(n)
    results[n] = piF
    extra = ""
    if sp.isprime(n):
        extra = f"  [p-1 = {sp.factorint(n-1)}, p+1 = {sp.factorint(n+1)}, pi_F = {sp.factorint(piF) if piF > 1 else {}}]"
    print(f"  n={n}: core {c}/{n**3}, cycles {lens}, pi_F({n}) = {piF}{extra}  [{time.time()-t0:.0f}s]", flush=True)
print(f"\n  HEADLINE: pi_F(10) = {results[10]}   (Fibonacci's clock: pi(10) = 60)", flush=True)
from math import lcm as _l
print(f"  CRT check: lcm(pi_F(2), pi_F(5)) = {_l(results[2], results[5])} vs pi_F(10) = {results[10]}", flush=True)
print("\nDONE.", flush=True)

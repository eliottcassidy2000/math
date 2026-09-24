#!/usr/bin/env python3
"""Orchestrator's independent audit of the arXiv:2502.20642v1 refutation (THM-4471), written without reading the lane's code.

Coefficients of the paper's Section 3 transcribed independently from the arXiv HTML.  Checks (each raises on failure):
  R1  the paper's Theorem 3.1 inequality holds at every pair x, y < 400 (its Collatz algebra is right);
  R2  condition (5) (pairwise reading, lambda = 0, A = 1/2, B = 2, M = 2) holds at every pair < 400;
  R3  the proof's needed alternative, alt1(p,Tp) or alt2(Tp,p) (Theorem 2.1 form, no B), fails exactly at the odd p >= 3 (p < 5000);
  R4  the claimed decay d(T^n x,T^(n+1)x)^2 <= (1/2)^n d(x,Tx)^2 fails at x = 5, n = 1;
  R5  counterexample 1: X = {0,1,2}, T = 3-cycle, M = 4 -- all hypotheses of Theorem 2.3(5) hold, no fixed point;
  R6  counterexample 2: T(x) = x+1 on (N,|.|) with the paper's constants (A=1/2, B=2, M=2) -- all hypotheses hold (x,y<300);
  R7  SHEET control: the paper's table, verbatim, satisfies all hypotheses for the 3n-1 shortcut map (x,y<400),
      whose orbit 5 -> 7 -> 10 -> 5 never reaches 1.
"""
from fractions import Fraction as Fr
def T(x):
    if x == 1: return 1
    return x // 2 if x % 2 == 0 else (3 * x + 1) // 2
def Tm(x): return x // 2 if x % 2 == 0 else (3 * x - 1) // 2
def cls(x): return 'one' if x == 1 else ('even' if x % 2 == 0 else 'odd')
TAB = {('one','one'):(1,0,0,0,-1,1), ('one','even'):(1,0,0,-1,0,1), ('one','odd'):(0,0,0,-2,1,2),
       ('even','one'):(1,0,0,-1,1,0), ('even','even'):(1,0,-1,0,-1,1), ('even','odd'):(0,0,-2,1,-2,2),
       ('odd','one'):(0,0,0,-2,2,1), ('odd','even'):(0,-2,0,1,2,-2)}
def coeffs(x, y):
    c = TAB.get((cls(x), cls(y)))
    if c: return c
    k, l = (x - 1) // 2, (y - 1) // 2
    b0 = -2 if k - l <= -2 else (2 if k - l >= 2 else k - l)
    d0 = -2 if ((k-l <= -2 and 11*k-10*l+1 <= 0) or (k-l >= 2 and -10*k+11*l+1 <= 0)) else -1
    e0 = 2 if (k-l <= -2 and 11*k-10*l+1 <= 0) else 0
    z0 = 2 if (k-l >= 2 and -10*k+11*l+1 <= 0) else 0
    return (2, b0, -b0, d0, e0, z0)
d = lambda a, b: abs(a - b)
def wgp(c, x, y, F):
    a, b, g, dl, e, z = c
    return a*d(F(x),F(y))**2 + b*d(x,F(y))**2 + g*d(F(x),y)**2 + dl*d(x,y)**2 + e*d(x,F(x))**2 + z*d(y,F(y))**2
def alt1(c, A, B):
    a, b, g, dl, e, z = c; s = a + z + 2*min(b,0); r = dl + e + 2*min(b,0)
    return s > 0 and -r <= A*s and a + b + z >= B
def alt2(c, A, B):
    a, b, g, dl, e, z = c; s = a + e + 2*min(g,0); r = dl + z + 2*min(g,0)
    return s > 0 and -r <= A*s and a + g + e >= B
A, B, M = Fr(1,2), 2, 2
N = 400
assert all(wgp(coeffs(x,y), x, y, T) <= 0 for x in range(1,N) for y in range(1,N)); print("R1  Theorem 3.1 holds for x,y < 400: ok")
assert all(alt1(coeffs(x,y),A,B) or alt2(coeffs(x,y),A,B) for x in range(1,N) for y in range(1,N)); print("R2  condition (5) holds pairwise for x,y < 400: ok")
unc = [p for p in range(1,5000) if T(p) != p and not (alt1(coeffs(p,T(p)),A,0) or alt2(coeffs(T(p),p),A,0))]   # Theorem 2.1 form: no B
assert unc == list(range(3,5000,2)), unc[:10]; print("R3  uncovered steps p < 5000 are exactly the odd p >= 3 (", len(unc), "): ok")
assert d(T(8),8)**2 == 16 and 16 > A*d(5,T(5))**2 and T(5) == 8; print("R4  decay claim fails at x = 5, n = 1 (16 > 9/2): ok")
# R5 3-cycle
Tc = {0:1, 1:2, 2:0}; F3 = lambda x: Tc[x]
def c3(x,y):
    if x == y: return (4,0,-4,0,0,0)
    if Tc[x] == y: return (1,-4,4,0,0,0)
    return (1,4,-4,0,0,0)
assert all(wgp(c3(x,y),x,y,F3) <= 0 and (alt1(c3(x,y),A,B) or alt2(c3(x,y),A,B)) and max(map(abs,c3(x,y))) <= 4 for x in Tc for y in Tc)
print("R5  3-point 3-cycle satisfies Theorem 2.3(5) with M = 4, has no fixed point: ok")
# R6 successor
S = lambda x: x + 1
cs = lambda x, y: (1,1,-2,0,0,0) if x >= y else (1,-2,1,0,0,0)
assert all(wgp(cs(x,y),x,y,S) <= 0 and (alt1(cs(x,y),A,B) or alt2(cs(x,y),A,B)) and max(map(abs,cs(x,y))) <= M for x in range(1,300) for y in range(1,300))
print("R6  x -> x+1 on (N,|x-y|) satisfies Theorem 2.3(5) with the paper's constants, has no fixed point: ok")
# R7 SHEET control
assert all(wgp(coeffs(x,y),x,y,Tm) <= 0 and (alt1(coeffs(x,y),A,B) or alt2(coeffs(x,y),A,B)) for x in range(1,N) for y in range(1,N))
o = [5]
for _ in range(6): o.append(Tm(o[-1]))
assert o == [5,7,10,5,7,10,5]
print("R7  the paper's table, verbatim, satisfies all hypotheses for 3n-1 (x,y < 400); 3n-1 orbit of 5 is the 3-cycle", o[:3], ": ok")
print("ALL CHECKS PASSED")

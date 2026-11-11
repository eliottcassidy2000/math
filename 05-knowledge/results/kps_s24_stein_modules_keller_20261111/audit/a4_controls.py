"""a4_controls.py - controls on C^2 with the independent engine tw.py (two primes, two weightings).

Expected from Dimca-Saito (dim H^2 = b_1 of the generic fibre):
   xy -> 1;  x^3 + y^3 -> 4;  x + x^2 y -> 1;  x^2 y -> 1;  f g for an automorphism -> 1;
additional (audit): x^2 y^2 -> 2 (two components C^*), x^3 y^2 -> 1, x^2 y^2 + x -> ?, (x^2 - 1) y -> 2 (x-line minus two points... see AUDIT.md)
"""
import sys
from tw import RingA, add, mul, scal, power, table, build, codim, P1, P2

A = RingA()
x, y, one = {(1, 0): 1}, {(0, 1): 1}, {(0, 0): 1}


def pw(a, k):
    return power(A, a, k, one)


f1 = add(x, pw(y, 2))
g1 = add(y, pw(f1, 3))
f2 = add(x, pw(y, 2))
g2 = add(y, pw(f2, 2))
U2 = add(f2, scal(one, 2), -1)
V2 = add(g2, scal(one, 3), -1)
cases = [
    ("W = x y                                   expected 1", mul(A, x, y), (2, 4, 6), (8, 12)),
    ("W = x^3 + y^3                             expected 4", add(pw(x, 3), pw(y, 3)), (3, 5, 7), (10, 14)),
    ("W = x + x^2 y  (no critical point)        expected 1", add(x, mul(A, pw(x, 2), y)), (2, 4, 6), (12, 18)),
    ("W = x^2 y      (critical line)            expected 1", mul(A, pw(x, 2), y), (2, 4, 6), (12, 16)),
    ("W = f g, (f,g) = (x+y^2, y+(x+y^2)^3)     expected 1", mul(A, f1, g1), (2, 4, 6), (16, 22)),
    ("W = (f-2)(g-3), (f,g) = (x+y^2, y+(x+y^2)^2)  expected 1", mul(A, U2, V2), (2, 4, 6, 10), (16, 22)),
    ("W = x^2 y^2  (gcd 2: two components C^*)  expected 2  [HYP-9668 (E) says 1 for x^p y^q]", mul(A, pw(x, 2), pw(y, 2)), (2, 4, 6), (12, 16)),
    ("W = x^3 y^3  (gcd 3)                      expected 3", mul(A, pw(x, 3), pw(y, 3)), (2, 4, 6), (14, 18)),
    ("W = x^3 y^2  (coprime)                    expected 1", mul(A, pw(x, 3), pw(y, 2)), (2, 4, 6), (14, 18)),
]
for p in (P1, P2):
    print("=== prime", p)
    for label, W, Ns, Nps in cases:
        table(A, W, Ns, Nps, p, label)
# a second weighting (x of weight 2, y of weight 3) for three of the cases
A2 = RingA(2, 3)
print("=== prime", P1, " weights (2, 3)")
for label, W, Ns, Nps in (cases[2], cases[3], cases[6]):
    table(A2, W, (6, 12, 18), (36, 48), P1, label)

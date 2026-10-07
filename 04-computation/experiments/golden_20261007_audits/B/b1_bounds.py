#!/usr/bin/env python3
"""THM-4591 / HYP-9230 arithmetic: continued fraction of log2 3, convergents, semiconvergents,
||f log2 3||, and the lemma's min-element bounds maximised over every (L, k) with L <= 301993."""
from mpmath import mp, mpf, log, floor, nint, ceil, power
mp.dps = 80
x = log(3) / log(2)
# continued fraction
a = []; y = x
for _ in range(25):
    q = int(floor(y)); a.append(q); y = 1 / (y - q)
print("CF log2 3 =", a)
p2, q2, p1, q1 = 0, 1, 1, 0
conv = []
for ai in a:
    p, q = ai * p1 + p2, ai * q1 + q2
    conv.append((p, q)); p2, q2, p1, q1 = p1, q1, p, q
print("convergents:", conv[:16])
# semiconvergents (all intermediate fractions j = 1..a_i)
semi = []
P2, Q2, P1, Q1 = 0, 1, 1, 0
for i, ai in enumerate(a[:16]):
    for j in range(1, ai + 1):
        semi.append((j * P1 + P2, j * Q1 + Q2))
    P2, Q2, P1, Q1 = P1, Q1, ai * P1 + P2, ai * Q1 + Q2
print("semiconvergent denominators <= 700:", sorted({q for p, q in semi if q <= 700}))
def dist(f):
    v = f * x
    return abs(v - nint(v))
for f in (5, 12, 17, 29, 41, 53, 306, 665, 15601):
    print(f"||{f} log2 3|| = {mp.nstr(dist(f), 6)}   nearest p = {int(nint(f * x))}")
th = dist(665)
for C in (280, 550, 551):
    print(f"1/(C theta_665^2) with C={C}: {mp.nstr(1 / (C * th ** 2), 6)} octaves")

# bounds: positive side  min x <= 1/(2^(L/k) - 3)   (2^L > 3^k)
#         negative side  min|x| <= 1/(3 - 2^(L/k))  (3^k > 2^L)
best_pos = (mpf(0), None); best_neg = (mpf(0), None)
LMAX = 301993
for L in range(1, LMAX + 1):
    kp = int(floor(L / x))          # largest k with 3^k < 2^L  (L/x never an integer)
    kn = kp + 1                     # smallest k with 3^k > 2^L
    if kp >= 1:
        bp = 1 / (power(2, mpf(L) / kp) - 3)
        if bp > best_pos[0]:
            best_pos = (bp, (L, kp))
    bn = 1 / (3 - power(2, mpf(L) / kn))
    if bn > best_neg[0]:
        best_neg = (bn, (L, kn))
print("max positive bound over L <= 301993:", mp.nstr(best_pos[0], 12), "at", best_pos[1])
print("max negative bound over L <= 301993:", mp.nstr(best_neg[0], 12), "at", best_neg[1])
print("scan used XB = 7216102493, YB = 10295871817:",
      "XB >= floor(pos):", 7216102493 >= int(floor(best_pos[0])),
      " YB >= floor(neg):", 10295871817 >= int(floor(best_neg[0])))
# bound at the next convergent
print("bound at 301994/190537:", mp.nstr(1 / (power(2, mpf(301994) / 190537) - 3), 6))
print("bound at 27/17:", mp.nstr(1 / (power(2, mpf(27) / 17) - 3), 6), " at 19/12:", mp.nstr(1 / (3 - power(2, mpf(19) / 12)), 6))

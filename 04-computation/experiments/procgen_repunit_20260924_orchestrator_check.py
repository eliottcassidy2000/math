#!/usr/bin/env python3
"""Orchestrator's independent check of THM-4473 (digit chains, rotation, repunits), written without reading the lane's code.
  D1  Haar law of the last digit of consecutive odd Syracuse terms U(n) = (3n+1)/2^v, computed exactly by
      enumerating n mod 10*2^K (K = 18, tail < 2^-16): 3 -> 5 with probability 1; every other entry is 2^j/15.
  D2  mod-9 stationary law (8,16,11,4,2,22)/63 on (1,2,4,5,7,8), from the exact 6-periodic transition rule; mod 3 i.i.d. (0,1/3,2/3).
  D3  rotation: for b = 10 and every k <= 6, the fixed points of x -> b x mod (b^k - 1) are exactly the repdigits c*R_k.
  D4  the parity-vector map Q (Bernstein-Lagarias): Q(-1) = -1, Q(1/3) = 1/3, {1, -1/3} and {-1/5, 5/7} are 2-cycles (exact rationals);
      among odd-numerator rationals p/q, |p| <= 121, odd q < 60, these are all fixed points and 2-cycles.
  D5  T^k(2^k - 1) = 3^k - 1 and T^(k+1)(2^k - 1) = (3^k - 1)/2 for k <= 200.
"""
from fractions import Fraction as F
from collections import Counter
import math
# D1
K = 18; M = 10 * 2**K
cnt = {d: Counter() for d in (1,3,5,7,9)}
for n in range(1, M, 2):
    m = 3*n + 1; k = (m & -m).bit_length() - 1
    if k >= K - 2: continue
    cnt[n % 10][(m >> k) % 10] += 1
assert set(cnt[3]) == {5}
for d in (1,5,7,9):
    tot = sum(cnt[d].values())
    got = sorted(round(c * 15 / tot, 3) for c in cnt[d].values())
    assert got == [1.0, 2.0, 4.0, 8.0], (d, got)
assert abs(cnt[9][9] / sum(cnt[9].values()) - 8/15) < 1e-4
print("D1  Syracuse last-digit law: 3 -> 5 certain; other rows are a permutation of (1,2,4,8)/15; 9 -> 9 = 8/15: ok")
# D2
inv2 = [pow(2, -k, 9) for k in range(1, 7)]            # 2^-k mod 9, k = 1..6
w = [F(2**(5-j), 63) for j in range(6)]                 # P(k = j+1 mod 6) = 2^-(j+1) * 64/63
units = [1,2,4,5,7,8]
P = {r: Counter() for r in units}
for r in units:
    for j in range(6): P[r][((3*r+1) * inv2[j]) % 9] += w[j]
pi = {r: F(1,6) for r in units}
for _ in range(200): pi = {s: sum(pi[r] * P[r][s] for r in units) for s in units}
target = dict(zip(units, [F(x,63) for x in (8,16,11,4,2,22)]))
assert all(abs(pi[s] - target[s]) < F(1,10**12) for s in units)
for r in units:
    m3 = Counter(); [m3.__setitem__(s % 3, m3[s % 3] + P[r][s]) for s in P[r]]
    assert m3[0] == 0 and m3[1] == F(1,3) and m3[2] == F(2,3)
print("D2  mod-9 stationary law (8,16,11,4,2,22)/63 and mod-3 rows (0,1/3,2/3): ok")
# D3
for k in range(1, 7):
    Nk = 10**k - 1; R = Nk // 9
    fixed = [x for x in range(Nk) if (10 * x) % Nk == x]
    assert fixed == [c * R for c in range(9)], (k, fixed[:12])
print("D3  fixed points of digit rotation (base 10, k <= 6) are exactly the repdigits c*R_k: ok")
# D4
def Tq(x): return x/2 if x.numerator % 2 == 0 else (3*x + 1)/2
def Q(x, maxlen=600):
    seen = {}; bits = []; y = x; j = 0
    while y not in seen:
        seen[y] = j; bits.append(y.numerator % 2); y = Tq(y); j += 1
        if j > maxlen: raise RuntimeError
    s = seen[y]; per = bits[s:]
    return sum(F(b) * 2**i for i, b in enumerate(bits[:s])) + F(2**s) * sum(F(b) * 2**i for i, b in enumerate(per)) / (1 - F(2)**len(per))
assert Q(F(-1)) == -1 and Q(F(1,3)) == F(1,3) and Q(F(1)) == F(-1,3) and Q(F(-1,3)) == 1 and Q(F(-1,5)) == F(5,7) and Q(F(5,7)) == F(-1,5)
fix = set(); two = set()
for q in range(1, 60, 2):
    for p in range(-121, 122, 2):                       # odd numerators
        if math.gcd(p, q) != 1: continue
        x = F(p, q)
        try:
            y = Q(x)
            if y == x: fix.add(x)
            elif y.numerator % 2 and Q(y) == x: two.add(frozenset((x, y)))
        except RuntimeError: pass
assert fix == {F(-1), F(1,3)}, fix
assert two == {frozenset((F(1), F(-1,3))), frozenset((F(-1,5), F(5,7)))}, two
print("D4  Q fixed points {-1, 1/3}; odd 2-cycles {1,-1/3} and {-1/5,5/7} (odd p, |p|<=121, odd q<60): ok")
# D5
def T(n): return n // 2 if n % 2 == 0 else (3*n + 1) // 2
for k in range(1, 201):
    x = 2**k - 1
    for _ in range(k): x = T(x)
    assert x == 3**k - 1 and T(x) == (3**k - 1) // 2
print("D5  T^k(2^k-1) = 3^k-1 and T^(k+1)(2^k-1) = (3^k-1)/2 for k <= 200: ok")
print("ALL CHECKS PASSED")

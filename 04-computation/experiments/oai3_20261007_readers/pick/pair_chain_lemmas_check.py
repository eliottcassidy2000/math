"""Brute-force check of the local lemmas used in the a.s.-absorption proof, over all states
-KMAX <= k <= KMAX and |num| <= NMAX (num = 3^max(0,-k) * e), both driving bits beta.
 L1 (fair coin c): c = beta if k >= 0, beta xor (e mod 2) if k < 0.
 L2 (scale bound): |f'| <= A(c)|f| + add, A(0)=1/2, A(1)=3/2, add <= 1/2, where |f| = |e| 3^-max(k,0);
     departures from k=0 satisfy |f'| <= |f|/2 + 1/6 for both coins.
 L3 (flip direction): for k != 0 and e odd: |k'| = |k|-1 iff c = 1.
 L4 (runs, e even): k = 0: both coins keep v2 -> v2-1, no additive term;
     k odd: exactly one coin value gives e' even; k even != 0: v2 dynamics with m = v2(3^|k| - 1):
     c=0: v-1; c=1: min(v,m)-1 if v != m, >= m if v == m; additive term only when c = 1.
"""
from fractions import Fraction as Fr
from pair_chain_verify import step

def v2(n):
    if n == 0: return 10**9
    n = abs(n); c = 0
    while n % 2 == 0: n //= 2; c += 1
    return c

def e_of(k, num):
    return Fr(num, 3 ** (-k if k < 0 else 0))

def fabs(k, num):
    return abs(e_of(k, num)) / (3 ** max(k, 0))

KMAX, NMAX = 9, 1500
bad = {"L1": 0, "L2": 0, "L3": 0, "L4": 0}; checked = 0
worst_add = Fr(0)
for k in range(-KMAX, KMAX + 1):
    for num in range(-NMAX, NMAX + 1):
        e = e_of(k, num); par = num & 1
        f = fabs(k, num)
        res = {}
        for beta in (0, 1):
            c = beta if k >= 0 else beta ^ par
            k2, n2, _ = step(k, num, beta)
            f2 = fabs(k2, n2)
            A = Fr(1, 2) if c == 0 else Fr(3, 2)
            add = f2 - A * f
            if k == 0 and par == 1:
                if f2 > f / 2 + Fr(1, 6): bad["L2"] += 1
            else:
                if add > Fr(1, 2): bad["L2"] += 1
                worst_add = max(worst_add, add)
            if k != 0 and par == 1:
                toward = abs(k2) == abs(k) - 1
                if toward != (c == 1): bad["L3"] += 1
            res[c] = (k2, n2, add, f2)
            checked += 1
        if par == 0 and not (k == 0 and num == 0):
            ev = {c: (res[c][1] & 1) == 0 for c in (0, 1)}
            if k == 0:
                for c in (0, 1):
                    if v2(res[c][1]) != v2(num) - 1 or res[c][3] != (Fr(1, 2) if c == 0 else Fr(3, 2)) * f: bad["L4"] += 1
            elif abs(k) % 2 == 1:
                if ev[0] == ev[1]: bad["L4"] += 1
            else:
                m = v2(3 ** abs(k) - 1); v = v2(num)
                if v2(res[0][1]) != v - 1 and v < 10**8: bad["L4"] += 1
                w1 = v2(res[1][1])
                if v != m and w1 != min(v, m) - 1: bad["L4"] += 1
                if v == m and w1 < m: bad["L4"] += 1
                if res[0][3] != f / 2: bad["L4"] += 1     # no additive term at c = 0
print(f"checked {checked} (state, bit) pairs; violations: {bad}; max additive term seen (non-departure): {float(worst_add):.4f}")

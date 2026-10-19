"""Observation: the sorted elements of the Christoffel (bulb) free cycle of 3x + (2^q - 3^p) have consecutive differences
that are 3-smooth numbers 2^a 3^b, forming consecutive 3-smooth numbers. Test for all coprime (q, p), q <= QMAX."""
from math import gcd, floor
import sys
QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
def christoffel(p, q): return [floor(k*p/q) - floor((k-1)*p/q) for k in range(1, q+1)]
def carry(w):
    P, Bc, Q = 1, 0, 1
    for b in w:
        if b: P, Bc, Q = 3*P, 3*Bc + Q, 2*Q
        else: Q = 2*Q
    return Bc
def smooth3(n):
    if n <= 0: return False
    while n % 2 == 0: n //= 2
    while n % 3 == 0: n //= 3
    return n == 1
def exps(n):
    a = b = 0
    while n % 2 == 0: n //= 2; a += 1
    while n % 3 == 0: n //= 3; b += 1
    return (a, b)
allsmooth = 0; consec = 0; total = 0; fails = []
for q in range(2, QMAX+1):
    for p in range(1, q):
        if gcd(p, q) != 1: continue
        w = christoffel(p, q)
        elems = sorted(carry(w[i:] + w[:i]) for i in range(q))
        D = 2**q - 3**p
        diffs = [elems[i+1] - elems[i] for i in range(q-1)]
        total += 1
        sm = all(smooth3(d) for d in diffs)
        allsmooth += sm
        # consecutive 3-smooth numbers?
        if sm and q > 1:
            lo, hi = min(diffs), max(diffs)
            S = []
            a = 0
            while 2**a <= hi:
                b = 0
                while 2**a * 3**b <= hi:
                    v = 2**a * 3**b
                    if v >= lo: S.append(v)
                    b += 1
                a += 1
            S.sort()
            if S == sorted(diffs): consec += 1
            elif q <= 12: fails.append((q, p, D, diffs))
        elif q <= 12: fails.append((q, p, D, diffs))
        if q <= 13:
            lo, hi = min(diffs), max(diffs)
            print(f"(q,p)=({q},{p}) D={D:7d}: min elem {elems[0]}, window [{lo}={exps(lo)}, {hi}={exps(hi)}], ratio {hi/lo:.4f}, #diffs {q-1}, sorted diffs {sorted(diffs)}")
print(f"\n{total} Christoffel cycles with q <= {QMAX}: all differences 3-smooth in {allsmooth}; consecutive 3-smooth in {consec}")
print("first failures (q<=12):", fails[:6])

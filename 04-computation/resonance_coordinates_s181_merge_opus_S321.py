# opus-2026-07-16-S321 -- HYP-6980: the S181 merge -- do the new coordinates
# (sawtooth profile, Y* spectrum, maxtree) separate S181's matched-energy twins?
# Sets: AP {1..13} (E=1469, TIGHT); 2AP-1 = {1,3,...,25} (E=1469, L=0.116);
# GAP 7x2 = {1..7} u {10..16} (E~1225, L~0.09); near-AP {21..32} (E=1245, L=0.033).
# Coordinates per set: E(S); M(S) exact; Sum-rho over pairs; sawtooth-sum
# (Sum [rho - 4/169]); maxtree; Y*-spectrum summary (min, median).
from fractions import Fraction
from math import gcd
from collections import defaultdict
import itertools, statistics

def T_of(a, b):
    s = a + b
    T, m = 0, 0
    while s - 13*m > 0:
        T += (min(2*a, s - 13*m)) * (1 if m == 0 else 2)
        m += 1
    return T

def rho(x1, x2):
    g = gcd(x1, x2)
    a, b = x1//g, x2//g
    if a > b: a, b = b, a
    return Fraction(T_of(a, b), 13*a*b)

def Ystar(x1, x2, Q=13):
    best = None
    for q in range(1, Q+1):
        p0 = round(q*x2/x1) if x1 <= x2 else round(q*x1/x2)
        lo, hi = (x1, x2) if x1 <= x2 else (x2, x1)
        for pp in (p0-1, p0, p0+1):
            if pp < 1: continue
            v = abs(q*hi - pp*lo)
            if best is None or v < best: best = v
    return best

def energy(S):
    cnt = defaultdict(int)
    for a in S:
        for b in S: cnt[a+b] += 1
    return sum(c*c for c in cnt.values())

def exact_M(S, qmax=350):
    best = Fraction(0)
    for q in range(2, qmax+1):
        for p in range(1, q):
            if gcd(p, q) != 1: continue
            t = Fraction(p, q)
            m = min(min((v*t) % 1, 1 - (v*t) % 1) for v in S)
            if m > best: best = m
    return best

def maxtree(xs):
    n = len(xs)
    W = {}
    for i, j in itertools.combinations(range(n), 2):
        W[(i, j)] = rho(xs[i], xs[j])
    intree = {0}; tot = Fraction(0)
    while len(intree) < n:
        best = None
        for i in intree:
            for j in range(n):
                if j not in intree:
                    w = W[(min(i,j), max(i,j))]
                    if best is None or w > best[0]: best = (w, j)
        tot += best[0]; intree.add(best[1])
    return tot

SETS = {
    'AP {1..13} (tight)': list(range(1, 14)),
    '2AP-1 {1,3,..,25}': list(range(1, 26, 2)),
    'GAP 7x2 {1..7}u{10..16}': list(range(1, 8)) + list(range(10, 17)),
    'near-AP {21..32}': list(range(21, 33)),
}
print(f"{'set':28s} {'|S|':3s} {'E':5s} {'M':10s} {'L=M-1/(n+1)':11s} "
      f"{'saw-sum':10s} {'maxtree':8s} {'Y*min':5s} {'Y*med':5s}")
for name, S in SETS.items():
    E = energy(S)
    M = exact_M(S)
    L = float(M - Fraction(1, len(S)+1))
    pairs = list(itertools.combinations(S, 2))
    saw = sum(rho(a, b) - Fraction(4, 169) for a, b in pairs)
    mt = maxtree(S)
    ys = sorted(Ystar(a, b) for a, b in pairs)
    print(f"{name:28s} {len(S):3d} {E:5d} {str(M):10s} {L:11.4f} "
          f"{float(saw):+10.4f} {float(mt):8.4f} {ys[0]:5d} "
          f"{statistics.median(ys):5.0f}")
print()
print("VERDICT: does any new coordinate order the matched-E twins correctly")
print("(AP tight vs 2AP-1 loose; near-AP tighter than GAP) where E cannot?")

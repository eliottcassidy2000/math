#!/usr/bin/env python3
"""Independent audit (written from scratch) of Section 2 of
05-knowledge/results/chessboard_weave_20261006.md:
Thm 2.1 (two-speed delta), Thm 2.2 (scaffold speed sets), Cor 2.3 (central block of 2m x 2m board),
AP lonely times, Thm 2.4 (first lonely time tau at level 1/3).
Exact arithmetic (Fractions / scaled integers).
Method for delta: g(t) = min_i ||t v_i|| is piecewise linear with every piece of nonzero slope,
breakpoints in {k/(2v_i)} u {k/(v_i+v_j)} u {k/|v_i-v_j|}; the max over [0,1] is attained at a breakpoint.
Method for 'meets a closed box' and tau: exact intersection of unions of closed intervals.
"""
from fractions import Fraction as Fr
from math import gcd, floor
import itertools

fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)
        print("FAIL:", msg)


def nd(x):  # ||x|| distance to nearest integer, x Fraction
    f = x - floor(x)
    return min(f, 1 - f)


def delta_exact(v):
    cands = set()
    for a in v:
        for k in range(0, 2 * a + 1):
            cands.add(Fr(k, 2 * a))
    for a, b in itertools.combinations(v, 2):
        for d in (a + b, abs(a - b)):
            if d:
                for k in range(0, d + 1):
                    cands.add(Fr(k, d))
    best, arg = Fr(-1), []
    for t in cands:
        g = min(nd(t * a) for a in v)
        if g > best:
            best, arg = g, [t]
        elif g == best:
            arg.append(t)
    return best, sorted(arg)


# ---------- Theorem 2.1 ----------
print("=== Thm 2.1 ===")
cnt40 = 0
for b in range(1, 61):
    for a in range(1, b + 1):
        if gcd(a, b) != 1:
            continue
        if a < b <= 40:
            cnt40 += 1
        d, _ = delta_exact([a, b])
        s = a + b
        check(d == Fr(s // 2, s), f"delta({a},{b})")
print("coprime pairs a<b<=40:", cnt40)
print("Thm 2.1 formula checked for all coprime 1<=a<=b<=60")
vals = sorted(set((delta_exact([a, b])[0], (a, b)) for b in range(1, 30) for a in range(1, b + 1) if gcd(a, b) == 1))
print("smallest deltas:", [(str(d), p) for d, p in vals[:4]])

# ---------- Theorem 2.2 ----------
print("=== Thm 2.2 ===")
for m in range(1, 6):
    odd = list(range(1, 2 * m, 2))
    even = list(range(2, 2 * m + 1, 2))
    ap = list(range(1, 2 * m + 1))
    d1, _ = delta_exact(odd)
    d2, _ = delta_exact(even)
    d3, _ = delta_exact(ap)
    check(d1 == Fr(1, 2), f"odd set m={m}")
    check(d2 == Fr(1, m + 1), f"even set m={m}")
    check(d3 == Fr(1, 2 * m + 1), f"AP 1..2m m={m}")
    print(f"m={m}: delta(odd)={d1}, delta(even)={d2}, delta(1..{2*m})={d3}")
print("delta(4,12,20,28) =", delta_exact([4, 12, 20, 28])[0])


# ---------- interval machinery ----------
def safe_intervals(a, lo, hi, T=Fr(1)):
    """closed intervals of t in [0,T] with frac(t a) in [lo,hi] (0<lo<=hi<1)."""
    out = []
    i = 0
    while True:
        L, R = (i + lo) / a, (i + hi) / a
        if L > T:
            break
        out.append((L, min(R, T)))
        i += 1
    return out


def intersect(A, B):
    out = []
    i = j = 0
    while i < len(A) and j < len(B):
        L = max(A[i][0], B[j][0])
        R = min(A[i][1], B[j][1])
        if L <= R:
            out.append((L, R))
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return out


# ---------- AP lonely times ----------
print("=== AP lonely times ===")
for n in range(1, 13):
    lo, hi = Fr(1, n + 1), Fr(n, n + 1)
    S = [(Fr(0), Fr(1))]
    for j in range(1, n + 1):
        S = intersect(S, safe_intervals(j, lo, hi))
    pts = sorted(set(L for (L, R) in S if L == R and L < 1))
    nonpoint = [I for I in S if I[0] != I[1]]
    expect = [Fr(k, n + 1) for k in range(1, n + 1) if gcd(k, n + 1) == 1]
    check(not nonpoint and pts == expect, f"AP lonely times n={n}")
    # permutation points
    for t in expect:
        coords = sorted((j * t - floor(j * t)) for j in range(1, n + 1))
        check(coords == [Fr(i, n + 1) for i in range(1, n + 1)], f"perm point n={n}")
print("AP lonely times = {k/(n+1): gcd(k,n+1)=1} and permutation points: checked n=1..12")

# ---------- Corollary 2.3 ----------
print("=== Cor 2.3 ===")


def meets_box(a, b, lo, hi):
    return len(intersect(safe_intervals(a, lo, hi), safe_intervals(b, lo, hi))) > 0


for m in range(2, 13):
    N = 2 * m
    lo, hi = Fr(m - 1, N), Fr(m + 1, N)
    exc_closed, exc_open = [], []
    for s in range(2, 4 * m + 3):
        for a in range(1, s):
            b = s - a
            if gcd(a, b) != 1 or a > b:
                continue
            # closed box test
            if not meets_box(a, b, lo, hi):
                exc_closed.append((a, b))
            # open box test: intersection with positive length or interior point
            I = intersect(safe_intervals(a, lo, hi), safe_intervals(b, lo, hi))
            # interior of box met iff some t with both coords strictly inside: shrink by epsilon test
            # exact: interior met iff intersection contains a nondegenerate interval
            if not any(R > L for (L, R) in I):
                exc_open.append((a, b))
            check(meets_box(a, b, lo, hi) == (delta_exact([a, b])[0] >= lo), f"box vs delta {a},{b},m={m}")
    pred_closed = sorted((a, s - a) for s in range(2, m) for a in range(1, s) if s % 2 == 1 and gcd(a, s - a) == 1 and a <= s - a)
    pred_open = sorted((a, s - a) for s in range(2, m + 1) for a in range(1, s) if s % 2 == 1 and gcd(a, s - a) == 1 and a <= s - a)
    check(sorted(exc_closed) == pred_closed, f"closed exceptions 2m={N}")
    check(sorted(exc_open) == pred_open, f"open exceptions 2m={N}")
    print(f"{N}x{N}: miss closed central 2x2 block: {sorted(exc_closed)}; miss open block: {sorted(exc_open)}")
# knight on 6x6: contact points with the closed block
I = intersect(safe_intervals(1, Fr(1, 3), Fr(2, 3)), safe_intervals(2, Fr(1, 3), Fr(2, 3)))
print("knight vs [1/3,2/3]^2 (6x6 central block / 3x3 centre): t-set", I,
      "points", [((L * 1) % 1, (L * 2) % 1) for (L, R) in I])
# odd boards with a central single square: smallest board with an exceptional rider
print("odd boards (central square [m/N,(m+1)/N]^2, N=2m+1), closed:")
for N in range(3, 12, 2):
    m = (N - 1) // 2
    lo, hi = Fr(m, N), Fr(m + 1, N)
    exc = [(a, s - a) for s in range(2, 4 * N) for a in range(1, s) if gcd(a, s - a) == 1 and a <= s - a and not meets_box(a, s - a, lo, hi)]
    print(f"   {N}x{N}: exceptions {exc}")

# ---------- Theorem 2.4 ----------
print("=== Thm 2.4 ===")


def tau(a, b):
    # earliest t>0 with frac(ta), frac(tb) in [1/3, 2/3]; work in units 1/(3ab)
    i = j = 0
    while True:
        A = ((3 * i + 1) * b, (3 * i + 2) * b)
        B = ((3 * j + 1) * a, (3 * j + 2) * a)
        L, R = max(A[0], B[0]), min(A[1], B[1])
        if L <= R:
            return Fr(L, 3 * a * b)
        if A[1] < B[1]:
            i += 1
        else:
            j += 1


LIM = 200
allv = {}
for b in range(2, LIM + 1):
    for a in range(1, b):
        if gcd(a, b) != 1:
            continue
        t = tau(a, b)
        allv[(a, b)] = t
        # brute check that t is safe and nothing earlier
        check(nd(t * a) >= Fr(1, 3) and nd(t * b) >= Fr(1, 3), f"tau safe {a},{b}")
        check(t >= Fr(1, 3 * a), f"tau lower bound {a},{b}")
        if a == 1:
            exp = Fr(1, 3) if b % 3 else Fr(1, 3) + Fr(1, 9 * (b // 3))
            check(t == exp, f"tau(1,{b})")
        elif b < 2 * a:
            check(t == Fr(1, 3 * a), f"tau a<b<2a {a},{b}")
        else:
            check(b > 2 * a and t <= Fr(2, 3 * a) <= Fr(1, 3), f"tau b>2a {a},{b}")
# independent brute-force cross-check of tau for small pairs with candidate times k/(3ab) grid
for b in range(2, 41):
    for a in range(1, b):
        if gcd(a, b) != 1:
            continue
        # all interval endpoints are multiples of 1/(3ab); the earliest safe time is one of them
        D = 3 * a * b
        first = None
        for k in range(1, D + 1):
            t = Fr(k, D)
            if nd(t * a) >= Fr(1, 3) and nd(t * b) >= Fr(1, 3):
                first = t
                break
        check(first == allv[(a, b)], f"tau brute {a},{b}")
top = sorted(allv.items(), key=lambda kv: -kv[1])[:8]
print("largest tau over coprime a<b<=%d:" % LIM, [(p, str(v)) for p, v in top])
print("pairs with tau>1/3:", sorted([p for p, v in allv.items() if v > Fr(1, 3)])[:12], "...")
check(all(p[0] == 1 and p[1] % 3 == 0 for p, v in allv.items() if v > Fr(1, 3)), "tau>1/3 only at (1,3k)")
check(max(allv.values()) == Fr(4, 9) and [p for p, v in allv.items() if v == Fr(4, 9)] == [(1, 3)], "sup tau = 4/9 only at (1,3)")
print("number of coprime pairs a<b<=%d checked: %d" % (LIM, len(allv)))
print()
print("TOTAL FAILURES:", len(fails))

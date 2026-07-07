#!/usr/bin/env python3
r"""
metagraph_pureblue_n7_kps_S66.py   (kps-S66, HYP-4967)

Verify the pure-blue count at n=7 cheaply, to test the conjecture
   pure-blue(n) = floor((n+1)/2) - [n even]   ( (n+1)/2 odd n, n/2-1 even n )
   n=3..7 predicted: 2,1,3,2,4.

Trick: pure-blue classes are SC classes whose ENTIRE tiling fiber is grid-symmetric
(blue-mult == tc). Only the blue sub-cube Fix(g) = Q_{e(n)} (2^e tilings, e=9 at n=7 =>
512, cheap) can contain blue tilings. Enumerate those, get their classes + blue-mult,
compute tc=H/|Aut| per class, and count classes with blue-mult == tc. That avoids the
full 2^15 enumeration. (For a class touched by the blue cube, blue-mult here counts ALL
its grid-sym tilings; tc = H/|Aut| is its total tiling count; pure-blue iff equal.)
"""
from itertools import permutations

def tiles_of(n):
    T = []
    for y in range(1, n - 1):
        for x in range(n, y + 1, -1):
            if x - y >= 2: T.append((x, y))
    return T

def grid_orbits(n, T):
    """orbits of the grid involution on tiles -> free coordinates for the blue sub-cube."""
    pos = {t: i for i, t in enumerate(T)}
    seen = set(); orbits = []
    for i, t in enumerate(T):
        if i in seen: continue
        x, y = t; j = pos[(n - y + 1, n - x + 1)]
        orb = {i, j}; seen |= orb; orbits.append(sorted(orb))
    return orbits

def tour_from_bits(n, T, bits):
    A = [[0]*(n+1) for _ in range(n+1)]
    for k in range(2, n+1): A[k][k-1] = 1
    for (x, y), b in zip(T, bits):
        if b == 0: A[x][y] = 1
        else: A[y][x] = 1
    return A

def canon(n, A):
    best = None
    for p in permutations(range(1, n+1)):
        key = 0
        for a in range(n):
            pa = p[a]; Aa = A[pa]
            for b in range(n):
                if a != b and Aa[p[b]]:
                    key |= 1 << (a*n + b)
        if best is None or key < best: best = key
    return best

def aut_size(n, A):
    cnt = 0
    for p in permutations(range(1, n+1)):
        ok = True
        for a in range(1, n+1):
            for b in range(1, n+1):
                if a != b and A[a][b] != A[p[a-1]][p[b-1]]: ok = False; break
            if not ok: break
        if ok: cnt += 1
    return cnt

def ham(n, A):
    dp = {}
    for i in range(1, n+1): dp[(1 << (i-1), i)] = 1
    full = (1 << n) - 1
    for mask in range(1 << n):
        for last in range(1, n+1):
            c = dp.get((mask, last), 0)
            if not c: continue
            for nx in range(1, n+1):
                if not (mask >> (nx-1)) & 1 and A[last][nx]:
                    k = (mask | (1 << (nx-1)), nx); dp[k] = dp.get(k, 0) + c
    return sum(dp.get((full, i), 0) for i in range(1, n+1))

for n in (5, 6, 7):
    T = tiles_of(n); m = len(T); orbits = grid_orbits(n, T)
    e = len(orbits)
    assert e == ((n-1)**2)//4, (e, ((n-1)**2)//4)
    # enumerate the 2^e blue tilings
    blue_mult = {}; reps = {}
    for code in range(1 << e):
        bits = [0]*m
        for oi, orb in enumerate(orbits):
            v = (code >> oi) & 1
            for idx in orb: bits[idx] = v
        A = tour_from_bits(n, T, bits)
        ck = canon(n, A)
        blue_mult[ck] = blue_mult.get(ck, 0) + 1
        if ck not in reps: reps[ck] = A
    pureblue = 0; details = []
    for ck, bm in blue_mult.items():
        A = reps[ck]; H = ham(n, A); au = aut_size(n, A); tc = H // au
        if bm == tc:
            pureblue += 1; details.append((H, au, tc, bm))
    pred = ((n+1)//2) - (1 if n % 2 == 0 else 0)
    print(f"n={n}: blue sub-cube 2^{e}={1<<e} tilings; SC-classes touched={len(blue_mult)}; "
          f"PURE-BLUE={pureblue}  (conjecture floor((n+1)/2)-[n even]={pred}: {pureblue==pred})")
    for H, au, tc, bm in sorted(details):
        print(f"    pure-blue: H={H} |Aut|={au} tc={tc}")
print("DONE.")

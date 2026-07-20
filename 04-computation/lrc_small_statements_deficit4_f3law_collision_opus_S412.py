"""
opus-2026-07-19-S412 (HYP-8060): three small statements, connection-checked first.

CONNECTION CAUGHT by the mandated grep: canon's s470/s577 thread (incomplete round
tournament + two-nearest digraph; antipodal-tie census "n=14 has 7") is an ANCESTOR of
the S407 regularity bridge -- credited in the write-up.  kind-pasteur is running a
parallel becoming-the-concepts session (HYP-8050); scope here = my own S407/S409/S411
follow-ons only.

S1  THE DEFICIT-4 LAW, extended + mechanism: at t = 1/(N+1), GW_N's positions are the
    uniform grid with slot N-1 DELETED and slot N-3 DOUBLED (since 2(N-1) == N-3 mod
    N+1 -- the collision pair).  Verify deficit = 4 for gated N = 7..43 and that the
    delete+double reading reproduces it.
S2  THE F_3 1/N LAW, extended: M(F_3(N)) = 1/N for odd N != 1 (mod 3), N = 21..41;
    maximizers at (N+-1)/(3N).
S3  THE COLLISION PAIR: GW_N contains the pair (N-3, 2(N-1)) differing by N+1 --
    positions collide at EVERY witness-grid time k/(N+1); verify, and verify the
    tight AP has NO such pair: the collision is the GW mechanism's signature (the
    same tie-degeneracy that killed the coarse winding tournament, HYP-3093 --
    one thread's noise is another's signal).
"""
from math import gcd, comb
from fractions import Fraction

def scan0(V, qmax):
    bg, bq, wit = 0, 1, None
    for q in range(2, qmax+1):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            m = q
            for v in V:
                r = (v*a) % q
                r = min(r, q-r)
                if r < m:
                    m = r
                    if m*bq < bg*q: break
            if m*bq > bg*q: bg, bq, wit = m, q, (a, q)
    return Fraction(bg, bq), wit

def c3_at(V, a, q):
    pos = [Fraction(0)] + [Fraction((v*a) % q, q) for v in V]
    n = len(pos); out = [0]*n
    for i in range(n):
        for j in range(n):
            if i == j: continue
            d = (pos[j]-pos[i]) % 1
            if (0 < d < Fraction(1,2)) or ((d == 0 or d == Fraction(1,2)) and i < j):
                out[i] += 1
    return comb(n,3) - sum(comb(o,2) for o in out)

print("S1: deficit-4, gated N = 7..43 (N == 1 mod 6):")
for N in range(7, 44, 6):
    APn = list(range(1, N+1)); GWn = list(range(1, N-1)) + [N, 2*(N-1)]
    n = N + 1
    mx = n*(n*n-4)//24 if n % 2 == 0 else (n**3-n)//24
    cA, cG = c3_at(APn, 1, N+1), c3_at(GWn, 1, N+1)
    # delete+double reading: positions of GW = grid \ {N-1} with N-3 doubled
    posmulti = sorted([(v % (N+1)) for v in GWn])
    reading = (posmulti.count(N-3) == 2 and (N-1) not in posmulti)
    print(f"  N={N}: c3(AP)={cA} (max {mx}), c3(GW)={cG}, deficit={mx-cG}; "
          f"delete({N-1})+double({N-3}) reading holds: {reading}")

print("\nS2: F_3 1/N law, odd N != 1 mod 3, N = 21..41:")
for N in range(21, 42, 2):
    if N % 3 == 1: continue
    V = list(range(1, N-1)) + [N, 3*(N-1)]
    M, (a, q) = scan0(V, 4*N + 30)
    print(f"  N={N}: M={M} (1/N: {M == Fraction(1, N)}); t*={a}/{q} "
          f"(q=3N: {q == 3*N}, a=N+-1: {a in (N-1, N+1)})")

print("\nS3: the collision pair signature:")
for N in (7, 13, 19, 25):
    GWn = list(range(1, N-1)) + [N, 2*(N-1)]
    APn = list(range(1, N+1))
    collG = [(x, y) for i, x in enumerate(GWn) for y in GWn[i+1:]
             if (y - x) % (N+1) == 0]
    collA = [(x, y) for i, x in enumerate(APn) for y in APn[i+1:]
             if (y - x) % (N+1) == 0]
    print(f"  N={N}: GW collision pairs mod {N+1}: {collG}; AP: {collA}")

"""
opus-2026-07-19-S411 (HYP-8025): five proof/hypothesis EXTENSIONS, one per inhabited
concept (DFS style).

E1  MULTI-POINT GHOST PACKING (inhabiting THM-1292): adding the speed 13 to the point
    set buys a second soft gap: for the ladder {1..11, 13, 12m}, points {0..13}t give
    12 hard gaps + 2 soft (pairs differing by 12: (0,12),(1,13)), each soft >= M/m
    (12m in V) => M <= m/(12m+2), sharpening THM-1292's m/(12m+1).  For GW (K=2):
    M <= 1/13, sharpening 2/25.  Verify ceilings + case analysis adversarially.
E2  TOWER GATES AS PHASE-DENIAL (inhabiting S409's mechanism): the odd non-members of
    the D=3 tower row fall to the SAME (3, 2(N-1)) escape at Q = 2N+1 whenever
    3a == +-2 is solvable (N != 1 mod 3); verify M(F_3(N)) = 2/(2N+1) there and
    denial at members.
E3  THE c-SHEET BRANCH LEMMA (inhabiting the sheet tower): 11 of 12 speeds divisible
    by c >= 2, one coprime exception => M >= min(1/12, 1/2 - 1/(2c)) = 1/12.
    (Proof: LRC(12) on the quotient + best-of-c-sheets for the exception: the c sheet
    offsets b/c are 1/c-spaced, one lies in [1/2 - 1/(2c), 1/2 + 1/(2c)] distance
    >= 1/2 - 1/(2c) >= 1/4 >= 1/12.)  Verify exactly, c = 2..5.
E4  CROSS-N c3 DEFICIT (inhabiting the regularity bridge): c3 of GW_N vs AP_N at
    t = 1/(N+1) for gated N = 7, 13, 19: is the GW deficit from max ALWAYS 4?
E5  THETA-INVARIANT ACTIVE PAIR (inhabiting the deformation): the AP13 active pair
    under v^theta||vt|| at the theta-maximizer, theta = p/8: always (1, 13)?
"""
import random
from math import gcd, comb
from fractions import Fraction
random.seed(8025)

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

print("E1: multi-point ghost ceilings:")
for m in range(3, 9):
    V = list(range(1,12)) + [13, 12*m]
    M, _ = scan0(V, 24*m + 40)
    old, new = Fraction(m,12*m+1), Fraction(m,12*m+2)
    print(f"  ladder m={m}: M={M}; old ceiling {old}, NEW {new}: "
          f"M<=new {M <= new}{'  (actual '+str(M)+' = m/(12m+5))' if M==Fraction(m,12*m+5) else ''}")
GW = list(range(1,12)) + [13,24]
Mgw, _ = scan0(GW, 200)
print(f"  GW: M={Mgw}; old ceiling 2/25, NEW 1/13: holds {Mgw <= Fraction(1,13)}")
viol = 0
for _ in range(600):   # adversarial: shapes {1..11, 13, 12K} + noise speed
    K = random.randint(2, 9)
    extra = random.choice([x for x in range(14, 160) if x % 12])
    V = list(range(1,12)) + [13, 12*K, extra]
    M, _ = scan0(V, 130)
    if M > Fraction(K, 12*K + 2): viol += 1; print("  E1 VIOLATION:", V, M)
print(f"  adversarial (with extra speed, ceiling still applies): violations {viol}/600")

print("\nE2: D=3 tower row -- phase-denial anatomy:")
for N in (19, 21, 23, 25, 27, 29, 31):
    V = list(range(1, N-1)) + [N, 3*(N-1)]
    M, (a, q) = scan0(V, 6*N + 40)
    Q2 = 2*N + 1
    denied = (Q2 % 3 == 0)
    pred2 = Fraction(2, Q2)
    print(f"  F_3({N}): M={M} at t*={a}/{q}; D=2 escape denied (3|{Q2}): {denied}; "
          f"M == 2/(2N+1): {M == pred2}")

print("\nE3: c-sheet branch lemma, exact verification:")
ok3 = True
for c in range(2, 6):
    for _ in range(80):
        U = random.sample(range(1, 120), 11)
        b = random.choice([x for x in range(1, 300) if x % c])
        V = [c*u for u in U] + [b]
        M, _ = scan0(V, 90)
        if M < Fraction(1, 12):
            ok3 = False; print(f"  VIOLATION c={c}: {V} M={M}")
print(f"  c = 2..5, 80 random each: all M >= 1/12: {ok3}")

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

print("\nE4: cross-N c3 deficit at t = 1/(N+1):")
for N in (7, 13, 19):
    APn = list(range(1, N+1))
    GWn = list(range(1, N-1)) + [N, 2*(N-1)]
    n = N + 2  # runners incl. observer... positions: observer + N speeds
    mx = (N+1)*((N+1)**2 - 4)//24 if (N+1) % 2 == 0 else ((N+1)**3 - (N+1))//24
    cA = c3_at(APn, 1, N+1)
    cG = c3_at(GWn, 1, N+1)
    print(f"  N={N}: c3(AP)={cA}, c3(GW)={cG}, max={mx}; GW deficit = {mx - cG}")

print("\nE5: theta-invariance of the AP13 active pair (theta = p/8):")
AP = list(range(1, 14))
for p in range(0, 9):
    best, wit = None, None
    for q in range(2, 121):
        for a in range(1, q):
            if gcd(a, q) != 1: continue
            mm = None
            for v in AP:
                d = (v*a) % q; d = min(d, q-d)
                kv = Fraction((v**p) * (d**8), q**8)
                if mm is None or kv < mm: mm = kv
            if best is None or mm > best: best, wit = mm, (a, q)
    a, q = wit
    vals = []
    for v in AP:
        d = (v*a) % q; d = min(d, q-d)
        vals.append((Fraction((v**p)*(d**8), q**8), v))
    mn = min(vals)[0]
    act = sorted(v for val, v in vals if val == mn)
    print(f"  theta={p}/8: maximizer {wit}; active speeds {act}")

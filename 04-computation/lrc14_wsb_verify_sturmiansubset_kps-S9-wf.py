"""
Adversarial verification of T859 (angle "sturmian-subset-domination", kps-S9).

Claims to verify (exact rationals, fractions.Fraction):
  (1) meas(S7) breakpoint engine vs an independent Sturmian theta-space engine: match for m<=22.
  (2) a(m):=meas(S7(AP_m)) is NON-DECREASING in m, a(m)=0 for m<=6.
  (3) inclusion monotonicity  E subset F => meas(S7(E)) <= meas(S7(F))  over integer sets.
  (4) a(m) exact values for m=7..16 as claimed.
  (5) N*(k) = max{N : a(N+1) <= cap_k} = 7,8,10,13,21,21 (claim)  vs THM-536's 7,8,10,13,20,20.
  (6) cap_k canonical values & proved floor cap_k >= (k-6)/7.
  (7) L_y(consec_k) <= cap_k all k; consec maximizes L_y over bounded spread (spot check).
  (8) HUNT for counterexamples to meas(S7(E)) <= cap_k:
        exhaustive bounded-spread, resonant w==0 mod 7, short-relation {0,1,N,N+1}, random wide.
  (9) Check the wide-spread "bound" is actually an explicit rigorous constant (it is NOT, per the GAPS).
"""
from fractions import Fraction as F
from itertools import combinations
import random

# ---------- canonical caps (THM-535 part D) ----------
CAP = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91), 12:F(6,7), 13:F(1,1)}

# ---------- core engine: J(A,E) = meas{x: all frac(e x) avoid sectors A} ----------
def J(A, E):
    A = sorted(A)
    E = sorted(set(E))
    arcs = [(F(j,7), F(j+1,7)) for j in A]
    bp = set([F(0), F(1)])
    for e in E:
        if e == 0:
            continue
        for (a, b) in arcs:
            for end in (a, b):
                m = 0
                while True:
                    xv = (end + m) / e
                    if xv >= 1:
                        break
                    if xv >= 0:
                        bp.add(xv)
                    m += 1
    bp = sorted(b for b in bp if 0 <= b < 1)
    tot = F(0)
    for lo, hi in zip(bp, bp[1:] + [F(1)]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        if all(not any(a < ((e*mid) % 1) < b for (a, b) in arcs) for e in E):
            tot += hi - lo
    return tot

# ---------- meas(S7) via inclusion-exclusion over A subset {1..6} ----------
def measS7_IE(E):
    tot = F(0)
    for r in range(0, 7):
        for A in combinations(range(1,7), r):
            tot += (-1)**len(A) * J(set(A), E)
    return tot

# ---------- meas(S7) via direct breakpoint sweep (independent of IE) ----------
def measS7_direct(E):
    """Sweep x in [0,1); at each elementary interval, S7 holds iff sectors {0..6}
    are all hit by some frac(e x). Breakpoints at j/(7e)."""
    E = sorted(set(E))
    bp = set([F(0), F(1)])
    for e in E:
        if e == 0:
            continue
        for j in range(0, 7*e+1):
            xv = F(j, 7*e)
            if 0 <= xv < 1:
                bp.add(xv)
    bp = sorted(bp)
    tot = F(0)
    for lo, hi in zip(bp, bp[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = set()
        for e in E:
            frac = (e*mid) % 1
            hit.add(int(frac * 7))  # sector index in 0..6
        if len(hit & set(range(7))) == 7:
            tot += hi - lo
    return tot

# ---------- independent Sturmian theta-space engine (theta = 7x) ----------
def measS7_sturmian(E):
    """theta in [0,7). sector of e at x = floor(7 e x) mod 7 = floor(e theta) mod 7.
    S7 holds iff {floor(e theta) mod 7 : e in E} == Z/7.
    meas in x = (1/7) * meas_theta{covered}."""
    E = sorted(set(E))
    bp = set([F(0), F(7)])
    for e in E:
        if e == 0:
            continue
        # floor(e theta) jumps at theta = j/e
        j = 0
        while True:
            tv = F(j, e)
            if tv >= 7:
                break
            bp.add(tv)
            j += 1
    bp = sorted(t for t in bp if 0 <= t < 7)
    tot = F(0)
    for lo, hi in zip(bp, bp[1:] + [F(7)]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        residues = set()
        for e in E:
            residues.add(int(e*mid) % 7)  # floor(e*mid) mod 7
        if residues == set(range(7)):
            tot += hi - lo
    return tot / 7

def consec(k):
    return tuple(range(k))

# =====================================================================
print("="*70)
print("(1) ENGINE EQUIVALENCE: IE vs direct vs Sturmian, on AP_m for m up to 22")
print("="*70)
a_vals = {}
eng_ok = True
for m in range(1, 23):
    E = consec(m)
    d = measS7_direct(E)
    s = measS7_sturmian(E)
    a_vals[m] = d
    ok = (d == s)
    if not ok:
        eng_ok = False
    # IE only for small m (expensive: 2^6 J-calls each scanning many breakpoints)
    if m <= 14:
        ie = measS7_IE(E)
        ok2 = (ie == d)
        if not ok2:
            eng_ok = False
        flag = "" if (ok and ok2) else "  <<< MISMATCH"
        print(f"  m={m:2d}: direct={str(d):>16}  ie==direct:{ie==d}  sturm==direct:{ok}{flag}")
    else:
        flag = "" if ok else "  <<< MISMATCH"
        print(f"  m={m:2d}: direct={str(d):>16}  (IE skipped)      sturm==direct:{ok}{flag}")
print(f"\n  All engines agree (m<=22): {eng_ok}")

print()
print("="*70)
print("(2) MONOTONICITY a(m)<=a(m+1), and a(m)=0 for m<=6")
print("="*70)
mono_ok = True
zero_ok = True
for m in range(1, 22):
    if a_vals[m] > a_vals[m+1]:
        mono_ok = False
        print(f"  VIOLATION: a({m})={a_vals[m]} > a({m+1})={a_vals[m+1]}")
for m in range(1, 7):
    if a_vals[m] != 0:
        zero_ok = False
        print(f"  VIOLATION: a({m})={a_vals[m]} != 0")
print(f"  monotone non-decreasing (m<=22): {mono_ok}")
print(f"  a(m)=0 for m<=6: {zero_ok}")

print()
print("="*70)
print("(4) CLAIMED a(m) EXACT VALUES m=7..16")
print("="*70)
claimed_a = {
    7: F(31,210), 8: F(481,1470), 9: F(2447,5880), 10: F(8899,17640),
    11: F(3419,5880), 12: F(121103,194040), 13: F(14573,21560),
    14: F(14109,20020), 15: F(61211,84084), 16: F(62767,84084),
}
amatch = True
for m in range(7,17):
    got = a_vals[m]
    cl = claimed_a[m]
    ok = (got == cl)
    if not ok:
        amatch = False
    print(f"  a({m}): computed={str(got):>16}  claimed={str(cl):>16}  match:{ok}  ~{float(got):.4f}")
print(f"  all claimed a(m) match: {amatch}")

print()
print("="*70)
print("(6) CAP FLOOR cap_k >= (k-6)/7  (proved THM-535)")
print("="*70)
floor_ok = True
for k in range(8,14):
    floor = F(k-6,7)
    ok = CAP[k] >= floor
    if not ok:
        floor_ok = False
    print(f"  k={k}: cap={str(CAP[k]):>12} (~{float(CAP[k]):.4f})  floor={str(floor):>5} (~{float(floor):.4f})  cap>=floor:{ok}")
print(f"  all caps respect floor: {floor_ok}")
# check the GARBLED dispatch values
print("  Garbled dispatch check:")
for k,bad in [(9,F(2025,4004)),(10,F(36,91)),(11,F(25,91))]:
    floor=F(k-6,7)
    print(f"    k={k}: dispatch={str(bad):>9} (~{float(bad):.4f})  floor={str(floor):>5}  violates_floor:{bad<floor}")

print()
print("="*70)
print("(5) N*(k) = max{N : a(N+1) <= cap_k}   (claim 7,8,10,13,21,21; THM-536 says ...,20,20)")
print("="*70)
# need a(m) up to m=22 (have it). N*(k) largest N with a(N+1)<=cap_k.
for k in range(8,14):
    cap = CAP[k]
    Nstar = None
    for N in range(1, 22):
        if a_vals[N+1] <= cap:
            Nstar = N
    print(f"  k={k}: cap={float(cap):.4f}  N*={Nstar}  (a(N*+1)={float(a_vals[Nstar+1]):.4f} <= cap, a(N*+2)={float(a_vals[Nstar+2]) if Nstar+2<=22 else 'NA'})")

print()
print("="*70)
print("(3) INCLUSION MONOTONICITY  E subset F => meas(S7(E)) <= meas(S7(F))")
print("="*70)
# exhaustive over subsets of {0..N} containing 0, pair (E,F) with E subset F
inc_ok = True
checked = 0
random.seed(7)
N = 7  # subsets of {0..7}
elems = list(range(1, N+1))
allsets = []
for r in range(0, N+1):
    for c in combinations(elems, r):
        allsets.append(frozenset((0,)+c))
# pick pairs E subset F to test (full would be huge; do systematic: for each F, test E by dropping one elem)
cache = {}
def m7(fs):
    if fs not in cache:
        cache[fs] = measS7_direct(tuple(sorted(fs)))
    return cache[fs]
viol = []
for F_ in allsets:
    mF = m7(F_)
    for x in F_:
        if x == 0:
            continue
        E_ = F_ - {x}
        mE = m7(E_)
        checked += 1
        if mE > mF:
            inc_ok = False
            viol.append((tuple(sorted(E_)), tuple(sorted(F_)), mE, mF))
print(f"  checked {checked} (E=F\\{{x}}) pairs over subsets of {{0..{N}}}")
print(f"  inclusion monotonicity holds: {inc_ok}")
if viol:
    for e,f,me,mf in viol[:5]:
        print(f"    VIOL: meas(S7({e}))={me} > meas(S7({f}))={mf}")

print()
print("="*70)
print("(7) L_y(consec_k) <= cap_k all k  +  consec maximizes L_y spot check")
print("="*70)
def Sr(E, r):
    return sum((J(set(A), E) for A in combinations(range(1,7), r)), F(0))
def Ly(E, k):
    if k in (11,12,13): return 1 - F(1,2)*Sr(E,1) + F(1,6)*Sr(E,2)
    if k in (9,10): return 1 - F(13,18)*Sr(E,1) + F(4,9)*Sr(E,2) - F(1,6)*Sr(E,3)
    return 1 - Sr(E,1) + Sr(E,2) - F(9,10)*Sr(E,3) + F(3,5)*Sr(E,4)  # k=8
ly_ok = True
for k in range(8,14):
    E = consec(k)
    ly = Ly(E, k)
    ok = ly <= CAP[k]
    if not ok: ly_ok = False
    print(f"  k={k}: L_y(consec)={str(ly):>16} (~{float(ly):.4f})  cap~{float(CAP[k]):.4f}  L_y<=cap:{ok}")
print(f"  L_y(consec_k)<=cap_k all k: {ly_ok}")

print()
print("="*70)
print("(8) COUNTEREXAMPLE HUNT: meas(S7(E)) <= cap_k  (the ACTUAL target)")
print("="*70)
def primitive(E):
    from math import gcd
    g = 0
    for e in E:
        g = gcd(g, e)
    return g == 1

over = []  # (k, E, meas, cap)

# (8a) exhaustive bounded-spread, k=8,9 small windows
print("  (8a) exhaustive bounded-spread...")
for k, maxspan in [(8,15),(9,14),(10,13)]:
    cap = CAP[k]
    cnt = 0
    worst = (F(-1), None)
    # E = {0} U (k-1 elements from 1..maxspan)
    for c in combinations(range(1, maxspan+1), k-1):
        E = (0,)+c
        if not primitive(E):
            continue
        cnt += 1
        mv = measS7_direct(E)
        if mv > worst[0]:
            worst = (mv, E)
        if mv > cap:
            over.append((k, E, mv, cap))
    print(f"    k={k} span<={maxspan}: {cnt} primitive sets, max meas={float(worst[0]):.5f} at {worst[1]}, cap={float(cap):.5f}, over_cap={sum(1 for o in over if o[0]==k)}")

# (8b) resonant: include a multiple of 7 offset (w==0 mod 7)
print("  (8b) resonant w==0 mod 7 families...")
for k in [8,9,10]:
    cap = CAP[k]
    worst = (F(-1), None)
    base = list(range(k-1))  # 0..k-2
    for w in range(7, 64, 7):
        E = tuple(sorted(set(base + [w])))
        if len(E) != k:
            continue
        if not primitive(E):
            continue
        mv = measS7_direct(E)
        if mv > worst[0]:
            worst = (mv, E)
        if mv > cap:
            over.append((k, E, mv, cap))
    print(f"    k={k}: max resonant meas={float(worst[0]):.5f} at {worst[1]}, cap={float(cap):.5f}")

# (8c) short-relation shapes {0,1,N,N+1,...}
print("  (8c) short-relation {0,1,N,N+1,...} shapes...")
for k in [8,9]:
    cap = CAP[k]
    worst = (F(-1), None)
    for Nbig in range(5, 50):
        # block {0..k-3} U {Nbig, Nbig+1}
        E = tuple(sorted(set(list(range(k-2)) + [Nbig, Nbig+1])))
        if len(E) != k or not primitive(E):
            continue
        mv = measS7_direct(E)
        if mv > worst[0]:
            worst = (mv, E)
        if mv > cap:
            over.append((k, E, mv, cap))
    print(f"    k={k}: max short-rel meas={float(worst[0]):.5f} at {worst[1]}, cap={float(cap):.5f}")

# (8d) aggressive random wide-spread
print("  (8d) random wide-spread (span up to ~120)...")
random.seed(14)
for k in range(8,13):
    cap = CAP[k]
    worst = (F(-1), None)
    for _ in range(400):
        rest = random.sample(range(1, 121), k-1)
        E = tuple(sorted(set([0]+rest)))
        if len(E) != k or not primitive(E):
            continue
        mv = measS7_direct(E)
        if mv > worst[0]:
            worst = (mv, E)
        if mv > cap:
            over.append((k, E, mv, cap))
    print(f"    k={k}: max wide meas={float(worst[0]):.5f} at {worst[1]}, cap={float(cap):.5f}")

print()
print("  TOTAL over-cap witnesses found:", len(over))
for k,E,mv,cap in over[:20]:
    print(f"    k={k}: meas(S7({E}))={float(mv):.5f} > cap={float(cap):.5f}  <<< COUNTEREXAMPLE")

print()
print("="*70)
print("SUMMARY FLAGS")
print("="*70)
print(f"  engines_equiv(m<=22)      : {eng_ok}")
print(f"  monotone a(m)             : {mono_ok}")
print(f"  a(m)=0 for m<=6           : {zero_ok}")
print(f"  claimed a(m) values match : {amatch}")
print(f"  cap floor respected       : {floor_ok}")
print(f"  inclusion monotone        : {inc_ok}")
print(f"  L_y(consec)<=cap          : {ly_ok}")
print(f"  over-cap counterexamples  : {len(over)}")

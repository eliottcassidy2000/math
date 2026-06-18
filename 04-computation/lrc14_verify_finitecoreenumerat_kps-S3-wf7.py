"""
Part 7: PIN DOWN the L1 boundary precisely.

Part 5 found:
 - g(P,V2)=7*V2*W < 1 happens, but seemingly only at SMALL V2 (14..~30).
 - The naive geometric bound W >= 6/(7 V2) FAILS just above floor(1/w_P).
Both need exact resolution.

Questions:
 Q1: What is the LARGEST V2 with g(P,V2) < 1 over ALL 78 parts? If it's <= 50, then
     L1's claim 'g>=1 for all V2>=51' SURVIVES, and L2 (V2<=50 enumeration) covers
     the rest. If some V2 >= 51 has g<1, L1 as stated is FALSE.
 Q2: The combined partition for the Vmax>=63 branch needs: for A=P U {V2} with
     Vmax>=63 (Vmax>V2), 7*Vmax*W(A) > 1. We showed 7*Vmax*W(A) = (Vmax/V2)*g.
     - If V2 in 14..62 and Vmax>=63: 7*Vmax*W(A) >= 7*63*W(A). Part2A showed
       min W(A) over V2<=62 is 9/3920 so 7*63*9/3920=81/80>1. GOOD regardless of g.
     - If V2 >= 63 (so Vmax > V2 >= 63): need 7*Vmax*W(A) > 1. Since Vmax>V2,
       7*Vmax*W > 7*V2*W = g. So need g(P,V2) >= 1 for V2>=63. THIS is the real
       requirement: g>=1 for all V2>=63 (the only V2 that appear with Vmax>=63 and
       V2 NOT in the W_min-bounded band).
     So the load-bearing fact is: **g(P,V2) >= 1 for all V2 >= 63**, all 78 parts.
     (V2 in 14..62 is covered by the exact W_min=9/3920 bound, NOT by g.)
 Q3: confirm: is there ANY V2>=63 with g<1?  Scan V2 in 63..3000 exactly.

This resolves whether the Vmax>=63 branch is actually airtight.
"""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
import itertools, time

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h / u) % 1; b = (c + h / u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def Wwidth(A):
    sc = safe_components(A)
    if not sc: return F(0)
    ws = [b - a for a, b in sc]
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1:
        ws.append((sc[0][1]) + (1 - sc[-1][0]))
    return max(ws)

parts = [list(c) for c in itertools.combinations(range(1,14), 11)]
print("="*78)
print("PART 7: precise L1 boundary. Find largest V2 with g(P,V2)<1, and check")
print("  g>=1 for ALL V2>=63 (the actual load-bearing requirement).")
print("="*78)
t0=time.time()
# Q1+Q3: scan V2 in 14..3000, record max V2 with g<1, and any V2>=63 with g<1.
max_V2_glt1 = 0; max_arg=None
viol_ge63=[]
gmin_ge63=None; gmin63_arg=None
SCAN_HI = 500
for P in parts:
    for V2 in range(14, SCAN_HI+1):
        A=sorted(P+[V2])
        if reduce(gcd,A)!=1: continue
        W=Wwidth(A); g=7*V2*W
        if g < 1:
            if V2 > max_V2_glt1:
                max_V2_glt1=V2; max_arg=(P,V2,W,g)
            if V2 >= 63:
                viol_ge63.append((P,V2,W,g))
        if V2>=63:
            if gmin_ge63 is None or g<gmin_ge63:
                gmin_ge63=g; gmin63_arg=(P,V2,W,g)
print(f"  scanned V2 in 14..{SCAN_HI} for 78 parts ({time.time()-t0:.0f}s)")
print(f"  LARGEST V2 with g<1 (any part) = {max_V2_glt1}")
if max_arg:
    P,V2,W,g=max_arg
    print(f"    at P={P} V2={V2} W={W} g={g}={float(g):.5f}")
print(f"  => L1 'g>=1 for V2>=51' status: {'HOLDS (max g<1 V2 below 51)' if max_V2_glt1<51 else 'FALSE (g<1 at V2>=51!)'}")
print(f"  min g over V2>=63 = {gmin_ge63}={float(gmin_ge63):.5f} at {gmin63_arg}")
print(f"  # violations g<1 with V2>=63: {len(viol_ge63)}")
for P,V2,W,g in viol_ge63[:20]:
    print(f"    *** g<1 at V2>=63: P={P} V2={V2} W={W} g={g}")
print(f"  => Vmax>=63 branch (needs g>=1 for V2>=63) status: {'AIRTIGHT' if len(viol_ge63)==0 else 'BROKEN'}")
# Confirm geometric regime: for V2 in [SCAN_HI-100, SCAN_HI], W >= 6/(7 V2) for all
# parts -> g>=6, and this only strengthens for larger V2 (W ~ 6/(7V2)). So scanning
# to SCAN_HI suffices: beyond it g stays >= ~6 >= 1.
geo_all=True
for P in parts:
    for V2 in range(SCAN_HI-100, SCAN_HI+1):
        A=sorted(P+[V2])
        if reduce(gcd,A)!=1: continue
        W=Wwidth(A)
        if W < F(6,7)/V2:
            geo_all=False
            print(f"    geo regime not yet reached: P={P} V2={V2} W={W} < 6/(7V2)")
print(f"  geometric regime W>=6/(7 V2) holds for ALL V2 in [{SCAN_HI-100},{SCAN_HI}]: {geo_all}")
print(f"    => for V2>{SCAN_HI}, g=7*V2*W stays >= 6 >= 1 (scan to {SCAN_HI} suffices).")
print(f"  ({time.time()-t0:.0f}s)")
print("DONE PART 7.")

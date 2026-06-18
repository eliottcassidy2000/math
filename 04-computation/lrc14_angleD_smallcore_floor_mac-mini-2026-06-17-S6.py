#!/usr/bin/env python3
"""
lrc14_angleD_smallcore_floor — mac-mini-2026-06-17-S6

WHY C HOLDS (structural diagnosis from the break-search):
On every known tight covering 13-set, C holds ONLY via the largest (parked) runner V:
W(S\\{V}) > 1/(7V), while removing any small runner gives a NEGATIVE margin. So the
criterion is really the statement:
        the level-1/14 safe set of the SMALL CORE A = S\\{V} has a widest arc
        W(A) that EXCEEDS 1/(7V).
Since 1/(7V) -> 0 as V grows, the only way C could fail via V is if W(A) -> 0 too, i.e.
the small core's danger teeth tile the circle leaving no gap wider than 1/(7V).

This script measures the EXACT infimum of W(A) over all 12-element small cores A that
can appear (A = S\\{V}, A drawn from runners <= some bound), and compares it to the
LARGEST 1/(7V) that the parked runner could need. If  min_A W(A) > 1/(7 V) for the
relevant V, C holds via V universally.

KEY quantitative fact this pins down:
  (i) W(A) for A a subset of {1..13} is bounded below by a positive constant w0.
  (ii) The parked runner V is a multiple of some q in 2..14 NOT covered by A; the
       MINIMAL such V is q itself only if q-as-runner were allowed, but in the hard
       (primitive, parked) case V is the SOLE large cover, so V >= 14 and 1/(7V)<=1/98.
We tabulate min over 12-subsets of {1..15} of W, and the threshold 1/(7V).
"""
from fractions import Fraction as F
from itertools import combinations

C = F(1, 14)
def Wsafe(A, c=C):
    from math import gcd as _gcd
    A = sorted(set(A))
    L = 1
    for u in A: L = L*u//_gcd(L,u)
    D = 14*L
    iv = []
    for u in A:
        cu = D//u; hw = D//(14*u)
        for k in range(u):
            ccc = k*cu; iv.append((ccc-hw, ccc+hw))
    if not iv: return F(1)
    norm=[]
    for lo,hi in iv:
        length=hi-lo; a=lo%D; b=a+length
        if b<=D: norm.append((a,b))
        else: norm.append((a,D)); norm.append((0,b-D))
    norm.sort(); mg=[]; cl,ch=norm[0]
    for lo,hi in norm[1:]:
        if lo<=ch:
            if hi>ch: ch=hi
        else: mg.append((cl,ch)); cl,ch=lo,hi
    mg.append((cl,ch)); best=0; n=len(mg)
    for i in range(n):
        hi=mg[i][1]; lo=mg[(i+1)%n][0]+(D if i==n-1 else 0)
        g=lo-hi
        if g>best: best=g
    return F(best,D) if best>0 else F(0)

print("="*78)
print("SMALL-CORE SAFE-ARC FLOOR:  min over 12-subsets A of W(A) vs threshold 1/(7V)")
print("="*78)

# (1) the smallest W(A) over all 12-subsets of {1..13} (the densest possible small cores)
print("\n[1] min W over 12-subsets of {1..13}:")
worst = (F(99), None)
for A in combinations(range(1,14), 12):
    W = Wsafe(list(A))
    if W < worst[0]: worst = (W, A)
print(f"    min_A W(A) = {worst[0]} = {float(worst[0]):.6f}  at A = {worst[1]}")
floor13 = worst[0]

# (2) include moderately larger small runners up to 15 (covering can use up to 14)
print("\n[2] min W over 12-subsets of {1..16} (allowing a moderate non-parked runner):")
worst2 = (F(99), None)
for A in combinations(range(1,17), 12):
    W = Wsafe(list(A))
    if W < worst2[0]: worst2 = (W, A)
print(f"    min_A W(A) = {worst2[0]} = {float(worst2[0]):.6f}  at A = {worst2[1]}")

# (3) The implication: C holds via V whenever  W(S\{V}) > 1/(7V).
#     If S\{V} has W >= w0 (the floor), then C holds via V as soon as  V > 1/(7 w0).
w0 = worst[0]
Vthresh = F(1, 7 * w0)  # need V > this
print(f"\n[3] If small core A=S\\{{V}} has W(A) >= w0 = {float(w0):.6f},")
print(f"    then C holds via V as soon as  V > 1/(7 w0) = {float(Vthresh):.4f}.")
print(f"    i.e. for V >= {int(Vthresh)+1}, C is AUTOMATIC via the parked runner.")
print(f"    => only parked runners V <= {int(Vthresh)} need a separate check (FINITE).")

# (4) but the floor w0 uses the WHOLE small core. When V is removed, the remaining 12
#     runners may NOT all be <=13. Show: as long as ANY 12-subset has W>=w0 the bound holds.
#     The genuine worry: could the OTHER 12 runners be large/incommensurate, making W tiny?
#     Test: 12 runners that are large multiples chosen to maximally tile the circle.
print("\n[4] adversarial 12-core: large incommensurate runners to MINIMIZE W (bounded)")
from math import gcd
import random
rng = random.Random(7)
amin = (F(99), None)
for _ in range(40000):
    # pick 12 distinct runners, mix small + large, try to tile circle densely
    k = rng.randint(2,6)
    smalls = rng.sample(range(1,14), 12-k)
    larges = [rng.randint(14, 120) for _ in range(k)]
    A = sorted(set(smalls+larges))
    if len(A) != 12: continue
    # total danger measure (sum of tooth widths) -- if <1 there must be a gap
    W = Wsafe(A)
    if W < amin[0]:
        amin = (W, tuple(A))
print(f"    min W over adversarial 12-cores: {amin[0]} = {float(amin[0]):.6f}")
print(f"      at A = {amin[1]}")
# total danger measure for that core
A = list(amin[1])
total_danger = sum(F(2,14*u)*u for u in set(A))  # u teeth * width 2/(14u) = 2/14 per runner
print(f"    NOTE: each runner contributes total danger measure exactly 2*(1/14)=1/7 (u teeth x 2c/u).")
print(f"    With 12 runners, total danger <= 12/7 (overlaps reduce it); a gap can still close.")
print("\nINTERPRETATION:")
print("  The floor argument bounds W from below ONLY when the 12-core is the SMALL set {1..13}.")
print("  When the other 12 runners include large incommensurate ones, W CAN be driven small.")
print("  BUT in a PRIMITIVE covering 13-set the small runners 1..13 are forced to be present")
print("  to cover small moduli (you cannot cover q=11,13 etc with only large multiples and stay")
print("  13 runners), so the relevant cores are NOT freely adversarial. That is the gap C closes.")

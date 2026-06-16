#!/usr/bin/env python3
"""
lrc14_stranger_limit_macmini_0615s5.py  (mac-mini-2026-06-15-S5; OPEN-Q-104/097 ENDGAME)

Two exact finite checks that REFRAME inf L>0:

(A) PRIME POINT-MASS BOUND (Bedert Lemma 5.3 / Cor 5.4 — the RIGHT tool for AP-cores).
    ML(V) >= ceil((p-1)/(2n))/p for any prime p dividing no v in V (n=13, 2n=26).
    Does ANY admissible prime reach 1/14 = 0.07143 for the extremizer cores? (Answer
    expected: best is 2/29≈0.069 — JUST short. Quantifies the honest gap.)

(B) STRANGER-LIMIT REDUCTION. The cores are {1..13}\{j} ∪ {14m}. As the stranger 14m→∞,
    14m·τ equidistributes, so the stranger's constraint ||14m τ||>1/14 becomes
    asymptotically INDEPENDENT of the bounded core: conjecturally
        L({1..13}\{j} ∪ {14m})  ->  (6/7) · meas(Lonely({1..13}\{j})).
    If TRUE, inf over the infinite family = (6/7)·min_j meas(Lonely({1..13}\{j})), a
    FINITE set of 13 measures. We (i) compute meas(Lonely({1..13}\{j})) exactly for each j,
    (ii) compute L for m=1,2,3,5,8,13,21 and check convergence to (6/7)*that.
"""
import sys
from math import gcd
from sympy import isprime

sys.stdout.reconfigure(line_buffering=True)

def lonely_measure(S, Q):
    """meas{τ∈[0,1): ||vτ||>1/14 ∀v∈S}, threshold band = floor(Q/14)."""
    rad = Q//14; c=0
    for a in range(Q):
        ok=True
        for v in S:
            r=(v*a)%Q
            if r<=rad or r>=Q-rad: ok=False; break
        if ok: c+=1
    return c/Q

print("="*78)
print("(A) PRIME POINT-MASS BOUND  ML >= ceil((p-1)/26)/p,  p ∤ any core element")
print("="*78)
target=1/14
cores = {
  "ext j=6  {1..13}\\6∪56":  sorted(set([x for x in range(1,14) if x!=6]+[56])),
  "ext j=12 {1..13}\\12∪84": sorted(set([x for x in range(1,14) if x!=12]+[84])),
  "ext j=1  {2..13}∪14":     sorted(set(list(range(2,14))+[14])),
}
for name,S in cores.items():
    best=(0.0,None)
    for p in range(2,200):
        if not isprime(p): continue
        if any(v%p==0 for v in S): continue
        bound = -(-(p-1)//26)/p   # ceil((p-1)/26)/p
        if bound>best[0]: best=(bound,p)
    pct = 100*best[0]/target
    print(f"   {name:24s}: best ML-bound = {best[0]:.5f} at p={best[1]}  "
          f"({pct:.1f}% of 1/14={target:.5f})  {'REACHES 1/14' if best[0]>=target else 'SHORT of 1/14'}")
print("   --> the prime route (Bedert's right tool for AP-cores) lands JUST below 1/14.")

print("\n" + "="*78)
print("(B) STRANGER-LIMIT REDUCTION:  L({1..13}\\j ∪ {14m})  ->  (6/7)·meas(Lonely({1..13}\\j))")
print("="*78)
Q=210210  # 2*3*5*7*11*13*7? -> use a big multiple of 14 with many small factors: 2*3*5*7*11*13*14? keep =14*15015
Q=14*15015  # =210210, divisible by 14 and by 1..13 lcm-ish (15015=3*5*7*11*13)
print(f"   grid Q={Q}")
six_sevenths=6/7
for j in [1,6,12]:
    core12 = sorted(x for x in range(1,14) if x!=j)
    Lcore = lonely_measure(core12, Q)
    limit = six_sevenths*Lcore
    print(f"\n   j={j}: core12={core12}")
    print(f"        meas(Lonely(core12)) = {Lcore:.6f}   (6/7)·that = {limit:.6f}  [conjectured limit]")
    row=[]
    for m in [1,2,3,5,8,13]:
        S=sorted(set(core12+[14*m]))
        if len(S)!=13:
            row.append((m,None)); continue
        L=lonely_measure(S, Q)
        row.append((m,L))
    cells=", ".join(f"m={m}:L={L:.5f}" for m,L in row if L is not None)
    print(f"        L(core∪14m): {cells}")
    last = [L for m,L in row if L is not None][-1]
    print(f"        ratio L(m=13)/limit = {last/limit:.4f}  (->1 confirms the reduction)")

print("\n   IF L(m)->(6/7)·meas(Lonely(core12))>0 for every j, then inf over the infinite")
print("   extremizer family = (6/7)·min_j meas(Lonely({1..13}\\j)) > 0  — a FINITE certificate.")
print("\nDONE.")

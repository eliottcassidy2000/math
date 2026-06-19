#!/usr/bin/env python3
"""
lrc14_sector_thm534check_macmini_0618s7b.py  (mac-mini-2026-06-18-S7, ANGLE B / adversarial)

ADVERSARIAL CHECK of THM-534's exact dual functionals. Use the REAL g(t) duals:
  k=11,12,13 (R=2): L = 1 - (1/2)S_1 + (1/6)S_2          g=(t-3)(t-4)/12
  k=9,10     (R=3): L = 1 - (13/18)S_1 + (4/9)S_2 - (1/6)S_3   g=-(t-2)(t-3)(t-6)/36
  k=8        (R=4): L = 1 - S_1 + S_2 - (9/10)S_3 + (3/5)S_4   g=(t-1)(t-2)(t-4)(t-5)/40

CHECK 1: verify each g(t)>=1[t=0] on t=0..6 (the proof of the per-E bound). PROVED if yes.
CHECK 2: verify meas(S7(E)) <= L_y(E) on a bank (the Bonferroni inequality).
CHECK 3 (the live gap): does consec MAXIMIZE the CORRECT-degree L_y for k=8,9,10? Re-run the
  bounded-spread sweep with the RIGHT functional. The deg-2-functional beater
  {0,2,4,5,6,7,8,10,12,14} at k=10 -- does it beat consec on the deg-3 functional? And is it
  still <= cap_10? Widen the spread to stress THM-534's "consec maximizes" claim.
CHECK 4: my Angle-B value-add: the L_y-beaters at k=10,11 are 2-DILATES of perforated APs.
  Since meas(S7) AND each S_r is SCALE-invariant, a 2-dilate has the SAME L_y as its primitive
  {0,1,2,...}/gcd. So {0,2,4,5,6,7,8,10,12,14} is NOT primitive? gcd=1 (has odd 5,7). Check.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd, comb
sys.stdout.reconfigure(line_buffering=True)

def factorial_moments(E, R=6):
    E=sorted(set(E)); Enz=[e for e in E if e!=0]
    bps=set([F(0),F(1)])
    for e in Enz:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); S=[F(0)]*(R+1)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; L=x1-x0
        hit=set(int(7*e*xm)%7 for e in E)
        free=6-len([s for s in hit if s!=0])
        for r in range(R+1): S[r]+=L*comb(free,r)
    return S
def measS7_mom(S): return sum((-1)**r*S[r] for r in range(len(S)))

duals = {
  8:  (lambda S: F(1)-S[1]+S[2]-F(9,10)*S[3]+F(3,5)*S[4],  lambda t: F((t-1)*(t-2)*(t-4)*(t-5),40)),
  9:  (lambda S: F(1)-F(13,18)*S[1]+F(4,9)*S[2]-F(1,6)*S[3], lambda t: F(-(t-2)*(t-3)*(t-6),36)),
  10: (lambda S: F(1)-F(13,18)*S[1]+F(4,9)*S[2]-F(1,6)*S[3], lambda t: F(-(t-2)*(t-3)*(t-6),36)),
  11: (lambda S: F(1)-F(1,2)*S[1]+F(1,6)*S[2], lambda t: F((t-3)*(t-4),12)),
  12: (lambda S: F(1)-F(1,2)*S[1]+F(1,6)*S[2], lambda t: F((t-3)*(t-4),12)),
  13: (lambda S: F(1)-F(1,2)*S[1]+F(1,6)*S[2], lambda t: F((t-3)*(t-4),12)),
}
cap = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91), 12:F(6,7), 13:F(1)}

print("CHECK 1: g(t)>=1[t=0] on t=0..6 for each dual (the PROOF of per-E bound)")
for k in [8,9,10,11]:
    Ly,g=duals[k]
    gs=[g(t) for t in range(7)]
    ok=all(gs[t]>=(1 if t==0 else 0) for t in range(7))
    print(f"  k={k}: g(t)={[str(x) for x in gs]}  valid(g>=1[t=0])={ok}")

print()
print("CHECK 3 (LIVE GAP): does consec maximize the CORRECT-degree L_y for k=8,9,10? widen spread.")
def gen(k,maxE):
    out=[]
    for rest in itertools.combinations(range(1,maxE+1),k-1):
        E=(0,)+rest; g=0
        for e in E: g=gcd(g,e)
        if g!=1: continue
        out.append(E)
    return out
for k in [8,9,10]:
    Ly,_=duals[k]; ck=cap[k]
    AP=tuple(range(k)); apLy=Ly(factorial_moments(AP,6))
    for maxE in [k+6, k+8]:
        shapes=gen(k,maxE)
        beat=[]; overcap=[]
        mx=(apLy,AP)
        for E in shapes:
            v=Ly(factorial_moments(E,6))
            if v>mx[0]: mx=(v,E)
            if v>apLy: beat.append((v,E))
            if v>ck: overcap.append((v,E))
        beat.sort(reverse=True)
        print(f"  k={k}, box maxE<={maxE} ({len(shapes)} sets): L_y(AP)={float(apLy):.5f} cap={float(ck):.4f}")
        print(f"     consec_is_max={mx[1]==AP}; #beaters={len(beat)}; #OVER-CAP={len(overcap)}"
              + (f"  TOP beater {float(beat[0][0]):.5f} at {beat[0][1]}" if beat else ""))
        if overcap:
            overcap.sort(reverse=True)
            print(f"     *** OVER-CAP shapes (L_y>cap, THM-534 BREACH): {[(float(v),E) for v,E in overcap[:3]]}")

print()
print("CHECK 4: scale-invariance of L_y. {0,2,4,5,6,7,8,10,12,14} primitive? its L_y vs primitive.")
E=(0,2,4,5,6,7,8,10,12,14)
g=0
for e in E: g=gcd(g,e)
print(f"  gcd(E)={g} (primitive={g==1})")
Ly,_=duals[10]
print(f"  L_y_deg3(E)={float(Ly(factorial_moments(E,6))):.5f}, cap_10={float(cap[10]):.4f}")
print("\nDONE.")

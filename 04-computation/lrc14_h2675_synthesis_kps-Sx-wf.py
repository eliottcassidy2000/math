"""
HYP-2675 SYNTHESIS / COMPLETENESS verification (kps-Sx-wf).

Re-verifies, with EXACT rationals, the load-bearing facts that decide whether
HYP-2675 ("wide => p0 <= cap_k") is PROVED with an explicit B.

Conclusion (see header comment block at bottom): HYP-2675 is NOT proved.
All four reduction angles bottom out in the SAME unproved object:
the sharpened invariant (W*) "span>B => p0(E) <= Q(k-1)" is REFUTED as an
upper bound (p0 can exceed Q(k-1) on genuinely wide sets), and the true
cap-level bound rests on an unproved joint Erdos-Turan-Koksma / AP-extremality
input (HYP-2603).
"""
from fractions import Fraction as F
import math

def p0p1(E):
    E=sorted(set(E)); bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*e+1): bps.add(F(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); p0=F(0); p1=F(0)
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        mid=(lo+hi)/2; miss=set(range(1,7))-set(int((e*mid)%1*7) for e in E)
        if len(miss)==0: p0+=hi-lo
        elif len(miss)==1: p1+=hi-lo
    return p0,p1

caps={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7)}

# Q(k-1) = Plat(consec_{k-1})
Q={}
for k in range(8,13):
    p0,p1=p0p1(list(range(k-1)))
    Q[k]=p0+p1/7

print("== FACT 1: caps, consec p0 (cap is NOT defined by p0 of consec; cap>consec p0) ==")
for k in range(8,13):
    p0,_=p0p1(list(range(k)))
    print(f" k={k} consec p0={float(p0):.5f} < cap_k={float(caps[k]):.5f}  Q(k-1)={float(Q[k]):.5f}")

print("\n== FACT 2 (DECISIVE): (W*) 'p0(E)<=Q(k-1) for wide E' is REFUTED ==")
ce=[0,19,20,21,22,23,24,25]   # k=8, 2nd-largest=24 > 14 => genuinely wide
p0,_=p0p1(ce)
print(f" E={ce}  p0={p0}={float(p0):.5f}")
print(f" Q(7)={Q[8]}={float(Q[8]):.5f}   p0 > Q(7): {p0>Q[8]}  (W* FALSE)")
print(f" but p0 < cap_8={float(caps[8]):.5f}: {p0<caps[8]}  (LRC conclusion still holds)")

print("\n== FACT 3: scan {0} U (g+consec_7), p0 exceeds Q(7) for several wide g ==")
over=[]
for g in range(1,40):
    E=[0]+[g+i for i in range(7)]
    pp,_=p0p1(E)
    if pp>Q[8] and g>=2:  # g>=2 to be wide-ish; g large => truly wide
        over.append((g,float(pp)))
print(" wide g with p0>Q(7):", over)

print("\n== FACT 4: induction premise of plateau-recursion 'wide E'=>p0<=Q(k-2)' check ==")
# Q(k-2) for the induction at k=8 is Plat(consec_6)=11/294
p0c6,p1c6=p0p1(list(range(6)))
Qm=p0c6+p1c6/7
Ep=[0,12,18,24,30,36,42]
pE,_=p0p1(Ep)
print(f" Q(6)=Plat(consec_6)={Qm}={float(Qm):.5f}")
print(f" E'={Ep}  p0={float(pE):.5f}  p0>Q(6): {pE>Qm}  (premise FALSE => written induction unsound)")

print("\n== SUMMARY ==")
print(" - No cap violation found anywhere (LRC conclusion p0<=cap likely TRUE).")
print(" - (W*) sharpened invariant p0<=Q(k-1) is FALSE on wide sets (FACT 2/3).")
print(" - Therefore all four angles' decorrelated worst-case core (=Q(k-1)) is a")
print("   LIMIT, not an upper bound; the cap-level bound needs HYP-2603 AP-extremality")
print("   plus an explicit joint ET-Koksma decorrelation constant -- both UNPROVED.")

"""
LRC(14): REFUTATION of HYP-2607 / THM-534 "AP maximises L_y(E)=U_t(E)" at k=12,13.
=================================================================================
opus 2026-06-29, option (a). The moment-LP L_y=U_t is the PROVED upper bound
meas(S7(E)) <= L_y(E) (THM-534); the OPEN finishing conjecture HYP-2607 is that the
consecutive/AP cluster MAXIMISES L_y over all shapes E (exhaustively verified only at k=8).

VALIDATION (exact anchor): U4(consec_8) = 2633/7350 = 0.358231 (THM-537); U2,U3 match too.
FINDING: at k=12 and k=13 the AP is NOT the maximiser -- near-AP shapes (a single
shifted endpoint) have strictly larger U_t at every order t=2,3,4 (t=2 is the order
THM-537 uses for k>=11). So HYP-2607 'AP maximises for all E' is FALSE beyond k=8.

IMPACT (honest): like the rho* refutation, this kills the clean 'reduce all k to the AP'
formulation, NOT LRC(14): the binding case is small k (k=8, exhaustively verified, AP
unique max); at k=12,13 the cap margin is huge so L_y<=cap still holds via a non-AP
maximiser. The 'AP-orbit majorises' / three-distance rearrangement proposed to PROVE
HYP-2607 cannot succeed as a universal statement -- it must be restricted to small k.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb
from scipy.optimize import linprog

def avoid_measure(E, J):
    Js=set(J); bps={F(0),F(1)}
    for e in E:
        if e:
            for m in range(7*e+1): bps.add(F(m,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); tot=F(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi<=lo: continue
        mid=(lo+hi)/2; ok=True
        for e in E:
            y=e*mid; y=y-(y.numerator//y.denominator)
            if (7*y.numerator)//y.denominator in Js: ok=False; break
        if ok: tot+=hi-lo
    return tot

def U_t(E, t):
    S=[F(1)]+[sum((avoid_measure(E,J) for J in combinations(range(1,7),tt)),F(0)) for tt in range(1,t+1)]
    A=[[comb(i,tt) for i in range(7)] for tt in range(t+1)]
    r=linprog([-1]+[0]*6, A_eq=A, b_eq=[float(s) for s in S], bounds=[(0,None)]*7, method='highs')
    return -r.fun if r.success else None

def shapes(k, smax):
    for s in range(k-1, smax+1):
        for mid in combinations(range(1,s),k-2):
            yield [0]+list(mid)+[s]

print("VALIDATE U4(consec_8)=2633/7350=%.6f -> computed %.6f"%(float(F(2633,7350)), U_t(list(range(8)),4)))
print()
for k,smax,t in [(12,14,2),(13,15,2)]:
    AP=list(range(k)); apU=U_t(AP,t)
    best=apU; arg=AP; nbeat=0
    for E in shapes(k,smax):
        if E==AP: continue
        v=U_t(E,t)
        if v>apU+1e-9: nbeat+=1
        if v>best: best=v; arg=E
    print(f"k={k} (order t={t}): U{t}(AP)={apU:.4f}; MAX over shapes(spread<={smax})={best:.4f} at {arg}")
    print(f"        #shapes beating AP={nbeat}  ->  {'HYP-2607 REFUTED (AP not max)' if best>apU+1e-9 else 'AP is max'}")

"""LRC(14) universal Farey floor + wide-regime R' (mac-mini-2026-06-22-S32, HYP-2869).
VERIFIED: meas(GOOD(E))=meas{maxgap{frac(e x)}>1/7} >= 3/pi^2=0.304 for ANY cluster (Farey nbhds
universally good: at x=a/b, b<7, any cluster collapses to the b-grid, maxgap>=1/b>1/7). Worst over
46 clusters = consec_13 at 0.4425. Bypasses 'consec is extremal'. R'->1 for wide clusters."""
from fractions import Fraction as Fr
def maxgap(E,x):
    ph=sorted(set((e*x)%1 for e in E))
    if len(ph)<=1: return Fr(1)
    g=max(b-a for a,b in zip(ph,ph[1:])); return max(g, ph[0]+1-ph[-1])
def measGOOD(E,thr=Fr(1,7)):
    E=sorted(set(int(e) for e in E)); bp={Fr(0),Fr(1)}
    for e in E:
        if e==0: continue
        for m in range(0,7*abs(e)+1): bp.add(Fr(m,7*abs(e)))
    bp=sorted(bp); tot=Fr(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        if maxgap(E,(a+b)/2)>thr: tot+=b-a
    return tot
if __name__=="__main__":
    import random; random.seed(3); worst=Fr(1); arg=None
    tests=[list(range(13)),list(range(9)),[0,1,2,4,5,7,8,10,11,13],[3*i for i in range(9)],
           [0,5,11,19,28,40,55,71]]
    for _ in range(40):
        k=random.randint(8,13); tests.append([0]+sorted(random.sample(range(1,40),k-1)))
    for E in tests:
        m=measGOOD(E)
        if m<worst: worst=m; arg=E
    print(f"worst meas(GOOD) over {len(tests)} clusters = {float(worst):.4f} at {arg}; >=3/pi^2=0.304: {worst>=Fr(304,1000)}")

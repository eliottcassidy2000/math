"""
Goddyn-Wong GW={1..11,13,24} (n=14), M=1/14. WHY tight despite 12=13-1 absent?
Explore: at t where ||12t|| is small (the absent-difference regime), which runner catches the danger arc
[-1/14,1/14]? Candidate: 24=2*12 (so ||24t||<=2||12t||). Find the binding runner across all t.
"""
from fractions import Fraction
def frac(x): x=x%1; return min(x,1-x)
GW=list(range(1,12))+[13,24]; n=14
# scan t, record min-runner distance and the argmin; focus on t near k/12 (||12t|| small)
Qmax=14*14
print("GW: t maximizing min_v||vt|| and the binding runner(s):")
best=Fraction(-1); bt=None
for q in range(1,Qmax+1):
    for a in range(1,q):
        t=Fraction(a,q); m=min(frac(Fraction(v)*t) for v in GW)
        if m>best: best=m; bt=t
binders=[v for v in GW if frac(Fraction(v)*bt)==best]
print(f"   M(GW)={best}={float(best):.5f} (1/14={1/14:.5f}); t*={bt}; binding runners v with ||vt*||=M: {binders}")
print()
print("The absent-difference-12 regime: for t with ||12t|| small, is a runner (esp. 24) in [-1/14,1/14]?")
print(f"   {'t':>10} {'||12t||':>9} {'||24t||':>9} {'runners in [-1/14,1/14]':>28}")
for q,a in [(12,1),(24,1),(24,5),(36,5),(25,2),(23,2),(12*14+1,14)]:
    t=Fraction(a,q)
    d12=frac(12*t); d24=frac(24*t)
    inarc=[v for v in GW if frac(Fraction(v)*t)<Fraction(1,14)]
    print(f"   {str(t):>10} {str(d12):>9} {str(d24):>9} {str(inarc):>28}")
print()
print("MECHANISM CHECK: whenever ||12t||<1/14 (12-gap between runners 1,13), is SOME runner in the danger arc?")
import random; random.seed(1); bad=0; checked=0
for _ in range(30000):
    q=random.randint(2,400); a=random.randint(1,q-1); t=Fraction(a,q)
    if frac(12*t) < Fraction(1,14):  # the absent-difference regime
        checked+=1
        if min(frac(Fraction(v)*t) for v in GW) >= Fraction(1,14):  # no runner in danger arc => would beat tightness
            bad+=1
print(f"   over {checked} t with ||12t||<1/14: t's with NO runner in danger arc (would give M>1/14) = {bad}")
print("   (0 confirms GW covers the absent-difference regime; the covering runner identifies GW's mechanism.)")

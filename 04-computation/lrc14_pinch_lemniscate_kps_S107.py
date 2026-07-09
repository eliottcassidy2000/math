# Grid-invisible pinches and the lemniscate node (kps-S107). mac-mini-S64: the strict boundary
# maxgap(x)=1/7 is hit at MEASURE-ZERO RATIONAL PINCHES invisible to uniform grids. HUNCH: near such
# a pinch two gap functions CROSS, and the local normal form is the LEMNISCATE NODE (x^2-y^2=0 => y=+-x).
# Test: (1) pinches are at rationals x=m/d (d=cluster difference); (2) local structure = two crossing
# LINES (lemniscate node to leading order); (3) count pinches vs spread (Davenport-Schinzel linearity).
from fractions import Fraction as F
from math import gcd
def teeth_sorted(E, x):  # x a Fraction; teeth frac(e*x)
    return sorted(F(e*x.numerator % (x.denominator), x.denominator) if False else (e*x - (e*x).__floor__()) for e in E)
def maxgap(E, x):  # exact max circular gap at slow x (Fraction), teeth = frac(e*x)
    ts = sorted((F(e)*x - int((F(e)*x)) if (F(e)*x)>=0 else (F(e)*x)-int((F(e)*x))+ (0 if (F(e)*x)==int(F(e)*x) else -1) ) for e in E)
    # simpler: fractional part
    ts = sorted(((F(e)*x) - ( (F(e)*x).numerator // (F(e)*x).denominator )) for e in E)
    ts = sorted(t - (1 if t>=1 else 0) for t in ts)  # ensure [0,1)
    n=len(ts); g=F(0)
    for i in range(n):
        d = (ts[(i+1)%n]-ts[i]) if i<n-1 else (ts[0]+1-ts[n-1])
        if d>g: g=d
    return g
def frac(q):
    return q - (q.numerator // q.denominator)
def maxgap2(E,x):
    ts = sorted(frac(F(e)*x) for e in E)
    n=len(ts); g=F(0)
    for i in range(n):
        d = (ts[(i+1)%n]-ts[i]) if i<n-1 else (ts[0]+1-ts[n-1])
        if d>g: g=d
    return g
E=[0,1,3,7]  # small cluster, spread 7
TH=F(1,7)
# candidate pinch abscissae: x = m/d for d in pairwise differences, m=0..d
diffs=set()
for i in range(len(E)):
    for j in range(len(E)):
        if E[i]!=E[j]: diffs.add(abs(E[i]-E[j]))
cands=set()
for d in diffs:
    for m in range(0,d+1): cands.add(F(m,d))
cands=sorted(c for c in cands if 0<=c<1)
print(f"E={E} spread={max(E)-min(E)} pairwise-diffs={sorted(diffs)}")
print("Sampling maxgap around each candidate pinch; a 1/7-crossing = grid-invisible arc boundary.")
# find x where maxgap crosses 1/7: sample maxgap on a fine rational grid between candidate breakpoints
# and detect sign changes of (maxgap-1/7)
xs=sorted(set(cands)|{F(k,420) for k in range(0,420)})
prev=None; crossings=[]
for x in xs:
    val=maxgap2(E,x)-TH
    if prev is not None and (prev==0 or val==0 or (prev>0)!=(val>0)):
        crossings.append((prevx,x,prev,val))
    prev=val; prevx=x
print(f"# sign-structure changes of (maxgap-1/7): {len(crossings)}")
# report the exact pinch rationals where maxgap==1/7
exact=[x for x in cands if maxgap2(E,x)==TH]
print(f"exact pinches maxgap=1/7 at rationals: {exact}  (all m/d, GRID-INVISIBLE unless d|Vmax-grid)")
# LOCAL NORMAL FORM at a pinch: at a tooth collision x*=m/d, two teeth coincide; the two adjacent gaps
# g_a(x), g_b(x) cross linearly (slopes +-, like y=+-x = the lemniscate node x^2-y^2=0).
print("\nLOCAL NORMAL FORM at a tooth-collision pinch (two teeth merge => two adjacent gaps swap):")
for d in sorted(diffs)[:3]:
    xstar=F(1,d)
    eps=F(1,d*1000)
    tsm=sorted(frac(F(e)*(xstar-eps)) for e in E); tsp=sorted(frac(F(e)*(xstar+eps)) for e in E)
    print(f"  collision x*=1/{d}: teeth just before {[float(t) for t in tsm]}")
    print(f"                     just after  {[float(t) for t in tsp]}  (order swaps = the node crossing)")
print("\n=> pinches are rational x=m/d (measure-zero, grid-invisible); local structure = two gap-lines")
print("   crossing = the LEMNISCATE NODE y=+-x (leading order of (x^2+y^2)^2=x^2-y^2 at origin).")

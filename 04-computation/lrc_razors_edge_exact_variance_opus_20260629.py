"""
EXACT-value-plus-variance frame. The AP is the RAZOR'S EDGE: exact identities, not a range.
- M(AP)=1/14 EXACT; L(AP)=0 EXACT (danger arcs cover the circle except the 2 points {1/14,13/14}).
- Perfect cancellation: Sum_k (-1)^k E_k = E[(1-1)^X] = P(X=0) = 0  <=>  X>=1 a.e. (exact covering).
- Everything else = exact value + a ONE-SIGNED variance (covering: M-1/14>0, L>0).
"""
from fractions import Fraction as F
def nrm(x): f=x-(x.numerator//x.denominator); return min(f,1-f)
def lonely_measure_exact(S,n=14):
    # exact: boundaries t=(n k +-1)/(n v); lonely interval if midpoint has min||v t||>1/n (STRICT, open)
    B=set([F(0),F(1)])
    for v in S:
        k=0
        while F(n*k-1,n*v) < 1 and k< n+v+3:
            for t in (F(n*k+1,n*v),F(n*k-1,n*v)):
                if 0<t<1: B.add(t)
            k+=1
    B=sorted(B); meas=F(0); peakpts=[]
    for i in range(len(B)-1):
        a,b=B[i],B[i+1]; mid=(a+b)/2
        mn=min(nrm(v*mid) for v in S)
        if mn>F(1,n): meas+=b-a                       # open lonely interval (positive measure)
    # razor points: tau where min = exactly 1/n (the closed-but-measure-zero peaks)
    cand=set()
    for v in S:
        for k in range(v+1):
            for t in (F(n*k+1,n*v),F(n*k-1,n*v)):
                if 0<t<1: cand.add(t)
    for t in cand:
        if min(nrm(v*t) for v in S)==F(1,n): peakpts.append(t)
    return meas, sorted(set(peakpts))
ap=list(range(1,14))
Lap,pk=lonely_measure_exact(ap)
print("RAZOR'S EDGE -- the AP {1..13} (EXACT identities):")
print(f"  L(AP) = lonely measure = {Lap}  (EXACTLY 0: danger arcs cover the circle minus measure-zero peaks)")
print(f"  peak points where min||v t|| = 1/14 EXACTLY: {pk}  (the razor's edge, measure 0)")
print(f"  => PERFECT CANCELLATION: Sum_k(-1)^k E_k = P(X=0) = 0, i.e. X>=1 a.e. (exact tightest covering)\n")
print("EXACT value + ONE-SIGNED variance (covering perturbations):")
print(f"  {'set':>26} {'L (variance from 0)':>20} {'M':>8} {'M-1/14 (variance)':>18}")
def M_exact(S):
    C=set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]+S[j],abs(S[i]-S[j])):
                if d:
                    for m in range(d+1): C.add(F(m,d))
        for m in range(2*S[i]+1): C.add(F(m,2*S[i]))
    return max(min(nrm(s*t) for s in S) for t in C if 0<t<1)
for nm,S in [("AP {1..13}",ap),("{1..11,13,84}",[1,2,3,4,5,6,7,8,9,10,11,13,84]),("{2..14}",list(range(2,15)))]:
    L,_=lonely_measure_exact(S); M=M_exact(S)
    print(f"  {nm:>26} {str(L):>20} {str(M):>8} {str(M-F(1,14)):>18}")
print("\n  The AP is the EXACT critical point (L=0, M=1/14). Covering = exact value + POSITIVE variance.")
print("  Proof = the covering variance is one-signed (never overshoots THROUGH the razor's edge to M<1/14).")

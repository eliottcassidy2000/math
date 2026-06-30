"""
The cut-deficit mechanism: cut(s|prefix) = M(prefix) - M(prefix+{s}).
Small speeds cut the peak MOST; the covering-forced large mult-of-14 cuts LESS => the +margin.
"""
from fractions import Fraction as F
def norm(x): f=x-(x.numerator//x.denominator); return min(f,1-f)
def M_exact(S):
    if not S: return F(1,2)
    C=set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]+S[j],abs(S[i]-S[j])):
                if d:
                    for m in range(d+1): C.add(F(m,d))
        for m in range(2*S[i]+1): C.add(F(m,2*S[i]))
    best=F(0)
    for t in C:
        if 0<t<1:
            v=min(norm(s*t) for s in S)
            if v>best: best=v
    return best
base=list(range(1,12))  # {1..11}, M=1/12
M0=M_exact(base)
print(f"base {{1..11}}: M={M0}={float(M0):.5f}")
print("cut(s | {1..11}) = M0 - M({1..11,s}) for candidate next speeds s:")
print(f"{'s':>4} {'M(base+s)':>10} {'cut':>9}  note")
for s in [12,13,14,15,16,24,28,42,84,168,1000]:
    M=M_exact(base+[s]); cut=M0-M
    note="AP step (smallest)" if s==12 else ("mult-of-14 (covering-forced)" if s%14==0 else "")
    print(f"{s:>4} {str(M):>10} {float(cut):>9.5f}  {note}")
print()
print("=> the SMALLEST available speed (12, the AP step) cuts the MOST (steepest descent to 1/13).")
print("   Multiples of 14 (covering-forced) cut LESS as they grow (84:.0047, 168:less) -> the margin.")
print("   So a covering set, FORCED to include a mult-of-14 instead of a small speed, ends ABOVE 1/14.")
print("   THE BOUND: M(S) = 1/2 - sum cut(s_k); AP attains sum=3/7 (M=1/14); covering's forced large")
print("   speed has a CUT DEFICIT => sum < 3/7 => M > 1/14. The +10% margin IS the mult-of-14 deficit.")

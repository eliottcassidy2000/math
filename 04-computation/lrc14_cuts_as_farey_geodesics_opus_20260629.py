"""
Cuts as rational geodesics on the additive Farey skeleton.
RESONANCE LAW: a speed s cuts the witness a/q  iff  s*a mod q is small (exact hit iff q | s*a, i.e. q|s).
The AP traces the Farey LEFT SPINE 1/2->1/3->...->1/14 (consecutive Farey neighbors, det=1),
speed k self-resonant with denominator k. The tight witness 1/14 SURVIVES iff no speed is a mult of 14.
A covering set is FORCED to include a mult-of-14 => it KILLS 1/14 => witness pushed to a deeper Farey
rational (denom 89) => the +margin.
"""
from fractions import Fraction as F
from math import gcd
def norm(x): f=x-(x.numerator//x.denominator); return min(f,1-f)
def M_arg(S):
    if not S: return F(1,2),F(1,2)
    C=set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]+S[j],abs(S[i]-S[j])):
                if d:
                    for m in range(d+1): C.add(F(m,d))
        for m in range(2*S[i]+1): C.add(F(m,2*S[i]))
    best=F(0);arg=F(0)
    for t in C:
        if 0<t<1:
            v=min(norm(s*t) for s in S)
            if v>best: best,arg=v,t
    return best,arg
def fareypath(S,label):
    print(f"\n{label}: witness path along the Farey skeleton")
    print(f"  {'prefix':>16} {'M':>7} {'witness a/q':>11} {'q':>4} {'next s':>6} {'q|s?':>5} {'cut':>8}")
    prevq=None
    for k in range(len(S)+1):
        M,t=M_arg(S[:k]); q=t.denominator
        nxt=S[k] if k<len(S) else None
        hit = (nxt is not None and (nxt*t.numerator)%q==0)
        if k<len(S):
            M2,_=M_arg(S[:k+1]); cut=M-M2
        else: cut=F(0)
        det=""
        if prevq is not None: det=""  # Farey-neighbor check below
        print(f"  {str(S[:k]) if k<4 else '{1..'+str(S[k-1])+'}':>16} {float(M):>7.4f} {str(t):>11} {q:>4} "
              f"{str(nxt):>6} {str(hit):>5} {float(cut):>8.5f}")
        prevq=q
ap=list(range(1,14))
cov=[1,2,3,4,5,6,7,8,9,10,11,13,84]
fareypath(ap,"AP {1..13} (additive spine)")
fareypath(cov,"COVERING {1..11,13,84}")
# Farey-neighbor check on AP spine 1/2,1/3,...,1/14
print("\nFarey-neighbor check on AP spine (|a_i b_{i+1} - a_{i+1} b_i| should be 1):")
spine=[F(1,q) for q in range(2,15)]
ok=all(abs(spine[i].numerator*spine[i+1].denominator-spine[i+1].numerator*spine[i].denominator)==1 for i in range(len(spine)-1))
print(f"  1/2,1/3,...,1/14 are consecutive Farey neighbors (geodesic): {ok}")
# Tight-witness killing: which speeds kill 1/14?
print("\nThe tight witness 1/14 (k/14): killed by speed s iff ||s/14||<1/14 iff 14|s.")
for s in [13,14,28,84,15,27]:
    killed = norm(F(s,14))<F(1,14)
    print(f"  s={s:>3}: ||{s}/14|| = {str(norm(F(s,14))):>5} {'< 1/14 -> KILLS 1/14' if killed else '>=1/14 -> 1/14 survives'}  (14|s: {s%14==0})")

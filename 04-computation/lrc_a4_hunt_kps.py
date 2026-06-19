import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
def fn(x):
    r=x-(x.numerator//x.denominator); return min(r,1-r)
def maxmin(S):
    S=sorted(set(S)); cand=set()
    for i in range(len(S)):
        for j in range(i,len(S)):
            for d in (S[i]+S[j],S[i]-S[j]):
                if d==0: continue
                d=abs(d)
                for m in range(1,d): cand.add(Fraction(m,d))
    best=Fraction(0); bt=Fraction(0)
    for t in cand:
        mv=min(fn(v*t) for v in S)
        if mv>best: best,bt=mv,t
    return best,bt
# a=4 target 4/(4k+3). Try structured families that kill t=1/(k-1) and t=1/(k-2).
print("target a=4: 4/(4k+3).  Testing structured 'remove 2, add 2 large' families across k (k%12 noted)")
fams = {
 "rm{k-1}+4(k-1)":        lambda k: list(range(1,k-1))+[k,4*(k-1)],
 "rm{k-1,k-2}+4(k-1),k":  lambda k: list(range(1,k-2))+[k-... ] if False else None,
}
def test(k):
    floor=Fraction(1,k+1); a4=Fraction(4,4*k+3); a3=Fraction(3,3*k+2)
    cands=[]
    # family 1: {1..k-2,k,4(k-1)}
    cands.append(("rm(k-1)+4(k-1)", sorted(set(list(range(1,k-1))+[k,4*(k-1)]))))
    # family 2: {1..k-3,k-1,k,4(k-2)}  (remove k-2, quadruple it)
    cands.append(("rm(k-2)+4(k-2)", sorted(set(list(range(1,k-2))+[k-1,k,4*(k-2)]))))
    # family 3: {1..k-3,k-1,k, 4(k-1)} 
    cands.append(("rm(k-2)+4(k-1)", sorted(set(list(range(1,k-2))+[k-1,k,4*(k-1)]))))
    # family 4: {1..k-3, k-1, 3(k-1), 4(k-1)}? two large
    cands.append(("rm(k-2)+3,4(k-1)", sorted(set(list(range(1,k-2))+[k-1,3*(k-1),4*(k-1)]))))
    out=[]
    for name,S in cands:
        if len(S)!=k: continue
        M,t=maxmin(tuple(S))
        a = M.numerator if (M.numerator*(k+1)-M.denominator)==1 else None
        if M<a3:  # below the a=3 level => a>=4 territory
            out.append((name,M,a,S))
    return out,a4,a3
for k in range(7,40):
    out,a4,a3=test(k)
    hits=[(n,str(M),a) for n,M,a,S in out]
    if hits:
        print(f"  k={k:2d}(%12={k%12},%6={k%6}): BELOW a=3! {hits}  [a=4 target {a4}]")

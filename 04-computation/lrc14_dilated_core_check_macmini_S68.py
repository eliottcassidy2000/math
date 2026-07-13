#!/usr/bin/env python3
"""mac-mini-S68 CRITICAL: does a DILATED core c*{1..12} + killer 182 break the 14/183 floor?
opus balance predicts M(c*{1..12},182) = 14/(182+c) < 14/183 for c>1. Test M EXACTLY.
If M < 14/183 => floor broken (major). If M >= 14/183 => balance is a loose lower bound; the
true M uses a better (non-balance) witness => floor holds & the large-s trade is a mirage."""
from fractions import Fraction as F
from math import gcd

def M_exact(S, Qmax):
    """exact M(S)=max_{q<=Qmax,a} (1/q) min_i ||a v_i||_q. Local maxima of min-dist are at such
    rationals; Qmax must exceed the true witness denominator. Report argmax too."""
    best=F(0); arg=None
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            mind=q
            for v in S:
                r=(a*v)%q; d=r if r<=q-r else q-r
                if d<mind: mind=d
                if d==0: break
            if mind>0 and F(mind,q)>best: best=F(mind,q); arg=(q,a,mind)
    return best,arg

def is_covering(S,n=14): return all(any(v%q==0 for v in S) for q in range(2,n+1))
def is_primitive(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1

target=F(14,183)
print(f"floor 14/183 = {float(target):.6f}\n")
print(f"{'S':52s} | prim | cover | balance pred | M_exact | verdict")
print("-"*112)
for c in [1,2,3,5,7,11]:
    core=[c*k for k in range(1,13)]
    S=sorted(core+[182])
    prim=is_primitive(S); cov=is_covering(S)
    pred=F(14,182+c)  # opus balance = (1/13)*182/(182+c)
    # need Qmax > witness denom; deep well is 183; dilated may be larger. Use generous Qmax.
    Qmax=max(2*max(S), 500)
    M,arg=M_exact(S,Qmax)
    verdict = "BELOW FLOOR!" if M<target else ("=floor" if M==target else "above")
    print(f"{str(S):52s} | {str(prim):5s}| {str(cov):5s}| {str(pred):>10}={float(pred):.5f} | {float(M):.5f}={M} | {verdict}  (q*={arg[0] if arg else '?'})")

print("\nInterpretation: if M_exact >= 14/183 despite balance predicting less, opus's balance is")
print("a LOWER bound realized only at the core-resonance witness; the TRUE M(S) uses a better")
print("base (larger q) => the large-s/dilated-core trade is a MIRAGE and 14/183 stands.")
print("Also check: is the dilated-core config even primitive+covering, or excluded upstream?")

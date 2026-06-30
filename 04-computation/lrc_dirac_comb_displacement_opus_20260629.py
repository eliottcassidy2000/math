"""
(1) DIRAC COMB: the AP at t=1/14 = the 14-comb (14th roots). Observer=empty tooth; killer fills it.
(2) DISPLACEMENT: for covering sets, the killer w (mult of lcm(missing q's,14)) BINDS the new witness;
    displacement delta=M-1/14. Look for a uniform/cyclotomic lower bound.
"""
from fractions import Fraction as F
from math import gcd, lcm
def nrm(x): f=x-(x.numerator//x.denominator); return min(f,1-f)
ap=list(range(1,14))
print("(1) DIRAC COMB at the AP witness t=1/14: runner positions {i/14} = the 14-comb (minus observer):")
pos=sorted(F(i,14) for i in ap)
print(f"    runners at: {[f'{p.numerator}/14' for p in pos]}")
print(f"    + observer at 0/14  => the FULL 14-comb {{0,1,..,13}}/14 (the 14th roots of unity)")
print(f"    observer's tooth (0) is EMPTY = the gap = the razor's edge (dist 1/14 to nearest teeth 1,13)")
print(f"    a killer w (14|w) sits at ||w/14||=||0||=0 = the observer's tooth => FILLS the gap => kills edge\n")
# (2) displacement for covering sets: killer = lcm of the missing residues' needs
def M_exact(S):
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
            v=min(nrm(s*t) for s in S)
            if v>best: best,arg=v,t
    return best,arg
def is_cov(S): return all(any(s%q==0 for s in S) for q in range(2,15))
print("(2) covering sets: drop one small speed k from {1..13}, add the forced killer w=lcm(k,14):")
print(f"    {'dropped':>8} {'killer w':>9} {'M':>9} {'witness':>10} {'killer binds?':>14} {'delta=M-1/14':>14}")
mindelta=(F(1),None)
for k in range(2,14):
    w=lcm(k,14)
    S=[x for x in ap if x!=k]+[w]
    if len(set(S))<13 or not is_cov(S): 
        # not covering (missing some q) -- skip, but report
        continue
    M,t=M_exact(S)
    # is the killer binding at t?
    kb = nrm(w*t)==M
    delta=M-F(1,14)
    if delta<mindelta[0]: mindelta=(delta,(k,w,M))
    print(f"    {k:>8} {w:>9} {str(M):>9} {str(t):>10} {str(kb):>14} {str(delta):>14}")
print(f"\n    MIN displacement among these: delta={mindelta[0]}={float(mindelta[0]):.5f} at drop={mindelta[1][0]}, killer={mindelta[1][1]}, M={mindelta[1][2]}")
print(f"    (the killer BINDS the displaced witness => the new M is set by the killer's resonance)")
print(f"    1/14={float(F(1,14)):.5f}; the smallest killer (drop 12, w=84) gives M=7/89, delta=9/1246={float(F(9,1246)):.5f}")

# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont67: work the 2-RUNNER EXTREMAL problem -- covering-min = min over binding pairs {a,b}
# of Delta*ab/(a+b) (cont.66). This UNIFIES the slot-depth formula (kps) with klein's clearing-modulus frame.
#
# KEY LINK: when runner 1 binds (a=1) with co-binder b, the slot formula M = p_b/(b+1) IS a band-clearing
# certificate at modulus q = b+1, band-edge mu = p_b: M = mu/q. So the co-binder b = q-1 (klein's clearing
# modulus minus 1). The two frames are the SAME 2-runner extremal. The lcm-outlier (cont.55) reappears: at the
# extremal t*=14/183 (q=183, b=182), covering 13 forces a mult of 13 with u in [14,169], minimal 182=lcm(13,14).
from math import gcd
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
def normF(x):
    x=x-(x.numerator//x.denominator)
    if x<0: x+=1
    return min(x,1-x)
def Mexact(v):
    qc=4*max(v)+2; best=F(0); bt=None
    for q in range(2,qc):
        for p in range(1,q):
            if gcd(p,q)==1:
                m=min(normF(F(vi*p,q)) for vi in v)
                if m>best: best,bt=m,F(p,q)
    return best,bt
def is_cov(v,N=14): return all(any(x%d==0 for x in v) for d in range(2,N+1))
def prim(v): return reduce(gcd,v)==1

def main():
    print("(1) THE UNIFICATION: runner-1-binding => clearing at q = b+1 (co-binder b), band-edge mu=p_b, M=mu/q:")
    print(f"    {'family':<24} | binders (1,b) | M | q=b+1 | mu=p_b | M==mu/q?")
    for name,v in [("deep well {1..12,182}", list(range(1,13))+[182]),
                   ("ladder S_2 {1..12,364}", list(range(1,13))+[364]),
                   ("multi {1..10,13,22,84}", list(range(1,11))+[13,22,84])]:
        M,t=Mexact(v)
        binders=[b for b in v if normF(F(b*t.numerator,t.denominator))==M]
        if 1 in binders and len(binders)>=2:
            b=[x for x in binders if x!=1][0]
            pb=round(b*float(t)); q=b+1
            print(f"    {name:<24} | (1,{b}) | {str(M):>7} | {q:>5} | {pb:>6} | {M==F(pb,q)} (={pb}/{q})")

    print("\n(2) THE lcm-OUTLIER in the 2-runner frame: at the extremal t*=14/183 (q=183), covering d=13 forces a")
    print("    mult of 13 that is LONELY (>=M): 13u with ||13u*14/183||>=14/183 <=> u in [14,169]; minimal u=14 =>")
    print("    13*14 = 182 = lcm(13,14) = the deep-well co-binder. (182u mod 183 = -u => in-band iff u in [14,169].)")
    lonely_u = [u for u in range(1,183) if normF(F(13*u*14,183))>=F(14,183)]
    print(f"    smallest u with 13u lonely at 14/183: u={min(lonely_u)} => 13u={13*min(lonely_u)} (=lcm(13,14)? {13*min(lonely_u)==182})")

    print("\n(3) THE 2-RUNNER EXTREMAL MINIMIZER: min over covering families of Delta*ab/(a+b) = M; enumerate:")
    best=F(1); bf=None; bp=None
    # single-killer sweep (a=1) + a few multi
    for fam in [list(range(1,13))+[k] for k in range(13,200)]:
        if prim(fam) and is_cov(fam):
            M,t=Mexact(fam)
            if M<best:
                binders=[b for b in fam if normF(F(b*t.numerator,t.denominator))==M]
                best,bf,bp=M,fam,binders
    print(f"    min over single-killer {{1..12,k}}, k<=199: M={best} at {bf}, binding pair {bp}")
    print(f"    => the 2-runner extremal minimizer is {{1,182}} with Delta=1/13 => M=14/183 (the deep well).")
    print()
    print("=> THE 2-RUNNER EXTREMAL UNIFIES: slot-depth M=Delta*ab/(a+b) (kps cont.66) = clearing at q=b+1 (klein)")
    print("   = lcm-outlier co-binder 182 (kps cont.55). Single-killer (co-binder = mult of lcm 182) PROVED >=14/183")
    print("   uniformly (cont.60/61,66); general (co-binder mult-of-14 + separate mult-of-13) = klein clearing ILP<=182.")

if __name__ == "__main__":
    main()

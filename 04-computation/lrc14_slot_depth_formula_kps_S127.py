# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont66: bound the backbone slot depth uniformly in m via the SLOT-DEPTH FORMULA.
#
# cont.65 backbone lens: at the covering-min t*, exactly 2 runners bind (equioscillation). The depth of the
# backbone slot = the 2-runner balance value. DERIVE it: binders a,b (speeds), arcs at p_a/a, p_b/b (p=round(a t*)),
# t* between them, one rising (slope a) one falling (slope b): a(t*-p_a/a)=M, b(p_b/b-t*)=M
#   => M(1/a+1/b) = p_b/b - p_a/a =: Delta   =>   M = Delta / (1/a + 1/b) = Delta * a*b/(a+b).
# UNIFORM IN m: for the RUNNER-1-contested (single-killer) slot, a=1, backbone b=14m at arc 1/13 (m mult of 13),
# Delta=1/13 => M = 14m/(13(14m+1)) >= 14/183 for all m>=13 (min at m=13, the deep well). For FILLER-contested
# (multi-killer), a=filler>=2; the slot is looser (deeper) but the uniform bound is the open analytic part.
from math import gcd
from fractions import Fraction as F
def normF(x):
    x = x - (x.numerator // x.denominator)
    if x < 0: x += 1
    return min(x, 1 - x)
def Mexact(v):
    qc = 4*max(v)+2; best=F(0); bt=None
    for q in range(2,qc):
        for p in range(1,q):
            if gcd(p,q)==1:
                m=min(normF(F(vi*p,q)) for vi in v)
                if m>best: best,bt=m,F(p,q)
    return best,bt

def slot_formula(v):
    M,t = Mexact(v)
    binders = [b for b in v if normF(F(b*t.numerator, t.denominator))==M]
    if len(binders) < 2: return M,t,binders,None,None
    a,b = binders[0], binders[1]
    pa = round(a*float(t)); pb = round(b*float(t))
    Delta = F(pb,b) - F(pa,a)
    D = abs(Delta) / (F(1,a)+F(1,b))
    return M,t,(a,b),(F(pa,a),F(pb,b),abs(Delta)),D

def main():
    print("SLOT-DEPTH FORMULA: M = Delta/(1/s_a+1/s_b), Delta = gap between the two binders' arcs. Verify M==formula:")
    print(f"{'family':<26} | M | binders (a,b) | arcs | Delta | formula D | D==M?")
    fams = {
        "deep well {1..12,182}": list(range(1,13))+[182],
        "ladder S_2 {1..12,364}": list(range(1,13))+[364],
        "multi {1..11,13,84}": list(range(1,12))+[13,84],
        "multi {1..10,13,22,84}": list(range(1,11))+[13,22,84],
    }
    for name,v in fams.items():
        M,t,binders,arcs,D = slot_formula(v)
        if D is None:
            print(f"{name:<26} | {str(M):>7} | {binders} (need 2)")
            continue
        pa,pb,Delta = arcs
        print(f"{name:<26} | {str(M):>7} | {binders} | {pa},{pb} | {str(Delta):>5} | {str(D):>7} | {D==M}")

    print("\nUNIFORM-IN-m bound, RUNNER-1-contested (single-killer) slot: a=1, backbone 14m at arc 1/13, Delta=1/13:")
    print(f"    {'m':>4} | backbone 14m | M = 14m/(13(14m+1)) | >= 14/183?")
    tgt = F(14,183)
    for m in [13, 26, 39, 130, 1300]:
        Mm = F(14*m, 13*(14*m+1))
        print(f"    {m:>4} | {14*m:>11} | {str(Mm):>18}={float(Mm):.6f} | {Mm>=tgt} (min at m=13: {Mm==tgt})")
    print("    => 14m/(13(14m+1)) >= 14/183 <=> 183m >= 13(14m+1)=182m+13 <=> m >= 13. PROVED uniformly in m>=13")
    print("       (single-killer needs m mult of 13 => m>=13; min at m=13 = deep well = 14/183). = the Lean ladder.")
    print()
    print("FILLER-contested (multi-killer): a=filler>=2, so M = Delta*a*b/(a+b) with a>=2. Examples above:")
    print("   {5,84} Delta=1/60 -> 7/89=0.0787;  {1,22} Delta -> 2/23=0.087 -- both LOOSER (deeper) than 14/183.")
    print("   Faster binder (a>=2 vs runner-1's a=1) => shallower-denominator balance => deeper slot. But the")
    print("   uniform bound needs the equioscillation constraint (all others >= M forces Delta*a*b/(a+b)>=14/183)")
    print("   -- the open analytic/finite-check content (opus Fourier). The formula LOCALIZES the crux to a 2-runner balance.")

if __name__ == "__main__":
    main()

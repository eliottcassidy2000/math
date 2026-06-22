"""
The sharpest single-speed disproof (kps-S31p round 2): {1..12, v}. One speed v can kill b=13
(v==0 mod13) OR b=14 (v==0 mod14) but to kill BOTH needs v==0 mod lcm(13,14)=182 -- huge.
Then {1..12,182} kills b=1..14 at EXACT points, smallest exact-survivor b=15 => M<=1/15<1/14?
PROOF response: 182 is huge => equidistributes => thin arc misses the b=13,14 NEIGHBORHOODS =>
they survive => M ~ 1/13. Verify the whole spectrum; min over v must be 1/14 (the AP, v=13).
"""
from fractions import Fraction as F
def nf(x):
    r=x%1; return min(r,1-r)
def M_wit(S):
    S=sorted(set(abs(s) for s in S if s!=0)); C=set()
    for i in range(len(S)):
        for j in range(i,len(S)):
            for comb in {S[i]+S[j], abs(S[i]-S[j])}:
                if comb:
                    for k in range(1,comb): C.add(F(k,comb))
    best=F(0); arg=None
    for t in C:
        if 0<t<1:
            m=min(nf(s*t) for s in S)
            if m>best: best=m; arg=t
    return best,arg
print("{1..12, v} SPECTRUM -- which resonances does v kill, and resulting M:  1/14 =",float(F(1,14)))
base=list(range(1,13))
for v in [13,14,26,28,39,42,91,98,182,364,2366]:
    S=base+[v]; M,t=M_wit(S)
    kills=[b for b in (13,14) if v%b==0]
    print(f"  v={v:5d} (kills b={kills if kills else 'neither'} of 13,14): M={str(M):>7s}={float(M):.4f} "
          f"{'<1/14!!' if M<F(1,14) else ''}  witness b={t.denominator}")
mn=min(M_wit(base+[v])[0] for v in range(13,400))
print(f"\n  MIN over v in [13,399] of M({{1..12,v}}) = {mn} = {float(mn):.4f}  (at v=13 = the AP)")
print(f"  => {'COUNTEREXAMPLE' if mn<F(1,14) else 'floor 1/14 HOLDS'}: one speed kills b=13 XOR b=14 cheaply;")
print("     killing BOTH needs v=182 (huge) which EQUIDISTRIBUTES and misses the neighborhoods (M~1/13).")
print("     13 speeds give 13 cheap killers -- exactly enough for b=1..13, never b=1..14. The AP is tight.")

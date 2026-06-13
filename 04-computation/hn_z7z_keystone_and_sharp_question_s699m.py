"""Extend S687. (A) SHARP QUESTION: are the depth-PGF roots at ±2π/3? Compute angles; clarify the
LRC cyclotomic/±2π/3 structure lives in the WITNESS clock points (n-th roots, THM-403) & the
worry-set round-tournament F(T,x) (S44), NOT the depth PGF. (B) z^7−z KEYSTONE: 7-coloring =
ℤ[ζ6]/(2+ζ6)≅𝔽₇ (Eisenstein prime of norm 7); W6 ℂ-roots incl ±2π/3 = the χ=3 gadget; unifies
χ≤7 (𝔽₇ palette) with the forbidden tournament value 7=Φ3(2)=N(2+ζ6)=|PG(2,2)|. (C) field/Heegner
tower. opus-2026-06-06-S699m."""
import cmath, math
from fractions import Fraction as F
try:
    import numpy as np; HAVE=True
except Exception: HAVE=False
def dist(x): x%=1; return min(x,1-x)
def depth_dist(V,delta):
    eps=set([F(0)])
    for v in V:
        for k in range(v+1):
            for s in(1,-1): eps.add((F(k)+s*delta)/v % 1)
    pts=sorted(eps); p={}; L=len(pts)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        d=sum(1 for v in V if dist(v*mid)<delta); p[d]=p.get(d,F(0))+ln
    return p
def Neis(a,b): return a*a+a*b+b*b   # norm in Z[zeta6], zeta6=e^{iπ/3}: N(a+bζ)=a²+ab+b²
def main():
    print("(A) SHARP QUESTION — depth-PGF roots of the worry-set: at ±2π/3 (=120°)?")
    for V,n in [((1,2,3,4),5),((1,2,3,4,5),6)]:
        p=depth_dist(V,F(1,n)); m=max(p); co=[float(p.get(k,F(0))) for k in range(m+1)]
        if HAVE:
            hi=co[::-1]
            while len(hi)>1 and abs(hi[0])<1e-12: hi=hi[1:]
            r=np.roots(hi)
            angs=sorted(round(math.degrees(cmath.phase(x)),1) for x in r if abs(x)>1e-9)
            print(f"   AP-ish V={V} (n={n}): nonzero-root angles(deg)={angs}  (±120°? {any(abs(abs(a)-120)<2 for a in angs)})")
    print("   => depth-PGF roots NOT at ±120°. The LRC cyclotomic/±2π/3 lives elsewhere:")
    print("   WITNESS clock points = n-th roots of unity (THM-403): for 3|n, t=1/3 ↔ angle 2π/3 IS a")
    print("   worry-set lonely time; and the worry-set ROUND tournament's F(T,x) Lee-Yang zeros cluster")
    print("   at ±2π/3 (opus-S44). So the Eisenstein angle is the WITNESS/tournament face, not the PGF.\n")
    print("(B) z^7−z KEYSTONE (the W6 / 7-coloring / Φ3(2) unification, extending S687):")
    print(f"   7 = N(2+ζ6) = a²+ab+b² at (2,1) = {Neis(2,1)}  ⟹ (2+ζ6) is an Eisenstein PRIME of norm 7")
    print(f"   ⟹ ℤ[ζ6]/(2+ζ6) ≅ 𝔽₇  = the HEXAGONAL 7-COLORING palette (Isbell χ≤7) — color = residue mod the norm-7 prime")
    print(f"   z^7−z = z(z^6−1): over ℂ the roots = {{0}}∪{{6th roots}} = W6 (center+hexagon, χ=3 gadget),")
    roots6=[round(math.degrees(cmath.phase(cmath.exp(2j*math.pi*k/6))),0) for k in range(6)]
    print(f"          hexagon root angles(deg)={sorted(roots6)} — INCLUDES ±120°=±2π/3 (the Eisenstein angle);")
    print(f"   over 𝔽₇ the roots = all of 𝔽₇ (Fermat) = the 7 colors. ONE polynomial, two readings.")
    print(f"   ⟹ the HN UPPER BOUND 7 = the FORBIDDEN tournament H-value 7 = Φ3(2) = N(2+ζ6) = |PG(2,2)| — the SAME 7.\n")
    print("(C) FIELD / HEEGNER TOWER (extending S687's conjecture):")
    print("   χ=3: ℚ(√−3) (Eisenstein ℤ[ζ6], roots of unity, Mahler measure 1)  — d=3 Heegner")
    print("   χ=4: +ℚ(√−11) (Moser spindle: cosθ=5/6, e^{iθ} root of 3z²−5z+3, disc −11, Mahler measure 3) — d=11 Heegner")
    heeg=[1,2,3,7,11,19,43,67,163]
    print(f"   Heegner numbers (class number 1): {heeg}; χ=3↦d=3, χ=4↦d=11 (both Heegner).")
    print("   CONJECTURE (sharpening S687): each chromatic step adjoins a new class-number-1 (Heegner)")
    print("   imaginary-quadratic rotation field; χ(ℝ²)=5 would adjoin √−19 (next Heegner). The class")
    print("   number 1 = the rigidity (unique factorization of rotations) needed to FORCE a new color.")
if __name__=='__main__': main()

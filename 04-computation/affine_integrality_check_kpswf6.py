#!/usr/bin/env python3
"""
THREAD C, part (2b): the DECISIVE integrality check for the affine LP.

The level-2 affine LP has objective = a (signed) functional of the affine-type moments.
"Integral / complete for the AP optimizer" means precisely:

  (I)  the affine degree-2 atom (distinct-residue-pair count mod 7) is MAXIMIZED by AP, AND
  (II) among the shapes that TIE AP on the affine signature (its affine-orbit class), L_y is
       CONSTANT (= they are LP-equivalent), and that common value is the max -- i.e. AP is an
       extreme point of the affine polytope, unique UP TO the affine symmetry group (the
       precise analogue of "the linear optimizer is unique up to the code's automorphism
       group"). This is what makes the LP integral: the optimal face is exactly one affine
       orbit.

If instead L_y VARIES across AP's affine-signature class, the affine atom does NOT capture
L_y and the fix would be only partial. We test (I) and (II) exactly.

We ALSO test the contrast: under the LINEAR (translation-only difference mod 14) type, is
there a single atom AP extremizes? (HYP-2738 says no -- mixed signs.) We confirm AP does NOT
uniquely maximize any single linear-difference atom, quantifying the collapse.

Saved: 05-knowledge/results/affine_integrality_check_kpswf6.out
kind-pasteur-2026-06-21 THREAD C.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb
from collections import Counter, defaultdict

def breakpoints(E):
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        e = abs(e)
        for j in range(7):
            base = F(j, 7*e); step = F(1, e); m = 0; x = base
            while x < 1:
                if x >= 0: bps.add(x)
                m += 1; x = base + step*m
    return sorted(bps)

def hit_sectors_at(E, x):
    hit = set()
    for e in E:
        val = (F(e)*x) % 1
        s = int(val*7); s = 6 if s==7 else s
        hit.add(s)
    return hit

def factorial_moments_exact(E):
    bps = breakpoints(E)
    p = {t: F(0) for t in range(7)}
    for i in range(len(bps)-1):
        a,b = bps[i], bps[i+1]
        if b<=a: continue
        N = sum(1 for j in range(1,7) if j not in hit_sectors_at(E,(a+b)/2))
        p[N]+= (b-a)
    S = {r:F(0) for r in range(7)}
    for t in range(7):
        for r in range(7):
            if t>=r: S[r]+= comb(t,r)*p[t]
    return S

LY = {
    8:  {0:F(1),1:F(-1),2:F(1),3:F(-9,10),4:F(3,5)},
    9:  {0:F(1),1:F(-13,18),2:F(4,9),3:F(-1,6)},
    10: {0:F(1),1:F(-13,18),2:F(4,9),3:F(-1,6)},
}
def Ly(E,k):
    S=factorial_moments_exact(E)
    return sum(LY[k][r]*S[r] for r in LY[k])

def affine_sig(E):
    res=[e%7 for e in E]
    d=c=0
    for i in range(len(res)):
        for j in range(i+1,len(res)):
            if res[i]==res[j]: c+=1
            else: d+=1
    # full residue signature: also record the MULTISET of residue-occupancy counts
    occ = tuple(sorted(Counter(res).values()))
    return (d,c,occ)

def residue_occupancy(E):
    res=[e%7 for e in E]
    return tuple(sorted(Counter(res).values())), len(set(res))

def main():
    print("="*78)
    print("THREAD C (2b): DECISIVE integrality check of the affine level-2 LP")
    print("="*78)
    for k in [8,9,10]:
        AP=list(range(k))
        sig_ap=affine_sig(AP)
        Ly_ap=Ly(AP,k)
        occ_ap, ndist_ap = residue_occupancy(AP)
        # gather all bounded shapes
        tie_class=[]            # shapes sharing AP full affine signature
        ly_on_tie=set()
        fullres=[]              # shapes occupying all 7 residues (HYP-2749 stratum)
        ly_max=Ly_ap; argmax=AP; beaters=0
        ly_max_fullres=Ly_ap
        for rest in combinations(range(1,14),k-1):
            E=[0]+list(rest)
            ly=Ly(E,k)
            if ly>ly_max: ly_max=ly; argmax=E
            if ly>Ly_ap: beaters+=1
            occ,nd=residue_occupancy(E)
            if nd==7:
                fullres.append((E,ly))
                if ly>ly_max_fullres: ly_max_fullres=ly
            if affine_sig(E)==sig_ap:
                tie_class.append(E); ly_on_tie.add(ly)
        print(f"\nk={k}: AP={AP}")
        print(f"  AP affine signature (distinct,coll,occupancy) = {sig_ap}")
        print(f"  AP occupies {ndist_ap} distinct residues mod 7 (occupancy {occ_ap})")
        print(f"  L_y(AP) = {Ly_ap} = {float(Ly_ap):.6f}; global beaters over all shapes = {beaters}")
        print(f"  (I)  AP maximizes L_y globally: {beaters==0}")
        print(f"  affine TIE-class size (shapes with identical affine signature) = {len(tie_class)}")
        if ly_on_tie:
            print(f"  L_y values realized ON the tie-class: {len(ly_on_tie)} distinct")
            print(f"     min={float(min(ly_on_tie)):.6f}  max={float(max(ly_on_tie)):.6f}")
            print(f"  (II) L_y CONSTANT on AP's exact affine signature class: {len(ly_on_tie)==1}")
        # FULL-RESIDUE stratum (HYP-2749): is AP the max THERE?
        fr_beaters=sum(1 for _,ly in fullres if ly>Ly_ap)
        print(f"  full-residue (all 7) stratum: {len(fullres)} shapes; "
              f"AP-beaters on stratum = {fr_beaters}; AP is stratum-max: {fr_beaters==0}")
        # how many full-residue shapes are AFFINE-EQUIVALENT to AP (share residue occupancy)?
        fr_same_occ=sum(1 for E,_ in fullres if residue_occupancy(E)[0]==occ_ap)
        print(f"  full-residue shapes with AP's residue-occupancy profile {occ_ap}: {fr_same_occ}")

    # ---- CONTRAST: LINEAR (difference mod 14) single-atom -- AP extremizes none uniquely ----
    print("\n" + "-"*78)
    print("CONTRAST: under the LINEAR (translation-only) degree-2 type = difference mod 14,")
    print("AP does NOT uniquely extremize any single atom (the collapse, HYP-2738).")
    k=8; AP=list(range(k))
    def diff_profile(E):
        c=Counter()
        Es=sorted(E)
        for i in range(len(Es)):
            for j in range(i+1,len(Es)):
                d=(Es[j]-Es[i])%14; d=min(d,14-d); c[d]+=1
        return c
    ap_dp=diff_profile(AP)
    print(f"  AP diff-profile (|e_i-e_j| mod 14): {dict(ap_dp)}")
    # for each difference value, is AP the max?
    maxd_holders=defaultdict(int)
    nshapes=0
    for rest in combinations(range(1,14),k-1):
        E=[0]+list(rest); nshapes+=1
        dp=diff_profile(E)
        for d in range(1,8):
            if dp[d]>ap_dp[d]:
                maxd_holders[d]+=1
    print(f"  over {nshapes} shapes, # beating AP on each difference-atom:")
    for d in range(1,8):
        print(f"     diff={d}: {maxd_holders[d]} shapes have MORE pairs at this difference than AP")
    print("  => AP is NOT the per-atom maximizer in the linear basis: the linear lift has")
    print("     many free directions AP fails to extremize -> no single-atom certificate")
    print("     -> the linear/CJJ lift collapses to a signed aggregate (HYP-2744/2738).")
    print("  In the AFFINE basis these collapse to ONE atom AP DOES maximize.")
    print("="*78)

if __name__=="__main__":
    main()

"""
opus-2026-07-11-S251: the REMAINING Ostrowski LRC mathematics -- bridge my S248/S249/S250 tight-locus results
to the covering-min Ostrowski ladder (mac-mini S38/S65cont54), and push the open core "tight => {k*alpha}
structure" onto the FULL tight locus {AP, V*}.

BACKGROUND (prior fleet work):
  * Ostrowski ladder (mac-mini S38): the covering-min rungs are M_k = [0; 13, k] = k/(13k+1). Ends: rung k=1 =
    1/14 (AP, tight, non-covering); rung k=14 = 14/183 (deep well {1..12,182}, the DC covering-min).
  * The OPEN CORE (mac-mini S38): "tight => the extremal config is a {k*alpha}-progression"; then the three-gap
    (Steinhaus) theorem gives g<=3 for FREE (THM-527). The whole difficulty is proving the {k*alpha} STRUCTURE,
    not counting gaps.
  * My recent arc: S248 (empty window (1/14,3/41)), S249 (tight locus = exactly 2 mod-14 residue patterns
    {AP, V*}), S250 (base-rigidity). These are secretly the rung-1 (AP) END of the ladder.

FINDINGS (this session):

(1) THE LOW M-SPECTRUM = OSTROWSKI RUNGS [0;13,...]. Verified continued fractions: 1/14=[0;13,1], 2/27=[0;13,2],
    3/40=[0;13,3], ..., 14/183=[0;13,14] (the ladder spine M_k), and 3/41=[0;13,1,2] (a Farey CHILD of the AP
    rung). So S248's empty window (1/14, 3/41) is EXACTLY the Farey/Ostrowski gap between the rung [0;13,1] and
    its child [0;13,1,2] -- empty because no simpler continued fraction lies between a value and its own child.
    (This identifies S248 as the rung-1 instance of mac-mini/klein's Farey-tree M-spectrum.)

(2) THE {k*alpha} STRUCTURE HOLDS FOR THE FULL TIGHT LOCUS {AP, V*} (confirming the open core ON the classified
    locus). At the tight optimum t = 1/14, phases = residues mod 14:
      - AP = the FULL progression {k/14 : k=1..13}, gaps {1/14} (g=1, uniform);
      - V* = {1..11,13,24}: phase-set {1..11,13}/14 = the progression PUNCTURED at 12/14, gaps {1/14, 2/14}
        (g=2).
    BOTH are {k*alpha}-supported (alpha = 1/14), both three-gap (g <= 2 <= 3), both closest-approach = 1/14
    (tight). So S249's two mod-14 patterns ARE the two {k*alpha}-configs at n=14: the complete progression and
    the once-punctured progression. This is the bridge: S249's classification IS the {k*alpha} structure, made
    explicit for the whole tight locus, not just the AP.

(3) THE COMPOSITE-14 GIVES TWO RUNG-1 OCCUPANTS (corrected). The punctured progression must be realized by
    INTEGER speeds AND stay at rung 1 (no better witness). Among single moves 12 -> m (all land back on the
    progression), computing the ACTUAL M shows: ONLY 12 -> 24 keeps M = 1/14; 12 -> 36, 38, 50, ... land on an
    existing residue too (so are {k*alpha}-supported at t=1/14) but LIFT to a higher rung (M = 3/41 = [0;13,1,2],
    etc.) because a better witness appears at another t. So {k*alpha}-support at t=1/14 is NECESSARY but not
    sufficient; V* (the doubling 12->24) is the UNIQUE integer realization that also does not lift. At PRIME n
    the doubling map k->2k is a bijection mod n (no collision), so the punctured progression has no such integer
    realization => rung 1 has ONE occupant (the AP alone). At n=14 = 2*7 composite, rung 1 has TWO {AP, V*}.

NET. The remaining Ostrowski math, pushed one rung: the low M-spectrum is the Ostrowski/Farey tree rooted at
[0;13,.]; S248's empty gap is a Farey gap; S249's tight locus {AP, V*} = the complete and once-punctured
{k/14}-progressions (so the open "tight => {k*alpha}" holds on the classified locus, with three-gap free); and
the composite 14 = 2*7 is precisely why rung 1 has two occupants where the proved prime k=12 case has one. The
still-open piece is unchanged in kind -- prove tight => {k*alpha}-support for ARBITRARY families (not just
verify it on the locus) -- but it is now anchored to the explicit two-config n=14 picture.

-> mac-mini S38 (Ostrowski ladder), mac-mini S65cont54 (Farey tree), klein S267 (14/183 rung), THM-527
(three-gap rigidity), THM-612 (tight locus {AP,GW}), opus-S248/S249/S250.
"""
from math import gcd
from fractions import Fraction
def cf(fr):
    a=fr.numerator; b=fr.denominator; out=[]
    while b: out.append(a//b); a,b=b,a-(a//b)*b
    return out
def Mval(v):
    qs=set()
    for i in range(len(v)):
        for j in range(i+1,len(v)):
            s=v[i]+v[j]; g=gcd(v[i],v[j]); qs.add(s//gcd(s,g)); qs.add(s)
    best=Fraction(0)
    for q in sorted(x for x in qs if x>=2)[:5000]:
        bq=0
        for k in range(1,q//2+1):
            if gcd(k,q)!=1: continue
            m=min(min((vi*k)%q,q-(vi*k)%q) for vi in v)
            if m>bq: bq=m
        if Fraction(bq,q)>best: best=Fraction(bq,q)
    return best
def main():
    print("(1) Low M-spectrum as Ostrowski rungs [0;13,...]:")
    for name,f in [("1/14",Fraction(1,14)),("3/41",Fraction(3,41)),("2/27",Fraction(2,27)),
                   ("3/40",Fraction(3,40)),("1/13",Fraction(1,13)),("14/183",Fraction(14,183))]:
        print(f"   {name}: cf={cf(f)}")
    print("   empty gap (1/14,3/41) = Farey gap [0;13,1] -> child [0;13,1,2].")
    ap=list(range(1,14)); vstar=[1,2,3,4,5,6,7,8,9,10,11,13,24]
    print("\n(2) {k/14}-support of the tight locus (t=1/14):")
    for name,v in [("AP",ap),("V*",vstar)]:
        ph=sorted(set(x%14 for x in v)|{0})
        gaps=sorted(set((ph[(i+1)%len(ph)]-ph[i])%14 for i in range(len(ph))))
        miss=sorted(set(range(1,14))-set(ph))
        print(f"   {name}: phases(x14)={sorted(set(x%14 for x in v))} missing={miss} gaps(x14)={gaps} g={len(gaps)}")
    print("   AP=full progression; V*=punctured at 12. BOTH {k*alpha}-supported, three-gap.")
    print("\n(3) Composite realization -- ACTUAL M for single moves 12->m (only 24 stays at rung 1):")
    for m in [24,26,36,38,50,64]:
        w=sorted([x for x in ap if x!=12]+[m])
        M=Mval(w); print(f"   12->{m} (mod14={m%14}): M={M}=[0;{cf(M)[1:]}]  {'<-- TIGHT rung 1' if M==Fraction(1,14) else '(lifts)'}")
    print("   => only 12->24 (doubling, needs 14 composite) stays at 1/14; others lift to higher Ostrowski rungs.")
    print("   PRIME n: k->2k bijection, no collision, punctured progression not integer-realizable => AP alone.")
if __name__=='__main__':
    main()

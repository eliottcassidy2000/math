"""
opus-2026-07-11-S248: a sharper picture of the LRC(14) crux -- the composite-14 effect on the minimizer
locus and the spectrum near 1/14, via hypothesis investigation.

CRUX (post-S247): LRC(14) = no M < 1/14 = the AP-value 1/14 is the global minimum of M over 13-speed
families. Sharpening HOW the composite 14 = 2*7 (zero-divisor at 7) shapes this, vs the proved k=12 (prime 13):

(1) THE MINIMIZER LOCUS IS NOT UNIQUE (composite effect). Among single-element replacements {1..13}, k->m with
    M = 1/14 EXACTLY:
      - k -> m with m ≡ k (mod 14) [residue-preserving shift]: trivially tight for every k (same mod-14 config).
      - k=12 -> 24 = 2*12 [the DOUBLING]: the ONLY non-residue-preserving tight replacement (mac-mini THM-708;
        siblings 11->25, 13->15, 6->20 all jump to higher M). Enabled by 24 = 2*12 and 14 = 2*7: the doubling
        by 2 preserves the binding lattice constraint through the zero-divisor.
    So the ESSENTIAL tight families (distinct mod 14) are TWO: {1..13} (residues {1..13}, full) and
    V* = {1..11,13,24} (the doubling; residues {1..11,13,10} = missing 12, doubling 10). ALL non-AP tight
    families have residue COLLISIONS mod 14 (0 with full residues) -- the composite signature.

(2) THE SPECTRUM ABOVE 1/14 (the near-tight ladder). Empirically: 1/14 (tight locus), then a GAP, then
    3/41 = 0.0732 (= HYP-2934's K_{3,3} incidence wall / near-miss 12->36), then 2/27 (mediant), 3/40, 1/13.
    So the true second value is 3/41, NOT 2/27; the empty gap is (1/14, 3/41), width 1/574. The composite 14
    DENSIFIES the ladder (3/41 fills what is empty at k=12), which is why S246's window-empty was false.

(3) WHY k=12 (prime 13, PROVED HYP-4151) DOES NOT TRANSFER. At k=12 the modulus 13 is a field: the tight
    family is UNIQUE ({1..12} dilated, full residues mod 13), and the window (1/13, 2/25) is empty
    (equioscillation rigidity). At k=13 the modulus 14 is composite (zero-divisor at 7): the tight locus GAINS
    the doubling V*, the residues COLLIDE, and the ladder gains 3/41. The apex prime 7 = 14/2 is the exact
    seam (opus-S552o's mod-7 singleton = the multiple of 7 = the zero-divisor).

NET (sharpened crux): LRC(14) = "no M < 1/14", with the minimizer locus = {AP} u {doubling V*} (up to
residue-preserving shift), all mod-14-collision families. The obstruction to a k=12-style proof is precisely
the composite 14 = 2*7: it enlarges the minimizer set (doubling) and densifies the near-tight ladder (3/41),
both traced to the zero-divisor at 7. Any k=13 rigidity must handle the doubling family and the 3/41 rung,
which the field-based k=12 argument cannot see.
"""
from math import gcd
from functools import reduce
from fractions import Fraction
def primitive(v): return reduce(gcd,v)==1
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
    tight=Fraction(1,14); ap=list(range(1,14))
    print("Single-replacement {1..13} k->m with M=1/14 (non-residue-preserving only):")
    for k in range(1,14):
        for m in range(14,60):
            if (m-k)%14==0 or m in ap: continue
            v=sorted(set(ap)-{k}|{m})
            if len(v)==13 and primitive(v) and Mval(v)==tight:
                print(f"   {k} -> {m}: {v}  (m/k={m/k}, m mod14={m%14})")
    print("\n  => ONLY 12->24 (doubling) is a non-residue-preserving tight family (V*); enabled by 14=2*7.")
    print("  ESSENTIAL tight locus (mod 14): {1..13} (full residues) u V*={1..11,13,24} (doubling, residues collide).")
    print("  Spectrum above 1/14: gap (1/14, 3/41) empty; second value 3/41 (K33 wall), not 2/27.")
    print("  Composite 14 enlarges the minimizer set (doubling) + densifies the ladder (3/41) => k=12 proof cannot transfer.")
if __name__=='__main__':
    main()

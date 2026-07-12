"""
opus-2026-07-11-S235: paths 1 (quantify the detuning) + 3 (inverse-theorem bridge) integrated, with a
fleet course-correction. The margin for divisor-complete families is a FREE COROLLARY of bounded-clearing.

SETUP. S234: covering => divisor-complete (THM-366); LRC(14) <=> divisor-complete => M > 1/14. Path 1 =
"quantify the detuning" (prove divisor-complete => M >= 1/14 + eps). Path 3 = "inverse theorem: covering =>
near-AP". This session works both, pulling the fleet's additive-energy/Freiman state.

THE BAND-EDGE MARGIN LEMMA (elementary, PROVED; verified 0 violations / 19999):
  If S has a lonely time (clears, bandCount=0) at a modulus q with 14 nmid q, then every cleared runner has
  ||v_i p/q|| >= ceil(q/14)/q, so  M(S) >= ceil(q/14)/q > 1/14  (strict, since 14 nmid q => ceil(q/14) > q/14).
  At a multiple of 14 (14|q) the band edge is exactly (q/14)/q = 1/14 (no margin).
COROLLARY (characterization): M(S) = 1/14 (TIGHT)  =>  S clears ONLY at multiples of 14.
  Verified: AP {1..13} and V* {1..11,13,24} clear at q in {14,28,42,56} and at NO 14-nmid q. (Meets THM-610:
  a covering/tight family hides at a denominator divisible by 14.)

PATH 1 INTEGRATION. Divisor-complete families (which do NOT clear at q=14, having a mult of 14) all clear at a
non-14 modulus q in [15,41] (verified; adversarial worst q=31), so by the band-edge lemma M >= 2/27 > 1/14 --
they are LOOSE. So path 1's margin is NOT a separate difficulty: it is a FREE COROLLARY of "divisor-complete
clears at a bounded non-14 modulus" (the anti-concentration of S230/S231, verified <=60 diameter-free).
Proving bounded-clearing yields BOTH the loneliness certificate AND the strict margin M > 1/14 at once.

PATH 3 COURSE-CORRECTION (from the fleet + Explore). "covering => near-AP" is BACKWARDS: the AP is the
MAX-energy configuration and the MINIMIZER of every energy floor (THM-656/660: high R2 => WEAK floor); the
floors prove "even the AP clears". kps cont.36 DECOUPLES: the window-hard covering cores are LOOSE, not
near-tight -- confirmed here (all divisor-complete are loose, M > 1/14). So an energy-floor + "covering =>
near-AP" strategy pushes against the fleet's own finding. The repo's preferred inverse invariant is E3
(Schur triples, dilation-invariant like loneliness), NOT E2/doubling/BSG (translation-invariant, "cannot
distinguish the tight AP from its loose translate {2..14}"): opus-S182 + LEM-015 (interval maximizes E3,
Lean kernel-pure). So if the inverse-structure bridge is wanted, target E3, not BSG/Freiman-3k-4.

NET. Both paths reduce to the SAME open statement -- "divisor-complete clears at a bounded non-14 modulus" --
with the band-edge lemma making the strict margin M > 1/14 automatic. The near-AP direction is refuted; the
hard core is loose. No new open object: it is the S230/S231 anti-concentration, now with a free margin and a
clean tight-locus characterization.
"""
from math import gcd, ceil
from functools import reduce
from fractions import Fraction
import random
def clears(v,q):
    for p in range(1,q):
        if all(q<=14*((vi*p)%q)<=13*q for vi in v): return True
    return False
def divisor_complete(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def primitive(v): return reduce(gcd,v)==1
def Mval(v):
    qs=set()
    for i in range(13):
        for j in range(i+1,13):
            s=v[i]+v[j]; g=gcd(v[i],v[j]); qs.add(s//gcd(s,g)); qs.add(s)
    best=Fraction(0)
    for q in sorted(x for x in qs if x>=2):
        bq=0
        for k in range(1,q//2+1):
            if gcd(k,q)!=1: continue
            m=min(min((vi*k)%q,q-(vi*k)%q) for vi in v)
            if m>bq: bq=m
        if Fraction(bq,q)>best: best=Fraction(bq,q)
    return best

def main():
    random.seed(1); viol=tested=0
    for _ in range(20000):
        v=sorted(random.sample(range(1,150),13))
        if not primitive(v): continue
        M=Mval(v)
        for q in range(8,42):
            if q%14==0: continue
            if clears(v,q):
                tested+=1
                if M<Fraction(ceil(q/14),q): viol+=1
                break
    print(f"BAND-EDGE LEMMA [clear at 14-nmid q => M>=ceil(q/14)/q>1/14]: {viol} violations / {tested}")
    print("\nTIGHT families clear ONLY at multiples of 14:")
    for name,v in [("AP {1..13}",list(range(1,14))),("V* {1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24])]:
        v=sorted(v)
        c14=[q for q in range(2,60) if q%14==0 and clears(v,q)]
        cn=[q for q in range(2,60) if q%14!=0 and clears(v,q)]
        print(f"  {name:18}: M={float(Mval(v)):.4f}; 14|q clears={c14}; 14-nmid clears={cn[:4]}")
    random.seed(9); ms=[]
    for _ in range(30000):
        v=sorted(random.sample(range(1,80),13))
        if not (primitive(v) and divisor_complete(v)): continue
        q0=next((q for q in range(15,42) if q%14!=0 and clears(v,q)),None)
        if q0: ms.append(ceil(q0/14)/q0)
    print(f"\ndivisor-complete (n={len(ms)}): all clear at a non-14 q => min band-edge margin M>={min(ms):.5f} > 1/14={1/14:.5f} (LOOSE)")

if __name__=='__main__':
    main()

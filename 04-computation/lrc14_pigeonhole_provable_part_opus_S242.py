"""
opus-2026-07-11-S242: the LAST PROVABLE PART of the residual -- the pigeonhole-forced sub-class. A fully
rigorous theorem (three proved lemmas, no anti-concentration) clearing ~44% of divisor-complete families.

THEOREM (pigeonhole clearing, PROVED). Let v be a 13-speed family and q in {15,21,25,27} a composite with
prime factors <= 13 (danger band {0,+-1}; danger residues {1,q-1} coprime to q). If v has NO multiple of q
and FEWER THAN phi(q)/2 speeds coprime to q, then v clears at q, hence M(v) >= 2/q > 1/14.

PROOF (three proved lemmas):
  [auto-safe, S241] every speed with gcd(v_i,q)>1 is inert: at a unit multiplier p, v_i p mod q shares a
     factor with q, so it is not in the danger band {0,+-1} (unless q|v_i, excluded).
  [pigeonhole, this session] the < phi(q)/2 coprime-to-q speeds occupy < phi(q)/2 unit fold-classes, so some
     unit +-fold-class r is EMPTY. At the unit p = r^{-1}, no coprime speed hits {+-1} and the structured
     speeds are inert => bandCount(v,q,p) = 0 (a lonely multiplier).
  [band-edge, S235] clearing at q with 14 nmid q => M(v) >= ceil(q/14)/q = 2/q > 1/14. QED.

VERIFIED END-TO-END: pigeonhole-forced => actually clears, 0 violations / 1749 forced (of 4000 divisor-complete
spread families). Coverage: 43.7% of divisor-complete are pigeonhole-forced (no anti-concentration needed).
Forcing modulus distribution: q=27 (54%), q=25 (28%), q=21 (16%), q=15 (2%). Drivers: >=5 speeds div by 3
(forces q=27), >=4 speeds div by 5 (forces q=25).

TOTAL ELEMENTARY-PROVABLE COVERAGE (of ALL 13-speed families):
  - not-divisor-complete (t=1/d ladder, THM-366): ~91.5%  [misses a mult of some d<=14 => t=1/d witness]
  - divisor-complete & pigeonhole-forced (this theorem): ~44% of the 8.5% = ~3.7%
  => ~95.2% of all families are provably lonely (M > 1/14) by ELEMENTARY means (no anti-concentration).
  The remaining ~4.8% -- divisor-complete families whose coprime sub-family is large enough to potentially
  cover the fold-classes at every bounded composite -- is the genuine ANTI-CONCENTRATION CORE (fold-classes
  must CLUSTER, not just be few); this is where the hypothetical covering counterexamples would live, and it
  is the sole remaining open part.
"""
from math import gcd, ceil
from functools import reduce
import random
def bandCount(v,q,p): return sum(1 for vi in v if not (q<=14*((vi*p)%q)<=13*q))
def clears(v,q): return any(bandCount(v,q,p)==0 for p in range(1,q))
def th(q): return sum(1 for r in range(1,q) if gcd(r,q)==1)//2
def coprime_size(v,q): return sum(1 for x in v if gcd(x,q)==1)
def divisor_complete(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def primitive(v): return reduce(gcd,v)==1
def longest_AP(v):
    s=set(v); best=1
    for a in v:
        for d in range(1,max(v)//2+1):
            L=1;x=a+d
            while x in s:L+=1;x+=d
            if L>best:best=L
    return best
BAND1=[15,21,25,27]
def pigeon_forced_q(v,q): return (not any(x%q==0 for x in v)) and coprime_size(v,q) < th(q)

def main():
    random.seed(1); pool=[]; tries=0
    while len(pool)<4000 and tries<800000:
        tries+=1
        v=sorted(random.sample(range(1,150),13))
        if primitive(v) and divisor_complete(v) and longest_AP(v)<=7: pool.append(v)
    viol=forced=0
    for v in pool:
        for q in BAND1:
            if pigeon_forced_q(v,q):
                forced+=1
                if not clears(v,q): viol+=1
                break
    prov=sum(1 for v in pool if any(pigeon_forced_q(v,q) for q in BAND1))
    print(f"THEOREM check: pigeonhole-forced => clears: {viol} violations / {forced} forced (rigorous chain)")
    print(f"PROVABLE sub-class: {prov}/{len(pool)} = {100*prov/len(pool):.1f}% of divisor-complete (no anti-concentration)")
    from collections import Counter
    wq=Counter(next(q for q in BAND1 if pigeon_forced_q(v,q)) for v in pool if any(pigeon_forced_q(v,q) for q in BAND1))
    print(f"  forcing modulus: {dict(wq)}")
    print(f"TOTAL elementary-provable: ~91.5% (t=1/d ladder) + ~{0.085*prov/len(pool)*100:.1f}% (pigeonhole) "
          f"= ~{91.5+0.085*prov/len(pool)*100:.1f}% of all families; remainder = anti-concentration core")

if __name__=='__main__':
    main()

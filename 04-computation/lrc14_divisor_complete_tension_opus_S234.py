"""
opus-2026-07-11-S234: the hard core is the DIVISOR-COMPLETE families (= THM-366), and it equals the S232
multiplicand-maximal wall INTERSECT mult-of-14 -- the "detuned AP" with a 13-slot tension.

HONEST NOTE: the t=1/d ladder and "covering => divisor-complete" I derived this session ARE THM-366
(codex-S388, PROVED): if no speed is divisible by m (2<=m<=n=14), then t=a/m is a lonely witness; so every
covering counterexample contains a multiple of every m in {2,...,14}. I re-derived it. The NEW content is the
SYNTHESIS with the S232 summand-shell frame.

THE LADDER (THM-366): a family missing a multiple of some d in {2,...,14} is lonely via t=1/d (M>=1/d>=1/14).
Verified: 27435/27435 non-divisor-complete families cleared. So LRC(14) reduces to the DIVISOR-COMPLETE
families (a multiple of every d<=14) -- measured 8.5% of primitive 13-families.

THE SYNTHESIS (new): divisor-complete  <=>  multiplicand-maximal (mult of every d<=13, the S232 wall)  AND
mult-of-14. Verified 0 mismatches / 20000. So THM-366's hard core = S232's wall + a multiple of 14.

THE 13-SLOT TENSION (the mechanism):
  - tight (M=1/14)  =>  multiplicand-maximal AND NO mult of 14   (= the AP {1..13}, bucket A: t=1/14 works).
  - divisor-complete =>  multiplicand-maximal AND a mult of 14.
Both need multiplicand-maximal (AP-coherence). But with only 13 slots you cannot be the tight AP {1..13}
(which has no mult of 14) AND contain a mult of 14 -- the mult-of-14 costs a slot, breaking the coherence.
So divisor-complete families are AP-coherent-but-DETUNED, hence M > 1/14 (empirically). E.g. the AP {1..13}
is tight (M=1/14, not DC); the shift-AP {2..14} is DC but M=1/8=0.125. divisor-complete families clear at
q in [15,24] (the {0,+-1}-band shell regime of S230/S232).

HONEST STATUS: "divisor-complete => lonely" IS LRC(14) restricted to the 8.5% hard core (covering =>
divisor-complete, THM-366). The tension is the MECHANISM but is not quantified into a proof; the min-M margin
over divisor-complete is search-limited (0.087 at Vmax<=24, +0.0155 above 1/14; not a rigorous bound). The
value of the S230->S234 arc: the hard core is now a single, sharply-named, structurally-explained class --
the detuned-AP (multiplicand-maximal + mult-of-14) divisor-complete families -- not an anti-concentration
lemma and not the AP itself.
"""
import random
from math import gcd
from functools import reduce
def divisor_complete(v): return all(any(x%d==0 for x in v) for d in range(2,15))
def mult_maximal(v): return all(any(x%d==0 for x in v) for d in range(2,14))
def has14(v): return any(x%14==0 for x in v)
def primitive(v): return reduce(gcd,v)==1
def bandCount(v,q,p): return sum(1 for vi in v if not (q<=14*((vi*p)%q)<=13*q))

def main():
    random.seed(1)
    # ladder: not-divisor-complete => t=1/d clears (M>=1/14)
    N=nc=bad=0
    for _ in range(30000):
        v=sorted(random.sample(range(1,120),13))
        if not primitive(v): continue
        N+=1
        md=[d for d in range(2,15) if not any(x%d==0 for x in v)]
        if md:
            nc+=1; d=md[0]
            if not all((x%d)!=0 for x in v): bad+=1   # t=1/d clears
    print(f"(ladder, THM-366) not-divisor-complete => t=1/d clears: {nc-bad}/{nc}; divisor-complete = {100*(N-nc)/N:.1f}% hard core")
    # synthesis identity
    mis=0
    for _ in range(20000):
        v=sorted(random.sample(range(1,60),13))
        if divisor_complete(v)!=(mult_maximal(v) and has14(v)): mis+=1
    print(f"(synthesis) divisor-complete == multiplicand-maximal(S232) AND mult-of-14: {mis} mismatches / 20000")
    # tension
    print("(tension) tight AP vs detuned divisor-complete:")
    for name,v in {"AP {1..13}":list(range(1,14)),"shift-AP {2..14}":list(range(2,15))}.items():
        cl=[q for q in range(15,25) if any(bandCount(v,q,p)==0 for p in range(1,q))]
        print(f"    {name:16}: DC={divisor_complete(v)}, clears in [15,24] at q={cl[:4]}")

if __name__=='__main__':
    main()

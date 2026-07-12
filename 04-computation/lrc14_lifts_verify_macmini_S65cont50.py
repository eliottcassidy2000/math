#!/usr/bin/env python3
"""cont.50: HONESTY CHECK on cont.49's '<=6 distinct lifts'. For GENERIC (random) primitive
large-diameter DC 13-families, compute (a) min distinct-lifts over candidate scales L, and
(b) the odd-runner count. If either exceeds 6, the '<=6' is FALSE for generic DC (cont.49's
count was construction-specific) => the large-diameter route is multi-scale (THM-688), not few-lifts."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random
def is_covering(v):
    return all(any(x%q==0 for x in v) for q in range(2,15))
def make_random_dc(seed, Vmax):
    random.seed(seed)
    for _ in range(2000):
        v=sorted(set(random.sample(range(1,Vmax),13)))
        if len(v)!=13: continue
        if reduce(gcd,v)!=1: continue
        if is_covering(v): return v
    return None
def min_distinct_lifts(v):
    # over candidate scales L, distinct lifts k_i = round(v_i/L); min over L in [Vmax/12 .. Vmax]
    best=99; bestL=None
    Vmax=max(v)
    for L in range(max(2,Vmax//13), Vmax+1):
        k=set(round(x/L) for x in v)
        if len(k)<best: best=len(k); bestL=L
    return best,bestL
def odd_count(v): return sum(1 for x in v if x%2==1)
print("HONESTY CHECK: is '<=6 distinct lifts / <=6 odd' true for GENERIC large-diameter DC?")
print(f"{'seed':>5s} {'Vmax':>6s} {'diam':>6s} {'min-lifts':>10s} {'#odd':>6s}")
random.seed(50)
maxlifts=0; maxodd=0; n=0
for Vmax in [100,300,800,2000]:
    for s in range(6):
        v=make_random_dc(s*17+Vmax, Vmax)
        if v is None: continue
        ml,L=min_distinct_lifts(v); oc=odd_count(v); n+=1
        maxlifts=max(maxlifts,ml); maxodd=max(maxodd,oc)
        flag = "" if ml<=6 and oc<=6 else " *** EXCEEDS 6"
        print(f"  {s:>5d} {Vmax:>6d} {max(v)-min(v):>6d} {ml:>10d} {oc:>6d}{flag}")
# adversarial: DC with FEW special runners + many free ODD (breaks <=6 odd)
print("\nADVERSARIAL: DC = 1 mult-of-840 (all even conds) + 1 odd-mult-of-45045 + 11 free odd:")
adv=sorted(set([840,45045]+[3*i+1 for i in [1,3,5,7,9,11,13,15,17,19,21] if (3*i+1)%2==1][:11]))
adv=[840,45045]+[x for x in [7,11,13,17,19,23,25,29,31,37,41] ]
adv=sorted(set(adv))[:13]
while len(adv)<13: adv.append(max(adv)+2)
print(f"  family {adv[:6]}...: covering={is_covering(adv)}, #odd={odd_count(adv)}, min-lifts={min_distinct_lifts(adv)[0]}")
print(f"\nGENERIC DC: max min-lifts = {maxlifts}, max #odd = {maxodd} over {n} families")
print(f"VERDICT: '<=6' {'HOLDS for generic (empirical)' if maxlifts<=6 and maxodd<=6 else 'FAILS -- not a theorem, route is multi-scale THM-688'}")

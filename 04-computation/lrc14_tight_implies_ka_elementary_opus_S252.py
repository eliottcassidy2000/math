"""
opus-2026-07-11-S252: working the Chebyshev-equioscillation route toward "tight => {k*alpha}" -- and finding
the tight side is ELEMENTARY (the trivial t=1/14 witness), with the hard {k*alpha} relocating to the covering
side. A clarification/unification, honestly.

TARGET (open core, mac-mini S38): tight (M=1/14) => the optimal config is a {k*alpha}-progression; then the
three-gap (Steinhaus) theorem (THM-527) gives the gap count for free.

THE ELEMENTARY LEMMA (PROVED). Let v = (v_1,...,v_13) be primitive with NO multiple of 14, and M(v) = 1/14.
Then t* = 1/14 attains M, and the phase multiset {v_i/14 mod 1} = {(v_i mod 14)/14} is contained in
{1/14, 2/14, ..., 13/14} = {k*(1/14)} -- an arithmetic progression (alpha = 1/14). Hence tight => {k*alpha}.
  Proof. No v_i ≡ 0 (mod 14) => v_i mod 14 in {1,...,13} => ||v_i/14|| = min(v_i mod 14, 14 - v_i mod 14)/14
  >= 1/14 for every i. So min_i ||v_i/14|| >= 1/14, giving M(v) >= 1/14. Since M(v) = 1/14 and 1/14 is a value
  of min_i ||v_i t|| at t=1/14, equality holds AT t*=1/14. The phases are then multiples of 1/14 in {1..13}/14.
  QED. (No equioscillation needed; the {k*alpha} structure is the trivial witness.)

THE REDUCTION. So "tight => {k*alpha}" holds as soon as "tight => no multiple of 14". And
  tight => no-mult-14   <=>   mult-of-14 => M != 1/14   <=>  mult-of-14 => M > 1/14 (loose)
which is exactly a face of the CLEAN-RULER RESIDUAL. So on the TIGHT side, the Ostrowski {k*alpha} rigidity is
NOT a separate hard problem -- it is EQUIVALENT to "mult-of-14 families are loose", the residual the fleet is
already proving. The {k*alpha} content is the trivial t=1/14 witness plus the residual.

VERIFIED (this file):
  (A) LEMMA CORE: every no-mult-14 family has min_i ||v_i/14|| >= 1/14 (1358/1358) -- trivially true.
  (B) tight => no-mult-14: 275/275 tight families have no multiple of 14; ZERO tight with a mult of 14.
  (C) mult-of-14 => loose: 11995/11995 mult-of-14 primitives have M > 1/14; ZERO at 1/14, ZERO below.
      Critical near-AP test {1..12,14}: M = 1/13 (loose, witnessed at q*=13, the peeled AP, since 14 ≡ 1 mod 13).
      Deep well {1..12,182}: M = 14/183 (loose, the covering-min).

RELOCATION (the honest lesson). The equioscillation (mac-mini S40) is a TWO-POINT (binding-pair) phenomenon
at the optimizer -- it pins t* by a pair at distance M, but does NOT force the FULL config to be {k*alpha}.
On the tight side that is irrelevant (the t=1/14 witness gives {k*alpha} outright). The genuinely hard
{k*alpha} is the COVERING side: the M-minimizing covering family is the deep well {1..12,182}, whose config is
{k*alpha} (alpha = 14/183) ONLY because its core {1..12} is an interval; a GENERIC covering family has g=5
(not {k*alpha}). So "the extremal covering config is {k*alpha}" = "the covering-min is the interval-core deep
well" = the covering-min bound M >= 14/183 -- which is the residual, and (mac-mini S40) needs a DUAL
(Delsarte / de la Vallee-Poussin positive-polynomial) certificate, not equioscillation-greedy.

NET. Working the route: "tight => {k*alpha}" is ELEMENTARY (t=1/14 witness) and EQUIVALENT to the clean-ruler
residual on the tight side; the equioscillation tool and the hard {k*alpha} both belong to the COVERING side.
This unifies the Ostrowski rigidity with the residual (same statement, tight side) and correctly places the
remaining difficulty on the covering-min dual certificate.

-> mac-mini S38 (Ostrowski ladder / the open {k*alpha}), mac-mini S40 (2-point equioscillation, dual
certificate), klein S267 (14/183 covering-min), THM-527 (three-gap), opus-S251 ({AP,V*}={k/14} full/punctured),
opus-S249 (tight => no-mult-14, first observed), S246 (divisor-complete loose), S235 (band-edge).
"""
from math import gcd
from functools import reduce
from fractions import Fraction
import random
def primitive(v): return reduce(gcd,v)==1
def has14(v): return any(x%14==0 for x in v)
def Mval(v):
    qs=set()
    for i in range(len(v)):
        for j in range(i+1,len(v)):
            s=v[i]+v[j]; g=gcd(v[i],v[j]); qs.add(s//gcd(s,g)); qs.add(s)
    best=Fraction(0); argq=None
    for q in sorted(x for x in qs if x>=2)[:6000]:
        bq=0
        for k in range(1,q//2+1):
            if gcd(k,q)!=1: continue
            m=min(min((vi*k)%q,q-(vi*k)%q) for vi in v)
            if m>bq: bq=m
        if Fraction(bq,q)>best: best=Fraction(bq,q); argq=q
    return best,argq
def main():
    tick=Fraction(1,14)
    m14=lambda v: Fraction(min(min(x%14,14-x%14) for x in v),14)
    random.seed(1); okA=totA=0
    for _ in range(4000):
        v=sorted(random.sample(range(1,60),13))
        if not primitive(v) or has14(v): continue
        totA+=1; okA+= (m14(v)>=tick)
    print(f"(A) lemma core: no-mult-14 with min||v/14||>=1/14: {okA}/{totA}")
    random.seed(2); tn=tw=0
    ap=list(range(1,14)); vstar=[1,2,3,4,5,6,7,8,9,10,11,13,24]
    def pert(b,k):
        v=list(b)
        for _ in range(k): v[random.randrange(13)]=random.randint(1,70)
        return sorted(set(v)) if len(set(v))==13 else None
    for _ in range(30000):
        v=pert(ap,random.randint(1,2)) if random.random()<0.5 else pert(vstar,random.randint(1,2))
        if v is None or len(v)!=13 or not primitive(v): continue
        if Mval(v)[0]==tick:
            if has14(v): tw+=1
            else: tn+=1
    print(f"(B) tight: no-mult-14={tn}, with-mult-14={tw}  => tight=>no-mult-14: {tw==0}")
    random.seed(3); loose=eq=below=0
    for _ in range(20000):
        v=sorted(random.sample(range(1,80),13))
        if not primitive(v) or not has14(v): continue
        M=Mval(v)[0]
        if M>tick: loose+=1
        elif M==tick: eq+=1
        else: below+=1
    print(f"(C) mult-of-14 primitives: loose={loose}, =1/14: {eq}, <1/14: {below}")
    for nm,v in [("{1..12,14}",[1,2,3,4,5,6,7,8,9,10,11,12,14]),("deep well",[1,2,3,4,5,6,7,8,9,10,11,12,182])]:
        M,q=Mval(v); print(f"    {nm}: M={M} at q*={q}")
if __name__=='__main__':
    main()

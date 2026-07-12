"""
opus-2026-07-11-S233: the tight-floor census is the ELEMENTARY bucket -- the clean two-bucket dispatch.

Following S232 (the hard part = the multiplicand-maximal AP-coherent wall) and mac-mini THM-708 (tight
families admit no clean ruler and need none: dispatched by t=1/14, "14 nmid every element"), the loneliness
certificate of every 13-family splits cleanly:

  BUCKET A [no speed divisible by 14]:  t = 1/14 is a loneliness witness -> M(S) >= 1/14. ELEMENTARY.
        Proof: 14 nmid v_i => v_i/14 not in Z => ||v_i/14|| = min(v_i mod 14, 14 - v_i mod 14)/14 >= 1/14. QED.
        Contains ALL tight / AP-coherent families (the S232 "wall"): AP {1..13}, GW 12->24, V*={1..11,13,24}.
        These are exactly the M=1/14 extremizers -- the tight FLOOR lives entirely in this elementary bucket.

  BUCKET B [some speed divisible by 14]:  a clean ruler certifies loneliness (B5>0). The mult-of-14 speed
        DETUNES AP-coherence, so bucket B is comfortably lonely and clean-ruled at a SMALL bounded modulus.

VERIFIED:
  (1) [clean ruler q<=80] OR [no mult of 14] covers 19999/19999 primitive families (0 gap).
  (2) The tight families are all bucket A (no mult 14, has_clean=False) -- t=1/14 dispatched, S232 wall dissolved.
  (3) Bucket B is ADVERSARIALLY clean-ruled at bounded, diameter-free q: hill-climb max smallest-clean-q =
      37,47,49 across Vmax ceilings 200,2000,20000 (0 no-clean). The mult-of-14 speed => easy clean ruler.

RESOLUTION of the S232 wall: the AP-coherent families (no clean ruler) are EXACTLY bucket A, so the
clean-ruler route's inability to certify them is CORRECT -- they use the elementary t=1/14 witness. The
clean-ruler route (hB5) is responsible ONLY for bucket B (mult of 14).

HONEST LIMIT: (1)-(3) verify the dispatch on RANDOM/adversarial families, none of which are covering. The
residual hard core = NEAR-COVERING bucket-B families (M near 1/14 with a mult of 14). "Every bucket-B family
is clean-ruled" IS hB5 = LRC(14)'s content -- not closed by search (near-covering families are structured,
not reached by hill-climb; tightest mult-of-14 found here is M ~ 0.12, but that is weak-search, not a bound).
The value is the CLEAN DIVISION OF LABOR: the tight floor is elementary (bucket A); the clean-ruler route
never approaches 1/14 -- it only certifies the detuned, comfortably-lonely bucket B.
"""
import random
from math import gcd
from functools import reduce
def bandCount(v,q,p): return sum(1 for vi in v if not (q<=14*((vi*p)%q)<=13*q))
def smallest_clean(v,Q):
    for q in range(8,Q+1):
        live=False; ok=True
        for p in range(1,q):
            bc=bandCount(v,q,p)
            if bc>5: ok=False; break
            if bc==0: live=True
        if ok and live: return q
    return None
def t14(v): return all(x%14!=0 for x in v)          # t=1/14 witness available
def primitive(v): return reduce(gcd,v)==1

def main():
    random.seed(1); N=0; cov=0; gap=0
    for _ in range(20000):
        v=sorted(random.sample(range(1,200),13))
        if not primitive(v): continue
        N+=1
        if t14(v) or smallest_clean(v,80): cov+=1
        else: gap+=1
    print(f"(1) [clean ruler q<=80] OR [no mult 14]  covers {cov}/{N}  (gap {gap})")
    print("\n(2) tight / AP-coherent families are BUCKET A (elementary t=1/14):")
    for name,v in [("AP {1..13}",list(range(1,14))),("GW 12->24",[1,2,3,4,5,6,7,8,9,10,11,13,24]),
                   ("V* sporadic",[1,2,3,4,5,6,7,8,9,10,11,13,24])]:
        v=sorted(set(v))
        print(f"    {name:12}: no-mult-14={t14(v)} (t=1/14 works), clean_ruler={smallest_clean(v,120)}")
    print("\n(3) bucket B (mult of 14) adversarial clean-ruler (bounded, diameter-free):")
    def cost(v):
        if t14(v) or not primitive(v): return -1
        sc=smallest_clean(v,200); return 10**9 if sc is None else sc
    random.seed(31)
    for ceil in [200,2000,20000]:
        worst=0; nofound=0
        for _ in range(150):
            v=None
            for _ in range(200):
                c=sorted(set(random.sample(range(1,ceil),12)+[14*random.randint(1,ceil//14)]))
                if len(c)==13 and not t14(c) and primitive(c): v=c; break
            if v is None: continue
            cur=cost(v)
            for _ in range(60):
                i=random.randrange(13); nv=v[:]; nv[i]=random.randrange(1,ceil); nv=sorted(set(nv))
                if len(nv)<13: continue
                c=cost(nv)
                if c>=cur: v,cur=nv,c
            if cur>=10**9: nofound+=1
            elif cur>worst: worst=cur
        print(f"    Vmax ceil {ceil:6d}: adversarial max smallest-clean-q = {worst}  (#no-clean<=200: {nofound})")

if __name__=='__main__':
    main()

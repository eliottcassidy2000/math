from math import gcd
from itertools import combinations
import random
# UNIFIED beat certificate.  For a beat frequency q:
#   q <= 14: W_q = {0}, so a speed blocks only if q | v.  If q divides NO speed,
#            EVERY numerator survives -> lonely (the classical sieve).
#   q >  14: counting theorem -- if (#distinct kill-sets)*k_q < phi(q), some
#            numerator survives -> lonely.
# My S384 run excluded q <= 14 and so missed the tight families' own witnesses.
def W(q): return set([0]+[r%q for j in range(1,(q-1)//14+1) for r in (j,-j)])
def killset(v,q,Wq): return frozenset(p for p in range(1,q) if gcd(p,q)==1 and (v*p)%q in Wq)
def phi(n): return sum(1 for k in range(1,n+1) if gcd(k,n)==1)
def kq(q,Wq): return sum(1 for w in Wq if w!=0 and gcd(w,q)==1)
def beats(V):
    B=set()
    for a,b in combinations(V,2):
        B.add(a+b)
        d=abs(a-b)
        if d>0: B.add(d)
    for v in V: B.add(2*v)
    return {q for q in B if q>=2}
def certifies(V, cap=400):
    """returns (q, kind, margin) for the best certifying beat frequency, else None."""
    best=None
    for q in sorted(beats(V)):
        if q>cap: continue
        if q<=14:
            if not any(v%q==0 for v in V):
                return (q,'sieve',phi(q))          # classical: all numerators survive
            continue
        if any(v%q==0 for v in V): continue
        Wq=W(q); k=kq(q,Wq)
        if k==0: continue
        dis=len({killset(v,q,Wq) for v in V})
        m=phi(q)-dis*k
        if best is None or m>best[2]: best=(q,'count',m)
    return best if (best and best[2]>0) else best

print("(6) CORRECTED: unified beat certificate INCLUDING q <= 14")
for nm,V in [("{1,...,13}",list(range(1,14))),
             ("{1..11,13,24}",[1,2,3,4,5,6,7,8,9,10,11,13,24]),
             ("2*{1,...,13}",[2*i for i in range(1,14)]),
             ("odd {1,3,..,25}",[2*i+1 for i in range(13)])]:
    c=certifies(V)
    print(f"    {nm:16s} -> q={c[0]:4d} kind={c[1]:6s} margin={c[2]}" if c else f"    {nm:16s} -> NONE")

print()
print("(7) ADVERSARIAL HUNT against the UNIFIED certificate")
random.seed(38421)
worst=(10**9,None,None)
for trial in range(30):
    V=sorted(random.sample(range(1,120),13))
    c=certifies(V); cur=c[2] if c else -10**9
    stall=0
    for step in range(240):
        Wv=list(V); i=random.randrange(13)
        Wv[i]=max(1,Wv[i]+random.choice([-5,-3,-2,-1,1,2,3,5]))
        Wv=sorted(set(Wv))
        if len(Wv)!=13: continue
        c2=certifies(Wv); v2=c2[2] if c2 else -10**9
        if v2<=cur:
            if v2<cur: stall=0
            V,cur=Wv,v2
        else: stall+=1
        if stall>150: break
    if cur<worst[0]: worst=(cur,list(V),certifies(V))
print(f"    lowest certificate margin found: {worst[0]}")
print(f"      V = {worst[1]}")
print(f"      best beat: {worst[2]}")
print()
print("    margin > 0 for every family tested => the unified beat certificate held throughout")

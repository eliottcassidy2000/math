from math import gcd
from functools import reduce
from itertools import combinations
import random
def W(q): return set([0]+[r%q for j in range(1,(q-1)//14+1) for r in (j,-j)])
def killset(v,q,Wq): return frozenset(p for p in range(1,q) if gcd(p,q)==1 and (v*p)%q in Wq)
def phi(n): return sum(1 for k in range(1,n+1) if gcd(k,n)==1)
def kq(q,Wq): return sum(1 for w in Wq if w!=0 and gcd(w,q)==1)
def pmclasses(V,q): return len({min(v%q,(-v)%q) for v in V})

print("(1c) CORRECTED: is distinct kill-sets <= |{+-v mod q}| (the safe direction)?")
random.seed(3872)
viol=0; n=0; strict=0
for _ in range(150):
    q=random.randint(15,140); Wq=W(q)
    V=sorted(random.sample(range(1,300),13))
    a=len({killset(v,q,Wq) for v in V}); b=pmclasses(V,q)
    n+=1
    if a>b: viol+=1
    if a<b: strict+=1
print(f"    {n} pairs: violations of distinct <= pmclasses = {viol}; strictly fewer = {strict}")
print("    => the pairing lemma direction holds; the converse fails, so pmclasses is a")
print("       CONSERVATIVE proxy (real blocking capacity is often lower still).")

def beats(V):
    B=set()
    for a,b in combinations(V,2):
        B.add(a+b); d=abs(a-b)
        if d>0: B.add(d)
    for v in V: B.add(2*v)
    return {q for q in B if q>=2}
def certifies(V,cap=500,conservative=True):
    best=None
    for q in sorted(beats(V)):
        if q>cap: continue
        if q<=14:
            if not any(v%q==0 for v in V): return (q,'sieve',phi(q))
            continue
        if any(v%q==0 for v in V): continue
        Wq=W(q); k=kq(q,Wq)
        if k==0: continue
        dis = pmclasses(V,q) if conservative else len({killset(v,q,Wq) for v in V})
        m=phi(q)-dis*k
        if best is None or m>best[2]: best=(q,'count',m)
    return best

print()
print("(3c) TARGETED ATTACK, corrected: speeds in a fixed class r mod m")
print("     (differences become multiples of m -- highly composite, low phi/k)")
random.seed(3873)
worst=(10**9,None,None)
tested=0
for m in [210,180,120,90,60,42,30]:
    for r in range(1,m):
        if gcd(r,m)!=1 and r!=m//2: continue
        V=sorted({r+m*i for i in range(0,13)})
        if len(V)!=13: continue
        if reduce(gcd,V)!=1: continue          # dilation: only primitive families count
        tested+=1
        c=certifies(V); mg=c[2] if c else -10**9
        if mg<worst[0]: worst=(mg,list(V),c)
print(f"     {tested} primitive fixed-class families tested")
print(f"     lowest certificate margin: {worst[0]}")
print(f"       V = {worst[1][:6]}... (first 6)")
print(f"       best beat = {worst[2]}")
print()
print("(4c) WHY THIS CONSTRUCTION CANNOT BE PUSHED FURTHER")
print("     To make DIFFERENCES highly composite, put all speeds in one class mod m.")
print("     Then sums are 2r + multiples of m.  To make SUMS composite too you need")
print("     2r = 0 mod m, i.e. r = 0 or m/2 -- and both force gcd(V) > 1, i.e. a")
print("     NON-PRIMITIVE family, which dilation invariance (THM-1050) reduces away.")
for m in [210,120,60]:
    for r in [0,m//2]:
        V=sorted({r+m*i for i in range(1,14)})
        if len(V)==13:
            print(f"       m={m:4d} r={r:4d}: gcd(V) = {reduce(gcd,V):4d}  "
                  f"{'PRIMITIVE' if reduce(gcd,V)==1 else 'non-primitive -> dilation reduces it'}")

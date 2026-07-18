from math import gcd
from functools import reduce
import random
def lonely_at(V,p,q):
    for v in V:
        r=(v*p)%q
        if min(r,q-r)*14 < q: return False
    return True
def min_den(V,Qmax=3000):
    for q in range(1,Qmax+1):
        for p in range(1,q):
            if gcd(p,q)!=1: continue
            if lonely_at(V,p,q): return q
    return None
def prim(V): return reduce(gcd,V)==1
print("(7) HARDER ADVERSARIAL HUNT -- is 32 a ceiling or just where I stopped?")
random.seed(37311)
best=(0,None); hist=[]
for trial in range(45):
    V=sorted(random.sample(range(1,600),13))
    if not prim(V): continue
    cur=min_den(V) or 9999
    stall=0
    for step in range(700):
        W=list(V); 
        for _ in range(random.choice([1,1,2])):
            i=random.randrange(13); W[i]=max(1,W[i]+random.choice([-8,-5,-3,-2,-1,1,2,3,5,8]))
        W=sorted(set(W))
        if len(W)!=13 or not prim(W): continue
        d=min_den(W) or 9999
        if d>=cur:
            if d>cur: stall=0
            V,cur=W,d
        else: stall+=1
        if stall>300: break
    hist.append(cur)
    if cur>best[0]: best=(cur,list(V))
hist.sort()
print(f"    best over 45 hill-climbs: min denominator = {best[0]}")
print(f"      V = {best[1]}")
print(f"    distribution of hill-climb optima: {hist}")
print()
print("(8) SANITY: verify the witness family really has no small lonely rational")
V=best[1]
if V:
    print(f"    checking all q <= {best[0]} exhaustively for V ...")
    found=[q for q in range(1,best[0]+1) for p in range(1,q) if gcd(p,q)==1 and lonely_at(V,p,q)]
    print(f"    denominators admitting a lonely point below the reported minimum: {found}")
    print(f"    gcd(V) = {reduce(gcd,V)}   (must be 1 for primitivity)")

from math import gcd
import random
def danger(q):
    lo=-(-q//14); hi=(13*q)//14
    return set(range(0,lo))|set(range(hi+1,q))
def clears_at(v,q):
    dg=danger(q); return any(all((x*p)%q not in dg for x in v) for p in range(1,q))
QS=[q for q in range(15,32) if q%14]  # window, 14-nmid
def clears_window(v): return any(clears_at(v,q) for q in QS)
def longest_AP(v):
    vs=sorted(set(v)); best=1
    for i in range(len(vs)):
        for d in range(1,vs[-1]):
            L=1; x=vs[i]
            S=set(vs)
            while x+d in S: L+=1; x+=d
            best=max(best,L)
    return best
random.seed(3)
# search for primitive 13-sets that FAIL to clear at every window q (the "tight" window-covers)
fails=[]
n=0
while n<40000:
    v=sorted(random.sample(range(1,60),13))
    g=0
    for x in v: g=gcd(g,x)
    if g!=1: continue
    n+=1
    if not clears_window(v): fails.append(tuple(v))
print(f"searched {n} primitive 13-sets in [1..59]: {len(fails)} fail to clear at EVERY window q")
for v in fails[:8]:
    print(f"  window-cover: {v}  longest-AP={longest_AP(list(v))}")
# check the AP {1..13} and near-APs
for name,v in [('AP {1..13}',tuple(range(1,14))), ('{2..14}',tuple(range(2,15))),
               ('AP+shift {1..12,14}',(1,2,3,4,5,6,7,8,9,10,11,12,14))]:
    print(f"  {name}: clears_window = {clears_window(list(v))}  longest-AP={longest_AP(list(v))}")

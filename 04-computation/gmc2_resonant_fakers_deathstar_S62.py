from fractions import Fraction as Fr
from math import factorial
from itertools import product
def pmul(p,q):
    r={}
    for (a,b),u in p.items():
        for (c,d),v in q.items():
            k=(a+c,b+d); r[k]=r.get(k,Fr(0))+u*v
    return r
def Epow(P,m):
    Pm={(0,0):Fr(1)}
    for _ in range(m): Pm=pmul(Pm,P)
    return sum((c*factorial(a) for (a,b),c in Pm.items() if a==b), Fr(0))
def charges(P): return sorted({a-b for (a,b),c in P.items() if c!=0})
def two_sided(P):
    ch=charges(P); return any(c<0 for c in ch) and any(c>0 for c in ch)
half=[Fr(x,2) for x in range(-4,5)]
mons=[(2,0),(1,0),(0,1),(0,2)]; names=["Z^2","Z","W","W^2"]
print("TWO-SIDED fakers on {+2,+1,-1,-2} (E1=E2=E3=0), with kill-m and E4:")
found=[]
for co in product(half,repeat=4):
    P={m:c for m,c in zip(mons,co) if c!=0}
    if not P or not two_sided(P): continue
    if any(Epow(P,m)!=0 for m in (1,2,3)): continue
    kill=next((m for m in range(1,10) if Epow(P,m)!=0),None)
    found.append((co,kill,Epow(P,4)))
for co,kill,E4 in found:
    s="+".join(f"{c}*{n}" for c,n in zip(co,names) if c!=0)
    print(f"  {s:38s} kill-m={kill}  E4={E4}")
print(f"total two-sided fakers: {len(found)}; all killed by m=4: {all(k==4 for _,k,_ in found)}")
# structural: are all these the (a,d)=(±1/2,∓1/2) parity family boxeph named?
print("all have (a,d)=(+-1/2,-+1/2):",
      all((co[0],co[3]) in {(Fr(1,2),Fr(-1,2)),(Fr(-1,2),Fr(1,2))} for co,_,_ in found))

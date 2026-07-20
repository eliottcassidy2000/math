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

# boxeph's exact parity fake on span {+2,+1,-1,-2}: P=(Z+W)(1+(Z-W)/2)
Pbox={(2,0):Fr(1,2),(1,0):Fr(1),(0,1):Fr(1),(0,2):Fr(-1,2)}
print("boxeph parity fake P=(1/2)Z^2+Z+W-(1/2)W^2, charges",charges(Pbox),"two-sided",two_sided(Pbox))
print("  E[P^m], m=1..6:", [str(Epow(Pbox,m)) for m in range(1,7)])
print("  => fakes null through m=3, killed at m=4:", all(Epow(Pbox,m)==0 for m in (1,2,3)) and Epow(Pbox,4)!=0)

# Systematic: all P = a Z^2 + b Z + c W + d W^2, coeffs in half-integers, find every two-sided
# faker (E1=E2=E3=0) and record worst kill-m.  Tests boxeph route-3 "dies by m=4".
half=[Fr(x,2) for x in range(-4,5)]  # -2,-1.5,...,2
mons=[(2,0),(1,0),(0,1),(0,2)]
fakers=0; twosided_fakers=0; worst=0; survivors=[]
for co in product(half,repeat=4):
    if all(c==0 for c in co): continue
    P={m:c for m,c in zip(mons,co) if c!=0}
    if any(Epow(P,m)!=0 for m in (1,2,3)): continue
    fakers+=1
    ts=two_sided(P)
    if ts: twosided_fakers+=1
    kill=next((m for m in range(1,10) if Epow(P,m)!=0),None)
    if kill is None: survivors.append((co,charges(P),ts))
    else: worst=max(worst,kill)
print(f"\nspan {{+2,+1,-1,-2}}, half-int box: fakers(E1=E2=E3=0)={fakers}, two-sided={twosided_fakers}")
print(f"  worst kill-m among fakers that die (m<=9): {worst}")
print(f"  fakers surviving E1..E9=0: {len(survivors)}; all one-sided: {all(not ts for _,_,ts in survivors)}")
for s in survivors[:8]: print("   surv:",[str(x) for x in s[0]],"charges",s[1])

import sys, itertools, cmath, math
from fractions import Fraction
from math import comb
sys.stdout.reconfigure(encoding="utf-8")
def measS7(E):
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]; mid=(lo+hi)/2; hit=set()
        for e in E:
            v=e*mid; v-=v.numerator//v.denominator; hit.add((v.numerator*7)//v.denominator)
        if all(j in hit for j in range(1,7)): tot+=hi-lo
    return tot
def M7(k): return sum(Fraction((-1)**t)*comb(6,t)*Fraction(7-t,7)**(k-1) for t in range(7))
def support6_count(E,H):
    Enz=[e for e in sorted(set(E)) if e!=0]; d=len(Enz); cnt=0
    for nv in itertools.product(range(-H,H+1),repeat=d):
        sup=sum(1 for x in nv if x!=0 and x%7!=0)
        if sup<6: continue
        if sum(n*e for n,e in zip(nv,Enz))!=0: continue
        cnt+=1
    return cnt
bank={"AP":list(range(8)),"near-AP":[0,1,2,3,4,5,6,8],"odd-AP":[0,1,3,5,7,9,11,13],
      "primes":[0,2,3,5,7,11,13,17],"squares":[0,1,4,9,16,25,36,49]}
print("name        corr        R6(H=2)  R6(H=3)")
data=[]
for name,E in bank.items():
    corr=float(measS7(E)-M7(len(E)))
    r2=support6_count(E,2); r3=support6_count(E,3)
    data.append((name,corr,r2,r3))
    print(f"{name:10s} {corr:+.6f}   {r2:6d}  {r3:6d}")
print("\nMonotone check (sorted by corr desc):")
for name,corr,r2,r3 in sorted(data,key=lambda t:-t[1]):
    print(f"  {name:10s} corr={corr:+.6f}  R6(2)={r2}  R6(3)={r3}")

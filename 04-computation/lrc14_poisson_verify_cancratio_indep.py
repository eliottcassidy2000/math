import sys, itertools
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
# P_T = meas of orbit avoiding sectors in T = exact via cells; iid = ((7-|T|)/7)^{k-1}
def P_T_measure(E,Tset):
    # measure of x s.t. NO point falls in any sector of Tset
    E=sorted(set(e for e in E if True)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]; mid=(lo+hi)/2; ok=True
        for e in E:
            v=e*mid; v-=v.numerator//v.denominator; sec=(v.numerator*7)//v.denominator
            if sec in Tset: ok=False; break
        if ok: tot+=hi-lo
    return tot
SUBS=[frozenset(c) for r in range(7) for c in itertools.combinations(range(1,7),r)]
bank={"AP":list(range(8)),"near-AP":[0,1,2,3,4,5,6,8],
      "primes":[0,2,3,5,7,11,13,17],"squares":[0,1,4,9,16,25,36,49]}
print("name      sum_T|P_T-iid|   |corr|     ratio")
for name,E in bank.items():
    k=len(E); corr=measS7(E)-M7(k)
    tot=Fraction(0)
    for T in SUBS:
        PT=P_T_measure(E,set(T)); iid=Fraction(7-len(T),7)**(k-1)
        tot+=abs(PT-iid)
    ratio=float(tot)/abs(float(corr))
    print(f"{name:9s} {float(tot):.6f}     {abs(float(corr)):.6f}   {ratio:.1f}")

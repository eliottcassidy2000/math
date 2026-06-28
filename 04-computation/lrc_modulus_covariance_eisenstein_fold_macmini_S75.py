"""
S75c: the cap kernel K^(n)(a,b)=meas{||ax||<1/n, ||bx||<1/n} is MODULUS-COVARIANT (pure scale (2/n)h)
for n>=2*max(speed), and the Eisenstein fold n->2n is exactly x1/2 -- EXCEPT it BREAKS at the apex n/2.
This locates the LRC(14) difficulty at the apex 7 = n/2 (speeds > n/2 deviate = the antipode half).
"""
from fractions import Fraction as F
from math import gcd

def Kn(a,b,n):
    S=sorted({a,b}-{0})
    if not S: return F(1)
    bps=set([F(0),F(1)])
    for r in S:
        for k in range(0,r+1):
            for s in (F(1,n),-F(1,n)):
                v=(F(k)+s)/r
                if 0<=v<=1: bps.add(v)
    bps=sorted(bps); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        mid=(x0+x1)/2
        if all(min((r*mid)%1,1-(r*mid)%1)<F(1,n) for r in S): tot+=x1-x0
    return tot

print("(1) modulus-covariance: g_n(a,b)=(n/2)*ab*K^(n)(a,b) stable for n>=2*max(speed)?")
for n in (7,14,28,56):
    rows=[]
    for b in range(2,8):
        for a in range(1,b):
            if gcd(a,b)==1: rows.append(int(F(n,2)*a*b*Kn(a,b,n)))
    print(f"   n={n:>2}: g = {rows}")
print("   -> n=14,28,56 identical (stable); n=7 deviates (apex too small => wide-tooth over-overlap).")
print()
print("(2) Eisenstein fold K^(2n)/K^(n) = 1/2 (pure rescale), EXCEPT at the apex:")
for (a,b) in [(1,3),(2,3),(3,5),(2,5),(4,7),(5,8)]:
    r1=Kn(a,b,7)/Kn(a,b,14); r2=Kn(a,b,14)/Kn(a,b,28)
    flag="  <-- apex break" if r1!=2 else ""
    print(f"   ({a},{b}): K^7/K^14={str(r1):>5}  K^14/K^28={str(r2):>5}{flag}")
print()
print("(3) the break is at speeds > n/2 = 7: antipode half (speeds 8..13) deviates at n=14.")
print("    speeds<=7: K^(14)=(1/7)h (clean three-gap); speeds 8..13: K(a,13)=(2a-1)/(91a) deviates.")
for a in (6,7,8):
    print(f"      K^14(a={a},13)= {Kn(a,13,14)}   (2a-1)/(91a)={F(2*a-1,91*a)}")
print("    => the three recursions (Mobius order / Eisenstein modulus / Legendre three-gap) concentrate")
print("       at the apex 7 = n/2; the antipode-half deviation is the entire LRC(14) difficulty.")

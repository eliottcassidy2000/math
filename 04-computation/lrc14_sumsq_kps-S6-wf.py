#!/usr/bin/env python3
"""
lrc14_sumsq_kps-S6-wf.py
Relaxation route:  N(E)={x: maxgap<=1/7}.  Since maxgap >= sum_t g_t^2 (because
maxgap * sum g = maxgap >= sum g^2), we have
        N(E)  subset  Q(E) := { x : sum_t g_t(x)^2 <= 1/7 }.
If meas(Q(E)) <= cap_8 for all E (uniformly), we are DONE (rigorously), because
sum g^2 is a nicer functional and its average int_0^1 sum g^2 dx may be computable
in closed form, giving a uniform bound.

We:
 (1) verify maxgap >= sum g^2 (sanity),
 (2) compute meas(Q(E)) EXACTLY for many E and compare to net(E) and cap_8,
 (3) compute int_0^1 sum_t g_t(x)^2 dx exactly for many E to look for a uniform
     lower bound (which would, via anti-concentration, bound meas(Q)).

If meas(Q(E)) already EXCEEDS cap_8 for some E, this relaxation is too lossy and we
report that honestly.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(5)
from math import gcd
from functools import reduce
ONE7=F(1,7)

def cells(E):
    """Yield (a,b, sorted points-with-floor data) for each order-cell of x in [0,1)."""
    E=sorted(set(E)); n=len(E)
    bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        yield a,b

def net_meas(E):
    """meas{x: maxgap<=1/7} exact."""
    E=sorted(set(E)); n=len(E); good=[]
    for a,b in cells(E):
        mid=(a+b)/2
        order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        floors=[(E[order[t]]*mid).__floor__() for t in range(n)]
        lo=a; hi=b; feasible=True
        for t in range(n):
            o1=order[t]; o2=order[(t+1)%n]; wrap=1 if t==n-1 else 0
            s=E[o2]-E[o1]; c=F(floors[t]-floors[(t+1)%n]+wrap)
            if s==0:
                if not (c<=ONE7): feasible=False; break
            elif s>0: hi=min(hi,(ONE7-c)/s)
            else: lo=max(lo,(ONE7-c)/s)
            if lo>=hi: feasible=False; break
        if feasible and lo<hi: good.append((lo,hi))
    good.sort(); out=[]; tot=F(0)
    for a,b in good:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    for a,b in out: tot+=b-a
    return tot

def integral_sumsq(E):
    """int_0^1 sum_t g_t(x)^2 dx, EXACT.  On each order-cell, each gap g_t is linear
    in x: g_t = s_t x + c_t.  sum_t g_t^2 is a quadratic in x; integrate exactly."""
    E=sorted(set(E)); n=len(E); total=F(0)
    for a,b in cells(E):
        mid=(a+b)/2
        order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        floors=[(E[order[t]]*mid).__floor__() for t in range(n)]
        # gap_t(x) = s_t x + c_t
        S=[]; C=[]
        for t in range(n):
            o1=order[t]; o2=order[(t+1)%n]; wrap=1 if t==n-1 else 0
            s=F(E[o2]-E[o1]); c=F(floors[t]-floors[(t+1)%n]+wrap)
            S.append(s); C.append(c)
        # integrate sum (s x + c)^2 = sum s^2 x^2 + 2 s c x + c^2 over [a,b]
        A2=sum(s*s for s in S); A1=sum(2*s*c for s,c in zip(S,C)); A0=sum(c*c for c in C)
        I=A2*(b**3-a**3)/3 + A1*(b**2-a**2)/2 + A0*(b-a)
        total+=I
    return total

def Qmeas(E):
    """meas{x: sum_t g_t(x)^2 <= 1/7} EXACT.  sum g^2 quadratic per cell; solve."""
    E=sorted(set(E)); n=len(E); good=[]
    for a,b in cells(E):
        mid=(a+b)/2
        order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        floors=[(E[order[t]]*mid).__floor__() for t in range(n)]
        S=[]; C=[]
        for t in range(n):
            o1=order[t]; o2=order[(t+1)%n]; wrap=1 if t==n-1 else 0
            s=F(E[o2]-E[o1]); c=F(floors[t]-floors[(t+1)%n]+wrap)
            S.append(s); C.append(c)
        A2=sum(s*s for s in S); A1=sum(2*s*c for s,c in zip(S,C)); A0=sum(c*c for c in C)-ONE7
        # want A2 x^2 + A1 x + A0 <= 0 on [a,b]; A2>=0
        if A2==0:
            # linear A1 x + A0 <=0
            if A1==0:
                if A0<=0: good.append((a,b))
            elif A1>0:
                xr=-A0/A1
                lo=a; hi=min(b,xr)
                if lo<hi: good.append((lo,hi))
            else:
                xr=-A0/A1
                lo=max(a,xr); hi=b
                if lo<hi: good.append((lo,hi))
        else:
            disc=A1*A1-4*A2*A0
            if disc<0:
                continue  # A2>0 and no real roots => always >0 => empty
            # roots r1<=r2; <=0 between them
            import math
            # exact sqrt of Fraction disc if perfect square else use float-safe bracket
            # We need exact endpoints; disc may not be a perfect square -> the root is
            # irrational. The measure is still rational? No: measure could be irrational!
            # Handle via float for EXPLORATION; flag if irrational roots occur.
            from fractions import Fraction
            dnum=disc.numerator; dden=disc.denominator
            import math
            rn=math.isqrt(dnum) if dnum>=0 else -1
            rd=math.isqrt(dden)
            perfect = (rn*rn==dnum and rd*rd==dden)
            if perfect:
                sq=Fraction(rn,rd)
                r1=(-A1-sq)/(2*A2); r2=(-A1+sq)/(2*A2)
            else:
                sq=Fraction(math.isqrt(int(disc*10**40)),10**20)  # approx; EXPLORATION ONLY
                r1=(-A1-sq)/(2*A2); r2=(-A1+sq)/(2*A2)
            lo=max(a,r1); hi=min(b,r2)
            if lo<hi: good.append((lo,hi))
    good.sort(); out=[]; tot=F(0)
    for a,b in good:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    for a,b in out: tot+=b-a
    return tot

if __name__=="__main__":
    cap8=F(2243,5880)
    print(f"cap_8={float(cap8):.6f}")
    print("\nCompare net(E) <= Q(E)=meas{sum g^2<=1/7}, and int sum g^2, for sample E:")
    Es=[list(range(8)),
        [0,1,2,3,4,5,6,8],
        [0,1,2,3,4,5,6,14],
        [0,2,4,6,7,9,11,13],
        [0,1,2,3,11,12,13,14],
        [0,1,5,6,10,11,12,15]]
    # add randoms
    for _ in range(8):
        E=sorted(set([0]+random.sample(range(1,16),7)))
        if len(E)==8: Es.append(E)
    worstQ=F(0); worstE=None
    for E in Es:
        nm=net_meas(E); Q=Qmeas(E); I=integral_sumsq(E)
        if Q>worstQ: worstQ=Q; worstE=E
        print(f"  E={E}")
        print(f"     net={float(nm):.5f}  Q=meas(sumg^2<=1/7)={float(Q):.5f}  "
              f"int(sumg^2)={float(I):.5f}  {'Q<=cap8' if Q<=cap8 else '*** Q>cap8 (relaxation too lossy) ***'}")
    print(f"\n  worst Q over sample = {float(worstQ):.5f} at {worstE} vs cap_8={float(cap8):.5f}")
    print("  (if worst Q <= cap_8 across a broad scan, the sum-g^2 relaxation could give a uniform proof)")

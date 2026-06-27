#!/usr/bin/env python3
"""Test the owner's resolvent/quintic hint on the miss-count PGF (Lee-Yang object).
miss-PGF P(z)=sum_t p_t z^t (t=empty inner sectors, deg 6); its derivative is a quintic.
Look for: dyadic structure (OCF weights 2,4,8,16), coefficients 120/320, Lee-Yang zeros.
kind-mendel-2026-06-27-S9."""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
def gl(xs): return reduce(lambda a,b:a*b//gcd(a,b),[x for x in xs if x],1)
def miss_pgf(E):
    "exact distribution p_t = meas{N_E=t empty inner sectors 1..6}, return [p_0..p_6]"
    pos=[e for e in E if e]; D=7*gl(pos); bset=set([0,D])
    for e in pos:
        step=D//(7*e); x=0
        while x<=D: bset.add(x); x+=step
    bps=sorted(bset); p=[F(0)]*7
    for a,b in zip(bps,bps[1:]):
        if b<=a: continue
        mid2=a+b; hit=set()
        for e in E: hit.add((7*((e*mid2)%(2*D)))//(2*D))
        N=len([j for j in range(1,7) if j not in hit])
        p[N]+=F(b-a,D)
    return p
def evalpoly(coeffs,z): return sum(c*z**i for i,c in enumerate(coeffs))

for k,lbl in [(7,'consec_7'),(8,'consec_8'),(13,'consec_13 (n=14 AP)')]:
    E=list(range(k)); p=miss_pgf(E)
    print(f"=== {lbl}: miss-PGF P(z)=sum p_t z^t ===")
    print("  p_t =", [str(x) for x in p], " (sum=%s)"%sum(p))
    # rational P, clear denominators to integer polynomial
    den=gl([x.denominator for x in p]); ip=[int(x*den) for x in p]
    print(f"  integer PGF (x{den}):", ip)
    # roots numerically (Lee-Yang)
    import cmath
    # companion via numpy-free Durand-Kerner
    deg=6
    while ip and ip[-1]==0: ip.pop()
    n=len(ip)-1
    if n>=1:
        co=[c/ip[-1] for c in ip]  # monic-ish (float)
        roots=[ (0.4+0.9j)**i for i in range(n)]
        for _ in range(200):
            new=[]
            for i in range(n):
                num=evalpoly(co,roots[i])
                den2=1
                for j in range(n):
                    if j!=i: den2*= (roots[i]-roots[j])
                new.append(roots[i]-num/den2)
            roots=new
        rr=sorted(roots, key=lambda z:abs(z))
        print("  PGF roots (|z|):", [f"{abs(z):.3f}" for z in rr], " min|z|=%.3f"%min(abs(z) for z in rr))
        print("  PGF roots:", [f"{z.real:.3f}{'+' if z.imag>=0 else ''}{z.imag:.3f}i" for z in rr])
print()
print("=== the owner's resolvent quartic x^4+10x^3-120x^2-320x+1024, roots 2,-4,8,-16 ===")
print("  check symmetric fns of {2,-4,8,-16}: e1=%d e2=%d e3=%d e4=%d"%(
      2-4+8-16, 2*-4+2*8+2*-16-4*8-4*-16+8*-16, 2*-4*8+2*-4*-16+2*8*-16-4*8*-16, 2*-4*8*-16))
print("  => roots are 2*(-2)^k (dyadic, ratio -2) = OCF weights 2^d with alternating sign; e2=-120,e3=320,e4=1024=2^10")
print("  120 = max ΔH(n=7)=5!; 320 = #(β3=1 at n=6) = #distinct-H(n=8) = THM-155 denom")

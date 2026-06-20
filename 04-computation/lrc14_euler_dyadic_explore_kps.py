import sys, itertools
from fractions import Fraction
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
# (A) Euler f=g, ground it.
N=20
f=[0]*(N+1); f[0]=1
for n in range(1,N+1):  # f=prod(1+x^n): partitions into DISTINCT parts
    for j in range(N,n-1,-1): f[j]+=f[j-n]
g=[0]*(N+1); g[0]=1
for n in range(1,N+1,2):  # g=prod 1/(1-x^odd): partitions into ODD parts
    for j in range(n,N+1): g[j]+=g[j-n]
print("(A) Euler: distinct-part counts == odd-part counts?", f==g, "  first:", f[:12])

# (B) doubling map on the 7 sectors mod 7: order of 2, the <2>-orbits.
orb={}; seen=set()
for s in range(7):
    if s in seen: continue
    o=[]; x=s
    while x not in o: o.append(x); x=(2*x)%7
    for y in o: seen.add(y)
    orb[s]=o
print("(B) doubling z->2z on sectors mod 7: <2>-orbits =", {min(o):o for o in {tuple(v) for v in orb.values()}})
print("    => 6 inner sectors split into TWO size-3 orbits {1,2,4} and {3,6,5}; 0 fixed. (2 has order 3 mod 7.)")

# (C) Glaisher 2-adic decomposition of consec {1..13}: each e = 2^a * odd.
def dyadic(E):
    d={}
    for e in E:
        b=e; a=0
        while b%2==0: b//=2; a+=1
        d.setdefault(b,[]).append(a)
    return d
print("(C) consec {1..13} dyadic (odd-part -> 2-power exponents):", dyadic(range(1,14)))

# (D) does the 7-sector cover decompose along the two <2>-orbits?
# meas(S7) vs meas{hit all of {1,2,4} (+sector0)} and meas{hit all of {3,5,6}}.
def hits(E,x):
    h=set()
    for e in E:
        v=e*x; v=v-(v.numerator//v.denominator); h.add((v.numerator*7)//v.denominator)
    return h
def meas_cover(E, target):
    """meas{x: all sectors in `target` are hit}."""
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        if target<=hits(E,(lo+hi)/2): tot+=hi-lo
    return tot
A={0,1,2,4}; B={0,3,5,6}; ALL=set(range(7))
for E in [[0,1,2,3,4,5,6,7,8], list(range(14)), [0,1,2,3,4,5,6,7,9]]:
    pA=meas_cover(E,A); pB=meas_cover(E,B); p7=meas_cover(E,ALL)
    # if independent-ish: p7 ?= pA*pB / something
    print(f"   E={E[:9]}{'...' if len(E)>9 else ''}: meas(hit<2>orbit{{1,2,4}})={float(pA):.4f} meas(hit{{3,5,6}})={float(pB):.4f} meas(S7)={float(p7):.4f}  pA*pB={float(pA*pB):.4f}  p7/(pA*pB)={float(p7/(pA*pB)) if pA*pB else 0:.3f}")

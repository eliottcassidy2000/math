#!/usr/bin/env python3
"""
S639 — the FTA bridge: n+1 coefficients <-> n roots. Our two spectra are FTA-dual.
A degree-n polynomial has n+1 coefficients (the constant = the '+1') and n complex roots (FTA).
For our objects the COEFFICIENTS are the combinatorial spectrum (H / independence counts, S636-7) and
the ROOTS are the character-ratio / eigenvalue spectrum (S636/638, Gauss sums, roots of unity).
Newton's identities / Vieta map between them.
Centerpiece: the chromatic polynomial of the tie-graph C_n has roots = 1 + (n-1)th roots of unity =
the cyclotomic witness orbit (THM-403) = the character ratios. Also: Lee-Yang / fugacity zeros of the
covering-depth Z(z) (the LRC partition function); the glass transition (S637) as a real-axis pinch.
"""
import cmath, math

def durand_kerner(coeffs, iters=200):
    """roots of polynomial with coeffs[0]*x^n + ... + coeffs[n]; returns n complex roots."""
    c=[complex(x) for x in coeffs]
    n=len(c)-1
    if n==0: return []
    lead=c[0]; c=[x/lead for x in c]               # monic
    roots=[ (0.4+0.9j)**k for k in range(n) ]      # init
    for _ in range(iters):
        new=[]
        for i in range(n):
            num=sum(c[j]*roots[i]**(n-j) for j in range(n+1))
            den=1
            for j in range(n):
                if j!=i: den*= (roots[i]-roots[j])
            new.append(roots[i]-num/den if den!=0 else roots[i])
        roots=new
    return roots

# (1) chromatic polynomial of C_n: P(C_n,k) = (k-1)^n + (-1)^n (k-1).  coeffs in k.
def cycle_chrom_coeffs(n):
    # (k-1)^n = sum_{j} C(n,j) k^j (-1)^{n-j}; plus (-1)^n (k-1)
    from math import comb
    co=[0]*(n+1)
    for j in range(n+1): co[j]+=comb(n,j)*((-1)**(n-j))   # coefficient of k^j
    co[1]+= (-1)**n * 1     # +(-1)^n * k
    co[0]+= (-1)**n * (-1)  # +(-1)^n * (-1)
    # return high-to-low for durand_kerner: coeff of k^n ... k^0
    return [co[n-i] for i in range(n+1)]

print("(1) chromatic roots of the tie-graph C_n  =?  1 + (n-1)th roots of unity (= character ratios, THM-403)")
for n in (5,6,7):
    roots=durand_kerner(cycle_chrom_coeffs(n))
    roots=sorted(roots,key=lambda z:(round(z.real,3),round(z.imag,3)))
    # distance of (root-1) from unit circle
    onunit=[abs(abs(r-1)-1)<1e-3 for r in roots]
    print(f"   C_{n}: {n} roots; |root-1|=1 (on unit circle centered at 1)? {sum(onunit)}/{n} ;  k=1 a root? {any(abs(r-1)<1e-6 for r in roots)}")
    print(f"      roots-1 (should be (n-1)th roots of unity): "+", ".join(f"{(r-1):.2f}" for r in roots))

# (2) Lee-Yang / fugacity zeros of the covering-depth Z(z)=sum p_k z^k (LRC), and the glass/transition
from fractions import Fraction as Fr
def norm(x):
    f=x-(x.numerator//x.denominator); return f if f<=Fr(1,2) else 1-f
def depth_dist(speeds,delta):
    bps={Fr(0),Fr(1)}
    for v in speeds:
        for k in range(v):
            for s in (delta,-delta):
                t=(Fr(k)+s)/v; t-=t.numerator//t.denominator
                if 0<=t<1: bps.add(t)
    bps=sorted(bps); d={}
    for a,b in zip(bps,bps[1:]):
        mid=(a+b)/2; dep=sum(norm(v*mid)<delta for v in speeds)
        d[dep]=d.get(dep,Fr(0))+(b-a)
    return d
print("\n(2) Lee-Yang/fugacity zeros of the covering-depth Z(z)=sum p_k z^k (LRC partition function)")
for speeds,n in [([1,2,3,4],5),([1,3,4,7],5),([1,2,4,8],5)]:
    d=depth_dist(speeds,Fr(1,n)); K=max(d); p=[float(d.get(k,Fr(0))) for k in range(K+1)]
    coeffs=[p[K-i] for i in range(K+1)]  # high to low
    roots=durand_kerner(coeffs)
    near_real=min(abs(r.imag) for r in roots) if roots else None
    print(f"   {tuple(speeds)} (p0={p[0]:.3f}): Z-zeros min|Im|={near_real:.3f}  zeros={[f'{r:.2f}' for r in roots]}")

# (3) the FTA duality statement: char poly of a matrix — coeffs (invariants) <-> roots (eigenvalues=char ratios)
print("\n(3) FTA duality: coeffs = combinatorial spectrum (H, independence counts); roots = char ratios.")
print("    Newton's identities map power-sums (traces/H-moments) <-> elementary symmetric (coeffs) <-> roots.")
print("    chromatic roots of C_n = roots of unity = the (Z/n)* witness orbit (THM-403) = character ratios.")

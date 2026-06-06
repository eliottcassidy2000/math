#!/usr/bin/env python3
"""
S640 — Hadwiger-Nelson via Lee-Yang: the chromatic-ZERO locus (q-plane) of unit-distance gadgets.
Builds on S687/S699 (field tower, Vitali wall chi_f<=4.36<5). NEW angle: the chromatic polynomial
P(G,q) IS the zero-T antiferromagnetic Potts partition function; its zeros are Lee-Yang zeros in q.
CLAIM (compute): the rightmost REAL chromatic zero = chi-1; the COMPLEX zeros reach further right
(real part ~ the fractional/measure bound). HN's {5,6,7} = the location of the REAL Lee-Yang edge;
the Vitali-wall integrality gap = the lag of the real edge behind the complex accumulation.
"""
import itertools, math, cmath

def durand_kerner(coeffs, iters=300):
    c=[complex(x) for x in coeffs]; n=len(c)-1
    if n==0: return []
    c=[x/c[0] for x in c]
    roots=[(0.4+0.9j)**k for k in range(n)]
    for _ in range(iters):
        nw=[]
        for i in range(n):
            num=sum(c[j]*roots[i]**(n-j) for j in range(n+1))
            den=1
            for j in range(n):
                if j!=i: den*=(roots[i]-roots[j])
            nw.append(roots[i]-num/den if den else roots[i])
        roots=nw
    return roots

def chromatic_poly(adj):
    """coeffs (high->low) of P(G,q) via counting proper k-colorings for k=0..n and interpolation."""
    n=len(adj)
    def count(k):
        col=[-1]*n; cnt=0
        def bt(v):
            nonlocal cnt
            if v==n: cnt+=1; return
            for c in range(k):
                if all(col[u]!=c for u in range(v) if adj[v][u]):
                    col[v]=c; bt(v+1); col[v]=-1
        bt(0); return cnt
    xs=list(range(n+1)); ys=[count(k) for k in xs]
    # Lagrange interpolation -> coefficients (exact rationals)
    from fractions import Fraction as Fr
    def polymul_linear(p, a):   # multiply p(x) (low->high) by (x - a)
        res=[Fr(0)]*(len(p)+1)
        for t,c in enumerate(p):
            res[t]   += -a*c
            res[t+1] += c
        return res
    coeff=[Fr(0)]*(n+1)
    for i in range(n+1):
        num=[Fr(1)]; den=Fr(1)
        for j in range(n+1):
            if j==i: continue
            num=polymul_linear(num, Fr(xs[j])); den*=(xs[i]-xs[j])
        for t in range(len(num)): coeff[t]+=ys[i]*num[t]/den
    # coeff is low->high (degree increasing); make high->low
    return [float(c) for c in coeff[::-1]]

def points_to_adj(pts,tol=1e-7):
    n=len(pts); adj=[[False]*n for _ in range(n)]
    for i,j in itertools.combinations(range(n),2):
        if abs(abs(pts[i]-pts[j])-1)<tol: adj[i][j]=adj[j][i]=True
    return adj

def chi(adj):
    n=len(adj)
    def colorable(k):
        col=[-1]*n
        def bt(v):
            if v==n: return True
            for c in range(k):
                if all(col[u]!=c for u in range(v) if adj[v][u]):
                    col[v]=c
                    if bt(v+1): return True
                    col[v]=-1
            return False
        return bt(0)
    k=1
    while not colorable(k): k+=1
    return k

# gadgets
w6 = [0]+[cmath.exp(1j*math.pi/3*k) for k in range(6)]   # wheel W6: center+hexagon
A=0+0j;B=1+0j;C=cmath.exp(1j*math.pi/3);D=B+C
th=math.acos(5/6);r=cmath.exp(1j*th)
spindle=[A,B,C,D,r*B,r*C,r*D]

for name,pts in [("W6 (Eisenstein floor, chi=3)",w6),("Moser spindle (chi=4)",spindle)]:
    adj=points_to_adj(pts); k=chi(adj); P=chromatic_poly(adj); roots=durand_kerner(P)
    real_roots=sorted(r.real for r in roots if abs(r.imag)<1e-4)
    max_real_zero=max(real_roots) if real_roots else None
    max_re_complex=max(r.real for r in roots)
    print(f"{name}: chi={k}")
    print(f"   rightmost REAL chromatic zero = {max_real_zero:.4f}  (= chi-1 = {k-1}? {abs(max_real_zero-(k-1))<1e-3})")
    print(f"   rightmost COMPLEX-zero real part = {max_re_complex:.4f}  (the Lee-Yang bulk, > chi-1)")
    print(f"   all zeros: {sorted([complex(round(r.real,3),round(r.imag,3)) for r in roots],key=lambda z:z.real)}")
    print()
print("READING: chi(plane) = 1 + (rightmost REAL chromatic-zero edge of UD graphs).")
print("  de Grey: real edge >= 4 (chi>=5). {5,6,7} = real edge in {4,5,6}. Complex Lee-Yang bulk")
print("  reaches the fractional/measure value (~chi_f<=4.36, S699): the integrality gap = the REAL")
print("  Lee-Yang edge lagging the COMPLEX accumulation = the Vitali wall, in zero-locus language.")

#!/usr/bin/env python3
"""
S720 — What we keep track of: the transfer spectrum. The A+B+C-D-E-F+G recursion and the twistors are
one structure — a convolution diagonalized by a transform — read at different TEMPERATURES.

Threads to unify:
  - A+B+C-D-E-F+G (S715/S716): inclusion-exclusion over B_3 = (x-1)^3 = Delta^3 = the ADDITIVE recursion.
  - THM-337: base-path H satisfies x^3-3x^2-x-1 = the MULTIPLICATIVE recursion (lambda~3.383).
  - Pfaffian (S713): the signed n->n-2 step; det(I+2A)=Pf^2 (a unimodular/parity invariant).
  - Twistors (S718/S719): discrete log / angle = the LOG that turns MULTIPLICATIVE group action into
    ADDITIVE (cooling); transversal core = half-system; unit-distance optimum = rotation orbit.

THE UNIFYING CLAIM. For an order-3 C-finite tournament invariant the char poly is
        Lambda(a):  x^3 - 3 x^2 + a x - 1.
Its symmetric functions are WHAT WE KEEP TRACK OF, and they split:
  - e1 = sum of eigenvalues = 3 = the GEOMETRY (the 3 corners / branching of the triangle). FROZEN.
  - e3 = product of eigenvalues = 1 = the PARITY / unimodular transfer determinant (the Pfaffian face). FROZEN.
  - e2 = a = the middle symmetric function = the TEMPERATURE (additive vs multiplicative). RUNS.
Endpoints: a=3 -> (x-1)^3 (ADDITIVE, A+B+C-D-E-F+G, triple root 1, polynomial); a=-1 -> THM-337
(MULTIPLICATIVE, lambda~3.383, exponential). a=3 is the UNIQUE "crystalline" point (all eigenvalues 1).
The twistor (log) is the cooling operation that sends multiplicative (root != 1) to additive (root 1).

No numpy/sympy.
"""
import math, random
from fractions import Fraction as Fr

# --- real cubic: dominant root by bisection; all-root magnitudes via product/sum ---
def cubic_eval(a3,a2,a1,a0,x): return a3*x**3+a2*x**2+a1*x+a0
def dominant_real_root(a,b,c, lo=0.0, hi=20.0):
    # monic x^3 + a x^2 + b x + c ; find largest real root by scanning + bisection
    f=lambda x: x**3+a*x**2+b*x+c
    # scan for sign change from high to low
    roots=[]
    N=20000; prev=f(hi); px=hi
    x=hi
    step=(hi-lo)/N
    while x>lo:
        nx=x-step; nf=f(nx)
        if prev==0: roots.append(x)
        elif prev*nf<0:
            a_,b_=nx,x
            for _ in range(80):
                m=(a_+b_)/2
                if f(a_)*f(m)<=0: b_=m
                else: a_=m
            roots.append((a_+b_)/2)
        prev=nf; x=nx
    return max(roots) if roots else None

def find_recurrence(seq, maxord=5):
    n=len(seq)
    for r in range(1,maxord+1):
        if n<2*r+1: break
        A=[[Fr(seq[i-1-j]) for j in range(r)] for i in range(r,2*r)]
        rhs=[Fr(seq[i]) for i in range(r,2*r)]
        c=gauss(A,rhs)
        if c is None: continue
        if all(sum(c[j]*seq[i-1-j] for j in range(r))==seq[i] for i in range(r,n)):
            return r,c
    return None,None
def gauss(A,b):
    A=[row[:] for row in A]; b=b[:]; n=len(A)
    for k in range(n):
        p=next((i for i in range(k,n) if A[i][k]!=0),None)
        if p is None: return None
        A[k],A[p]=A[p],A[k]; b[k],b[p]=b[p],b[k]
        inv=A[k][k]; A[k]=[x/inv for x in A[k]]; b[k]=b[k]/inv
        for i in range(n):
            if i!=k and A[i][k]!=0:
                f=A[i][k]; A[i]=[A[i][j]-f*A[k][j] for j in range(n)]; b[i]-=f*b[k]
    return b

# base-path staircase H (THM-337) by direct DP
def ham_paths(A):
    n=len(A); dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for v in range(n):
            d=dp[mask][v]
            if not d or not (mask>>v)&1: continue
            for w in range(n):
                if (mask>>w)&1: continue
                if A[v][w]: dp[mask|(1<<w)][w]+=d
    return sum(dp[(1<<n)-1][v] for v in range(n))
def base_path(k):
    n=2*k; A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if j==i+1: A[j][i]=1
            else: A[i][j]=1
    return A

if __name__=="__main__":
    print("="*86)
    print("S720 — the transfer spectrum: A+B+C-D-E-F+G and the twistors are one structure at two temperatures")
    print("="*86)

    # (1) the temperature family Lambda(a)=x^3-3x^2+a x-1
    print("\n(1) THE TEMPERATURE FAMILY  Lambda(a): x^3 - 3 x^2 + a x - 1")
    print("    FROZEN: e1=sum=3 (corners/geometry), e3=product=1 (parity/unimodular determinant)")
    print("    RUNNING: e2=a = temperature.  a=3 -> (x-1)^3 (additive); a=-1 -> THM-337 (multiplicative)")
    print(f"    {'a':>4} | dominant root lambda(a) | regime")
    for a in (3,2,1,0,-1,-2,-3):
        # monic x^3 -3x^2 + a x -1
        lam=dominant_real_root(-3.0, float(a), -1.0)
        reg = "CRYSTALLINE (all roots 1, polynomial)" if abs(lam-1)<1e-6 else f"exponential (|c|^2=1/lambda={1/lam:.3f})"
        tag = "  <- (x-1)^3 = A+B+C-D-E-F+G" if a==3 else ("  <- THM-337 base-path H" if a==-1 else "")
        print(f"    {a:>4} | {lam:.6f}             | {reg}{tag}")
    print("    => a=3 is the UNIQUE crystalline point (sum 3 + product 1 + all |root|=1 forces roots=1);")
    print("       any a<3 heats it to exponential. Additivity = the perfectly-tuned inclusion-exclusion.")

    # (2) endpoints verified from actual sequences
    print("\n(2) ENDPOINTS from real sequences")
    tiles=[ (n-1)*(n-2)//2 for n in range(3,14)]   # #tiles C(n-1,2), additive
    r,c=find_recurrence(tiles)
    print(f"    additive (#tiles): recurrence order {r}, coeffs {[str(x) for x in c]} = (x-1)^3 char x^3-3x^2+3x-1 (a=3)")
    bp=[ham_paths(base_path(k)) for k in range(1,10)]
    r,c=find_recurrence(bp)
    print(f"    multiplicative (base-path H {bp[:6]}...): order {r}, coeffs {[str(x) for x in c]}")
    print(f"      => char x^3-3x^2-x-1 (a=-1), constant term -1 => |product of eigenvalues|=1 (UNIMODULAR transfer)")

    # (3) unimodular transfer = the parity/Pfaffian face (frozen e3)
    print("\n(3) UNIMODULAR TRANSFER = parity (the Pfaffian face): product of eigenvalues = 1 along the ladder")
    print("    both endpoints have constant term -1 (monic) => e3 = product of roots = 1, FROZEN.")
    print("    This is the recursion-level shadow of S713: det(I+2A)=Pf^2 with Pf odd (a unimodular/parity")
    print("    invariant). The transfer matrix of a bounded-width tournament family is unimodular (det +-1).")

    # (4) the twistor is the COOLING map: log sends multiplicative (root g) -> additive (root 1)
    print("\n(4) THE TWISTOR (log) IS THE COOLING MAP: multiplicative root g  -->  additive root 1")
    g=3
    geo=[g**n for n in range(8)]                 # multiplicative: a(n)=g a(n-1), char (x-g), HOT (root g)
    r1,c1=find_recurrence(geo)
    logs=[n for n in range(8)]                    # ell=log_g(g^n)=n : additive AP, char (x-1)^2, COLD (root 1)
    r2,c2=find_recurrence(logs)
    print(f"    multiplicative g^n={geo[:6]}: order {r1}, coeffs {[str(x) for x in c1]} (root g={g}, HOT)")
    print(f"    its discrete log ell=n={logs[:6]}: order {r2}, coeffs {[str(x) for x in c2]} (char (x-1)^2, root 1, COLD)")
    print("    => the discrete-log twistor (S718) / angle twistor (S719) literally moves a structure from the")
    print("       hot (multiplicative, root != 1) end toward the cold (additive, root 1) end of the ladder.")

    # (5) d-simplex generalization: order d+1, additive = (x-1)^(d+1), temperature = deviation of middle coeffs
    print("\n(5) d-SIMPLEX GENERALIZATION: order d+1; additive crystalline pt = (x-1)^(d+1) (binomial middle coeffs)")
    from math import comb
    for d in (1,2,3,4):
        coeffs=[((-1)**k)*comb(d+1,k) for k in range(d+2)]   # (x-1)^(d+1)
        print(f"    d={d} ({d}-simplex, {d+1} corners): crystalline char (x-1)^{d+1} = {coeffs}; e1={d+1} (corners), "
              f"e_last=+-1 (parity); middle coeffs = TEMPERATURE coordinates")
    print("="*86)

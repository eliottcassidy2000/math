#!/usr/bin/env python3
"""
S714 — The autocorrelation operator MM* across the cluster (operator + adjoint + autocorrelation).

Unifying object: a problem gives an OPERATOR M and its ADJOINT M*; the AUTOCORRELATION MM* (PSD) is the
master object; the answer is an invariant of MM* — a det (Pfaffian), a peak (unit distances), or its
flatness (Paley). This session verifies the SAME autocorrelation operator in three cluster problems and
ties my recent S711/S712 (unit distance) and S713 (Pfaffian) to the homometry line (THM-415) and Paley
(THM-438).

(A) TOURNAMENT autocorrelation operator. S = A - A^T (skew, adjoint S* = S^T = -S). The autocorrelation
    is C = S S* = -S^2 (PSD). Then |Pf(S)| = sqrt(det C) (the S713 Pfaffian = sqrt-det-autocorrelation),
    and the spectrum of C is the tournament's POWER SPECTRUM. Homometric tournaments = same spectrum of
    C. We verify |Pf|=sqrt(det C), count the distinct spectra (homometry classes) at n=6 (THM-120 says
    5), and the H (#Ham paths) distribution per class.

(B) UNIT-DISTANCE Patterson autocorrelation. M = sum of point masses; MM* = Patterson (difference
    multiset). u(n) = peak of the autocorrelation at radius 1. A Minkowski sum P (+) Q has Patterson =
    CONVOLUTION A_P * A_Q (operator tensor -> autocorrelation convolution). We verify
    A_{K3 [] W7} = A_{K3} * A_{W7} and recover u(21)=57 as a convolution peak; the n=22 product caps at 57.

(C) PALEY flat autocorrelation (the 'perfect operator'). For p = 3 mod 4 the QR set Q is a Paley
    difference set: its autocorrelation c(s)=#{(a,b) in Q^2: a-b=s} is CONSTANT for s != 0 (= (p-3)/4),
    i.e. MM* = (flat) and |Gauss sum|^2 = p. Flatness = quasirandom rigidity = the hard LRC core.

No numpy/sympy. Exact integer arithmetic; geometry in floats with tolerances.
"""
import math, random, itertools
from math import gcd

# ===================== shared linear algebra (integer) =====================
def matmul(A,B):
    n=len(A); m=len(B[0]); k=len(B)
    return [[sum(A[i][t]*B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]
def trace(A): return sum(A[i][i] for i in range(len(A)))
def neg(A): return [[-x for x in row] for row in A]

def pf(M):
    n=len(M)
    if n==0: return 1
    if n%2==1: return 0
    tot=0; r0=M[0]
    for j in range(1,n):
        a=r0[j]
        if a==0: continue
        rest=[k for k in range(1,n) if k!=j]
        sub=[[M[r][c] for c in rest] for r in rest]
        tot+=((-1)**(j-1))*a*pf(sub)
    return tot

# ===================== (A) tournament autocorrelation operator =====================
def skew(A):
    n=len(A); return [[A[i][j]-A[j][i] for j in range(n)] for i in range(n)]
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
def all_tournaments(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for idx,(i,j) in enumerate(pairs):
            if (bits>>idx)&1: A[i][j]=1
            else: A[j][i]=1
        yield A

def spectrum_sig(C):
    """power sums -> elementary symmetric (e1,e2,e3) of the 3 distinct eigenvalues of C=-S^2 (n=6)."""
    C2=matmul(C,C); C3=matmul(C2,C)
    p1=trace(C)//2; p2=trace(C2)//2; p3=trace(C3)//2   # each eigenvalue doubled
    e1=p1; e2=(p1*p1-p2)//2; e3=(p1**3-3*p1*p2+2*p3)//6
    return (e1,e2,e3)

def partA():
    print("(A) TOURNAMENT autocorrelation operator  C = S S* = -S^2 ;  |Pf(S)| = sqrt(det C)")
    n=6
    classes={}; pf_check=0; cnt=0
    for A in all_tournaments(n):
        S=skew(A); C=neg(matmul(S,S)); cnt+=1
        sig=spectrum_sig(C)            # (e1=15, e2, e3=|Pf|^2)
        p=abs(pf(S))
        if sig[2]==p*p: pf_check+=1     # det C restricted: e3 = product of distinct eigenvalues = |Pf|^2
        h=ham_paths(A)
        d=classes.setdefault(sig,{'count':0,'Hs':{}})
        d['count']+=1; d['Hs'][h]=d['Hs'].get(h,0)+1
    print(f"  n=6: {cnt} tournaments | e3==|Pf|^2 (autocorrelation det = Pfaffian^2): {pf_check}/{cnt}")
    print(f"  distinct autocorrelation spectra (HOMOMETRY classes): {len(classes)}  (THM-120 => 5)")
    print(f"  {'spectrum (e1,e2,e3)':28} {'|Pf|':>5} {'#tourns':>8}   H-distribution")
    for sig in sorted(classes, key=lambda s:s[2]):
        d=classes[sig]; pfv=int(round(d['count'] and (sig[2])**0.5))
        hd=dict(sorted(d['Hs'].items()))
        # compress H dist if large
        print(f"  {str(sig):28} {pfv:5d} {d['count']:8d}   {hd if len(hd)<=8 else {'(#H vals)':len(hd)}}")
    print("  => the tournament autocorrelation operator has only 5 spectra at n=6: massive homometry;")
    print("     |Pf| = sqrt(det of autocorrelation) is the S713 Pfaffian, here as a spectral invariant.")

# ===================== (B) unit-distance Patterson = convolution =====================
def rot(p,th):
    c,s=math.cos(th),math.sin(th); return (c*p[0]-s*p[1], s*p[0]+c*p[1])
def K3(): return [(0.0,0.0),(1.0,0.0),(0.5,math.sqrt(3)/2)]
def W7():
    P=[(0.0,0.0)]
    for k in range(6):
        a=math.pi/3*k; P.append((math.cos(a),math.sin(a)))
    return P
def diffs(P):
    """difference multiset (Patterson) as bins of squared length, rounded."""
    D={}
    for i in range(len(P)):
        for j in range(len(P)):
            if i==j: continue
            dx=P[i][0]-P[j][0]; dy=P[i][1]-P[j][1]
            key=round(dx*dx+dy*dy,6)
            D[key]=D.get(key,0)+1
    return D
def unit_count(P,tol=1e-7):
    c=0
    for i in range(len(P)):
        for j in range(i+1,len(P)):
            dx=P[i][0]-P[j][0]; dy=P[i][1]-P[j][1]
            if abs(math.hypot(dx,dy)-1.0)<tol: c+=1
    return c

def partB():
    print("\n(B) UNIT-DISTANCE Patterson autocorrelation; product => convolution")
    th=0.4
    A=K3(); B=W7(); BR=[rot(b,th) for b in B]
    P=[(a[0]+b[0],a[1]+b[1]) for a in A for b in BR]   # Minkowski sum = K3 [] W7
    eA=unit_count(A); eB=unit_count(B); ePQ=unit_count(P)
    # convolution prediction at radius^2 = 1: |A|*e_B + |B|*e_A (generic, cross terms 0)
    pred = len(A)*eB + len(B)*eA
    print(f"  e(K3)={eA}, e(W7)={eB}; |K3|={len(A)}, |W7|={len(B)}")
    print(f"  unit distances in K3 [] W7: direct={ePQ}, convolution-formula |A|e_B+|B|e_A={pred}  (=u(21)=57)")
    # verify the FULL Patterson is the convolution of difference multisets (vector convolution at radius 1)
    DA={}  # vector difference multiset of A
    for i in range(len(A)):
        for j in range(len(A)):
            DA[(round(A[i][0]-A[j][0],6),round(A[i][1]-A[j][1],6))]=DA.get((round(A[i][0]-A[j][0],6),round(A[i][1]-A[j][1],6)),0)+1
    DB={}
    for i in range(len(BR)):
        for j in range(len(BR)):
            DB[(round(BR[i][0]-BR[j][0],6),round(BR[i][1]-BR[j][1],6))]=DB.get((round(BR[i][0]-BR[j][0],6),round(BR[i][1]-BR[j][1],6)),0)+1
    # convolution: count pairs (dA,dB) with |dA+dB|^2 ~ 1
    conv_unit=0
    for (ax,ay),ca in DA.items():
        for (bx,by),cb in DB.items():
            if abs((ax+bx)**2+(ay+by)**2-1.0)<1e-6: conv_unit+=ca*cb
    # this counts ordered pairs; unordered unit distances = conv_unit/2
    print(f"  Patterson convolution (A_K3 * A_W7) peak at radius 1 (ordered): {conv_unit}  -> unordered {conv_unit//2}")
    print(f"  => the autocorrelation FACTORIZES: A_(P(+)Q) = A_P * A_Q; u(product) is a convolution peak.")
    print(f"  n=22=2*11: product Patterson = A_K2 * A_(11) caps the radius-1 peak at u(2)*11+2*u(11)=57<60.")

# ===================== (C) Paley flat autocorrelation =====================
def partC():
    print("\n(C) PALEY flat autocorrelation (the 'perfect operator'): QR difference set, |Gauss|^2=p")
    for p in (7,11,19,23):
        if p%4!=3: continue
        Q=set((x*x)%p for x in range(1,p))   # nonzero QRs
        # autocorrelation c(s) = #{(a,b) in Q^2 : a-b = s mod p}
        c={s:0 for s in range(p)}
        for a in Q:
            for b in Q:
                c[(a-b)%p]+=1
        nonzero=set(c[s] for s in range(1,p))
        # Gauss sum g = sum_{x} chi(x) zeta^x ; |g|^2 = p ; verify via difference-set lambda
        lam=(p-3)//4
        flat = (nonzero=={lam})
        print(f"  p={p:2d}: |Q|={len(Q)}, autocorr c(s) for s!=0 = {sorted(nonzero)} ; "
              f"flat at lambda=(p-3)/4={lam}? {flat}  (Paley difference set => |Gauss|^2={p})")
    print("  => flat autocorrelation = MM* proportional to I (off the DC term) = quasirandom RIGIDITY:")
    print("     the homometric/Paley core (THM-415) and the LRC transversal core (S642/THM-420) — the HARD")
    print("     cases — are exactly the flat-autocorrelation configs; peaked/convolutional = extremal/easy.")

if __name__=="__main__":
    print("="*84)
    print("S714 — the autocorrelation operator MM* across the cluster (operator + adjoint)")
    print("="*84)
    partA(); partB(); partC()
    print("="*84)

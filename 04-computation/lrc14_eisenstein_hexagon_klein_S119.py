#!/usr/bin/env python3
"""
klein-2026-07-03-S119 - THE HEXAGONAL GEOMETRY OF THE OPEN CORE (LRC covering-min).

Thesis: the LRC(n) covering-min certificate is a HEXAGONAL-LATTICE object.
  * Phi6(n) = n^2 - n + 1 = N(n - omega), the norm of the Eisenstein integer (n - omega),
    omega = e^{i pi/3} a primitive 6th root of unity (omega^2 = omega - 1).
  * Ring iso  Z[omega]/(n - omega)  ==  Z/Phi6(n),  omega |-> n.
    Hence MULTIPLICATION BY n mod Phi6(n) IS multiplication by omega = the 60-degree
    rotation of the triangular (Eisenstein) lattice.
  * Deep well DW(n) = {1,...,n-2} U {n^2-n}  (n-1 speeds; n^2-n = Phi6(n)-1 = -1 mod Phi6).
    At t* = n/Phi6(n): M(DW) = min_i ||v_i t*|| = n/Phi6(n) > 1/n, margin (n-1)/(n*Phi6(n)),
    and the phases are the omega-rotation orbit -- a hexagonal star.

Everything exact (Fraction / integer).  Verifies claims 1-4.
"""
from fractions import Fraction as F
from math import gcd
import itertools

def phi6(n): return n*n - n + 1

def cdist_q(a, q):
    r = a % q
    return F(min(r, q - r), q)

def order_mod(a, m):
    if gcd(a, m) != 1: return None
    x, k = a % m, 1
    while x != 1:
        x = (x * a) % m; k += 1
        if k > m + 2: return None
    return k

def factor(m):
    f=[]; d=2
    while d*d<=m:
        while m%d==0: f.append(d); m//=d
        d+=1
    if m>1: f.append(m)
    return f

print("="*80)
print("CLAIM 1: n is a primitive 6th root of unity mod Phi6(n); gcd(n,Phi6(n))=1.")
print("="*80)
print(f"{'n':>4} {'Phi6(n)':>9} {'factorization':>18} {'gcd':>4} {'ord(n mod Phi6)':>16}")
all_ord6 = True
for n in range(2, 25):
    q = phi6(n); o = order_mod(n, q)
    if n>=2 and o!=6 and not (n==2 and o==6): all_ord6 = all_ord6 and (o==6)
    print(f"{n:>4} {q:>9} {'*'.join(map(str,factor(q))):>18} {gcd(n,q):>4} {str(o):>16}")
print(f"  --> ord(n)=6 for all n>=2 (n a primitive 6th root mod Phi6(n)):  "
      f"{all(order_mod(n,phi6(n))==6 for n in range(2,25))}")

print()
print("="*80)
print("CLAIM 2: mult-by-n mod Phi6(n) == mult-by-omega == 60-degree hexagonal rotation.")
print("  In Z[omega], omega^2=omega-1, so M_omega (basis {1,omega}) = [[0,-1],[1,1]].")
print("="*80)
def matmul(A,B):
    return [[A[0][0]*B[0][0]+A[0][1]*B[1][0], A[0][0]*B[0][1]+A[0][1]*B[1][1]],
            [A[1][0]*B[0][0]+A[1][1]*B[1][0], A[1][0]*B[0][1]+A[1][1]*B[1][1]]]
Mw=[[0,-1],[1,1]]; P=[[1,0],[0,1]]
powers=[]
for k in range(1,7):
    P=matmul(P,Mw); powers.append([row[:] for row in P])
tr, det = Mw[0][0]+Mw[1][1], Mw[0][0]*Mw[1][1]-Mw[0][1]*Mw[1][0]
print(f"  M_omega=[[0,-1],[1,1]]  trace={tr}  det={det}")
print(f"  characteristic poly = x^2 - ({tr})x + ({det}) = x^2 - x + 1 = Phi6(x)  "
      f"==> eigenvalues are the primitive 6th roots of unity")
print(f"  trace = 1 = 2*cos(60deg)  ==> rotation by exactly 60 degrees")
print(f"  M_omega^6 = {powers[5]} = I ?  {powers[5]==[[1,0],[0,1]]}   (order exactly 6: "
      f"{[powers[k]==[[1,0],[0,1]] for k in range(6)]})")

print()
print("="*80)
print("CLAIM 3: DW(n)={1..n-2, n^2-n}; at t*=n/Phi6(n): M(DW)=n/Phi6(n)>1/n,")
print("         margin (n-1)/(n*Phi6(n)); t* is the maximizer (best time).")
print("="*80)
print(f"{'n':>4} {'M(DW)@t*':>13} {'=n/Phi6?':>9} {'1/n':>11} {'margin':>16} {'=(n-1)/nPhi6?':>13} {'t*optimal':>9}")
for n in range(4, 21):
    q=phi6(n); DW=list(range(1,n-1))+[n*n-n]; tstar=F(n,q)
    M_ts = min(cdist_q(v*n, q) for v in DW)
    target=F(n,q); margin=M_ts-F(1,n); margin_formula=F(n-1, n*q)
    # t* optimal? scan a/Q, Q up to 2*Phi6 (contains n/q and rules out nearby better times)
    best=F(0)
    for Q in range(2, 2*q+1):
        for a in range(1, Q//2+1):
            if gcd(a,Q)!=1: continue
            m=min(cdist_q(v*a, Q) for v in DW)
            if m>best: best=m
        if best>target: break
    print(f"{n:>4} {str(M_ts):>13} {str(M_ts==target):>9} {str(F(1,n)):>11} "
          f"{str(margin):>16} {str(margin==margin_formula):>13} {str(best==target):>9}")

print()
print("="*80)
print("CLAIM 3b: the phases {v*n mod Phi6(n)} for v in DW(n) form the omega-rotation ORBIT")
print("  (a hexagonal star).  Show for n=14: the sorted phase gaps and the 6-fold structure.")
print("="*80)
n=14; q=phi6(n); DW=list(range(1,n-1))+[n*n-n]
phases=sorted(set((v*n)%q for v in DW))
print(f"  n={n}, Phi6={q}, DW={DW}")
print(f"  phases (v*14 mod 183), sorted: {phases}")
gaps=[phases[(i+1)%len(phases)]-phases[i] + (q if i==len(phases)-1 else 0) for i in range(len(phases))]
print(f"  cyclic gaps: {gaps}   max gap={max(gaps)}/{q}={float(max(gaps))/q:.4f}  (2h=1/7={1/7:.4f})")
print(f"  min ||phase||: {min(min(p,q-p) for p in phases)}/{q}  (= M(DW)*Phi6 = numerator of M)")
# hexagonal reps: lift each residue r to the Eisenstein integer of least norm with x+14y=r mod 183
def eisen_norm(x,y): return x*x - x*y + y*y   # N(x+y omega), omega=exp(i pi/3)
def short_rep(r,q,n):
    best=None
    B=20
    for y in range(-B,B+1):
        # x + n y == r (mod q)  => x == r - n y (mod q)
        x0=(r - n*y)%q
        for x in (x0, x0-q):
            nm=eisen_norm(x,y)
            if best is None or nm<best[0]: best=(nm,x,y)
    return best
print("  Eisenstein short reps (norm, x, y) of each phase residue  [norm ~ (dist*Phi6)^2 scale]:")
for r in phases:
    nm,x,y=short_rep(r,q,n)
    print(f"    r={r:>3}: N={nm:>5}  ({x:>3}+{y:>3} omega)")

def M_family(S, Qcap=None):
    best=F(0); B=Qcap if Qcap else 2*(max(S)+1)
    for Q in range(2,B+1):
        for a in range(1,Q//2+1):
            if gcd(a,Q)!=1: continue
            m=min(cdist_q(v*a,Q) for v in S)
            if m>best: best=m
    return best
def is_covering(S,n):
    # mac-mini THM-610 definition: every q in 2..n divides at least one runner
    # (equivalently: no witness at any shallow modulus q<=n, since a q-divisible runner
    #  sits at 0 for all t=a/q).  This EXCLUDES the AP {1..n-1} (no multiple of n).
    for qq in range(2,n+1):
        if not any(v % qq == 0 for v in S): return False
    return True

print()
print("="*80)
print("CLAIM 4a (NEW, provable): covering FORCES the deep-well shape.")
print("  For shape {1..n-2, X}: covering <=> lcm(n-1,n)=n(n-1) | X, since 1..n-2 already")
print("  cover 2..n-2 but NEVER n-1 or n, so ONLY X can, requiring (n-1)|X and n|X.")
print("  Min such X = n(n-1) (the PRONIC) = the deep well's defect speed. Verify + M-min.")
print("="*80)
print(f"{'n':>4} {'covering X in [1,n^2] are exactly mult of n(n-1)?':>48} {'min covering X':>14} {'=n(n-1)?':>9}")
for n in range(4,13):
    pron=n*(n-1)
    cov_X=[X for X in range(1, n*n+1) if is_covering(list(range(1,n-1))+[X], n)]
    exactly_mult = cov_X == [X for X in range(1,n*n+1) if X % pron == 0]
    print(f"{n:>4} {str(exactly_mult):>48} {str(min(cov_X) if cov_X else None):>14} "
          f"{str((min(cov_X)==pron) if cov_X else 'NA'):>9}")

print()
print("CLAIM 4b: shape-restricted covering-min: over {1..n-2, k*n(n-1)}, M is minimized")
print("  at k=1 (the pronic deep well) and equals n/Phi6(n).  (M via exact optimum.)")
print(f"{'n':>4} {'n/Phi6':>12} {'M at k=1..4':>44} {'argmin k':>8} {'shape-min=n/Phi6?':>16}")
for n in range(4,11):
    q=phi6(n); target=F(n,q); pron=n*(n-1)
    Ms=[]
    for k in range(1,5):
        S=list(range(1,n-1))+[k*pron]
        Ms.append(M_family(S, Qcap=2*q+2))
    kmin=1+min(range(4), key=lambda i:Ms[i])
    print(f"{n:>4} {str(target):>12} {str([str(m) for m in Ms]):>44} {kmin:>8} {str(min(Ms)==target):>16}")

print()
print("CLAIM 4c (sanity): FULL brute covering-min over ALL primitive covering families")
print("  (speeds <= n(n-1)) for n=4,5 -- confirms shape restriction loses nothing.")
for n in (4,5):
    q=phi6(n); target=F(n,q); DW=tuple(list(range(1,n-1))+[n*(n-1)])
    Bmax=n*(n-1); best=None;bestS=None;cov=0
    for S in itertools.combinations(range(1,Bmax+1), n-1):
        g=0
        for v in S: g=gcd(g,v)
        if g!=1: continue
        if not is_covering(S,n): continue
        cov+=1
        M=M_family(S, Qcap=2*q+2)
        if best is None or M<best: best,bestS=M,S
    print(f"  n={n}: covering-min={best}(~{float(best):.5f}); n/Phi6={target}; "
          f"min==n/Phi6? {best==target}; argmin={bestS} (DW={DW}); #cov={cov}")
print()
print("DONE")

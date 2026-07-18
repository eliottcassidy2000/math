# opus-2026-07-17-S368 -- HYP-7510: THE K-FOLD FOLDED IDENTITY.
#
# THE REFRAME.  Write h = indicator of (-lam, lam) on the circle.  Then
#     mu(n_i D_{a_i}) = int prod_i h(a_i t) dt
#                     = sum over {n in Z^k : sum n_i a_i = 0} of prod_i hhat(n_i)
# with hhat(0) = 2lam and hhat(n) = sin(2 pi n lam)/(pi n).
# So the k-fold measure is a sum over the RESONANCE LATTICE
#     Lam(a) = ker(Z^k -> Z, n |-> sum n_i a_i),   rank k-1.
#
# WHY k=2 IS SPECIAL (and S367's degradation explained): rank 1.  A rank-1
# lattice is cyclic, so the sum is a single Bernoulli series and closes in
# terms of the tent function fold_M(r) = r(M-r).  That IS THM-965.  For k>=3
# the lattice has rank >= 2 and no single Bernoulli closes it.
#
# THE DECOMPOSITION.  Split lattice vectors by their SUPPORT:
#     mu = sum over S subset [k] of (2lam)^(k-|S|) * delta(S),
# where delta(S) = sum over FULL-SUPPORT vectors of Lam(a_S) of prod hhat.
# delta(empty) = 1 and delta({i}) = 0 (since n a_i = 0 forces n = 0).
# So the first correction is the PAIRS -- delta({i,j}) is exactly the THM-965
# fold term -- and delta(S) for |S| >= 3 is genuinely new joint structure.
# This is a CUMULANT expansion: it separates true k-body alignment from
# products of lower-body alignments.  That is the "joint alignment" object
# S367 showed was needed.
from fractions import Fraction as F
from math import sin, pi, gcd
import itertools, random
LAM = F(1,14)

def hhat(n):
    if n == 0: return float(2*LAM)
    return sin(2*pi*n*float(LAM))/(pi*n)

# ---------- exact measure, for ground truth ----------
def teeth01(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def inter(u,v):
    out,i,j=[],0,0
    while i<len(u) and j<len(v):
        a,b=max(u[i][0],v[j][0]), min(u[i][1],v[j][1])
        if a<b: out.append((a,b))
        if u[i][1]<v[j][1]: i+=1
        else: j+=1
    return out
def exact(S):
    cur=teeth01(S[0])
    for x in S[1:]:
        cur=inter(cur,teeth01(x))
        if not cur: break
    return sum(b-a for a,b in cur)

# ---------- the resonance sum, by brute enumeration ----------
def resonance(A, N):
    """sum over {n in Z^k, |n_i|<=N, sum n_i a_i = 0} of prod hhat(n_i)."""
    k=len(A); tot=0.0
    for n in itertools.product(range(-N,N+1), repeat=k):
        if sum(ni*ai for ni,ai in zip(n,A))==0:
            p=1.0
            for ni in n: p*=hhat(ni)
            tot+=p
    return tot

def delta(A, N):
    """full-support resonance sum for the family A."""
    k=len(A); tot=0.0
    for n in itertools.product([x for x in range(-N,N+1) if x!=0], repeat=k):
        if sum(ni*ai for ni,ai in zip(n,A))==0:
            p=1.0
            for ni in n: p*=hhat(ni)
            tot+=p
    return tot

print("(1) THE FOURIER FORMULA vs EXACT MEASURE")
print("    k  family            exact       resonance(N)   N")
random.seed(368)
for A in [(2,3),(3,5),(5,7),(2,3,5),(3,4,5),(2,5,7)]:
    ex=float(exact(list(A))); k=len(A)
    N = 200 if k==2 else 60
    rs=resonance(A,N)
    print(f"    {k}  {str(A):16s} {ex:.6f}   {rs:.6f}    {N}")

print()
print("(2) THE SUPPORT DECOMPOSITION (the k-fold folded identity)")
print("    mu = sum_S (2lam)^(k-|S|) delta(S),  delta(empty)=1, delta(single)=0")
print("    k=3:  mu = (2lam)^3 + 2lam*[d12+d13+d23] + d123")
tl=float(2*LAM)
for A in [(2,3,5),(3,4,5),(2,5,7),(3,5,8)]:
    ex=float(exact(list(A)))
    d2=[delta((A[i],A[j]),120) for i,j in [(0,1),(0,2),(1,2)]]
    d3=delta(A,60)
    rec = tl**3 + tl*sum(d2) + d3
    print(f"    {str(A):12s} exact {ex:.6f}  rebuilt {rec:.6f}   "
          f"[indep {tl**3:.6f}, pairs {tl*sum(d2):+.6f}, TRIPLE {d3:+.6f}]")
print()
print("  => the genuinely 3-body term d123 is separated from the pair terms.")
print("     Pairs are THM-965 (closed form); d123 is the new object.")

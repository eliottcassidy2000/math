#!/usr/bin/env python3
"""
death-star-2026-07-20-S59s (HYP-8165) -- dissect the owner's triangle, build the
n*2^x+1 unification table, identify the rational series.  All exact.
"""
from fractions import Fraction as Fr
from itertools import product

T = [
 [1],
 [2,1],
 [3,3,1],
 [4,6,5,1],
 [5,10,14,9,1],
 [6,15,30,37,17,1],
 [7,21,55,101,99,33,1],
]
N = len(T)
def get(r,c):  # 1-indexed row, 1-indexed col
    if 1<=r<=N and 1<=c<=len(T[r-1]): return T[r-1][c-1]
    return None

print("=== COLUMNS (0-indexed j = col-1) vs power sums Sum_{k=1}^{m} k^j ===")
def powersum(m,j): return sum(k**j for k in range(1,m+1))
for j in range(6):
    col = [get(r, j+1) for r in range(j+1, N+1)]
    ps  = [powersum(m, j) for m in range(1, N-j)]
    match = col == ps
    print(f"  col j={j}: {col}")
    print(f"    Sum k^{j}: {ps}   MATCH={match}" + ("" if match else "  <-- DEVIATES (Moser-type)"))

print("\n=== SUBDIAGONAL (col=row-1) and last diagonal ===")
sub = [get(r, r-1) for r in range(2, N+1)]
print(f"  subdiagonal T(r,r-1): {sub}   vs 2^x+1:", [2**x+1 for x in range(1,N)])
print(f"  main diagonal T(r,r): {[get(r,r) for r in range(1,N+1)]} (all 1)")

print("\n=== ROW SUMS and ANTIDIAGONAL (shallow) SUMS ===")
rs = [sum(T[r-1]) for r in range(1,N+1)]
print(f"  row sums: {rs}")
# shallow diagonals (Fibonacci-from-Pascal style): sum T(n-k, k+1)
def shallow(d):
    s=0
    for k in range(0, d+1):
        v=get(d-k, k+1)
        if v is not None: s+=v
    return s
print(f"  shallow-diagonal sums: {[shallow(d) for d in range(1,2*N)]}")

print("\n=== RECURRENCE SEARCH: T(r,c) = a*T(r-1,c) + b*T(r-1,c-1) + e*T(r-2,c-1) + g*T(r-2,c-2)? ===")
# collect equations for interior entries, solve least integer relation by brute rational fit
rows_eq = []
targets = []
feats = []
for r in range(3, N+1):
    for c in range(2, len(T[r-1])):  # interior-ish
        row = [get(r-1,c), get(r-1,c-1), get(r-2,c-1) or 0, get(r-2,c-2) or 0,
               get(r-1,c-2) or 0, get(r-2,c) or 0]
        if get(r-1,c) is None: continue
        feats.append((r,c,row, get(r,c)))
# try to solve for coefficients using first several equations (6 unknowns)
import itertools
A=[]; b=[]
for (r,c,row,tgt) in feats:
    A.append([Fr(x) for x in row]); b.append(Fr(tgt))
# gaussian solve least-squares-ish: take first 6 independent
def solve(A,b):
    m=len(A); n=len(A[0])
    M=[A[i][:]+[b[i]] for i in range(m)]
    piv=[]; r=0
    for col in range(n):
        p=next((i for i in range(r,m) if M[i][col]!=0),None)
        if p is None: continue
        M[r],M[p]=M[p],M[r]
        pv=M[r][col]; M[r]=[x/pv for x in M[r]]
        for i in range(m):
            if i!=r and M[i][col]!=0:
                f=M[i][col]; M[i]=[a-f*bb for a,bb in zip(M[i],M[r])]
        piv.append(col); r+=1
    # check consistency
    for i in range(r,m):
        if M[i][-1]!=0: return None
    sol=[Fr(0)]*n
    for i,col in enumerate(piv): sol[col]=M[i][-1]
    return sol, piv
res=solve(A,b)
if res:
    sol,piv=res
    labels=["T(r-1,c)","T(r-1,c-1)","T(r-2,c-1)","T(r-2,c-2)","T(r-1,c-2)","T(r-2,c)"]
    print("  consistent linear recurrence found (coeffs):")
    for l,s0 in zip(labels,sol):
        if s0!=0: print(f"    {l}: {s0}")
    # verify on ALL entries
    ok=all(sum(f*x for f,x in zip(sol,row))==tgt for (r,c,row,tgt) in feats)
    print("  verified on all interior entries:", ok)
else:
    print("  no consistent 6-term linear recurrence (of this shape) -- try OEIS")

print("\n=== THE n*2^x+1 TABLE (x=0..5 rows, n=0..6 cols) ===")
def f(x,n):
    if x==0: return n   # owner's boundary convention
    return n*2**x+1
print("       n=0  n=1  n=2  n=3  n=4  n=5  n=6")
for x in range(0,6):
    row=[f(x,n) for n in range(0,7)]
    tag = "  (x=0: naturals)" if x==0 else ("  <- 2n+1 (observer/Redei)" if x==1 else "")
    print(f"  x={x}: {row}"+tag)
print("  column n=1 (down):", [f(x,1) for x in range(0,6)], " <- 2^x+1 (hypotenuse/Fermat)")
print("  ODD check: n*2^x+1 for x>=1 is ODD for all n (n*2^x even) -> the family is ODD numbers.")

print("\n=== RATIONAL SERIES ===")
H=[sum(Fr(1,k) for k in range(1,m+1)) for m in range(1,7)]
print("  harmonic H_n:", [str(h) for h in H], " (matches 1,3/2,11/6,25/12,137/60)")
# identify series 2 candidates: 1,5/2,29/6,103/12,... = sum (2^k-1)/k ? and 1,5/2,29/3,109/12,1079/60 (owner)
s2a=[]; acc=Fr(0)
for k in range(1,7):
    acc+=Fr(2**k-1,k); s2a.append(acc)
print("  sum (2^k-1)/k:", [str(x) for x in s2a])
s2b=[]; acc=Fr(0)
for k in range(1,7):
    acc+=Fr(2**(k-1),k); s2b.append(acc)
print("  sum 2^(k-1)/k:", [str(x) for x in s2b])
# owner series 2: 1,5/2,29/3,109/12,1079/60 -> compute what fits denom 1,2,3,12,60 (=k!/(k-1)?? no)
owner2=[Fr(1),Fr(5,2),Fr(29,3),Fr(109,12),Fr(1079,60)]
print("  owner series 2:", [str(x) for x in owner2], " differences from H_n:", [str(owner2[i]-H[i]) for i in range(5)])
# try sum of k/(k) ... try sum_{k=1}^n (2^k+... )? try H_n weighted 2^{k}
s2c=[]; acc=Fr(0)
for k in range(1,7):
    acc+=Fr(2**k,k*(k+1)) if False else Fr(0)
# try sum 1/k * 2... report ratios owner2[i]/H[i]
print("  owner2 / H_n:", [str(owner2[i]/H[i]) for i in range(5)])

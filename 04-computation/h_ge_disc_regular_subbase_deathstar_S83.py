#!/usr/bin/env python3
"""h_ge_disc_regular_subbase_deathstar_S83.py (HYP-8697 cont.)
Work the REGULAR SUB-BASE H(reg)>=n*disc(reg) and the PFAFFIAN INJECTION.
Chain (n>=7): H(reg) >= n!/2^{n-1} [Szele avg] >= n*(n+1)^{(n-1)/2}/2^{n-1} >= n*disc(reg).
Test the crux (i) H(reg)>=avg; prove-check (ii) (n-1)!>=(n+1)^{(n-1)/2}; (iii) AM-GM disc bound.
"""
from math import comb, factorial
from itertools import combinations
import numpy as np, random
random.seed(3)

def scores(A,n): return [sum(A[i]) for i in range(n)]
def is_regular(A,n): return all(s==(n-1)//2 for s in scores(A,n))
def ham(A,n):
    full=(1<<n)-1; out=[sum(1<<j for j in range(n) if A[i][j]) for i in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for m in range(1<<n):
        for v in range(n):
            c=dp[m][v]
            if c:
                for w in range(n):
                    if not(m>>w&1) and out[v]>>w&1: dp[m|1<<w][w]+=c
    return sum(dp[full][v] for v in range(n))
def disc(A,n):
    K=np.array(A,float)-np.array(A,float).T
    return round(abs(np.linalg.det(np.eye(n)+K)))/2**(n-1)
def paley(p):
    QR=set((x*x)%p for x in range(1,p))
    return [[1 if ((j-i)%p) in QR else 0 for j in range(p)] for i in range(p)]
def rotational(n,k):  # beats next k
    return [[1 if (0<((j-i)%n)<=k) else 0 for j in range(n)] for i in range(n)]
def rand_regular(n,tries=2000):
    # start rotational, random arc-swaps preserving regularity (reverse a 3-cycle)
    A=[row[:] for row in rotational(n,(n-1)//2)]
    for _ in range(tries):
        a,b,c=random.sample(range(n),3)
        if A[a][b] and A[b][c] and A[c][a]:
            A[a][b]=A[b][c]=A[c][a]=0; A[b][a]=A[c][b]=A[a][c]=1
    return A

print("="*66,"\n(A/B) REGULAR SUB-BASE H(reg) >= n*disc(reg); H vs Szele avg n!/2^{n-1}\n","="*66)
# n=5 exhaustive regular
def all_t(n):
    P=list(combinations(range(n),2))
    for b in range(1<<len(P)):
        A=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(P):
            if b>>k&1: A[i][j]=1
            else: A[j][i]=1
        yield A
for n in (3,5):
    regs=[A for A in all_t(n) if is_regular(A,n)]
    avg=factorial(n)/2**(n-1)
    data=[(ham(A,n),disc(A,n)) for A in regs]
    minr=min(H/(n*d) if d>0 else 9 for H,d in data)
    minH=min(H for H,d in data); maxd=max(d for H,d in data)
    print(f"  n={n}: {len(regs)} regular labeled; H in [{minH},{max(H for H,d in data)}], "
          f"disc max={maxd}; min H/(n*disc)={minr:.3f}; Szele avg={avg:.2f}; min H>=avg? {minH>=avg}")
# n=7,9 representatives + random regulars
for n in (7,9):
    reps={'Paley' if n%4==3 else 'rot': paley(n) if n%4==3 else rotational(n,(n-1)//2),
          'rot': rotational(n,(n-1)//2)}
    samples=[rand_regular(n) for _ in range(6)]
    avg=factorial(n)/2**(n-1)
    allH=[]; 
    for name,A in reps.items():
        assert is_regular(A,n), name
        H=ham(A,n); d=disc(A,n); allH.append(H)
        print(f"  n={n} {name}: H={H}, disc={d:.1f}, n*disc={n*d:.1f}, ratio={H/(n*d):.3f}, H>=avg({avg:.1f})? {H>=avg}")
    sh=[ham(A,n) for A in samples]
    print(f"  n={n} random regulars: H min={min(sh)} (avg={avg:.1f}); all >= avg? {min(sh)>=avg}")

print("\n(ii) ARITHMETIC (n-1)! >= (n+1)^{(n-1)/2} :")
for n in range(3,16,2):
    lhs=factorial(n-1); rhs=(n+1)**((n-1)//2) if (n-1)%2==0 else 0
    print(f"   n={n}: (n-1)!={lhs} vs (n+1)^{{(n-1)/2}}={rhs}  -> holds? {lhs>=rhs}")

print("\n(iii) AM-GM: disc(reg) <= (1+n)^{(n-1)/2}/2^{n-1} (max at doubly-regular):")
for n in (3,7,11):
    A=paley(n) if n%4==3 else rotational(n,(n-1)//2)
    if is_regular(A,n):
        bound=(1+n)**((n-1)//2)/2**(n-1)
        print(f"   n={n}: disc={disc(A,n):.2f} <= AM-GM bound {bound:.2f}? {disc(A,n)<=bound+1e-9}")

print("\n(PF) PFAFFIAN INJECTION slack 2^{n-1}H - det(I+K) >= 0 over strong tournaments:")
def is_strong(A,n):
    r=[[A[i][j] for j in range(n)] for i in range(n)]
    for i in range(n): r[i][i]=1
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if r[i][k] and r[k][j]: r[i][j]=1
    return all(r[i][j] for i in range(n) for j in range(n))
for n in (5,6):
    mn=1e18; cnt=0
    for A in all_t(n):
        if not is_strong(A,n): continue
        cnt+=1
        slack=2**(n-1)*ham(A,n)-round(abs(np.linalg.det(np.eye(n)+(np.array(A,float)-np.array(A,float).T))))
        mn=min(mn,slack)
    print(f"   n={n}: {cnt} strong; min(2^{{n-1}}H - det(I+K))={mn} (>=0 confirms 2^{{n-1}}H>=det=disc*2^{{n-1}})")

#!/usr/bin/env python3
"""paley_hampath_ratio_deathstar_S85.py
Is H(Paley_p) = Theta(average n!/2^{n-1})? Compute H(Paley)/avg across primes and
compare to typical random tournaments (is Paley above the random-typical H?).
Resolves the 'quasirandom => H~avg' -- really H=Theta(avg), constant>1, which still
gives H(Paley) >= n*disc since n*disc/avg -> 0 super-exp. Non-overlapping with the
archaeology agents (this is COMPUTE, they READ)."""
from math import factorial
import random
random.seed(11)
def ham(A,n):
    full=(1<<n)-1; out=[sum(1<<j for j in range(n) if A[i][j]) for i in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for m in range(1<<n):
        row=dp[m]
        for v in range(n):
            c=row[v]
            if c:
                ov=out[v]
                for w in range(n):
                    if not(m>>w&1) and ov>>w&1: dp[m|1<<w][w]+=c
    return sum(dp[full][v] for v in range(n))
def paley(p):
    QR=set((x*x)%p for x in range(1,p))
    return [[1 if((j-i)%p)in QR else 0 for j in range(p)]for i in range(p)]
def rand_t(n):
    A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if random.random()<.5: A[i][j]=1
            else: A[j][i]=1
    return A
print(f"{'p':>3} {'H(Paley)':>12} {'avg=p!/2^{p-1}':>16} {'H/avg':>7} {'rand-typical H/avg (median of 15)':>32}")
for p in (3,7,11,19):
    avg=factorial(p)/2**(p-1)
    Hp=ham(paley(p),p)
    rr=sorted(ham(rand_t(p),p)/avg for _ in range(15))
    med=rr[len(rr)//2]
    print(f"{p:>3} {Hp:>12} {avg:>16.1f} {Hp/avg:>7.3f} {med:>32.3f}")
print("""
READING: if H(Paley)/avg is stable across p, H(Paley)=Theta(avg) (constant>1);
Paley/avg vs random-median/avg shows whether Paley is ABOVE the random-typical.
Either way, since n*disc(Paley)/avg = n(n+1)^{(n-1)/2}/n! -> 0 super-exponentially,
H(Paley) >= c*avg >> n*disc for large p -- the crux holds from Theta(avg), which is
exactly the quasirandom Ham-path order.""")

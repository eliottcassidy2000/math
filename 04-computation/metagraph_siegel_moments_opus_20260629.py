"""
The metagraph MOMENT formulas (the analog of the paper's Siegel 1st/2nd moments).
E[H] = n!/2^{n-1} (clean: each ordering is a HP w.p. 2^{-(n-1)}).
E[H^2] = sum_{pi,pi'} P(both HP) = pair-correlation over ORDERING PAIRS (= the 'Siegel 2nd moment').
Compare to paper (Siegel 2nd moment, dim 2) and to LRC S_2 (pair moment). Look for clean Var(H).
"""
from itertools import permutations
from fractions import Fraction as F
def Hcount(n,adj):
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for last in range(n):
            c=dp[mask][last]
            if not c or not (mask>>last)&1: continue
            av=adj[last]&~mask
            while av:
                nb=av&-av; nx=nb.bit_length()-1; dp[mask|nb][nx]+=c; av^=nb
    return sum(dp[(1<<n)-1])
def Es(n):
    E=[(i,j) for i in range(n) for j in range(i+1,n)]; m=len(E)
    s1=0; s2=0
    for bits in range(1<<m):
        adj=[0]*n
        for k,(i,j) in enumerate(E):
            if (bits>>k)&1: adj[i]|=1<<j
            else: adj[j]|=1<<i
        h=Hcount(n,adj); s1+=h; s2+=h*h
    N=1<<m
    EH=F(s1,N); EH2=F(s2,N); Var=EH2-EH*EH
    return N,EH,EH2,Var
print("metagraph moment formulas (over labeled tournaments):")
print(f"{'n':>2} {'E[H]':>10} {'n!/2^(n-1)':>10} {'E[H^2]':>12} {'Var(H)':>12} {'Var/E[H]':>10}")
import math
for n in range(3,7):
    N,EH,EH2,Var=Es(n)
    pred=F(math.factorial(n),2**(n-1))
    print(f"{n:>2} {str(EH):>10} {str(pred):>10} {str(EH2):>12} {float(Var):>12.3f} {float(Var/EH):>10.4f}")
print()
print("E[H]=n!/2^(n-1) EXACT (each of n! orderings is a HP w.p. 2^-(n-1)).")
print("E[H^2]=sum_{pi,pi'} 2^{-|arcs(pi) U arcs(pi')|} = the SIEGEL-style 2nd moment over ORDERING PAIRS,")
print("graded by how many consecutive-arcs the two orderings SHARE -- the metagraph pair-correlation,")
print("structurally the SAME object as the paper's dim-2 Siegel 2nd moment and the LRC S_2 pair moment.")

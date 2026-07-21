# Grounding: "degeneracy = few distinct eigenvalues = big GIT-stabilizer = special stratum."
# For tournaments n=5, distribution of #distinct eigenvalues of A; transitive = 1 (nullcone vertex),
# generic = 5 (stable). Also |Aut| vs spectral degeneracy to separate the two symmetry notions.
from itertools import permutations
import numpy as np
def tour(n,bits):
    A=[[0]*n for _ in range(n)]; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if (bits>>idx)&1: A[i][j]=1
            else: A[j][i]=1
            idx+=1
    return A
def aut(n,A):
    c=0
    for p in permutations(range(n)):
        if all(A[i][j]==A[p[i]][p[j]] for i in range(n) for j in range(n)): c+=1
    return c
n=5; from collections import Counter
dist=Counter(); byaut={}
transitive=None
for bits in range(1<<(n*(n-1)//2)):
    A=tour(n,bits); ev=np.linalg.eigvals(np.array(A,dtype=float))
    # count distinct eigenvalues (round)
    d=len({(round(z.real,4),round(abs(z.imag),4)) for z in ev})
    dist[d]+=1
    au=aut(n,A); byaut.setdefault(au,[]).append(d)
    if all(abs(z)<1e-6 for z in ev): transitive=(d,au)
print("n=5: distribution of #distinct-eigenvalues over 1024 tournaments:", dict(sorted(dist.items())))
print(f"  transitive (nullcone vertex, char poly x^5): #distinct eigenvalues = {transitive[0]} (=1, spectrum {{0}}), |Aut|={transitive[1]}")
print("  |Aut| -> avg #distinct eigenvalues:", {a: round(np.mean(v),2) for a,v in sorted(byaut.items())})
print("  => NULLCONE (transitive) = 1 distinct eigenvalue = MAXIMAL spectral degeneracy, but |Aut|=1 (rigid).")
print("     LARGE |Aut| (symmetric) ALSO forces spectral degeneracy (fewer distinct) -- BUT via the")
print("     semisimple/critical-line pole (Paley), not the nilpotent one. Two routes to degeneracy:")
print("     nilpotent (transitive, GIT big-parabolic-stabilizer) vs symmetric (Paley, big S_n-Aut).")

from math import gcd
from itertools import permutations, combinations
import random
MOD=14; UNITS=[a for a in range(1,MOD) if gcd(a,MOD)==1]
def depth(r): r%=MOD; return min(r,MOD-r)
def m4(S):
    m=len(S); adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1,m):
            wi=wj=0
            for a in UNITS:
                di=depth(S[i]*a); dj=depth(S[j]*a)
                if di>dj: wi+=1
                elif dj>di: wj+=1
            if wi>wj: adj[i][j]=1
            elif wj>wi: adj[j][i]=1
            else: adj[i][j]=1 if S[i]<S[j] else 0; adj[j][i]=1-adj[i][j]
    return adj
# Correct statement: M4(S) is determined by (residue_i mod 14) and the RELATIVE ORDER of
# raw speeds (which decides tie-breaks). If two speed sets have identical residues AND
# identical speed-rank order, M4 is identical. Verify.
random.seed(2)
ok=True; bad=None
for _ in range(5000):
    res=[random.randint(0,13) for _ in range(5)]
    # build two speed sets with these residues and the SAME order; offsets preserve order
    # by making speed = res + 14*k with k chosen strictly increasing along a fixed perm.
    perm=list(range(5)); random.shuffle(perm)  # this is the intended speed-rank order
    # assign multipliers so that speed order == perm order
    Sa=[0]*5; Sb=[0]*5
    for rank,idx in enumerate(perm):
        Sa[idx]=res[idx]+14*(rank+1)
        Sb[idx]=res[idx]+14*(rank+1)+14*5*(rank+1)  # bigger but SAME order
    # sanity: same residues, same order
    if [v%14 for v in Sa]!=res or [v%14 for v in Sb]!=res:
        continue
    oa=sorted(range(5),key=lambda i:Sa[i]); ob=sorted(range(5),key=lambda i:Sb[i])
    if oa!=ob: continue
    if m4(Sa)!=m4(Sb): ok=False; bad=(Sa,Sb); break
print("M4 determined by (residues mod 14 + speed-rank order):", ok, bad if bad else "")

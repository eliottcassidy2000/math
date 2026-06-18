# Decisive: can M4 (residue-orbit, mod 14) EVER produce a non-realized class?
# M4 depends ONLY on residues mod 14 + speed tie-break, so a large-speed search must
# match the residue-exhaustive forbidden set exactly. Confirm regular never appears.
from math import gcd
from itertools import combinations, permutations
MOD=14; U=[a for a in range(1,MOD) if gcd(a,MOD)==1]
def m4(S):
    m=len(S); adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1,m):
            wi=wj=0
            for a in U:
                ri=(S[i]*a)%MOD; rj=(S[j]*a)%MOD
                di=min(ri,MOD-ri); dj=min(rj,MOD-rj)
                if di>dj: wi+=1
                elif dj>di: wj+=1
            if wi>wj: adj[i][j]=1
            elif wj>wi: adj[j][i]=1
            else: adj[i][j]=1 if S[i]<S[j] else 0; adj[j][i]=1-adj[i][j]
    return adj
def score(adj,m): return tuple(sorted(sum(adj[i][j] for j in range(m) if j!=i) for i in range(m)))
found=0; total=0
seen_scores=set()
for S in combinations(range(1,41),5):
    g=0
    for v in S: g=gcd(g,v)
    if g!=1: continue
    total+=1
    sc=score(m4(S),5); seen_scores.add(sc)
    if sc==(2,2,2,2,2):
        found+=1
        if found<=3: print("REGULAR via M4:",S,flush=True)
print(f"M4 over {total} primitive 5-sets in 1..40: regular found = {found}",flush=True)
print(f"score sequences seen: {sorted(seen_scores)}",flush=True)

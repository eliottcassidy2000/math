import itertools, numpy as np
from collections import defaultdict
print("(A) WHAT DOES |Pf| = 3 MEAN FOR A 4-VERTEX TOURNAMENT?")
P4=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
tab=defaultdict(lambda: defaultdict(int))
for m in range(64):
    S=np.zeros((4,4),dtype=int)
    for b,(i,j) in enumerate(P4):
        v=1 if (m>>b)&1 else -1
        S[i,j]=v; S[j,i]=-v
    Pf = S[0,1]*S[2,3] - S[0,2]*S[1,3] + S[0,3]*S[1,2]
    cyc = sum(1 for (a,b,c) in itertools.combinations(range(4),3)
              if not(S[a,b]==S[a,c]==1 or S[b,a]==S[b,c]==1 or S[c,a]==S[c,b]==1))
    sc = tuple(sorted(sum(1 for j in range(4) if S[i,j]==1) for i in range(4)))
    tab[abs(Pf)][(cyc,sc)] += 1
for pf in sorted(tab):
    print(f"  |Pf| = {pf}:")
    for (cyc,sc),ct in sorted(tab[pf].items()):
        print(f"      cyclic-triples={cyc}  scores={sc}  count={ct}")
print()
print("(B) SO k COUNTS 4-SUBSETS OF A PARTICULAR TYPE. Confirm on 7-tournaments,")
print("    and check the evenness / the missing k=2.")
n=7
pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
base={(i,i+1) for i in range(n-1)}
free=[p for p in pairs if p not in base]
quads=list(itertools.combinations(range(n),4))
def kcount(S):
    k=0
    for (a,b,c,d) in quads:
        Pf = S[a,b]*S[c,d] - S[a,c]*S[b,d] + S[a,d]*S[b,c]
        if abs(Pf)==3: k+=1
    return k
seen=defaultdict(int)
for mask in range(1<<len(free)):
    S=np.zeros((n,n),dtype=np.int64)
    for (i,j) in base: S[i,j]=1; S[j,i]=-1
    for b,(i,j) in enumerate(free):
        v=1 if (mask>>b)&1 else -1
        S[i,j]=v; S[j,i]=-v
    seen[kcount(S)]+=1
print(f"    k values observed: {sorted(seen)}")
print(f"    k=2 present? {2 in seen}   all even? {all(v%2==0 for v in seen)}")
for k in sorted(seen): print(f"      k={k:3d}  ->  c3 = {35+8*k:4d}   ({seen[k]} switching classes)")

"""
Fourth pass: pin down the level at which P lives.
P is NOT iso-invariant (verify3). So 'P refines c3' is a statement about the tiling cube Q_F
with the FIXED staircase labeling, NOT about iso classes / the metagraph G_n.

Check: within a single ISO CLASS, does P take multiple values? (If yes, P is not an iso invariant
and cannot label nodes of G_n.) Group all 2^F tilings by iso class (canonical form via min over
relabelings) and report P-spread per class. Also confirm c3 IS constant per iso class (sanity).
"""
from itertools import combinations, product, permutations
from collections import defaultdict

def tiles(n): return [(a,b) for a in range(3,n+1) for b in range(1,a-1)]
def adj(n,bits,T):
    A=[[0]*(n+1) for _ in range(n+1)]
    for k in range(n,1,-1): A[k][k-1]=1
    for (a,b),bit in zip(T,bits):
        if bit==0: A[a][b]=1
        else: A[b][a]=1
    return A
def c3(A,n):
    t=0
    for i in range(1,n+1):
        for j in range(i+1,n+1):
            for k in range(j+1,n+1):
                if (A[i][j]+A[i][k],A[j][i]+A[j][k],A[k][i]+A[k][j])==(1,1,1):t+=1
    return t
def Pinv(A,n):
    P=0
    for w,x,y,z in combinations(range(1,n+1),4):
        if A[w][y]==A[x][z]: P+=1
    return P
def canon(A,n):
    best=None
    for perm in permutations(range(1,n+1)):
        bits=[]
        for i in range(1,n+1):
            for j in range(i+1,n+1):
                bits.append(A[perm[i-1]][perm[j-1]])
        t=tuple(bits)
        if best is None or t<best: best=t
    return best

out=[]
def log(*a):
    s=" ".join(str(x) for x in a); print(s); out.append(s)

def per_isoclass(n):
    T=tiles(n); F=len(T)
    cls_P=defaultdict(set); cls_c3=defaultdict(set)
    for bits in product([0,1],repeat=F):
        A=adj(n,list(bits),T)
        c=canon(A,n)
        cls_P[c].add(Pinv(A,n)); cls_c3[c].add(c3(A,n))
    nclasses=len(cls_P)
    P_constant_per_class = all(len(v)==1 for v in cls_P.values())
    c3_constant_per_class = all(len(v)==1 for v in cls_c3.values())
    multivalued=sum(1 for v in cls_P.values() if len(v)>1)
    return nclasses,P_constant_per_class,c3_constant_per_class,multivalued

log("=== P and c3 constancy within iso classes ===")
for n in (4,5,6):
    nc,pc,cc,mv=per_isoclass(n)
    log(f"n={n}: #iso classes={nc} ; c3 constant per class={cc} ; "
        f"P constant per class={pc} ; #classes with multi-valued P={mv}")

with open("05-knowledge/results/conn_verify4_crossingroad_kps-Sx-wf.out","w") as f:
    f.write("\n".join(out))

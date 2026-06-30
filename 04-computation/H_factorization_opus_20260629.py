"""
Practical recursive compression payoff: H (and more) FACTORIZES over the condensation.
  H(X (+) Y) = H(X) * H(Y)   (X dominates Y: once a Ham path enters Y it can't return).
So H(tournament) = PRODUCT of H over its strongly-connected primes. The primes are the
'H-atoms'; computing H reduces to the irreducibles. Verify, and list the prime H-atoms.
"""
from itertools import permutations, combinations

def pairs(n): return list(combinations(range(n),2))
def adj(n,bits):
    pr=pairs(n); A=[[False]*n for _ in range(n)]
    for k,(i,j) in enumerate(pr):
        if bits[k]: A[i][j]=True
        else: A[j][i]=True
    return A
def H_count(n,bits):
    A=adj(n,bits)
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for last in range(n):
            c=dp[mask][last]
            if not c or not (mask>>last)&1: continue
            for nxt in range(n):
                if (mask>>nxt)&1 or not A[last][nxt]: continue
                dp[mask|1<<nxt][nxt]+=c
    return sum(dp[(1<<n)-1])
def canon(n,bits):
    pr=pairs(n); idx={p:k for k,p in enumerate(pr)}; best=None
    for perm in permutations(range(n)):
        out=[bits[idx[(perm[a],perm[b])]] if perm[a]<perm[b] else 1-bits[idx[(perm[b],perm[a])]] for (a,b) in pr]
        t=tuple(out)
        if best is None or t<best: best=t
    return best
def sccs(n,bits):
    A=adj(n,bits); reach=[[A[i][j] or i==j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    if reach[k][j]: reach[i][j]=True
    comp=[-1]*n; c=0
    for i in range(n):
        if comp[i]<0:
            comp[i]=c
            for j in range(i+1,n):
                if comp[j]<0 and reach[i][j] and reach[j][i]: comp[j]=c
            c+=1
    groups={}
    for i in range(n): groups.setdefault(comp[i],[]).append(i)
    return list(groups.values())
def subcanon(n,bits,grp):
    pr=pairs(n); idx={p:k for k,p in enumerate(pr)}; m=len(grp)
    sub=[grp[t] for t in range(m)]
    sb=[]
    for (a,b) in pairs(m):
        va,vb=sub[a],sub[b]
        sb.append(bits[idx[(va,vb)]] if va<vb else 1-bits[idx[(vb,va)]])
    return m,canon(m,tuple(sb))

def all_iso(n):
    pr=pairs(n); m=len(pr); seen={}
    for x in range(1<<m):
        bits=tuple((x>>k)&1 for k in range(m))
        c=canon(n,bits)
        if c not in seen: seen[c]=bits
    return seen

# H-atoms: H of each SC prime
print("H-atoms (H of strongly-connected primes) per size:")
Hatom={}
for n in range(1,7):
    iso=all_iso(n)
    for c,bits in iso.items():
        if len(sccs(n,bits))==1:
            Hatom[(n,c)]=H_count(n,bits)
    vals=sorted(set(v for (m,_),v in Hatom.items() if m==n))
    print(f"  size {n}: SC primes have H in {vals}")

# verify H(T)=prod of H over condensation primes, for n=5,6
print("\nverify H multiplicative over condensation (n=5,6):")
for n in [5,6]:
    iso=all_iso(n); bad=0
    for c,bits in iso.items():
        prod=1
        for grp in sccs(n,bits):
            m,cc=subcanon(n,bits,grp)
            prod*= Hatom[(m,cc)]
        if prod!=H_count(n,bits): bad+=1
    print(f"  n={n}: H == product over primes for ALL {len(iso)} classes: {bad==0}")

print("\n=> H(tournament) = PRODUCT of H over its strongly-connected primes (the recursive")
print("   compression: factor into SC atoms, look up each atom's H). The atoms are the")
print("   incompressible cores; everything reducible multiplies cleanly.")

"""
NEW small useful object: the METAGRAPH ZETA  zeta_G(s) = sum_{iso-classes C} H(C)^{-s}.
Also the orbit identity sum_C H/|Aut| = 2^{C(n-1,2)} (tilings), and the paper's Gamma_0(N)
index analog. Computed by brute-force canonicalization, n=4,5,6.
"""
from itertools import permutations
def Hcount(n,adj):  # adj[i] = bitmask of out-neighbors
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
def edges(n):
    e=[]
    for i in range(n):
        for j in range(i+1,n): e.append((i,j))
    return e
def metagraph(n):
    E=edges(n); m=len(E); perms=list(permutations(range(n)))
    seen={}; 
    for bits in range(1<<m):
        # adjacency: edge (i,j): bit set => i->j else j->i
        adj=[0]*n
        for k,(i,j) in enumerate(E):
            if (bits>>k)&1: adj[i]|=1<<j
            else: adj[j]|=1<<i
        # canonical form: min over perms of the out-adjacency tuple
        canon=None; autc=0
        base=tuple(adj)
        for p in perms:
            padj=[0]*n
            for v in range(n):
                a=adj[v]; pv=p[v]
                while a:
                    b=a&-a; w=b.bit_length()-1; padj[pv]|=1<<p[w]; a^=b
            t=tuple(padj)
            if canon is None or t<canon: canon=t
        if canon in seen: continue
        # count Aut for the canonical rep
        cadj=list(canon); aut=0
        for p in perms:
            padj=[0]*n; ok=True
            for v in range(n):
                a=cadj[v]; pv=p[v]
                while a:
                    b=a&-a; w=b.bit_length()-1; padj[pv]|=1<<p[w]; a^=b
            if tuple(padj)==canon: aut+=1
        seen[canon]=(Hcount(n,cadj),aut)
    return seen
for n in [4,5,6]:
    cls=metagraph(n)
    Hs=sorted(h for h,a in cls.values())
    V=len(cls); 
    sumH=sum(h for h,a in cls.values())
    sH_over_aut=sum(h/a for h,a in cls.values())
    z1=sum(h**-1 for h,a in cls.values()); z2=sum(h**-2 for h,a in cls.values())
    m=n*(n-1)//2; mm=(n-1)*(n-2)//2
    print(f"n={n}: V={V} classes, H-spectrum={Hs}")
    print(f"      sum_C H/|Aut| = {sH_over_aut:.0f}  (should = 2^C(n-1,2) = {2**mm})  [orbit identity check]")
    print(f"      zeta_G(1)=sum 1/H = {z1:.5f}   zeta_G(2)=sum 1/H^2 = {z2:.5f}   zeta_G(-1)=sum H = {sumH}")

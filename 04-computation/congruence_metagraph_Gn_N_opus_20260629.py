"""
The CONGRUENCE METAGRAPH G_n(N): refine iso-classes by H mod N; track class-count and MASS (H/|Aut|).
Gamma_0(N) analog. Check: does the mass equidistribute mod N? Does H mod 2^k = the OCF 2-adic digits
(THM-466)? Which residues are forbidden (the H=7 / apex face)?
"""
from itertools import permutations
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
def classes(n):
    E=[(i,j) for i in range(n) for j in range(i+1,n)]; m=len(E); perms=list(permutations(range(n)))
    seen={}
    for bits in range(1<<m):
        adj=[0]*n
        for k,(i,j) in enumerate(E):
            if (bits>>k)&1: adj[i]|=1<<j
            else: adj[j]|=1<<i
        canon=None
        for p in perms:
            padj=[0]*n
            for v in range(n):
                a=adj[v]; pv=p[v]
                while a:
                    b=a&-a; w=b.bit_length()-1; padj[pv]|=1<<p[w]; a^=b
            t=tuple(padj)
            if canon is None or t<canon: canon=t
        if canon in seen: continue
        cadj=list(canon); aut=sum(1 for p in perms if tuple(
            (lambda pa: [pa[v] for v in range(n)])([ (lambda: 0)() for _ in range(n)]) )==canon) if False else 0
        # count aut properly
        aut=0
        for p in perms:
            padj=[0]*n
            for v in range(n):
                a=cadj[v]; pv=p[v]
                while a:
                    b=a&-a; w=b.bit_length()-1; padj[pv]|=1<<p[w]; a^=b
            if tuple(padj)==canon: aut+=1
        seen[canon]=(Hcount(n,cadj),aut)
    return list(seen.values())
from fractions import Fraction as F
for n in [5,6]:
    cl=classes(n)
    print(f"\n=== n={n}: {len(cl)} classes, total mass sum H/|Aut| = {sum(F(h,a) for h,a in cl)} (=2^{(n-1)*(n-2)//2}) ===")
    for N in [3,4,7,8]:
        cnt={}; mass={}
        for h,a in cl:
            r=h%N; cnt[r]=cnt.get(r,0)+1; mass[r]=mass.get(r,F(0))+F(h,a)
        tot=sum(mass.values())
        line=" ".join(f"r{r}:#{cnt.get(r,0)},m={float(mass.get(r,0)):.0f}" for r in range(N))
        print(f"  mod {N}: {line}")
        # equidistribution of mass?
        eq=all(abs(float(mass.get(r,0))/float(tot)-1/N)<0.02 for r in range(N) if True) if N in(3,7) else None

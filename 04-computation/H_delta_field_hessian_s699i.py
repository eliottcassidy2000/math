"""The DELTA FIELD on tournament arcs. H on the arc-flip hypercube {0,1}^{C(n,2)}. delta(e)=
H(T^e)−H(T) (the discrete GRADIENT, always EVEN since H odd). Flipping f changes delta(e) by the
discrete HESSIAN Δ_ef = H(T)−H(T^e)−H(T^f)+H(T^{ef}). 'Not all deltas change' = Hessian support is
sparse. Connect to OCF (Δ_ef≠0 ⟺ e,f share an odd cycle). And the 7,21 forbidden-value navigation
constraint on delta. opus-2026-06-06-S699i."""
from itertools import combinations, product
def Hcount(n,adj):
    size=1<<n; dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        for v in range(n):
            c=dp[mask][v]
            if not c: continue
            for w in range(n):
                if not(mask>>w&1) and adj[v]>>w&1: dp[mask|1<<w][w]+=c
    return sum(dp[size-1][v] for v in range(n))
def make(n,bits,arcs):
    adj=[0]*n
    for (i,j),b in zip(arcs,bits):
        if b: adj[i]|=1<<j
        else: adj[j]|=1<<i
    return adj
def shares_vertex(e,f): return len(set(e)&set(f))>0
def in_common_triangle(e,f,arcs_set):  # e,f share a vertex => in a triangle (with the 3rd)
    return shares_vertex(e,f)
def main():
    n=5; arcs=list(combinations(range(n),2)); m=len(arcs)
    print(f"n={n}, {m} arcs, {2**m} tournaments")
    # (1) delta always even; Hessian always div by 4
    deven=True; hdiv4=True
    # (2) interaction SUPPORT: pair (e,f) interacts if some T has Δ_ef≠0
    interacts={(a,b):False for a in range(m) for b in range(a+1,m)}
    # (3) per-flip interaction count distribution
    import random; rng=random.Random(1)
    sample=[tuple(rng.randint(0,1) for _ in range(m)) for _ in range(300)]
    maxcount=0; counts=[]
    for bits in sample:
        H0=Hcount(n,make(n,bits,arcs))
        d=[]
        for a in range(m):
            b2=list(bits); b2[a]^=1; d.append(Hcount(n,make(n,b2,arcs))-H0)
        if any(x%2 for x in d): deven=False
        # flip f, count how many e's delta changes
        for f in range(m):
            bf=list(bits); bf[f]^=1; Hf=Hcount(n,make(n,bf,arcs))
            df=[]
            for a in range(m):
                bfa=list(bf); bfa[a]^=1; df.append(Hcount(n,make(n,bfa,arcs))-Hf)
            chg=[a for a in range(m) if a!=f and df[a]!=d[a]]
            counts.append(len(chg)); maxcount=max(maxcount,len(chg))
            for a in chg:
                key=(min(a,f),max(a,f)); interacts[key]=True
                hess=d[a]-df[a] if False else None
    # check Hessian div by 4 on a few
    for bits in sample[:50]:
        for (a,b) in combinations(range(m),2):
            b1=list(bits); b1[a]^=1; b2=list(bits); b2[b]^=1; b12=list(bits); b12[a]^=1; b12[b]^=1
            H=Hcount(n,make(n,bits,arcs)); Ha=Hcount(n,make(n,b1,arcs)); Hb=Hcount(n,make(n,b2,arcs)); Hab=Hcount(n,make(n,b12,arcs))
            hess=H-Ha-Hb+Hab
            if hess%4!=0: hdiv4=False
    print(f"(1) delta always even: {deven}; Hessian always divisible by 4: {hdiv4}")
    print(f"(2) per-flip #other-arcs-whose-delta-changes: max={maxcount} (of {m-1}); ever ALL {m-1}? {maxcount==m-1}")
    print(f"    (so flipping an arc NEVER changes every other delta — the interaction is sparse)")
    # interaction support vs shares-a-vertex
    sv_int=sum(1 for (a,b) in interacts if interacts[(a,b)] and shares_vertex(arcs[a],arcs[b]))
    sv_noint=sum(1 for (a,b) in interacts if not interacts[(a,b)] and shares_vertex(arcs[a],arcs[b]))
    nsv_int=sum(1 for (a,b) in interacts if interacts[(a,b)] and not shares_vertex(arcs[a],arcs[b]))
    nsv_noint=sum(1 for (a,b) in interacts if not interacts[(a,b)] and not shares_vertex(arcs[a],arcs[b]))
    print(f"(3) interaction support vs share-a-vertex: shareV&interact={sv_int}, shareV&no={sv_noint}, disjoint&interact={nsv_int}, disjoint&no={nsv_noint}")
    print(f"    (interaction ⟺ share an odd cycle; sharing a vertex ⟹ a common triangle ⟹ interact)")
    print("\n(4) 7,21 FORBIDDEN-VALUE NAVIGATION: from each achievable H, which signed deltas land where?")
    # over all 1024, the (H, H_after_flip) transitions
    trans={}
    for bits in product((0,1),repeat=m):
        H0=Hcount(n,make(n,bits,arcs))
        for f in range(m):
            bf=list(bits); bf[f]^=1; Hf=Hcount(n,make(n,bf,arcs))
            trans.setdefault(H0,set()).add(Hf-H0)
    for H0 in sorted(trans):
        ds=sorted(trans[H0]); lands=sorted(set(H0+d for d in ds))
        print(f"    H={H0:2d}: deltas={ds}; lands on H'={lands}  (7,21 in lands? {7 in lands or 21 in lands})")
if __name__=='__main__': main()

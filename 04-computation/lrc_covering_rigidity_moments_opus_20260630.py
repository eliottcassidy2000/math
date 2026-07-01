"""
COVERING RIGIDITY. m(t) = #{v in S: ||vt||<=1/n} (multiplicity). Covering (tight) <=> m>=1 a.e.
Conserved: int m = (n-1)*2/n = 2-2/n, so overlap int(m-1) = 1-2/n. Explore moments + the m-histogram + the
LOWNESS (how the lonely set grows as we remove elements) -- hunting the invariant that forces the AP-skeleton.
"""
def frac(x): x=x%1.0; return min(x,1-x)
def mult(S,t,n): return sum(1 for v in S if frac(v*t)<=1.0/n+1e-12)
def profile(S,n,G):
    ms=[mult(S,i/G,n) for i in range(G)]
    lonely=sum(1 for x in ms if x==0)/G
    overlap=sum(x-1 for x in ms if x>=1)/G  # int(m-1) over covered part
    m2=sum(x*x for x in ms)/G
    from collections import Counter
    hist=Counter(ms)
    return lonely,overlap,m2,hist
for n in [8,14]:
    G=200*n
    AP=list(range(1,n))
    lon,ov,m2,hist=profile(AP,n,G)
    print(f"n={n} AP: lonely={lon:.4f} (should be ~0), overlap int(m-1)={ov:.4f} (=1-2/n={1-2/n:.4f}), int m^2={m2:.3f}")
    print(f"   m-histogram (fraction of [0,1) with multiplicity m): "+", ".join(f"m={k}:{v/G:.3f}" for k,v in sorted(hist.items())))
print()
print("LOWNESS: remove R elements from AP (no patch) -> lonely-set measure (the hole):")
for n in [14]:
    G=200*n; AP=list(range(1,n))
    import itertools
    print(f"  n={n}: single removals -> lonely measure per k:")
    singles=[(k,profile([v for v in AP if v!=k],n,G)[0]) for k in AP]
    for k,l in singles: print(f"     remove {k}: lonely={l:.4f}", end="")
    print()
    # which single k gives SMALLEST hole (most patchable-looking)?
    kbest=min(singles,key=lambda x:x[1]); print(f"     smallest hole at k={kbest[0]} (lonely {kbest[1]:.4f}) -- the patchable one (=12, GW)? {kbest[0]==12}")
    # double removals: does the hole roughly add up (disjoint) or overlap?
    print(f"  double removals (k1,k2) -> lonely; compare to sum of singles (disjoint if ~equal):")
    for (k1,k2) in [(12,2),(12,6),(2,3),(11,13),(6,7)]:
        l12=profile([v for v in AP if v not in (k1,k2)],n,G)[0]
        s1=dict(singles)[k1] if k1 in AP else 0; s2=dict(singles)[k2] if k2 in AP else 0
        print(f"     remove {{{k1},{k2}}}: lonely={l12:.4f}; sum singles={s1+s2:.4f}; disjoint~{abs(l12-(s1+s2))<0.005}")

#!/usr/bin/env python3
"""
death-star-2026-07-20-S61 (HYP-8265) -- the base-path-INDEPENDENT star invariant:
compute cap Gamma(P) over ALL Hamiltonian base paths P, in the common edge space
F_2^{C(n,2)} of K_n. kp-S128c109 flagged this as the candidate that might DESCEND
to iso classes where the single-path star group Gamma(P) does not.

Edge space: bit per edge (i,j), i<j. delta(v) = star of v (edges incident to v).
Star move for base path P: S_v^P = delta(v) XOR (path-edges incident to v).
Gamma(P) = span{S_v^P}. Compute cap_P Gamma(P) for n=4,5,6; interpret.
"""
from itertools import permutations, combinations

def run(n):
    edges=list(combinations(range(n),2)); E={e:k for k,e in enumerate(edges)}; m=len(edges)
    def ebit(i,j): return 1<<E[(min(i,j),max(i,j))]
    def delta(v): return sum(ebit(v,u) for u in range(n) if u!=v)
    # Hamiltonian paths (as vertex sequences), dedup reverse
    paths=set()
    for perm in permutations(range(n)):
        if perm[0]<perm[-1]: paths.add(perm)
    def gamma_of(path):
        pe=[(path[i],path[i+1]) for i in range(n-1)]
        pv={v:0 for v in range(n)}
        for a,b in pe: pv[a]^=ebit(a,b); pv[b]^=ebit(a,b)
        gens=[delta(v)^pv[v] for v in range(n)]
        return rref(gens,m)
    def rref(vecs,m):
        basis=[]
        for v in vecs:
            for b in basis:
                v=min(v,v^b)
            if v: basis.append(v); basis.sort(reverse=True)
        return basis
    def dim(basis): return len(basis)
    def intersect(A,B,m):
        # Zassenhaus over F_2: stack [a|a] and [b|0], rref, read bottom (0|x) rows
        rows=[(a<<m)|a for a in A]+[(b<<m)|0 for b in B]
        # rref on 2m-bit rows
        basis=[]
        for v in rows:
            for b in basis: v=min(v,v^b)
            if v: basis.append(v); basis.sort(reverse=True)
        low=(1<<m)-1
        inter=[]
        for b in basis:
            if (b>>m)==0 and (b&low):  # top half zero => intersection element
                inter.append(b&low)
        return rref(inter,m)
    Gs=[gamma_of(p) for p in paths]
    inter=Gs[0]
    for G in Gs[1:]:
        inter=intersect(inter,G,m)
        if not inter: break
    print(f"  n={n}: {len(paths)} Ham base paths; each Gamma(P) dim={dim(Gs[0])}; "
          f"cap Gamma dim = {dim(inter)}")
    if inter:
        print(f"    basis (edge-bitmasks): {inter}")
        # what edges does each basis vector flip?
        for b in inter:
            es=[edges[k] for k in range(m) if (b>>k)&1]
            print(f"      move flips edges {es}")
    else:
        print("    cap Gamma = {0}: NO nonzero base-path-independent star move.")
    return inter, edges

print("=== cap Gamma(P) over all Hamiltonian base paths (kp-S128c109 open computation) ===")
for n in [4,5,6]:
    run(n)
print("\n=== interpretation ===")
print("  If cap Gamma = 0: the base-path-independent star subgroup is trivial => confirms")
print("  kp's diagnosis that cut-space/star constructions are IRREDUCIBLY base-path-relative")
print("  (the base path is exactly what the iso quotient forgets). The right base-path-free")
print("  object is then the FULL cut space = SEIDEL SWITCHING (delta(v) with path edges),")
print("  whose 2^{C(n-1,2)} switching classes ARE the tiling hypercube (base-path picks the rep).")

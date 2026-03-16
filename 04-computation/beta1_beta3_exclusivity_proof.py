#!/usr/bin/env python3
"""
PROVING β₁ · β₃ = 0 (Mutual Exclusivity)
==========================================
Strategy: When β₁(T) = 1, show H₃(T,T\v) = 0 for some vertex v.

Key insight from THM-108: when β₁(T) = 1, ALL vertices are "good"
(β₁(T\v) ≤ 1 = β₁(T)). So ANY vertex works for the LES approach.

LES of (T, T\v):
  H₃(T\v) → H₃(T) → H₃(T,T\v) → H₂(T\v) = 0

Since H₂(T\v) = 0 (THM-108): H₃(T,T\v) → 0 exact.
And H₃(T\v) → H₃(T) → H₃(T,T\v) → 0 exact.

So β₃(T) ≤ β₃(T\v) + dim H₃(T,T\v).

Plan:
1. For β₁(T) = 1 tournaments at n=6, compute H₃(T,T\v) for all v
2. If always 0, then β₃(T) ≤ β₃(T\v), and by induction β₃ = 0

opus-2026-03-15-S72d
"""
import numpy as np
from collections import Counter
import sys, time

# Inline path homology
def enum_paths(A, n, p):
    if p == 0: return [(v,) for v in range(n)]
    if p < 0: return []
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j]: adj[i].append(j)
    paths = []; stack = []
    for s in range(n):
        stack.append(([s], 1<<s))
        while stack:
            path, vis = stack.pop()
            if len(path)==p+1: paths.append(tuple(path)); continue
            for u in adj[path[-1]]:
                if not (vis & (1<<u)):
                    stack.append((path+[u], vis|(1<<u)))
    return paths

def bcoeffs(path):
    return [((-1)**i, path[:i]+path[i+1:]) for i in range(len(path))]

def bdmat(ap, apm):
    if not ap or not apm: return np.zeros((max(len(apm),0),max(len(ap),0)))
    idx={p:i for i,p in enumerate(apm)}; M=np.zeros((len(apm),len(ap)))
    for j,p in enumerate(ap):
        for s,f in bcoeffs(p):
            if f in idx: M[idx[f],j]+=s
    return M

def omega_basis(A, n, p, ap, apm):
    d=len(ap)
    if d==0: return np.zeros((0,0))
    if p==0: return np.eye(d)
    aset=set(apm); naf={}; c=0
    for j,path in enumerate(ap):
        for s,f in bcoeffs(path):
            if len(set(f))==len(f) and f not in aset:
                if f not in naf: naf[f]=c; c+=1
    if c==0: return np.eye(d)
    P=np.zeros((c,d))
    for j,path in enumerate(ap):
        for s,f in bcoeffs(path):
            if f in naf: P[naf[f],j]+=s
    U,S,Vt=np.linalg.svd(P,full_matrices=True)
    r=sum(s>1e-10 for s in S); ns=Vt[r:].T
    return ns if ns.shape[1]>0 else np.zeros((d,0))

def betti_and_dims(A, n, max_dim):
    """Return betti numbers AND dimensions of each Ω_p."""
    al={p:([] if p<0 else enum_paths(A,n,p)) for p in range(-1, max_dim+2)}
    om={p:omega_basis(A,n,p,al[p],al[p-1]) for p in range(max_dim+2)}

    betti = []; dims = {}
    for p in range(max_dim+1):
        dp = om[p].shape[1] if om[p].ndim==2 else 0
        dims[p] = dp
        if dp==0: betti.append(0); continue
        bd = bdmat(al[p],al[p-1])@om[p]
        rp = sum(s>1e-8 for s in np.linalg.svd(bd,compute_uv=False)) if min(bd.shape)>0 else 0
        kd = dp-rp
        dp1 = om[p+1].shape[1] if om[p+1].ndim==2 else 0
        if dp1>0:
            bd1=bdmat(al[p+1],al[p])@om[p+1]
            im=sum(s>1e-8 for s in np.linalg.svd(bd1,compute_uv=False)) if min(bd1.shape)>0 else 0
        else: im=0
        betti.append(max(0,kd-im))
    return betti, dims

def tourn(n,bits):
    A=[[0]*n for _ in range(n)]; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if bits&(1<<idx): A[i][j]=1
            else: A[j][i]=1
            idx+=1
    return A

def delete_vertex(A, n, v):
    """Delete vertex v from tournament A."""
    B = []
    for i in range(n):
        if i == v: continue
        row = []
        for j in range(n):
            if j == v: continue
            row.append(A[i][j])
        B.append(row)
    return B

# === MAIN ===
print("=" * 70)
print("TESTING H₃(T, T\\v) = 0 WHEN β₁(T) = 1")
print("=" * 70)

# n=6: exhaustive check of all β₁=1 tournaments
n = 6; E = n*(n-1)//2
t0 = time.time()
b1_count = 0
h3_rel_nonzero = 0
max_rel_h3 = 0

for bits in range(1 << E):
    A = tourn(n, bits)
    b, dims = betti_and_dims(A, n, 3)

    if b[1] != 1: continue
    b1_count += 1

    # For each vertex v, compute β₃(T\v) and check the LES
    for v in range(n):
        Av = delete_vertex(A, n, v)
        bv, dv = betti_and_dims(Av, n-1, 3)

        # H₃(T,T\v) dimension:
        # From LES: 0 → H₃(T,T\v) → H₂(T\v) = 0
        # Wait, this is wrong. Let me reconsider.
        # LES: H₃(T\v) → H₃(T) → H₃(T,T\v) → H₂(T\v) = 0
        # Since H₂(T\v) = 0: im(δ) = 0, so H₃(T,T\v) → 0 exact
        # And ker(δ) = H₃(T,T\v)
        # From previous: im(H₃(T) → H₃(T,T\v)) = H₃(T,T\v) (surjective)
        # So dim H₃(T,T\v) = β₃(T) - dim(ker(H₃(T) → H₃(T,T\v)))
        #                   = β₃(T) - dim(im(H₃(T\v) → H₃(T)))

        # But we know β₃(T) ≤ 0 (because of exhaustive data at n=6),
        # so H₃(T,T\v) = 0 trivially.
        # This doesn't help for the PROOF — we need to show it algebraically.

    if b1_count % 1000 == 0:
        print(f"  Checked {b1_count} β₁=1 tournaments...", flush=True)

print(f"\nn=6: {b1_count} tournaments with β₁=1, time={time.time()-t0:.1f}s")

# What we really need: compute dim H₃(T,T\v) DIRECTLY from the relative complex
# Let me do this for a few cases
print("\n--- Direct relative H₃ computation at n=6 ---")

count = 0
for bits in range(1 << E):
    A = tourn(n, bits)
    b, _ = betti_and_dims(A, n, 3)
    if b[1] != 1: continue
    count += 1
    if count > 20: break  # just check 20 cases

    for v in range(n):
        # Relative chain complex: R_p = Ω_p(T) / Ω_p(T\v)
        # Paths in T that use vertex v
        # Compute paths in T and T\v
        al_T = {p: enum_paths(A, n, p) for p in range(5)}
        Av = delete_vertex(A, n, v)
        # Reindex vertices: {0,...,n-1}\{v} → {0,...,n-2}
        relabel = {}; idx = 0
        for i in range(n):
            if i != v: relabel[i] = idx; idx += 1
        al_Tv = {p: enum_paths(Av, n-1, p) for p in range(5)}

        # Paths in T that DON'T use v = paths that correspond to T\v paths
        # Paths in T that DO use v = relative generators
        for p in [3]:
            paths_with_v = [path for path in al_T[p] if v in path]
            paths_without_v = [path for path in al_T[p] if v not in path]
            # Check: paths_without_v should biject with al_Tv[p]
            # print(f"    p={p}: |A_p(T)|={len(al_T[p])}, with_v={len(paths_with_v)}, without_v={len(paths_without_v)}, |A_p(T\\v)|={len(al_Tv[p])}")

    if count == 1:
        # Detailed for first tournament
        print(f"  First β₁=1 tournament (bits={bits}):")
        print(f"    β = {b}")
        for v in range(n):
            Av = delete_vertex(A, n, v)
            bv, _ = betti_and_dims(Av, n-1, 3)
            print(f"    T\\{v}: β = {bv}")

# === ALTERNATIVE APPROACH ===
# Instead of relative homology, use Euler characteristic to constrain β₃
print("\n" + "=" * 70)
print("EULER CHARACTERISTIC ANALYSIS")
print("=" * 70)

# Compute χ for ALL n=6 tournaments
print("\nn=6 exhaustive Euler characteristic:")
chi_dist = Counter()
chi_by_type = {}
for bits in range(1 << E):
    A = tourn(n, bits)
    b, _ = betti_and_dims(A, n, 3)
    bt = tuple(b)
    chi = sum((-1)**p * v for p, v in enumerate(b))
    chi_dist[chi] += 1
    if bt not in chi_by_type:
        chi_by_type[bt] = chi
print(f"  χ distribution: {dict(chi_dist)}")
print(f"  By Betti type:")
for bt, chi in sorted(chi_by_type.items()):
    print(f"    β={list(bt)}: χ={chi}")

# Compute Euler char from chain dimensions (alternative formula)
# χ = Σ (-1)^p dim(Ω_p) should equal Σ (-1)^p β_p
print("\n--- Verify χ = Σ(-1)^p dim(Ω_p) ---")
for bits in [0, 1, 100, 500, 1000, 4800]:  # sample
    A = tourn(n, bits)
    b, dims = betti_and_dims(A, n, 3)
    chi_betti = sum((-1)**p * v for p, v in enumerate(b))
    chi_omega = sum((-1)**p * dims.get(p, 0) for p in range(4))
    print(f"  bits={bits}: β={b}, dims={dims}, χ(β)={chi_betti}, χ(Ω)={chi_omega}")

print("\nDone.")

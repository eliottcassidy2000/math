#!/usr/bin/env python3
"""
RELATIVE H₃(T, T\v) COMPUTATION
==================================
For β₁=1 tournaments at n=6, compute H₃(T, T\v) for v not in any free cycle.

If H₃(T, T\v) = 0 for such v, the induction proof of β₁·β₃=0 works:
  β₁(T)=1 → pick v not in free cycle → β₁(T\v)=1 → β₃(T\v)=0 (induction)
  → LES gives β₃(T) = dim H₃(T,T\v) = 0. QED.

The relative complex: R_p = chains in Ω_p(T) that "use" vertex v.
More precisely: the quotient Ω_p(T) / (image of Ω_p(T\v) under inclusion).

opus-2026-03-15-S72d
"""
import numpy as np
from collections import Counter
import time

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

def path_betti(A, n, max_dim):
    al={p:([] if p<0 else enum_paths(A,n,p)) for p in range(-1, max_dim+2)}
    om={p:omega_basis(A,n,p,al[p],al[p-1]) for p in range(max_dim+2)}
    betti = []
    for p in range(max_dim+1):
        dp = om[p].shape[1] if om[p].ndim==2 else 0
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
    return betti

def tourn(n,bits):
    A=[[0]*n for _ in range(n)]; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if bits&(1<<idx): A[i][j]=1
            else: A[j][i]=1
            idx+=1
    return A

def delete_vertex(A, n, v):
    B = []
    for i in range(n):
        if i == v: continue
        row = []
        for j in range(n):
            if j == v: continue
            row.append(A[i][j])
        B.append(row)
    return B

def is_3cycle(A, i, j, k):
    w = [A[i][j]+A[i][k], A[j][i]+A[j][k], A[k][i]+A[k][j]]
    return max(w) < 2

def is_free_cycle(A, n, i, j, k):
    """Check if 3-cycle {i,j,k} is free (not dominated)."""
    if not is_3cycle(A, i, j, k): return False
    for v in range(n):
        if v in (i,j,k): continue
        if A[v][i] and A[v][j] and A[v][k]: return False
        if A[i][v] and A[j][v] and A[k][v]: return False
    return True

def find_free_cycles(A, n):
    """Find all free 3-cycles."""
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if is_free_cycle(A, n, i, j, k):
                    cycles.append((i,j,k))
    return cycles

def compute_relative_h3(A, n, v, max_p=4):
    """
    Compute H₃(T, T\\v) via the LES.

    Using: H₃(T\\v) -> H₃(T) -> H₃(T,T\\v) -> H₂(T\\v) = 0
    Since H₂(T\\v) = 0: H₃(T) -> H₃(T,T\\v) surjective.
    ker = im(H₃(T\\v) -> H₃(T)).
    dim H₃(T,T\\v) = β₃(T) - dim(im(H₃(T\\v) -> H₃(T)))

    But we can compute it directly from the long exact sequence dimensions.
    """
    # Compute β₃(T) and β₃(T\\v) first
    b_T = path_betti(A, n, max_p)
    Av = delete_vertex(A, n, v)
    b_Tv = path_betti(Av, n-1, max_p)

    # From LES with β₂(T) = β₂(T\\v) = 0:
    # H₃(T\\v) --f--> H₃(T) --g--> H₃(T,T\\v) --h--> H₂(T\\v) = 0
    # Exact: im(f) = ker(g), im(g) = ker(h) = H₃(T,T\\v) (since h=0)
    # So g is surjective: dim H₃(T,T\\v) = β₃(T) - dim(ker g) = β₃(T) - dim(im f)
    # dim(im f) ≤ min(β₃(T\\v), β₃(T))

    # Continue LES upward:
    # H₄(T,T\\v) --δ--> H₃(T\\v) --f--> H₃(T) --g--> H₃(T,T\\v) --h--> 0
    # im(δ) = ker(f), so dim(im f) = β₃(T\\v) - dim(ker f) = β₃(T\\v) - dim(im δ)

    # We need to compute the inclusion-induced maps. For the FULL LES we need the
    # actual chain maps. Let me compute the connecting homomorphism dimension instead.

    # ALTERNATIVE: Direct computation via quotient complex.
    # R_p = Ω_p(T) restricted to paths using v
    # = {c ∈ Ω_p(T) : c uses v} (or more precisely, quotient by T\\v paths)

    # For now, just use the dimension formula:
    # dim H₃(T,T\\v) = β₃(T) - (β₃(T\\v) - dim(im δ: H₄(T,T\\v) -> H₃(T\\v)))
    # which simplifies when β₃(T\\v) = 0: dim H₃(T,T\\v) = β₃(T)

    return b_T, b_Tv

# === MAIN ===
print("=" * 70)
print("RELATIVE H₃ FOR β₁=1 TOURNAMENTS AT n=6")
print("=" * 70)

n = 6; E = n*(n-1)//2
t0 = time.time()

results = []
count = 0

for bits in range(1 << E):
    A = tourn(n, bits)
    b = path_betti(A, n, 3)

    if b[1] != 1: continue
    count += 1
    if count > 200: break  # sample

    free_cycles = find_free_cycles(A, n)
    free_verts = set()
    for c in free_cycles:
        free_verts.update(c)

    non_free_verts = [v for v in range(n) if v not in free_verts]

    for v in non_free_verts[:2]:  # check 2 non-free vertices per tournament
        Av = delete_vertex(A, n, v)
        bv = path_betti(Av, n-1, 3)

        results.append({
            'bits': bits,
            'v': v,
            'beta_T': b,
            'beta_Tv': bv,
            'free_verts': sorted(free_verts),
            'v_in_free': v in free_verts,
        })

    if count % 50 == 0:
        print(f"  {count} β₁=1 tournaments checked, t={time.time()-t0:.0f}s", flush=True)

print(f"\nProcessed {count} β₁=1 tournaments, {len(results)} (T,v) pairs, t={time.time()-t0:.1f}s")

# Analyze
print("\n--- Results for v NOT in free cycle ---")
b1_Tv = Counter()
b3_Tv = Counter()
rel_h3 = Counter()

for r in results:
    if r['v_in_free']: continue
    bv = r['beta_Tv']
    b1_Tv[bv[1]] += 1
    b3_Tv[bv[3] if len(bv) > 3 else 0] += 1
    # H₃(T,T\v) = β₃(T) since β₃(T\v)=0 when β₁(T\v)=1
    # Check: does β₁(T\v) = 1?
    if bv[1] == 1:
        h3_rel = 0 if len(bv) <= 3 else bv[3]  # if β₃(T\v)=0, then H₃(T,T\v) = β₃(T) = 0
        rel_h3[h3_rel] += 1

print(f"  β₁(T\\v) distribution: {dict(sorted(b1_Tv.items()))}")
print(f"  β₃(T\\v) distribution: {dict(sorted(b3_Tv.items()))}")

# KEY CHECK: When v not in free cycle, is β₁(T\v) always 1?
all_b1_1 = all(r['beta_Tv'][1] == 1 for r in results if not r['v_in_free'])
print(f"\n  *** β₁(T\\v) = 1 for ALL v not in free cycle? {all_b1_1} ***")

if all_b1_1:
    print("  This confirms the induction argument:")
    print("  β₁(T)=1 → pick v not in free cycle → β₁(T\\v)=1")
    print("  → β₃(T\\v)=0 by induction")
    print("  → β₃(T) = dim H₃(T,T\\v)")
    print("  Still need: H₃(T,T\\v) = 0")

# Show some examples
print("\n--- Sample (T, T\\v) pairs ---")
for r in results[:10]:
    if r['v_in_free']: continue
    print(f"  bits={r['bits']}, v={r['v']}: β(T)={r['beta_T']}, β(T\\v)={r['beta_Tv']}, free_verts={r['free_verts']}")

# Also check: for v IN free cycle, what happens?
print("\n--- Results for v IN free cycle ---")
in_free_results = [r for r in results if r['v_in_free']]
if not in_free_results:
    # Need to add these
    print("  (No in-free results in sample, extending...)")
    extra_count = 0
    for bits in range(1 << E):
        A = tourn(n, bits)
        b = path_betti(A, n, 3)
        if b[1] != 1: continue
        free_cycles = find_free_cycles(A, n)
        free_verts = set()
        for c in free_cycles:
            free_verts.update(c)
        for v in range(n):
            if v in free_verts:
                Av = delete_vertex(A, n, v)
                bv = path_betti(Av, n-1, 3)
                in_free_results.append({'bits': bits, 'v': v, 'beta_T': b, 'beta_Tv': bv})
                extra_count += 1
                if extra_count >= 20: break
        if extra_count >= 20: break

b1_Tv_in = Counter()
b3_Tv_in = Counter()
for r in in_free_results:
    bv = r['beta_Tv']
    b1_Tv_in[bv[1]] += 1
    b3_Tv_in[bv[3] if len(bv) > 3 else 0] += 1

print(f"  β₁(T\\v) for v in free cycle: {dict(sorted(b1_Tv_in.items()))}")
print(f"  β₃(T\\v) for v in free cycle: {dict(sorted(b3_Tv_in.items()))}")

print("\nDone.")

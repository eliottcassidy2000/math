#!/usr/bin/env python3
"""
MECHANISM ANALYSIS: Why β₁·β₃ = 0
===================================
At n=6: analyze the chain complex dimensions for β₁=1 vs β₃=1 tournaments.
What structural invariant separates these classes?

Key dimensions to track:
- dim(Ω_p) for p=0,...,4
- dim(ker ∂_p) and dim(im ∂_{p+1}) for each p
- The gap: β_p = ker(∂_p) - im(∂_{p+1})

Also track tournament invariants:
- Score sequence, 3-cycle count, domination structure
- Number of "free" 3-cycles (not dominated by external vertex)

opus-2026-03-15-S72d
"""
import numpy as np
from collections import Counter
import time

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

def full_analysis(A, n, max_dim=4):
    """Return detailed chain complex info."""
    al={p:([] if p<0 else enum_paths(A,n,p)) for p in range(-1, max_dim+2)}
    om={p:omega_basis(A,n,p,al[p],al[p-1]) for p in range(max_dim+2)}

    info = {'A_dims': {}, 'Omega_dims': {}, 'ker_dims': {}, 'im_dims': {}, 'betti': []}

    for p in range(max_dim+1):
        info['A_dims'][p] = len(al[p])
        dp = om[p].shape[1] if om[p].ndim==2 else 0
        info['Omega_dims'][p] = dp

        if dp==0:
            info['ker_dims'][p] = 0
            info['im_dims'][p] = 0
            info['betti'].append(0)
            continue

        bd = bdmat(al[p],al[p-1])@om[p]
        rp = sum(s>1e-8 for s in np.linalg.svd(bd,compute_uv=False)) if min(bd.shape)>0 else 0
        kd = dp-rp
        info['ker_dims'][p] = kd

        dp1 = om[p+1].shape[1] if om[p+1].ndim==2 else 0
        if dp1>0:
            bd1=bdmat(al[p+1],al[p])@om[p+1]
            im=sum(s>1e-8 for s in np.linalg.svd(bd1,compute_uv=False)) if min(bd1.shape)>0 else 0
        else: im=0
        info['im_dims'][p] = im
        info['betti'].append(max(0,kd-im))

    return info

def tourn(n,bits):
    A=[[0]*n for _ in range(n)]; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if bits&(1<<idx): A[i][j]=1
            else: A[j][i]=1
            idx+=1
    return A

def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                w = [A[i][j]+A[i][k], A[j][i]+A[j][k], A[k][i]+A[k][j]]
                if max(w) < 2: c3 += 1
    return c3

def count_free_3cycles(A, n):
    """Count 3-cycles not dominated by any external vertex."""
    free = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                w = [A[i][j]+A[i][k], A[j][i]+A[j][k], A[k][i]+A[k][j]]
                if max(w) >= 2: continue  # not a 3-cycle
                # Check if any external vertex dominates all of {i,j,k}
                dominated = False
                for v in range(n):
                    if v in (i,j,k): continue
                    if A[v][i] and A[v][j] and A[v][k]:
                        dominated = True
                        break
                if not dominated:
                    free += 1
    return free

def score_seq(A, n):
    return tuple(sorted([sum(A[i]) for i in range(n)]))

# === MAIN ===
print("=" * 70)
print("MECHANISM: WHY β₁·β₃ = 0")
print("=" * 70)

n = 6; E = n*(n-1)//2
t0 = time.time()

# Collect data by class
b1_data = []  # β₁=1 tournaments
b3_data = []  # β₃=1 tournaments
b0_data = []  # contractible

for bits in range(1 << E):
    A = tourn(n, bits)
    info = full_analysis(A, n, max_dim=4)
    b = info['betti']

    record = {
        'bits': bits,
        'betti': b,
        'Omega': info['Omega_dims'],
        'ker': info['ker_dims'],
        'im': info['im_dims'],
        'c3': count_3cycles(A, n),
        'free3': count_free_3cycles(A, n),
        'score': score_seq(A, n),
    }

    if b[1] == 1: b1_data.append(record)
    elif len(b) > 3 and b[3] == 1: b3_data.append(record)
    else: b0_data.append(record)

    if (bits+1) % 10000 == 0:
        print(f"  {bits+1}/{1<<E}, t={time.time()-t0:.0f}s, β₁=1:{len(b1_data)}, β₃=1:{len(b3_data)}", flush=True)

print(f"\nDone in {time.time()-t0:.1f}s")
print(f"  Contractible: {len(b0_data)}")
print(f"  β₁=1: {len(b1_data)}")
print(f"  β₃=1: {len(b3_data)}")

# === COMPARE CHAIN COMPLEX DIMENSIONS ===
print("\n" + "=" * 70)
print("CHAIN COMPLEX DIMENSIONS BY CLASS")
print("=" * 70)

for label, data in [("β₁=1", b1_data), ("β₃=1", b3_data)]:
    print(f"\n--- {label} ({len(data)} tournaments) ---")

    # Aggregate Ω dimensions
    omega_dist = Counter()
    ker_dist = Counter()
    im_dist = Counter()

    for r in data:
        omega_dist[tuple(r['Omega'].get(p,0) for p in range(5))] += 1
        ker_dist[tuple(r['ker'].get(p,0) for p in range(5))] += 1
        im_dist[tuple(r['im'].get(p,0) for p in range(5))] += 1

    print(f"  Ω dims distribution:")
    for dims, cnt in sorted(omega_dist.items(), key=lambda x: -x[1])[:10]:
        print(f"    Ω = {list(dims)}: {cnt}")

    print(f"  ker(∂) distribution:")
    for dims, cnt in sorted(ker_dist.items(), key=lambda x: -x[1])[:10]:
        print(f"    ker = {list(dims)}: {cnt}")

    print(f"  im(∂) distribution:")
    for dims, cnt in sorted(im_dist.items(), key=lambda x: -x[1])[:10]:
        print(f"    im = {list(dims)}: {cnt}")

# === KEY: COMPARE Ω₃ and ker/im at dimension 3 ===
print("\n" + "=" * 70)
print("DIMENSION 3 COMPARISON: THE KEY LEVEL")
print("=" * 70)

for label, data in [("β₁=1", b1_data), ("β₃=1", b3_data), ("contractible", b0_data[:100])]:
    print(f"\n--- {label} ---")
    om3_dist = Counter()
    ker3_dist = Counter()
    im3_dist = Counter()
    gap3_dist = Counter()

    for r in data:
        om3 = r['Omega'].get(3, 0)
        k3 = r['ker'].get(3, 0)
        i3 = r['im'].get(3, 0)
        om3_dist[om3] += 1
        ker3_dist[k3] += 1
        im3_dist[i3] += 1
        gap3_dist[k3 - i3] += 1

    print(f"  dim(Ω₃): {dict(sorted(om3_dist.items()))}")
    print(f"  dim(ker ∂₃): {dict(sorted(ker3_dist.items()))}")
    print(f"  dim(im ∂₄→Ω₃): {dict(sorted(im3_dist.items()))}")
    print(f"  β₃ = ker-im: {dict(sorted(gap3_dist.items()))}")

# === 3-CYCLE STRUCTURE ===
print("\n" + "=" * 70)
print("3-CYCLE AND FREE-CYCLE STRUCTURE")
print("=" * 70)

for label, data in [("β₁=1", b1_data), ("β₃=1", b3_data)]:
    print(f"\n--- {label} ---")
    c3_dist = Counter()
    free3_dist = Counter()
    score_dist = Counter()

    for r in data:
        c3_dist[r['c3']] += 1
        free3_dist[r['free3']] += 1
        score_dist[r['score']] += 1

    print(f"  c₃ (total 3-cycles): {dict(sorted(c3_dist.items()))}")
    print(f"  Free 3-cycles: {dict(sorted(free3_dist.items()))}")
    print(f"  Score sequences (top 5):")
    for sc, cnt in sorted(score_dist.items(), key=lambda x: -x[1])[:5]:
        print(f"    {list(sc)}: {cnt}")

# === DOMINATION TEST ===
print("\n" + "=" * 70)
print("DOMINATION TEST: DO ALL β₃=1 HAVE 0 FREE 3-CYCLES?")
print("=" * 70)

all_dominated = all(r['free3'] == 0 for r in b3_data)
print(f"  All β₃=1 tournaments have free3=0: {all_dominated}")

all_free = all(r['free3'] > 0 for r in b1_data)
print(f"  All β₁=1 tournaments have free3>0: {all_free}")

if all_dominated and all_free:
    print("\n  *** STRUCTURAL DICHOTOMY CONFIRMED ***")
    print("  β₁=1 ⟹ free 3-cycles exist")
    print("  β₃=1 ⟹ all 3-cycles dominated")
    print("  These are LOGICALLY INCOMPATIBLE ⟹ β₁·β₃ = 0")

# === EULER CHARACTERISTIC ===
print("\n" + "=" * 70)
print("EULER CHARACTERISTIC")
print("=" * 70)

chi_all = Counter()
for data in [b0_data, b1_data, b3_data]:
    for r in data:
        chi = sum((-1)**p * v for p, v in enumerate(r['betti']))
        chi_all[chi] += 1

print(f"  χ distribution (all n=6): {dict(sorted(chi_all.items()))}")

print("\nDone.")

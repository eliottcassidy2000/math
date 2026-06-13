#!/usr/bin/env python3
"""
FIX: Correct domination check for 3-cycles.
A 3-cycle {i,j,k} is DOMINATED if there exists external vertex v such that:
  v beats all of {i,j,k}, OR all of {i,j,k} beat v.
A "free" cycle has no such dominator.

opus-2026-03-15-S72d
"""
import numpy as np
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

def is_3cycle(A, i, j, k):
    """Check if {i,j,k} form a directed 3-cycle."""
    w = [A[i][j]+A[i][k], A[j][i]+A[j][k], A[k][i]+A[k][j]]
    return max(w) < 2  # each vertex has out-degree 1 within {i,j,k}

def count_free_3cycles_CORRECT(A, n):
    """Count 3-cycles not dominated by any external vertex.
    Dominated = ∃v ∉ {i,j,k}: v beats all three OR all three beat v."""
    free = 0
    total_cycles = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if not is_3cycle(A, i, j, k): continue
                total_cycles += 1
                dominated = False
                for v in range(n):
                    if v in (i,j,k): continue
                    # v beats all three?
                    if A[v][i] and A[v][j] and A[v][k]:
                        dominated = True; break
                    # all three beat v?
                    if A[i][v] and A[j][v] and A[k][v]:
                        dominated = True; break
                if not dominated:
                    free += 1
    return free, total_cycles

# === MAIN ===
print("=" * 70)
print("CORRECT DOMINATION CHECK FOR β₃=1 TOURNAMENTS")
print("=" * 70)

n = 6; E = n*(n-1)//2
t0 = time.time()

b1_free = []
b3_free = []
neither_free = []

for bits in range(1 << E):
    A = tourn(n, bits)
    b = path_betti(A, n, 3)

    free, total = count_free_3cycles_CORRECT(A, n)

    if b[1] == 1:
        b1_free.append((bits, free, total))
    elif len(b) > 3 and b[3] == 1:
        b3_free.append((bits, free, total))
    else:
        neither_free.append((bits, free, total))

    if (bits+1) % 10000 == 0:
        print(f"  {bits+1}/{1<<E}, t={time.time()-t0:.0f}s", flush=True)

print(f"\nDone in {time.time()-t0:.1f}s")

from collections import Counter

print(f"\n--- β₁=1 ({len(b1_free)} tournaments) ---")
fc_dist = Counter(f for _, f, _ in b1_free)
print(f"  Free 3-cycles: {dict(sorted(fc_dist.items()))}")
print(f"  All have free3 > 0: {all(f > 0 for _, f, _ in b1_free)}")

print(f"\n--- β₃=1 ({len(b3_free)} tournaments) ---")
fc_dist = Counter(f for _, f, _ in b3_free)
print(f"  Free 3-cycles: {dict(sorted(fc_dist.items()))}")
print(f"  All have free3 = 0: {all(f == 0 for _, f, _ in b3_free)}")

# Show some β₃=1 examples
print(f"\n  β₃=1 examples:")
for bits, free, total in b3_free[:5]:
    A = tourn(n, bits)
    sc = tuple(sorted([sum(A[i]) for i in range(n)]))
    print(f"    bits={bits}: score={list(sc)}, c₃={total}, free={free}")
    if free > 0:
        # Show the free cycles
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    if not is_3cycle(A, i, j, k): continue
                    dominated = False
                    for v in range(n):
                        if v in (i,j,k): continue
                        if (A[v][i] and A[v][j] and A[v][k]) or (A[i][v] and A[j][v] and A[k][v]):
                            dominated = True; break
                    if not dominated:
                        print(f"      FREE: {{{i},{j},{k}}}")

print(f"\n--- Contractible ({len(neither_free)} tournaments) ---")
fc_dist = Counter(f for _, f, _ in neither_free)
print(f"  Free 3-cycles: {dict(sorted(fc_dist.items()))}")

print(f"\n=== STRUCTURAL DICHOTOMY ===")
b3_all_dominated = all(f == 0 for _, f, _ in b3_free)
b1_all_free = all(f > 0 for _, f, _ in b1_free)
print(f"  β₃=1 => all 3-cycles dominated: {b3_all_dominated}")
print(f"  β₁=1 => some free 3-cycle exists: {b1_all_free}")
if b3_all_dominated and b1_all_free:
    print(f"  *** STRUCTURAL DICHOTOMY CONFIRMED ***")
    print(f"  β₃=1 and β₁=1 are logically incompatible")
else:
    print(f"  Structural dichotomy NOT confirmed. Need different mechanism.")

print("\nDone.")

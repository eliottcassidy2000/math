#!/usr/bin/env python3
"""
SEESAW IDENTITY VERIFICATION
==============================
From the chain complex with β₂ = 0:
  β₁ + β₃ = ker(d₁) + dim(Ω₃) - dim(Ω₂) - im(d₄)

This sum does NOT depend on im(d₂)! Verify this identity and check
whether β₁ + β₃ ≤ 1 always holds (which would prove β₁·β₃ = 0
since β₁ ∈ {0,1}).

Also check: is rank(d₂|Ω₂) = C(n,2) - n + 1 - β₁ always?
(This is the formula from the memory: "rank(∂_2|Ω_2) = C(n,2) - n + 1 - β_1")

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

def chain_complex_dims(A, n, max_dim=4):
    """Return all chain complex dimensions needed for the seesaw identity."""
    al={p:([] if p<0 else enum_paths(A,n,p)) for p in range(-1, max_dim+2)}
    om={p:omega_basis(A,n,p,al[p],al[p-1]) for p in range(max_dim+2)}

    result = {}
    for p in range(max_dim+1):
        dp = om[p].shape[1] if om[p].ndim==2 else 0
        result[f'Omega_{p}'] = dp

        if dp==0:
            result[f'ker_d{p}'] = 0
            result[f'im_d{p}'] = 0
            result[f'beta_{p}'] = 0
            continue

        bd = bdmat(al[p],al[p-1])@om[p]
        rp = sum(s>1e-8 for s in np.linalg.svd(bd,compute_uv=False)) if min(bd.shape)>0 else 0
        kd = dp-rp
        result[f'ker_d{p}'] = kd
        result[f'im_d{p}'] = rp  # rank of d_p: Omega_p -> Omega_{p-1}

        dp1 = om[p+1].shape[1] if om[p+1].ndim==2 else 0
        if dp1>0:
            bd1=bdmat(al[p+1],al[p])@om[p+1]
            im_from_above=sum(s>1e-8 for s in np.linalg.svd(bd1,compute_uv=False)) if min(bd1.shape)>0 else 0
        else: im_from_above=0
        result[f'im_d{p+1}_into_{p}'] = im_from_above
        result[f'beta_{p}'] = max(0, kd - im_from_above)

    return result

def tourn(n,bits):
    A=[[0]*n for _ in range(n)]; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if bits&(1<<idx): A[i][j]=1
            else: A[j][i]=1
            idx+=1
    return A

# === MAIN ===
print("=" * 70)
print("SEESAW IDENTITY VERIFICATION")
print("=" * 70)

for n in [5, 6]:
    E = n*(n-1)//2
    total = 1 << E
    print(f"\n--- n={n} EXHAUSTIVE ({total} tournaments) ---")
    t0 = time.time()

    kerd1_vals = Counter()
    imd2_vals = Counter()  # rank of d₂: Ω₂ → Ω₁ (this is im(d₂))
    imd3_into2_vals = Counter()  # rank of d₃: Ω₃ → Ω₂ restricted
    identity_check = Counter()
    b1_plus_b3 = Counter()

    # Track by β class
    imd2_by_class = {}
    omega_by_class = {}
    imd4_by_class = {}

    for bits in range(total):
        A = tourn(n, bits)
        r = chain_complex_dims(A, n, max_dim=4)

        b1 = r['beta_1']
        b3 = r['beta_3']
        b_class = 'b1' if b1 > 0 else ('b3' if b3 > 0 else 'neither')

        # im(d₂) = rank of d₂ restricted to Ω₂ = im_d2 (rank of boundary map)
        # But careful: im_d{p} as computed is the rank of d_p: Ω_p → Ω_{p-1}
        # So im_d2 = rank of d_2 map = dim(Ω_2) - ker(d_2)
        rank_d2 = r['im_d2']  # This is the rank of d₂|Ω₂ (how much of Ω₂ maps into Ω₁)

        # Wait - I need to be more careful. Let me re-examine:
        # im_d{p} = rp = rank of (bd_p @ omega_p), where bd_p = bdmat(al[p], al[p-1])
        # This is the dimension of the image of d_p in A_{p-1}
        # ker_d{p} = dim(Ω_p) - im_d{p}
        # im_d{p+1}_into_{p} = rank of d_{p+1} from Ω_{p+1} into Ω_p
        # beta_p = ker_d{p} - im_d{p+1}_into_{p}

        kerd1 = r['ker_d1']
        kerd1_vals[kerd1] += 1

        # im(d₂) into Ω₁: this is im_d2 (the rank of d₂: Ω₂ → Ω₁)
        im_d2 = r['im_d2']  # = rank(d₂|Ω₂)
        imd2_vals[(b_class, im_d2)] += 1

        # The image of d₃ into Ω₂
        im_d3 = r['im_d3_into_2'] if 'im_d3_into_2' in r else 0

        # Check: β₂ = ker_d2 - im_d3_into_2
        # ker_d2 = Omega_2 - im_d2
        ker_d2 = r['Omega_2'] - im_d2
        beta2_check = ker_d2 - im_d3

        # im(d₄) into Ω₃
        im_d4 = r.get('im_d4_into_3', 0)

        # Seesaw identity: β₁ + β₃ = ker(d₁) + dim(Ω₃) - dim(Ω₂) - im(d₄)
        # where im(d₄) = image of d₄ into Ω₃
        lhs = b1 + b3
        rhs = kerd1 + r['Omega_3'] - r['Omega_2'] - im_d4
        identity_check[(lhs, rhs, lhs == rhs)] += 1
        b1_plus_b3[(b1, b3, lhs)] += 1

        if b_class not in imd2_by_class:
            imd2_by_class[b_class] = Counter()
            omega_by_class[b_class] = Counter()
            imd4_by_class[b_class] = Counter()
        imd2_by_class[b_class][im_d2] += 1
        omega_by_class[b_class][(r['Omega_2'], r['Omega_3'])] += 1
        imd4_by_class[b_class][im_d4] += 1

        if (bits+1) % 10000 == 0:
            print(f"  {bits+1}/{total}, t={time.time()-t0:.0f}s", flush=True)

    print(f"  Done in {time.time()-t0:.1f}s")

    print(f"\n  ker(d₁) values: {dict(kerd1_vals)}")
    print(f"    Expected: C({n},2) - {n} + 1 = {n*(n-1)//2 - n + 1}")

    print(f"\n  im(d₂) = rank(d₂|Ω₂) by β class:")
    for k, v in sorted(imd2_by_class.items()):
        print(f"    {k}: {dict(sorted(v.items()))}")

    print(f"\n  (Ω₂, Ω₃) by β class:")
    for k, v in sorted(omega_by_class.items()):
        print(f"    {k}: {dict(sorted(v.items()))}")

    print(f"\n  im(d₄) by β class:")
    for k, v in sorted(imd4_by_class.items()):
        print(f"    {k}: {dict(sorted(v.items()))}")

    print(f"\n  Seesaw identity β₁+β₃ = ker(d₁)+Ω₃-Ω₂-im(d₄):")
    for (l, r_val, eq), cnt in sorted(identity_check.items()):
        print(f"    LHS={l}, RHS={r_val}, equal={eq}: {cnt}")

    print(f"\n  β₁ + β₃ distribution:")
    for (b1, b3, s), cnt in sorted(b1_plus_b3.items()):
        if cnt > 0:
            print(f"    (β₁={b1}, β₃={b3}), sum={s}: {cnt}")

print("\n" + "=" * 70)
print("KEY QUESTION: Can β₁+β₃ > 1?")
print("If β₁+β₃ ≤ 1 always, then since β₁ ∈ {0,1}, we get β₁·β₃=0.")
print("=" * 70)

print("\nDone.")

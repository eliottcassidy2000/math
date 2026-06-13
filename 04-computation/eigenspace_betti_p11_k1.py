#!/usr/bin/env python3
"""
Verify k=1 eigenspace Betti for P_11 at lower degrees (d≤7).
Just need to confirm β_d^(1)=0 for d≤5 and β_6^(1)=1.

opus-2026-03-13-S71c
"""
import sys, time
import numpy as np
sys.stdout.reconfigure(line_buffering=True)

def qr(p):
    return sorted(set((a*a)%p for a in range(1,p)))

def find_pth_root(p, prime):
    exp = (prime - 1) // p
    for g in range(2, prime):
        omega = pow(g, exp, prime)
        if omega != 1 and all(pow(omega, k, prime) != 1 for k in range(1, p)):
            return omega
    return None

def gauss_rank_mod(M, prime):
    A = M.astype(np.int64) % prime
    rows, cols = A.shape
    pivot_row = rank = 0
    for col in range(cols):
        found = -1
        for row in range(pivot_row, rows):
            if A[row, col] % prime != 0:
                found = row; break
        if found < 0: continue
        A[[pivot_row, found]] = A[[found, pivot_row]]
        inv = pow(int(A[pivot_row, col]), prime - 2, prime)
        A[pivot_row] = (A[pivot_row] * inv) % prime
        for row in range(rows):
            if row != pivot_row and A[row, col] != 0:
                A[row] = (A[row] - int(A[row, col]) * A[pivot_row]) % prime
        rank += 1; pivot_row += 1
    return rank

p = 11
S = qr(p)
m = len(S)
q = 463
omega = find_pth_root(p, q)
max_d = 7
print(f"P_{p}: m={m}, QR={S}, GF({q}), ω={omega}")
t0 = time.time()

# Enumerate
prev_seqs = [()]
prev_ps = {(): frozenset([0])}
prev_last = {(): 0}
all_seqs = {0: [()]}
for d in range(1, max_d + 1):
    curr = []; curr_ps = {}; curr_last = {}
    for seq in prev_seqs:
        ps_set = prev_ps[seq]; last = prev_last[seq]
        for s in S:
            ns = (last + s) % p
            if ns not in ps_set:
                nseq = seq + (s,)
                curr.append(nseq); curr_ps[nseq] = ps_set | {ns}; curr_last[nseq] = ns
    if not curr: break
    all_seqs[d] = curr; prev_seqs = curr; prev_ps = curr_ps; prev_last = curr_last
    print(f"  d={d}: A={len(curr)} ({time.time()-t0:.1f}s)")

allowed = {d: set(all_seqs.get(d, [])) for d in range(max_d + 1)}

# Compute for k=0 and k=1
for k in [1]:  # Skip k=0 (already computed in k0_betti_verify.out)
    omega_k = pow(omega, k, q)
    omega_dims = [1]
    boundary_ranks = {0: 0}

    for d in range(1, max_d + 1):
        A_d = all_seqs.get(d, [])
        A_dm1 = all_seqs.get(d-1, [])
        n_Ad = len(A_d); n_dm1 = len(A_dm1)
        if n_Ad == 0:
            omega_dims.append(0); boundary_ranks[d] = 0; continue

        junk_set = set(); face_junk = []; face_allowed = []
        for D in A_d:
            jf = []; af = []
            for fi in range(d + 1):
                if fi == 0:
                    fd = D[1:]; phase = pow(omega_k, D[0], q) if k != 0 else 1
                elif fi == d:
                    fd = D[:d-1]; phase = 1
                else:
                    fd = D[:fi-1] + ((D[fi-1]+D[fi])%p,) + D[fi+1:]; phase = 1
                sign = 1 if fi % 2 == 0 else q - 1
                coeff = (sign * phase) % q
                if fd in allowed[d-1]: af.append((fd, coeff))
                else: junk_set.add(fd); jf.append((fd, coeff))
            face_junk.append(jf); face_allowed.append(af)

        junk_list = sorted(junk_set)
        n_junk = len(junk_list)
        junk_idx = {j: i for i, j in enumerate(junk_list)}
        idx_dm1 = {s: i for i, s in enumerate(A_dm1)}

        C = np.zeros((n_junk, n_Ad), dtype=np.int64)
        for j, jf in enumerate(face_junk):
            for fd, coeff in jf:
                C[junk_idx[fd], j] = (C[junk_idx[fd], j] + coeff) % q
        rank_c = gauss_rank_mod(C, q) if n_junk > 0 else 0
        omega_dims.append(n_Ad - rank_c)

        CB = np.zeros((n_junk + n_dm1, n_Ad), dtype=np.int64)
        CB[:n_junk, :] = C
        for j, af in enumerate(face_allowed):
            for fd, coeff in af:
                row = n_junk + idx_dm1[fd]
                CB[row, j] = (CB[row, j] + coeff) % q
        rank_cb = gauss_rank_mod(CB, q) if (n_junk + n_dm1) > 0 else 0
        boundary_ranks[d] = rank_cb - rank_c
        print(f"  k={k} d={d}: Ω={omega_dims[-1]}, rank(∂)={boundary_ranks[d]} ({time.time()-t0:.1f}s)")

    betti = []
    for d in range(max_d + 1):
        od = omega_dims[d]; rd = boundary_ranks.get(d, 0); rd1 = boundary_ranks.get(d+1, 0)
        betti.append(od - rd - rd1)

    print(f"\n  k={k}: Ω={omega_dims}")
    print(f"  k={k}: ∂r={[boundary_ranks.get(d,0) for d in range(max_d+1)]}")
    print(f"  k={k}: β={betti}  (last entry is upper bound)")

    if k == 0:
        # Predictions for intermediate vanishing
        print(f"\n  k=0 intermediate vanishing check (β_d=0 for 1≤d≤m-1={m-1}):")
        for d in range(1, m):
            print(f"    β_{d} = {betti[d]} {'✓' if betti[d]==0 else '✗'}")
    if k == 1:
        print(f"\n  k=1 conjecture check (β_d=0 for d≠m+1={m+1}):")
        for d in range(max_d + 1):
            if d == m + 1:
                print(f"    β_{d} = {betti[d]} (expected 1) {'✓' if betti[d]==1 else '✗ (but last may be upper bound)'}")
            elif betti[d] != 0:
                print(f"    β_{d} = {betti[d]} (expected 0) {'✗' if d < max_d else '(upper bound, OK)'}")

print(f"\nTotal time: {time.time()-t0:.1f}s")

#!/usr/bin/env python3
"""
Verify the per-eigenspace Betti number decomposition for Paley tournaments.

CONJECTURE (opus-2026-03-13-S71c):
For Paley T_p (p≡3 mod 4), m=(p-1)/2:

  k=0 eigenspace:  β_d^(0) = { 1 if d=0,  m(m-3)/2 if d∈{m,m+1},  0 otherwise }
  k≠0 eigenspace:  β_d^(k) = { 1 if d=m+1,  0 otherwise }

Total:
  β_m = m(m-3)/2     (all from k=0)
  β_{m+1} = C(m+1,2) (= m(m-3)/2 + (p-1))
  χ = p

This script verifies by computing all eigenspace boundary ranks for P_7 and P_11.
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
    pivot_row = 0
    rank = 0
    for col in range(cols):
        found = -1
        for row in range(pivot_row, rows):
            if A[row, col] % prime != 0:
                found = row
                break
        if found < 0:
            continue
        A[[pivot_row, found]] = A[[found, pivot_row]]
        inv = pow(int(A[pivot_row, col]), prime - 2, prime)
        A[pivot_row] = (A[pivot_row] * inv) % prime
        for row in range(rows):
            if row != pivot_row and A[row, col] != 0:
                A[row] = (A[row] - int(A[row, col]) * A[pivot_row]) % prime
        rank += 1
        pivot_row += 1
    return rank

def compute_eigenspace_betti(p, q):
    """Compute Betti numbers for each eigenspace of P_p over GF(q)."""
    S = qr(p)
    m = len(S)
    omega = find_pth_root(p, q)
    max_d = 2 * m

    print(f"\nP_{p}: m={m}, QR={S}, GF({q}), ω={omega}")
    t0 = time.time()

    # Enumerate diff-sequences
    prev_seqs = [()]
    prev_ps = {(): frozenset([0])}
    prev_last = {(): 0}
    all_seqs = {0: [()]}

    for d in range(1, max_d + 1):
        curr = []
        curr_ps = {}
        curr_last = {}
        for seq in prev_seqs:
            ps_set = prev_ps[seq]
            last = prev_last[seq]
            for s in S:
                ns = (last + s) % p
                if ns not in ps_set:
                    nseq = seq + (s,)
                    curr.append(nseq)
                    curr_ps[nseq] = ps_set | {ns}
                    curr_last[nseq] = ns
        if not curr:
            break
        all_seqs[d] = curr
        prev_seqs = curr
        prev_ps = curr_ps
        prev_last = curr_last

    actual_max_d = max(all_seqs.keys())
    max_d = min(max_d, actual_max_d)
    allowed = {d: set(all_seqs.get(d, [])) for d in range(max_d + 1)}

    print(f"  Enumerated through d={max_d} ({time.time()-t0:.1f}s)")
    print(f"  A_d = {[len(all_seqs.get(d, [])) for d in range(max_d + 1)]}")

    # For each eigenspace k, compute Omega dims and boundary ranks
    results = {}
    for k in range(p):
        omega_k = pow(omega, k, q)
        omega_dims = [1]
        boundary_ranks = {0: 0}

        for d in range(1, max_d + 1):
            A_d = all_seqs.get(d, [])
            A_dm1 = all_seqs.get(d-1, [])
            n_Ad = len(A_d)
            n_dm1 = len(A_dm1)

            if n_Ad == 0:
                omega_dims.append(0)
                boundary_ranks[d] = 0
                continue

            # Compute faces with phase factors
            junk_set = set()
            face_data_junk = []
            face_data_allowed = []

            for D in A_d:
                jf = []
                af = []
                for fi in range(d + 1):
                    if fi == 0:
                        fd = D[1:]
                        phase = pow(omega_k, D[0], q) if k != 0 else 1
                    elif fi == d:
                        fd = D[:d-1]
                        phase = 1
                    else:
                        merged = (D[fi-1] + D[fi]) % p
                        fd = D[:fi-1] + (merged,) + D[fi+1:]
                        phase = 1
                    sign = (1 if fi % 2 == 0 else q - 1)
                    coeff = (sign * phase) % q
                    if fd in allowed[d-1]:
                        af.append((fd, coeff))
                    else:
                        junk_set.add(fd)
                        jf.append((fd, coeff))
                face_data_junk.append(jf)
                face_data_allowed.append(af)

            junk_list = sorted(junk_set)
            n_junk = len(junk_list)
            junk_idx = {j: i for i, j in enumerate(junk_list)}
            idx_dm1 = {s: i for i, s in enumerate(A_dm1)}

            # Constraint matrix
            C = np.zeros((n_junk, n_Ad), dtype=np.int64)
            for j, jf in enumerate(face_data_junk):
                for fd, coeff in jf:
                    C[junk_idx[fd], j] = (C[junk_idx[fd], j] + coeff) % q

            rank_c = gauss_rank_mod(C, q) if n_junk > 0 else 0
            omega_dims.append(n_Ad - rank_c)

            # Combined matrix
            CB = np.zeros((n_junk + n_dm1, n_Ad), dtype=np.int64)
            CB[:n_junk, :] = C
            for j, af in enumerate(face_data_allowed):
                for fd, coeff in af:
                    row = n_junk + idx_dm1[fd]
                    CB[row, j] = (CB[row, j] + coeff) % q
            rank_cb = gauss_rank_mod(CB, q) if (n_junk + n_dm1) > 0 else 0
            boundary_ranks[d] = rank_cb - rank_c

        # Betti
        betti = []
        for d in range(max_d + 1):
            od = omega_dims[d]
            rd = boundary_ranks.get(d, 0)
            rd1 = boundary_ranks.get(d + 1, 0)
            betti.append(od - rd - rd1)

        results[k] = {
            'omega': omega_dims,
            'ranks': [boundary_ranks.get(d, 0) for d in range(max_d + 1)],
            'betti': betti
        }

        if k == 0 or k == 1 or k == min(set(range(1, p)) - set(qr(p))):
            label = "k=0" if k == 0 else f"k={k} ({'QR' if k in set(S) else 'NQR'})"
            print(f"  {label}: Ω={omega_dims}, ∂r={results[k]['ranks']}, β={betti}")

    # Verify eigenspace identity
    all_equal = all(results[k]['omega'] == results[0]['omega'] for k in range(1, p))
    print(f"\n  Eigenspace identity (all Ω equal): {all_equal}")

    # Verify Betti decomposition
    print(f"\n  PER-EIGENSPACE BETTI:")
    for k in range(p):
        b = results[k]['betti']
        nz = [(d, b[d]) for d in range(len(b)) if b[d] != 0]
        label = "k=0" if k == 0 else f"k={k}"
        print(f"    {label}: nonzero at {nz}")

    # Total Betti
    total_betti = [sum(results[k]['betti'][d] for k in range(p)) for d in range(max_d + 1)]
    chi = sum((-1)**d * total_betti[d] for d in range(max_d + 1))
    print(f"\n  TOTAL β = {total_betti}")
    print(f"  χ = {chi} (expected: p={p})")

    # Verify conjecture
    print(f"\n  CONJECTURE VERIFICATION:")
    print(f"    m = {m}")
    print(f"    m(m-3)/2 = {m*(m-3)//2}")
    print(f"    C(m+1,2) = {(m+1)*m//2}")

    # k=0 predictions
    pred_k0 = [0] * (max_d + 1)
    pred_k0[0] = 1
    if m <= max_d:
        pred_k0[m] = m*(m-3)//2
    if m+1 <= max_d:
        pred_k0[m+1] = m*(m-3)//2

    actual_k0 = results[0]['betti']
    k0_match = (pred_k0 == actual_k0)
    print(f"    k=0: predicted {pred_k0}")
    print(f"    k=0: actual    {actual_k0}")
    print(f"    k=0 match: {k0_match}")

    # k≠0 predictions (each should have β_{m+1}=1 only)
    pred_kne0 = [0] * (max_d + 1)
    if m+1 <= max_d:
        pred_kne0[m+1] = 1

    all_kne0_match = True
    for k in range(1, p):
        if results[k]['betti'] != pred_kne0:
            all_kne0_match = False
            print(f"    k={k} MISMATCH: {results[k]['betti']}")

    print(f"    k≠0 all match β_{{m+1}}=1: {all_kne0_match}")

    return results, total_betti

print("=" * 70)
print("EIGENSPACE BETTI DECOMPOSITION FOR PALEY TOURNAMENTS")
print("=" * 70)

# Use a prime q with p | (q-1) for both p=7 and p=11
# Need q ≡ 1 mod 7 and q ≡ 1 mod 11
# q = 463 works: 463 = 66*7 + 1 = 42*11 + 1
q = 463

compute_eigenspace_betti(7, q)
compute_eigenspace_betti(11, q)

print("\n" + "=" * 70)

#!/usr/bin/env python3
"""Compute P_31 k=0 Omega_5 via sparse elimination. opus-2026-03-13-S71c"""
import sys, time, gc
sys.stdout.reconfigure(line_buffering=True)

PRIME = 104729

def qr(p):
    return sorted(set((a*a)%p for a in range(1,p)))

p = 31
S = qr(p)
m = len(S)
print(f"P_{p}: m={m}, QR={S}")

t0 = time.time()
prev_seqs = []
prev_ps = {}
prev_last = {}
allowed_prev = set()

# d=0
prev_seqs = [()]
prev_ps = {(): frozenset([0])}
prev_last = {(): 0}
allowed_prev = set([()])

for d in range(1, 6):
    curr_seqs = []
    curr_ps = {}
    curr_last = {}
    for seq in prev_seqs:
        ps = prev_ps[seq]
        last = prev_last[seq]
        for s in S:
            ns = (last + s) % p
            if ns not in ps:
                nseq = seq + (s,)
                curr_seqs.append(nseq)
                curr_ps[nseq] = ps | {ns}
                curr_last[nseq] = ns

    print(f"  d={d}: A={len(curr_seqs)} ({time.time()-t0:.1f}s)")

    if d < 5:
        # Only need to track allowed for previous degree
        if d == 4:
            allowed_d4 = set(curr_seqs)
        prev_seqs = curr_seqs
        prev_ps = curr_ps
        prev_last = curr_last
    else:
        # d=5: compute Omega
        n_Ad = len(curr_seqs)

        # Free memory
        del prev_ps, prev_last
        gc.collect()

        junk_set = set()
        face_data = []
        for D in curr_seqs:
            faces = []
            for fi in range(d + 1):
                if fi == 0: fd = D[1:]
                elif fi == d: fd = D[:d-1]
                else:
                    merged = (D[fi-1] + D[fi]) % p
                    fd = D[:fi-1] + (merged,) + D[fi+1:]
                is_allowed = (fd in allowed_d4)
                if not is_allowed:
                    junk_set.add(fd)
                    faces.append((fd, 1 if fi % 2 == 0 else -1))
            face_data.append(faces)

        junk_list = sorted(junk_set)
        n_junk = len(junk_list)
        junk_idx = {j: i for i, j in enumerate(junk_list)}
        print(f"  Junk rows: {n_junk} ({time.time()-t0:.1f}s)")

        del junk_set, allowed_d4, curr_seqs, curr_ps, curr_last
        gc.collect()

        # Sparse Gaussian elimination
        pivots = {}
        rank = 0
        for col_idx in range(n_Ad):
            if col_idx % 50000 == 0 and col_idx > 0:
                print(f"    col {col_idx}/{n_Ad}, rank={rank} ({time.time()-t0:.1f}s)")

            col = {}
            for fd, sign in face_data[col_idx]:
                row = junk_idx[fd]
                col[row] = (col.get(row, 0) + sign) % PRIME
                if col[row] == 0:
                    del col[row]

            for r in sorted(col.keys()):
                if r in pivots and r in col:
                    piv_data = pivots[r]
                    factor = col[r]
                    for pr, pv in piv_data.items():
                        col[pr] = (col.get(pr, 0) - factor * pv) % PRIME
                        if col[pr] == 0:
                            del col[pr]

            if not col:
                continue

            pivot_row = min(col.keys())
            inv = pow(col[pivot_row], PRIME - 2, PRIME)
            normalized = {r: (v * inv) % PRIME for r, v in col.items()}
            pivots[pivot_row] = normalized
            rank += 1

        omega_5 = n_Ad - rank
        dt = time.time() - t0
        print(f"\n  rank={rank}, Omega_5={omega_5} ({dt:.1f}s)")
        print(f"  P_31 k=0 Omega = [1, 15, 210, 2835, 36765, {omega_5}]")
        print(f"  Omega_5 mod 15 = {omega_5 % 15}")

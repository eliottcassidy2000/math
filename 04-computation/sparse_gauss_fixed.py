#!/usr/bin/env python3
"""
Re-verify Omega_5 values with FIXED sparse Gaussian elimination.

BUG FIX: The original sparse elimination iterated over sorted(col.keys())
once and subtracted pivots in order. But pivot subtraction introduces fill-in
at rows not in the original column, which then need to be eliminated against
their pivots. The fix is to restart the row scan after each pivot subtraction.

opus-2026-03-13-S71c
"""
import sys, time
sys.stdout.reconfigure(line_buffering=True)

PRIME = 104729

def qr(p):
    return sorted(set((a*a)%p for a in range(1,p)))

def sparse_rank_fixed(face_data, junk_idx, n_Ad, prime):
    """Sparse Gaussian elimination with CORRECT fill-in handling."""
    pivots = {}
    rank = 0
    for col_idx in range(n_Ad):
        if col_idx % 50000 == 0 and col_idx > 0:
            print(f"    col {col_idx}/{n_Ad}, rank={rank}", flush=True)

        col = {}
        for fd, sign in face_data[col_idx]:
            row = junk_idx[fd]
            col[row] = (col.get(row, 0) + sign) % prime
            if col[row] == 0:
                del col[row]

        # Eliminate with correct fill-in handling
        changed = True
        while changed:
            changed = False
            for r in sorted(col.keys()):
                if r in pivots and r in col:
                    piv_data = pivots[r]
                    factor = col[r]
                    for pr, pv in piv_data.items():
                        col[pr] = (col.get(pr, 0) - factor * pv) % prime
                        if col[pr] == 0:
                            del col[pr]
                    changed = True
                    break  # restart

        if not col:
            continue

        pivot_row = min(col.keys())
        inv = pow(col[pivot_row], prime - 2, prime)
        normalized = {r: (v * inv) % prime for r, v in col.items()}
        pivots[pivot_row] = normalized
        rank += 1

    return rank

def compute_omega5(p):
    S = qr(p)
    m = len(S)
    print(f"\nP_{p}: m={m}")

    t0 = time.time()
    prev_seqs = [()]
    prev_ps = {(): frozenset([0])}
    prev_last = {(): 0}

    for d in range(1, 5):
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
        if d < 4:
            prev_seqs = curr_seqs
            prev_ps = curr_ps
            prev_last = curr_last
        else:
            allowed_d4 = set(curr_seqs)
            prev_seqs = curr_seqs
            prev_ps = curr_ps
            prev_last = curr_last

    # d=5
    d = 5
    curr_seqs = []
    for seq in prev_seqs:
        ps = prev_ps[seq]
        last = prev_last[seq]
        for s in S:
            ns = (last + s) % p
            if ns not in ps:
                nseq = seq + (s,)
                curr_seqs.append(nseq)

    n_Ad = len(curr_seqs)
    print(f"  d={d}: A={n_Ad} ({time.time()-t0:.1f}s)")

    import gc
    del prev_ps, prev_last, prev_seqs
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

    del junk_set, allowed_d4, curr_seqs
    gc.collect()

    rank = sparse_rank_fixed(face_data, junk_idx, n_Ad, PRIME)
    omega_5 = n_Ad - rank
    dt = time.time() - t0
    print(f"  rank={rank}, Omega_5={omega_5} ({dt:.1f}s)")
    return omega_5

# Verify P_7 (known: Omega_5 = 6)
o7 = compute_omega5(7)
print(f"  P_7 expected: 6, got: {o7}, MATCH={o7==6}")

# Verify P_11 (known: Omega_5 = 460)
o11 = compute_omega5(11)
print(f"  P_11 expected: 460, got: {o11}, MATCH={o11==460}")

# Re-verify P_19
o19 = compute_omega5(19)
print(f"  P_19 (was 12602): {o19}")

# Re-verify P_23
o23 = compute_omega5(23)
print(f"  P_23 (was 50715): {o23}")

# Re-verify P_31
o31 = compute_omega5(31)
print(f"  P_31 (was 252065): {o31}")

# Re-verify P_43
o43 = compute_omega5(43)
print(f"  P_43 (was 1429652): {o43}")

print("\n" + "=" * 60)
print("SUMMARY:")
print(f"  P_7:  Omega_5 = {o7}")
print(f"  P_11: Omega_5 = {o11}")
print(f"  P_19: Omega_5 = {o19}")
print(f"  P_23: Omega_5 = {o23}")
print(f"  P_31: Omega_5 = {o31}")
print(f"  P_43: Omega_5 = {o43}")

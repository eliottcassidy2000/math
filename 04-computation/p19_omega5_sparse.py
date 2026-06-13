#!/usr/bin/env python3
"""
Compute P_19 k=0 eigenspace Omega_5 using sparse modular elimination.

The constraint matrix at d=5 has ~40k columns but is very sparse
(at most 6 nonzeros per column). Use dict-of-dicts sparse format.

opus-2026-03-13-S71c
"""
import sys, time
sys.stdout.reconfigure(line_buffering=True)

PRIME = 104729

def qr(p):
    return sorted(set((a*a)%p for a in range(1,p)))

def sparse_gauss_rank(cols_data, n_rows, prime):
    """Sparse Gaussian elimination mod prime.
    cols_data: list of dicts {row_idx: value} for each column.
    Returns rank.
    """
    # Pivot info: pivot_col[row] = column index of pivot in that row
    # pivot_val[row] = pivot value (should be 1 after normalization)
    pivots = {}  # row -> (col_idx, normalized column dict)

    rank = 0
    for col_idx, col in enumerate(cols_data):
        if col_idx % 5000 == 0 and col_idx > 0:
            print(f"    processed {col_idx}/{len(cols_data)} columns, rank={rank}", flush=True)

        # Make a working copy
        work = dict(col)

        # Eliminate using existing pivots
        rows_to_check = sorted(work.keys())
        for r in rows_to_check:
            if r in pivots and r in work:
                piv_col, piv_data = pivots[r]
                factor = work[r]  # since pivot is normalized to 1
                # Subtract factor * pivot column
                for pr, pv in piv_data.items():
                    work[pr] = (work.get(pr, 0) - factor * pv) % prime
                    if work[pr] == 0:
                        del work[pr]

        # Find new pivot (smallest row with nonzero entry)
        if not work:
            continue

        pivot_row = min(work.keys())
        pivot_val = work[pivot_row]

        # Normalize
        inv = pow(pivot_val, prime - 2, prime)
        normalized = {}
        for r, v in work.items():
            normalized[r] = (v * inv) % prime

        pivots[pivot_row] = (col_idx, normalized)
        rank += 1

    return rank

p = 19
S = qr(p)
m = len(S)
print(f"P_{p}: m={m}, QR={S}")
print(f"Computing Omega_5...")

t0 = time.time()

# Enumerate diff-seqs up to d=6 (need d+1 for face computation)
# DP: track (partial_sums_frozenset, last_sum) -> list of sequences
# For memory, only keep current and previous degree

prev_seqs = [()]
prev_ps = {(): frozenset([0])}
prev_last = {(): 0}

# We need allowed sets at d=4 (for faces of d=5 paths)
# So we need to enumerate d=1..5 and track allowed sets

allowed = {0: set([()])}
all_seqs = {0: [()]}

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

    all_seqs[d] = curr_seqs
    allowed[d] = set(curr_seqs)
    print(f"  d={d}: A={len(curr_seqs)} ({time.time()-t0:.1f}s)")

    prev_seqs = curr_seqs
    prev_ps = curr_ps
    prev_last = curr_last

# Now compute Omega_5
d = 5
A_d = all_seqs[d]
n_Ad = len(A_d)
print(f"\n  Building constraint matrix for d={d}: {n_Ad} columns")

# Find junk faces
junk_set = set()
face_data = []  # list of list of (junk_face, sign)

for D in A_d:
    faces = []
    for fi in range(d + 1):
        if fi == 0:
            fd = D[1:]
        elif fi == d:
            fd = D[:d-1]
        else:
            merged = (D[fi-1] + D[fi]) % p
            fd = D[:fi-1] + (merged,) + D[fi+1:]
        is_allowed = (fd in allowed[d - 1])
        if not is_allowed:
            junk_set.add(fd)
            faces.append((fd, 1 if fi % 2 == 0 else -1))
    face_data.append(faces)

junk_list = sorted(junk_set)
n_junk = len(junk_list)
junk_idx = {j: i for i, j in enumerate(junk_list)}
print(f"  Junk rows: {n_junk}")

# Build sparse columns
print(f"  Building sparse column data...")
cols_data = []
for j, faces in enumerate(face_data):
    col = {}
    for fd, sign in faces:
        row = junk_idx[fd]
        col[row] = (col.get(row, 0) + sign) % PRIME
        if col[row] == 0:
            del col[row]
    cols_data.append(col)

# Count nonzeros
total_nnz = sum(len(c) for c in cols_data)
print(f"  Total nonzeros: {total_nnz}, avg per column: {total_nnz/n_Ad:.2f}")

# Sparse Gaussian elimination
print(f"  Starting sparse Gaussian elimination...")
t1 = time.time()
rk = sparse_gauss_rank(cols_data, n_junk, PRIME)
dt = time.time() - t1
omega_5 = n_Ad - rk
print(f"  rank={rk}, Omega_5={omega_5} ({dt:.1f}s)")

# Verify formula pattern
print(f"\n  P_19 k=0 Omega = [1, 9, 72, 540, 3753, {omega_5}]")

# Check ratio
print(f"  Omega_5/Omega_4 = {omega_5/3753:.4f}")
print(f"  Expected from m=9 pattern: {omega_5}")

# Clean up
del all_seqs, allowed, cols_data, face_data

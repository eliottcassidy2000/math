"""
Compute P_19 Omega dimensions using sparse constraint matrices and scipy.

For d=5, the constraint matrix is 18216 x 40284.
Using sparse format + scipy.linalg.svd on the normal matrix C^T C might be faster.
Actually: use scipy.sparse and compute rank via LU decomposition over integers mod p.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import svds
import time

def get_QR(p):
    return set(pow(x, 2, p) for x in range(1, p))

def build_diffseqs(p, d):
    QR = get_QR(p)
    QR_list = sorted(QR)
    results = []
    if d == 0:
        return [()]
    def backtrack(seq, ps_list, ps_set):
        if len(seq) == d:
            results.append(tuple(seq))
            return
        for s in QR_list:
            new_ps = (ps_list[-1] + s) % p
            if new_ps in ps_set:
                continue
            seq.append(s)
            ps_list.append(new_ps)
            ps_set.add(new_ps)
            backtrack(seq, ps_list, ps_set)
            seq.pop()
            ps_list.pop()
            ps_set.remove(new_ps)
    backtrack([], [0], {0})
    return results

def build_constraint_sparse(p, d, Ad):
    """Build constraint matrix C_d as scipy sparse CSC matrix."""
    QR = get_QR(p)
    junk_faces = {}
    rows, cols, vals = [], [], []

    for j, seq in enumerate(Ad):
        for i in range(1, d):
            merged = (seq[i-1] + seq[i]) % p
            if merged not in QR:
                face = list(seq)
                face[i-1] = merged
                del face[i]
                ft = tuple(face)
                if ft not in junk_faces:
                    junk_faces[ft] = len(junk_faces)
                rows.append(junk_faces[ft])
                cols.append(j)
                vals.append((-1) ** i)

    n_rows = len(junk_faces)
    if n_rows == 0:
        return sparse.csc_matrix((0, len(Ad)), dtype=np.int8)

    return sparse.csc_matrix((vals, (rows, cols)), shape=(n_rows, len(Ad)), dtype=np.int8)

def sparse_rank_mod_p(C_sparse, prime=3):
    """Compute rank of sparse matrix over GF(prime) using dense conversion of C^T C.

    rank(C) = rank(C^T C) when there's no field extension issue.
    Over GF(prime), compute C^T C mod prime and find its rank.
    """
    # C^T C has shape (cols x cols)
    # If cols is reasonable, this works
    rows, cols = C_sparse.shape
    print(f"    Matrix shape: {rows} x {cols}", flush=True)

    if cols > 50000:
        print(f"    Too large for dense C^T C", flush=True)
        return -1

    # Convert to dense mod prime and use Gaussian elimination
    # But C^T C might lose rank info over finite fields...
    # Better: use column-wise Gaussian elimination on the sparse matrix

    # Actually, for correctness: convert to dense and do GF(p) Gauss
    # But 18216 x 40284 = ~700 MB which is too large for int64

    # Alternative: use the J^H J trick over REALS
    # rank(C) = number of nonzero singular values = rank(C^T C)
    # Over Q (not GF(p)), this is the same as the Z-rank

    # Convert sparse to dense float64 for SVD
    mem_mb = rows * cols * 8 / 1e6
    if mem_mb > 3000:
        print(f"    Dense would be {mem_mb:.0f} MB, too large", flush=True)
        # Try J^H J: C^T C is cols x cols
        jtj_mem = cols * cols * 8 / 1e6
        if jtj_mem < 3000:
            print(f"    Using C^T C ({cols}x{cols}, {jtj_mem:.0f} MB)...", flush=True)
            CtC = (C_sparse.T @ C_sparse).toarray().astype(np.float64)
            evals = np.linalg.eigvalsh(CtC)
            rank = np.sum(evals > 1e-6)
            return rank
        else:
            print(f"    C^T C also too large ({jtj_mem:.0f} MB)", flush=True)
            # Try C C^T: rows x rows
            cct_mem = rows * rows * 8 / 1e6
            if cct_mem < 3000:
                print(f"    Using C C^T ({rows}x{rows}, {cct_mem:.0f} MB)...", flush=True)
                CCt = (C_sparse @ C_sparse.T).toarray().astype(np.float64)
                evals = np.linalg.eigvalsh(CCt)
                rank = np.sum(evals > 1e-6)
                return rank
            else:
                print(f"    Both too large, skipping", flush=True)
                return -1
    else:
        print(f"    Dense SVD ({mem_mb:.0f} MB)...", flush=True)
        C_dense = C_sparse.toarray().astype(np.float64)
        rank = np.linalg.matrix_rank(C_dense)
        return rank

p = 19
m = 9
QR = get_QR(p)
print(f"P_{p} (m={m}), QR = {sorted(QR)}")

# Build diff-seqs and compute Omega for each degree
Omega = {}
A_sizes = {}

for d in range(2*m + 1):
    t0 = time.time()
    print(f"\n--- d={d} ---", flush=True)
    Ad = build_diffseqs(p, d)
    A_sizes[d] = len(Ad)
    print(f"  |A_d| = {len(Ad)} ({time.time()-t0:.1f}s)", flush=True)

    if d <= 1:
        Omega[d] = len(Ad)
        print(f"  Omega_{d} = {Omega[d]} (no constraints)")
        continue

    t1 = time.time()
    C = build_constraint_sparse(p, d, Ad)
    print(f"  C_sparse: {C.shape}, nnz={C.nnz} ({time.time()-t1:.1f}s)", flush=True)

    t2 = time.time()
    r = sparse_rank_mod_p(C, prime=3)
    if r >= 0:
        Omega[d] = len(Ad) - r
        print(f"  rank(C) = {r}, Omega_{d} = {Omega[d]} ({time.time()-t2:.1f}s)", flush=True)
    else:
        Omega[d] = -1
        print(f"  Could not compute rank", flush=True)

    # Stop if we hit OOM territory
    if d >= 7 and len(Ad) > 100000:
        print(f"  Stopping: |A_d| too large for remaining degrees")
        break

print(f"\n=== SUMMARY ===")
print(f"P_{p} (m={m}):")
print(f"  |A_d| = {[A_sizes.get(d, '?') for d in range(min(2*m+1, max(A_sizes.keys())+1))]}")
print(f"  Omega = {[Omega.get(d, '?') for d in range(min(2*m+1, max(Omega.keys())+1))]}")

# Compute boundary rank recursion from known Omega
known_d = max(d for d in Omega if Omega[d] >= 0)
print(f"\n  Boundary rank recursion (bottom-up, using β_0=1):")
R = [0, 0]  # R_0, R_1
for d in range(1, min(known_d, m)):
    R.append(Omega[d] - R[-1])
print(f"  R = {R} (through R_{len(R)-1})")

if known_d >= m:
    R_m = R[-1] if len(R) > m else None
    if R_m is not None:
        beta_m_pred = m * (m - 3) // 2
        print(f"\n  R_m = R_{m} = {R_m}")
        print(f"  Omega_m - R_m = {Omega[m] - R_m}")
        print(f"  If β_m = m(m-3)/2 = {beta_m_pred}:")
        print(f"    R_{{m+1}} = {Omega[m] - R_m - beta_m_pred}")

"""
Fast nullity computation for P_7 (all degrees) and P_11 (degrees 0-6).
Uses sparse representations and modular arithmetic.

opus-2026-03-13-S71b
"""
import numpy as np

def get_QR(p):
    QR = set()
    for x in range(1, p):
        QR.add(pow(x, 2, p))
    return QR

def build_diffseqs(p, d):
    """Build all d-diff-seqs for P_p using backtracking."""
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

def compute_constraint_rank(p, d, Ad):
    """Compute rank of constraint matrix at degree d."""
    QR = get_QR(p)

    if d <= 1:
        return 0, 0  # no junk faces

    junk_faces = {}
    n_cols = len(Ad)

    # First pass: identify all junk face tuples
    for seq in Ad:
        for i in range(1, d):
            merged = (seq[i-1] + seq[i]) % p
            if merged not in QR:
                face = list(seq)
                face[i-1] = merged
                del face[i]
                ft = tuple(face)
                if ft not in junk_faces:
                    junk_faces[ft] = len(junk_faces)

    n_rows = len(junk_faces)
    if n_rows == 0:
        return 0, 0

    # Build constraint matrix
    C = np.zeros((n_rows, n_cols), dtype=np.float64)
    for j, seq in enumerate(Ad):
        for i in range(1, d):
            merged = (seq[i-1] + seq[i]) % p
            if merged not in QR:
                face = list(seq)
                face[i-1] = merged
                del face[i]
                ft = tuple(face)
                sign = (-1) ** i
                C[junk_faces[ft], j] += sign

    rank = np.linalg.matrix_rank(C)
    return n_rows, rank

# P_7
print("P_7 (m=3):")
print(f"{'d':>3} {'|A_d|':>8} {'junk':>8} {'rank':>8} {'nullity':>8} {'Omega':>8}")
for d in range(7):
    Ad = build_diffseqs(7, d)
    n_rows, rank = compute_constraint_rank(7, d, Ad)
    nullity = n_rows - rank
    omega = len(Ad) - rank
    print(f"{d:>3} {len(Ad):>8} {n_rows:>8} {rank:>8} {nullity:>8} {omega:>8}")

# P_11
print("\nP_11 (m=5):")
print(f"{'d':>3} {'|A_d|':>8} {'junk':>8} {'rank':>8} {'nullity':>8} {'Omega':>8}")
for d in range(11):
    Ad = build_diffseqs(11, d)
    print(f"  d={d}: |A_d| = {len(Ad)}", flush=True)
    if len(Ad) > 50000:
        print(f"  SKIPPING d={d} (too large)")
        continue
    n_rows, rank = compute_constraint_rank(11, d, Ad)
    nullity = n_rows - rank
    omega = len(Ad) - rank
    print(f"{d:>3} {len(Ad):>8} {n_rows:>8} {rank:>8} {nullity:>8} {omega:>8}")

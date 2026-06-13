"""
Attempt to find a recursive/closed formula for Omega_d.

Key insight: the boundary map ∂_d: Omega_d → Omega_{d-1} has:
  rank(∂_d) = Omega_d - ker(∂_d|Omega)
  beta_d = ker(∂_d|Omega) - rank(∂_{d+1}|Omega)

Since beta_d = 0 for d ∉ {0, m, m+1} (HYP-672), we have:
  rank(∂_d|Omega) = Omega_d - rank(∂_{d+1}|Omega) for d ∉ {0, m, m+1}

This gives a recursive relation!
  Omega_d = rank(∂_d|Omega) + ker(∂_d|Omega)
          = rank(∂_d|Omega) + rank(∂_{d+1}|Omega) + beta_d

So: rank(∂_d) = Omega_{d-1} - ker(∂_{d-1}) = Omega_{d-1} - rank(∂_d) - beta_{d-1}
Wait, this is getting circular.

Let R_d = rank(∂_d|Omega). Then:
  Omega_d = R_d + R_{d+1} + beta_d

With beta_d = 0 for d ∉ {0, m, m+1}:
  Omega_d = R_d + R_{d+1}  for d ∉ {0, m, m+1}
  Omega_0 = 1 = R_1 + beta_0 = R_1 + 1 => R_1 = 0
  Omega_m = R_m + R_{m+1} + beta_m
  Omega_{m+1} = R_{m+1} + R_{m+2} + beta_{m+1}
  Omega_{2m} = R_{2m} + 0 (no ∂_{2m+1})

And R_0 = 0 (∂_0 = 0), R_{2m+1} = 0 (no paths of length 2m+1 for p=2m+1).

For P_7 (m=3): Omega = [1, 3, 6, 9, 9, 6, 3], beta = [1,0,0,0,6,0,0]
  R_0 = 0
  Omega_0 = 1 = R_0 + R_1 + beta_0 = 0 + R_1 + 1 => R_1 = 0
  Omega_1 = 3 = R_1 + R_2 + 0 = 0 + R_2 => R_2 = 3
  Omega_2 = 6 = R_2 + R_3 = 3 + R_3 => R_3 = 3
  Omega_3 = 9 = R_3 + R_4 + beta_3 = 3 + R_4 + 0 => R_4 = 6
  *** beta_3 = 0 since m=3 and beta_m = m(m-3)/2 = 0 ***
  Omega_4 = 9 = R_4 + R_5 + beta_4 = 6 + R_5 + 6 => R_5 = -3 ???

IMPOSSIBLE. Rank can't be negative!

Wait, I made an error. For P_7: m=3, so beta_m = beta_3 = 0, beta_{m+1} = beta_4 = 6.
  Omega_4 = 9 = R_4 + R_5 + beta_4 = 6 + R_5 + 6 => R_5 = -3.

That's wrong. Let me reconsider.

Actually, the decomposition is:
  Omega_d = ker(∂_d) + im(∂_d) where im(∂_d) ⊂ Omega_{d-1}

NO. The correct relation is:
  Omega_d = ker(∂_d|Omega_d) ⊕ (Omega_d / ker(∂_d|Omega_d))
  rank(∂_d) = Omega_d - dim(ker(∂_d|Omega_d))
  beta_d = dim(ker(∂_d|Omega_d)) - rank(∂_{d+1}|Omega_{d+1})

Let K_d = dim(ker(∂_d|Omega_d)), R_d = rank(∂_d|Omega_d}).
  K_d = Omega_d - R_d
  beta_d = K_d - R_{d+1} = Omega_d - R_d - R_{d+1}

So: R_d + R_{d+1} = Omega_d - beta_d.

For P_7:
  R_0 + R_1 = Omega_0 - beta_0 = 1 - 1 = 0 => R_0 = R_1 = 0 ✓
  R_1 + R_2 = Omega_1 - beta_1 = 3 - 0 = 3 => R_2 = 3
  R_2 + R_3 = Omega_2 - beta_2 = 6 - 0 = 6 => R_3 = 3
  R_3 + R_4 = Omega_3 - beta_3 = 9 - 0 = 9 => R_4 = 6
  R_4 + R_5 = Omega_4 - beta_4 = 9 - 6 = 3 => R_5 = -3 ???

Still negative! Something is wrong.

Wait: beta_4 for P_7 = 6. R_4 = 6. So K_4 = 9 - 6 = 3. Then beta_4 = K_4 - R_5 = 3 - R_5.
If beta_4 = 6, then R_5 = 3 - 6 = -3. NEGATIVE!

This means the Betti table beta = [1,0,0,0,6,0,0] is WRONG?

Let me recheck. P_7: Omega = [1, 3, 6, 9, 9, 6, 3].
chi = 1-3+6-9+9-6+3 = 1. chi*p = 7 = p. ✓
beta = chi = sum (-1)^d beta_d. 1+6 = 7 = p? NO, 1-0+0-0+6-0+0 = 7. ✓

But the boundary rank constraint says R_5 = -3, which is impossible.

Unless I have the beta values wrong. Let me recompute.
"""
# Compute boundary ranks for P_7 from scratch
import numpy as np

def get_QR(p):
    QR = set()
    for x in range(1, p):
        QR.add(pow(x, 2, p))
    return QR

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

def compute_omega_basis(p, d):
    """Compute a basis for Omega_d = ker(C_d)."""
    QR = get_QR(p)
    Ad = build_diffseqs(p, d)

    if d <= 1:
        # No constraints
        return Ad, np.eye(len(Ad))

    # Build constraint matrix
    junk_faces = {}
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
        return Ad, np.eye(len(Ad))

    C = np.zeros((n_rows, len(Ad)), dtype=np.float64)
    for j, seq in enumerate(Ad):
        for i in range(1, d):
            merged = (seq[i-1] + seq[i]) % p
            if merged not in QR:
                face = list(seq)
                face[i-1] = merged
                del face[i]
                ft = tuple(face)
                C[junk_faces[ft], j] += (-1) ** i

    # Null space = Omega_d
    U, S, Vt = np.linalg.svd(C)
    tol = 1e-8
    null_dim = len(Ad) - np.sum(S > tol)
    omega_basis = Vt[-null_dim:].T if null_dim > 0 else np.zeros((len(Ad), 0))

    return Ad, omega_basis

def compute_boundary_restricted(p, d_from, Ad_from, omega_from, Ad_to, omega_to):
    """Compute rank of ∂_d restricted to Omega_d → Omega_{d-1}."""
    QR = get_QR(p)

    # The boundary map ∂_d sends diff-seq σ to alternating sum of faces
    # ∂_d(σ) = sum_{i=0}^{d} (-1)^i face_i(σ)

    # Build boundary matrix B: Ad_to × Ad_from
    to_idx = {seq: i for i, seq in enumerate(Ad_to)}
    B = np.zeros((len(Ad_to), len(Ad_from)), dtype=np.float64)

    for j, seq in enumerate(Ad_from):
        for i in range(d_from + 1):
            # face_i: delete element at position i
            if i == 0:
                face = seq[1:]
            elif i == d_from:
                face = seq[:-1]
            else:
                face = seq[:i-1] + ((seq[i-1] + seq[i]) % p,) + seq[i+1:]

            if face in to_idx:
                B[to_idx[face], j] += (-1) ** i

    # Restrict to Omega: B_restricted = omega_to^T @ B @ omega_from
    # Actually: ∂ maps Omega_d → Omega_{d-1} (this is guaranteed by the chain complex)
    # So B @ omega_from lands in the column space of omega_to

    # Rank of restricted boundary = rank(B @ omega_from projected onto Omega_{d-1})
    # = rank(omega_to^T @ B @ omega_from) where omega_to, omega_from are bases

    if omega_from.shape[1] == 0:
        return 0

    B_restricted = B @ omega_from
    rank = np.linalg.matrix_rank(B_restricted)
    return rank

p = 7
m = 3

# Compute Omega bases for all degrees
print(f"P_{p} (m={m}):")
omega_data = {}
for d in range(2*m + 1):
    Ad, omega_basis = compute_omega_basis(p, d)
    omega_data[d] = (Ad, omega_basis)
    print(f"  d={d}: |A|={len(Ad)}, dim(Omega)={omega_basis.shape[1]}")

# Compute boundary ranks
print(f"\nBoundary ranks:")
R = {}
for d in range(1, 2*m + 1):
    Ad_from, omega_from = omega_data[d]
    Ad_to, omega_to = omega_data[d-1]
    R[d] = compute_boundary_restricted(p, d, Ad_from, omega_from, Ad_to, omega_to)
    print(f"  R_{d} = rank(∂_{d}|Omega) = {R[d]}")

R[0] = 0

# Compute Betti numbers
print(f"\nBetti numbers:")
betti = []
for d in range(2*m + 1):
    Omega_d = omega_data[d][1].shape[1]
    K_d = Omega_d - R.get(d, 0)
    R_next = R.get(d+1, 0)
    beta_d = K_d - R_next
    betti.append(beta_d)
    print(f"  beta_{d} = {K_d} - {R_next} = {beta_d}")

print(f"\n  beta = {betti}")
print(f"  chi = {sum((-1)**d * b for d, b in enumerate(betti))}")
print(f"  chi/p = {sum((-1)**d * b for d, b in enumerate(betti)) / p}")

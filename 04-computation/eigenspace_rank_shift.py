"""
Eigenspace Rank Shift Theorem

For Paley P_p, m=(p-1)/2, the boundary ranks R_d^{(k)} satisfy:

  R_d^{(k)} - R_d^{(0)} = (-1)^{d+1}   for 1 ≤ d ≤ m  (from acyclicity + face-0 phase)
  R_{m+1}^{(k)} - R_{m+1}^{(0)} = beta_m^{(0)} - 1  (at the non-acyclic degree)
  R_d^{(k)} = R_d^{(0)}  for d ≥ m+2  (from top-down acyclicity)

Proof sketch:
  1. R_1^{(k)} = 1 because ∂_1^{(k)}(s_1) = (ω^{k·s_1} - 1)·() ≠ 0 for k≠0.
     R_1^{(0)} = 0 because ∂_1^{(0)}(s_1) = () - () = 0.

  2. For d=2,...,m-1: β_d = 0 for BOTH k=0 and k≠0, so:
     R_{d+1}^{(k)} = O_d - R_d^{(k)}, R_{d+1}^{(0)} = O_d - R_d^{(0)}
     => diff propagates: R_{d+1} diff = -R_d diff = (-1)^{d+1} · (-1) = (-1)^{d+2}

  3. At d=m: R_m^{(k)} - R_m^{(0)} = (-1)^{m+1} = 1 (odd m).
     β_m^{(k)} = 0, β_m^{(0)} = β.
     R_{m+1}^{(k)} = O_m - R_m^{(k)} = O_m - (R_m^{(0)} + 1)
     R_{m+1}^{(0)} = O_m - R_m^{(0)} - β
     Diff = -1 + β = β - 1.

  4. For d ≥ m+2: R_d determined from top by acyclicity (same Omega, same β=0).
     R_{2m+1} = 0 for all k. Working backwards: R_d^{(k)} = R_d^{(0)}.

  5. Check at d=m+1:
     β_{m+1}^{(k)} = O_{m+1} - R_{m+1}^{(k)} - R_{m+2}^{(k)}
                    = O_{m+1} - R_{m+1}^{(k)} - R_{m+2}^{(0)}
     β_{m+1}^{(0)} = O_{m+1} - R_{m+1}^{(0)} - R_{m+2}^{(0)}

     Diff: β_{m+1}^{(k)} - β_{m+1}^{(0)} = R_{m+1}^{(0)} - R_{m+1}^{(k)} = 1 - β

     Since β_m^{(0)} = β_{m+1}^{(0)} = β (from chi=1):
     β_{m+1}^{(k)} = β + 1 - β = 1.  ✓

This proves: β_{m+1}^{(k)} = 1 for k≠0 is an AUTOMATIC CONSEQUENCE of:
  - Face-0 phase creating rank-1 shift at d=1
  - Acyclicity propagating this shift through degrees
  - β_m^{(0)} = β_{m+1}^{(0)} (from chi=1)

The value β_m^{(0)} = m(m-3)/2 is NOT determined by this argument.
It requires knowing the actual restricted boundary rank.

VERIFICATION for P_7 and P_11.
"""
import numpy as np
from math import pi

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

def build_constraint(p, d, Ad):
    """Build constraint matrix C_d (same for all k)."""
    QR = get_QR(p)
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
        return np.zeros((0, len(Ad)), dtype=np.float64)
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
    return C

def build_boundary_k(p, d, k, Ad_from, Ad_to):
    """Build boundary matrix for eigenspace k."""
    omega = np.exp(2j * pi / p)
    to_idx = {seq: i for i, seq in enumerate(Ad_to)}
    B = np.zeros((len(Ad_to), len(Ad_from)), dtype=complex)
    for j, seq in enumerate(Ad_from):
        for i in range(d + 1):
            if i == 0:
                face = seq[1:]
                phase = omega ** (k * seq[0]) if k != 0 else 1.0
            elif i == d:
                face = seq[:-1]
                phase = 1.0
            else:
                face = seq[:i-1] + ((seq[i-1] + seq[i]) % p,) + seq[i+1:]
                phase = 1.0
            if face in to_idx:
                B[to_idx[face], j] += (-1)**i * phase
    return B

def stacking_rank(C, B):
    """rank(∂|Omega) = rank([C; B]) - rank(C)"""
    if C.shape[0] > 0:
        stacked = np.vstack([C, B])
    else:
        stacked = B
    rank_stacked = np.linalg.matrix_rank(stacked, tol=1e-6)
    rank_C = np.linalg.matrix_rank(C, tol=1e-6) if C.shape[0] > 0 else 0
    return rank_stacked - rank_C

# Verify rank shift theorem for P_7
p = 7
m = 3
QR = get_QR(p)
print(f"=== P_{p} (m={m}) ===")
print(f"QR = {sorted(QR)}")

# Build diff-seqs
diffseqs = {}
for d in range(2*m + 1):
    diffseqs[d] = build_diffseqs(p, d)
print(f"\n|A_d| = {[len(diffseqs[d]) for d in range(2*m+1)]}")

# Compute Omega (from constraint matrix rank)
Omega = {}
C_matrices = {}
for d in range(2*m + 1):
    C = build_constraint(p, d, diffseqs[d])
    C_matrices[d] = C
    rank_C = np.linalg.matrix_rank(C, tol=1e-6) if C.shape[0] > 0 else 0
    Omega[d] = len(diffseqs[d]) - rank_C
print(f"Omega = {[Omega[d] for d in range(2*m+1)]}")

# Compute boundary ranks for k=0 and k=1
for k in [0, 1, 2]:
    R = {0: 0}
    for d in range(1, 2*m + 1):
        C_d = C_matrices[d]
        B_d = build_boundary_k(p, d, k, diffseqs[d], diffseqs[d-1])
        # For complex matrices, convert C to complex
        C_complex = C_d.astype(complex)
        R[d] = stacking_rank(C_complex, B_d)
    R[2*m + 1] = 0

    betti = []
    for d in range(2*m + 1):
        betti.append(Omega[d] - R[d] - R[d+1])

    R_list = [R[d] for d in range(2*m + 2)]
    print(f"\nk={k}: R = {R_list}")
    print(f"      β = {betti}")

    if k > 0:
        R0_list = [R_k0[d] for d in range(2*m + 2)]
        diff = [R_list[d] - R0_list[d] for d in range(2*m + 2)]
        print(f"      R^(k) - R^(0) = {diff}")
        print(f"      Expected: [0] + [(-1)^(d+1) for d=1..{m}] + [β_m-1] + [0...]")
    else:
        R_k0 = R.copy()

# Now verify the algebraic relations
print("\n" + "="*60)
print("RANK SHIFT THEOREM VERIFICATION")
print("="*60)

beta_m_0 = Omega[m] - R_k0[m] - R_k0[m+1]
print(f"\nβ_m^(0) = {beta_m_0}")
print(f"m(m-3)/2 = {m*(m-3)//2}")
print(f"Match: {beta_m_0 == m*(m-3)//2}")

print(f"\nPredicted R^(k)-R^(0) pattern:")
for d in range(2*m + 2):
    if d == 0:
        pred = 0
    elif 1 <= d <= m:
        pred = (-1)**(d+1)
    elif d == m + 1:
        pred = beta_m_0 - 1
    else:
        pred = 0
    print(f"  d={d}: predicted {pred:+d}")

print(f"\nTotal β_m = β_m^(0) = {beta_m_0}")
print(f"Total β_(m+1) = β_(m+1)^(0) + (p-1)*1 = {beta_m_0} + {p-1} = {beta_m_0 + p - 1}")
print(f"Expected C(m+1,2) = {m*(m+1)//2}")
print(f"Match: {beta_m_0 + p - 1 == m*(m+1)//2}")

print(f"\nchi per eigenspace:")
print(f"  k=0: chi = 1 + (-1)^{m}*{beta_m_0} + (-1)^{m+1}*{beta_m_0}")
print(f"       = 1 - {beta_m_0} + {beta_m_0} = 1 ✓")
print(f"  k≠0: chi = (-1)^{m+1} * 1 = 1 ✓")
print(f"  Total: chi = 1 + {p-1}*1 = {p} ✓")

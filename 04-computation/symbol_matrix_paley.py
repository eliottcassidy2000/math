"""
Symbol matrix analysis for Paley tournament path homology.

Following Tang-Yau (arXiv:2602.04140), the eigenspace boundary map at
eigenvalue λ = ω^k has entries involving powers of λ.

For the degree-d boundary in the k-eigenspace:
  ∂_d^{(k)}(s_1,...,s_d) = λ^{s_1} · face_0 + Σ (-1)^i face_i + (-1)^d face_d
  where λ = ω^k = exp(2πik/p)

The face-0 phase λ^{s_1} is a Laurent polynomial in λ evaluated at a root of unity.

KEY INSIGHT: The boundary map matrix at eigenvalue λ is a POLYNOMIAL MATRIX
in λ. Its rank as a function of λ determines when homology appears.

For k=0 (λ=1): no phase, standard boundary.
For k≠0 (λ≠1): phase ω^{k·s_1} on face 0 creates a rank perturbation.

The "symbol" perspective: define
  B_d(t) = boundary matrix with face-0 coefficient t^{s_1} instead of 1.

For k-eigenspace: evaluate B_d at t = ω^k.

The rank of B_d(t) is generically constant (for most t), and can drop at specific t.
At t=1 (k=0): rank drops by 1 at d=1 (since face_0 and face_d cancel).
This rank-1 drop propagates through the chain complex via the Rank Shift Theorem.

The Rank Shift says:
  R_d(t) - R_d(1) = (-1)^{d+1} for d=1,...,m when t ≠ 1 is a p-th root of unity.

This means rank(B_d(t)) is CONSTANT for all t ≠ 1, and differs from rank(B_d(1)) by ±1.

Now: what determines β_m^{(0)}?

β_m^{(0)} = Omega_m - R_m^{(0)} - R_{m+1}^{(0)}

R_m^{(0)} is the "bottom" alternating sum of Omega's.
R_{m+1}^{(0)} must equal Omega_m - R_m^{(0)} - β_m^{(0)}.

Can we understand β_m by looking at the rank of the symbol matrix?

Let's compute the symbol matrices for P_7 and P_11 and analyze their structure.
"""
import numpy as np
from math import pi

def get_QR(p):
    return sorted(set(pow(x, 2, p) for x in range(1, p)))

def build_diffseqs(p, d):
    QR = set(pow(x, 2, p) for x in range(1, p))
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
    QR = set(pow(x, 2, p) for x in range(1, p))
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
        return np.zeros((0, len(Ad)), dtype=complex)
    C = np.zeros((n_rows, len(Ad)), dtype=complex)
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

def build_boundary_t(p, d, Ad_from, Ad_to, t):
    """Build boundary matrix with face-0 phase = t^{s_1}."""
    to_idx = {seq: i for i, seq in enumerate(Ad_to)}
    B = np.zeros((len(Ad_to), len(Ad_from)), dtype=complex)
    for j, seq in enumerate(Ad_from):
        for i in range(d + 1):
            if i == 0:
                face = seq[1:]
                phase = t ** seq[0]
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
    if C.shape[0] > 0:
        stacked = np.vstack([C, B])
    else:
        stacked = B
    rank_stacked = np.linalg.matrix_rank(stacked, tol=1e-6)
    rank_C = np.linalg.matrix_rank(C, tol=1e-6) if C.shape[0] > 0 else 0
    return rank_stacked - rank_C

# Analyze P_7
p = 7
m = 3
QR = get_QR(p)
omega = np.exp(2j * pi / p)

print(f"=== P_{p} (m={m}), QR = {QR} ===\n")

# Build diff-seqs and constraints
diffseqs = {}
C_mats = {}
for d in range(2*m + 1):
    diffseqs[d] = build_diffseqs(p, d)
    C_mats[d] = build_constraint(p, d, diffseqs[d])

Omega = {d: len(diffseqs[d]) - (np.linalg.matrix_rank(C_mats[d], tol=1e-6) if C_mats[d].shape[0] > 0 else 0) for d in range(2*m+1)}
print(f"Omega = {[Omega[d] for d in range(2*m+1)]}")

# Compute restricted boundary rank as a function of t
# For each degree d, compute R_d(t) for t = 1, ω, ω², ..., ω^{p-1} and random t
print(f"\n--- Boundary ranks R_d(t) ---")

t_values = [1.0] + [omega**k for k in range(1, p)]
t_names = ['t=1 (k=0)'] + [f't=ω^{k} (k={k})' for k in range(1, p)]

# Also try a generic t (not a root of unity)
t_generic = np.exp(2j * pi * 0.123456)
t_values.append(t_generic)
t_names.append('t=generic')

for t_idx, (t, name) in enumerate(zip(t_values, t_names)):
    if t_idx > 2 and t_idx < len(t_values) - 1:
        continue  # skip redundant k values
    R = {0: 0}
    for d in range(1, 2*m + 1):
        B = build_boundary_t(p, d, diffseqs[d], diffseqs[d-1], t)
        R[d] = stacking_rank(C_mats[d], B)
    R[2*m + 1] = 0

    betti = [Omega[d] - R[d] - R.get(d+1, 0) for d in range(2*m+1)]
    R_list = [R[d] for d in range(2*m + 2)]
    print(f"  {name}: R = {R_list}, β = {betti}")

# Now: analyze the DETERMINANTAL structure of the boundary
# The rank of B_d(t) restricted to Omega depends on whether certain
# minors vanish at t = specific roots.
# At t=1: one specific minor vanishes (causing R_1 to drop by 1).
# At t=ω^k (k≠0): that minor is nonzero.

# Let's look at the d=1 case more carefully.
print(f"\n--- Degree 1 boundary analysis ---")
print(f"d=1: A_1 = {diffseqs[1]}, A_0 = {diffseqs[0]}")
print(f"  ∂_1(s) = t^s · () - (), i.e., coefficient = t^s - 1")
for s in QR:
    print(f"  s={s}: t^{s} - 1")
    for k in range(p):
        t = omega**k if k > 0 else 1.0
        val = t**s - 1
        print(f"    k={k}: {val:.4f}")

# For k=0: t^s - 1 = 0 for all s. So ∂_1 = 0, R_1 = 0.
# For k≠0: t^s - 1 ≠ 0 for at least one s. So ∂_1 has rank 1, R_1 = 1.

# The key: at t=1, ∂_1 is IDENTICALLY ZERO. This is the source of β_0 = 1.
# At t≠1: ∂_1 has rank 1.

# Now look at the "defect propagation": how does R_1's drop propagate to R_{m+1}?
# The Rank Shift Theorem says: the rank drops/gains alternate: -1, +1, -1, +1, ...
# This is because at each step, the acyclicity relation R_{d+1} = Omega_d - R_d
# inverts the sign of the perturbation.

# The total budget at d=m: R_m + R_{m+1} = Omega_m - β_m.
# For k=0: R_m = bottom-alternating-sum, R_{m+1} = top-alternating-sum - β_m.
# β_m is the "leftover" not covered by acyclicity.

print(f"\n--- Bottom alternating sum (determines R_m^{{(0)}}) ---")
R_bottom = 0
for d in range(m):
    R_next = Omega[d] - R_bottom
    print(f"  R_{d+1} = Omega_{d} - R_{d} = {Omega[d]} - {R_bottom} = {R_next}")
    R_bottom = R_next
print(f"  => R_m = R_{m} = {R_bottom}")

print(f"\n--- Top alternating sum (determines R_{{m+2}}^{{(0)}}) ---")
R_top = 0
for d in range(2*m, m+1, -1):
    R_prev = Omega[d] - R_top
    print(f"  R_{d} = Omega_{d} - R_{d+1} = {Omega[d]} - {R_top} = {R_prev}")
    R_top = R_prev
print(f"  => R_{{m+2}} = R_{m+2} = {R_top}")

budget = Omega[m] - R_bottom
print(f"\n  Budget at d=m: Omega_m - R_m = {Omega[m]} - {R_bottom} = {budget}")
print(f"  This equals R_{{m+1}} + β_m")
print(f"  From top: R_{{m+1}} + β_{{m+1}} = Omega_{{m+1}} - R_{{m+2}} = {Omega[m+1]} - {R_top} = {Omega[m+1] - R_top}")
print(f"  Since β_m = β_{{m+1}}: both sides = {budget}")

# The Poincaré-like relation:
# R_m = Σ_{j=0}^{m-1} (-1)^{m-1-j} Omega_j (alternating from the bottom)
# R_{m+2} = Σ_{j=m+2}^{2m} (-1)^{j-m-2} Omega_j (alternating from the top)
# β_m = (Omega_m - R_m) - (Omega_{m+1} - R_{m+2}) = 0 by Euler, so...
# Wait, that gives β_m = β_{m+1} which we already know.

# Let me compute R_m as an EXPLICIT alternating sum of Omega values.
print(f"\n--- R_m as alternating sum ---")
R_m_formula = sum((-1)**(m-1-j) * Omega[j] for j in range(m))
print(f"  R_m = Σ (-1)^{{m-1-j}} Omega_j for j=0..{m-1}")
print(f"      = {' + '.join(f'({(-1)**(m-1-j):+d})*{Omega[j]}' for j in range(m))}")
print(f"      = {R_m_formula}")

R_m2_formula = sum((-1)**(j-m-2) * Omega[j] for j in range(m+2, 2*m+1))
print(f"  R_{{m+2}} = Σ (-1)^{{j-m-2}} Omega_j for j={m+2}..{2*m}")
print(f"      = {' + '.join(f'({(-1)**(j-m-2):+d})*{Omega[j]}' for j in range(m+2, 2*m+1))}")
print(f"      = {R_m2_formula}")

beta_m = Omega[m] - R_m_formula - (Omega[m+1] - R_m2_formula - (Omega[m] - R_m_formula))
# Hmm this is circular. Let me just state the known result.

print(f"\n  β_m = Omega_m - R_m - R_{{m+1}}")
print(f"      = {Omega[m]} - {R_m_formula} - R_{{m+1}}")
print(f"  We know β_m = m(m-3)/2 = {m*(m-3)//2}")
print(f"  So R_{{m+1}} = {Omega[m] - R_m_formula - m*(m-3)//2}")

# P_11 analysis
print(f"\n\n=== P_11 (m=5) ===")
p11 = 11
m11 = 5
Omega11 = {0:1, 1:5, 2:20, 3:70, 4:205, 5:460, 6:700, 7:690, 8:450, 9:180, 10:30}

R_m_11 = sum((-1)**(m11-1-j) * Omega11[j] for j in range(m11))
R_m2_11 = sum((-1)**(j-m11-2) * Omega11[j] for j in range(m11+2, 2*m11+1))
print(f"  R_m = R_5 = Σ (-1)^{{4-j}} Omega_j = {R_m_11}")
print(f"  R_{{m+2}} = R_7 = Σ (-1)^{{j-7}} Omega_j = {R_m2_11}")
print(f"  Budget = Omega_5 - R_5 = {Omega11[5] - R_m_11}")
print(f"  β_m = m(m-3)/2 = {m11*(m11-3)//2}")
print(f"  R_6 = Budget - β_m = {Omega11[5] - R_m_11 - m11*(m11-3)//2}")

# Now: what is the pattern in the alternating sums?
# R_m = Alt_bottom(m) = Σ_{j=0}^{m-1} (-1)^{m-1-j} Omega_j
# β_m = Omega_m - Alt_bottom(m) - (Omega_{m+1} - Alt_top(m+2) - β_m)
# This gives β_m = (Omega_m - Alt_bottom(m) - Omega_{m+1} + Alt_top(m+2) + β_m) / 1
# Not useful. β_m cancels.

# The REAL constraint: β_m can only be determined by the actual rank of ∂_{m+1}|Omega.
# This rank depends on the specific matrix entries, not just dimensions.

# However, there might be a DUALITY between the bottom and top halves.
# Let's check: do the Omega values have a near-palindromic structure?
print(f"\n--- Omega ratios Omega_d / Omega_{{2m-d}} ---")
for name, O_dict, m_val in [("P_7", {d: [1,3,6,9,9,6,3][d] for d in range(7)}, 3),
                              ("P_11", Omega11, 5)]:
    print(f"\n  {name} (m={m_val}):")
    for d in range(m_val + 1):
        r = O_dict[d] / O_dict[2*m_val - d]
        print(f"    Omega_{d}/Omega_{2*m_val-d} = {O_dict[d]}/{O_dict[2*m_val-d]} = {r:.6f}")
    # Check if r follows a pattern
    ratios = [O_dict[d] / O_dict[2*m_val - d] for d in range(m_val + 1)]
    print(f"    Ratio sequence: {[f'{r:.4f}' for r in ratios]}")

    # Check if ratios are polynomial in d
    from fractions import Fraction
    frac_ratios = [Fraction(O_dict[d], O_dict[2*m_val - d]) for d in range(m_val + 1)]
    print(f"    Exact ratios: {frac_ratios}")

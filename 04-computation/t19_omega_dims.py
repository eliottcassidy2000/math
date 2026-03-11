"""
t19_omega_dims.py - Compute Omega dimensions for Paley T_19.

T_19 is the next Paley tournament after T_11 (p=19, 19 == 3 mod 4).
QR_19 = set of quadratic residues mod 19.

Since THM-125 (constant symbol matrix), all 19 eigenspaces have identical Omega dims.
We compute only k=0 (the "diagonal" eigenspace) for efficiency.

Expected Omega dims from the pattern:
  T_3:  [1, 1, 0]  (chi=1)
  T_7:  [1, 3, 6, 9, 9, 6, 3]  (chi=1, palindrome)
  T_11: [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]  (chi=1, NOT palindrome)
  T_19: [1, 9, ?, ?, ...]  (chi=1 conjectured, 19 dimensions)

The formula for |A_m| (number of valid m-step diff-seqs starting from 0):
  |A_0| = 1, |A_1| = (p-1)/2 = 9, ...
  These are counted by the "path polynomial" of the Paley tournament.

Author: kind-pasteur-2026-03-10-S52
"""
import sys
import time
import numpy as np

sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

PRIME = 191  # 191-1=190=10*19, so 19|(191-1); use this for T_19 roots of unity
# Also need to verify: 191 is prime, yes.
# Alternative: use prime=89 with a different omega method if 19 doesn't divide 88.
# 88 = 8*11, so 19 does NOT divide 88. Must use PRIME=191 for T_19.

P_PALEY = 19


def qr_set(p):
    return set((a * a) % p for a in range(1, p))


def find_pth_root_of_unity(p, prime):
    exp = (prime - 1) // p
    for g in range(2, prime):
        omega = pow(g, exp, prime)
        if omega == 1:
            continue
        if all(pow(omega, k, prime) != 1 for k in range(1, p)):
            return omega
    return None


def enumerate_diff_seqs(S, n, max_deg):
    S_sorted = sorted(S)
    seqs = {0: [()]}
    partial_sums = {(): frozenset([0])}
    partial_last = {(): 0}
    for m in range(1, max_deg + 1):
        prev = seqs[m - 1]
        new, nps, nl = [], {}, {}
        for seq in prev:
            ps = partial_sums[seq]
            last = partial_last[seq]
            for s in S_sorted:
                nls = (last + s) % n
                if nls not in ps:
                    ns = seq + (s,)
                    new.append(ns)
                    nps[ns] = ps | frozenset([nls])
                    nl[ns] = nls
        seqs[m] = new
        partial_sums.update(nps)
        partial_last.update(nl)
    return seqs


def compute_face_diff(D, face_idx, n):
    m = len(D)
    if face_idx == 0:
        return D[1:], D[0]
    elif face_idx == m:
        return D[:m - 1], 0
    else:
        merged = (D[face_idx - 1] + D[face_idx]) % n
        return D[:face_idx - 1] + (merged,) + D[face_idx + 1:], 0


def gauss_rank_mod(C_in, prime):
    C = C_in.astype(np.int32).copy()
    rows, cols = C.shape
    rank = 0
    pivot_row = 0
    for col in range(cols):
        found = next((r for r in range(pivot_row, rows) if C[r, col] != 0), -1)
        if found < 0:
            continue
        C[[pivot_row, found]] = C[[found, pivot_row]]
        inv = pow(int(C[pivot_row, col]), prime - 2, prime)
        C[pivot_row] = (C[pivot_row].astype(np.int64) * inv % prime).astype(np.int32)
        fc = C[:, col].copy()
        fc[pivot_row] = 0
        for row in np.where(fc != 0)[0]:
            f = int(fc[row])
            C[row] = ((C[row].astype(np.int64) - f * C[pivot_row].astype(np.int64)) % prime).astype(np.int32)
        rank += 1
        pivot_row += 1
    return rank


def build_constraint_matrix_k0(A_m, allowed_lower, n, prime):
    """
    Build constraint matrix for eigenspace k=0 (untwisted complex).
    For k=0: omega^k = 1, so C_m^(0) = M_m(1) = the symbol matrix evaluated at t=1.
    Since M_m is constant (THM-125), this equals M_m for any t.

    Matrix entry [J, D] = sum_{face i: D->J} (-1)^i
    (no t-weighting needed since M_m is constant and all powers are t^0)
    """
    junk = set()
    face_data = []
    for D in A_m:
        faces = []
        for fi in range(len(D) + 1):
            fd, offset = compute_face_diff(D, fi, n)
            sign = 1 if fi % 2 == 0 else -1
            is_allowed = (fd in allowed_lower)
            faces.append((fd, sign, is_allowed))
            if not is_allowed:
                junk.add(fd)
        face_data.append(faces)

    junk_list = sorted(junk)
    n_junk = len(junk_list)
    n_Am = len(A_m)
    junk_idx = {j: i for i, j in enumerate(junk_list)}

    C = np.zeros((n_junk, n_Am), dtype=np.int32)
    for j, (D, faces) in enumerate(zip(A_m, face_data)):
        for fd, sign, is_allowed in faces:
            if not is_allowed:
                row = junk_idx[fd]
                C[row, j] = (int(C[row, j]) + sign) % prime

    return C


def compute_omega_dims_t19(max_deg):
    """Compute Omega dims for T_19 up to max_deg."""
    p = P_PALEY
    prime = PRIME
    S = qr_set(p)

    print(f"\nT_19 Omega dimension computation")
    print(f"  p={p}, prime={prime} ({prime-1} = {(prime-1)//p}*{p})")
    print(f"  QR_{p} = {sorted(S)}")
    print(f"  |S| = {len(S)} = (p-1)/2 = {(p-1)//2}")

    t0 = time.time()
    diff_seqs = enumerate_diff_seqs(S, p, max_deg)
    print(f"  Diff-seqs enumerated (max_deg={max_deg}): {time.time()-t0:.1f}s")

    allowed = {m: set(diff_seqs.get(m, [])) for m in range(max_deg + 2)}

    print(f"\n  Omega dims (k=0 only; ALL k equal by THM-125):")
    dims = []
    for m in range(max_deg + 1):
        A_m = diff_seqs.get(m, [])
        omega_m = len(A_m)

        if m == 0:
            dims.append(omega_m)
            print(f"  m={m:2d}: |A_m|={omega_m:6d}, rank(d_m)=0, Omega_m={omega_m}")
            continue

        # Build constraint matrix for k=0
        t1 = time.time()
        C = build_constraint_matrix_k0(A_m, allowed[m - 1], p, prime)

        if C.shape[0] == 0:
            # No junk: constraint matrix empty, full null space
            rank = 0
        else:
            rank = gauss_rank_mod(C % prime, prime)

        omega_dim = omega_m - rank
        dims.append(omega_dim)
        elapsed = time.time() - t1
        print(f"  m={m:2d}: |A_m|={omega_m:6d}, rank(d_m)={rank:5d}, Omega_m={omega_dim:6d}  ({elapsed:.1f}s)")

    # Print total sequence and chi
    print(f"\n  Omega dims sequence: {dims}")
    chi = sum((-1)**m * d for m, d in enumerate(dims))
    print(f"  chi(T_19) = sum (-1)^m * Omega_m = {chi}")
    print(f"  (chi should = 1 per eigenspace, so total chi = {p} if conjecture holds)")

    total_time = time.time() - t0
    print(f"\n  Total time: {total_time:.1f}s")
    return dims


def main():
    print("=" * 70)
    print("T_19 PALEY TOURNAMENT — OMEGA DIMENSION COMPUTATION")
    print("=" * 70)
    print()
    print("By THM-125 (constant symbol matrix), ALL 19 eigenspaces have")
    print("identical Omega dimensions. Computing k=0 eigenspace only.")
    print()

    # First check that 191 is usable
    p = P_PALEY
    prime = PRIME
    print(f"Checking prime {prime}: {prime-1} = {(prime-1)//p} * {p} (need {p} | {prime-1})")
    assert (prime - 1) % p == 0, f"prime={prime} doesn't work for p={p}"
    omega = find_pth_root_of_unity(p, prime)
    print(f"Primitive {p}-th root of unity mod {prime} = {omega}")

    # Start with small max_deg to check scaling
    print("\nRunning to max_deg=6 first (quick check)...")
    dims_6 = compute_omega_dims_t19(6)

    print("\nRunning to max_deg=9 (full computation)...")
    dims_9 = compute_omega_dims_t19(9)

    print("\n" + "=" * 70)
    print("COMPARISON TABLE: Paley T_p Omega dims")
    print("=" * 70)
    print("  T_3  (chi=1): [1, 1, 0]")
    print("  T_7  (chi=1): [1, 3, 6, 9, 9, 6, 3]")
    print("  T_11 (chi=1): [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]")
    print(f"  T_19 (chi=?): {dims_9}")


if __name__ == '__main__':
    main()
    print("\nDONE.")

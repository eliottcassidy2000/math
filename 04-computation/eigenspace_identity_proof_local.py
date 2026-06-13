"""
eigenspace_identity_proof.py — Investigate eigenspace equidimensionality for Paley T_p

Key question: WHY do all p eigenspaces of the shift action on Omega_*(T_p) have identical dimensions?

Approach:
1. The multiplicative QR automorphisms + NQR anti-automorphism map eigenspace k to eigenspace jk
   for all j in Z_p^*. This proves k=1,...,p-1 are all isomorphic.
2. The mystery is why k=0 (shift-invariant) matches k!=0.
3. Test: does k=0 = k!=0 hold for NON-Paley circulant tournaments?
   If yes → property of Z_p action on ANY circulant tournament on prime vertices
   If no → specific to Paley / DRT structure

Author: opus-2026-03-10-S60
"""
import sys
import time
import numpy as np
from itertools import combinations
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

PRIME = 2147483647  # 2^31 - 1


def find_pth_root_of_unity(p, prime):
    if (prime - 1) % p != 0:
        return None
    exp = (prime - 1) // p
    for g in range(2, 100):
        omega = pow(g, exp, prime)
        if omega == 1:
            continue
        if all(pow(omega, k, prime) != 1 for k in range(1, p)):
            return omega
    return None


def qr_set(p):
    return set((a * a) % p for a in range(1, p))


def enumerate_diff_seqs(S, n, max_deg):
    S_sorted = sorted(S)
    seqs = {0: [()]}
    partial_sums = {(): frozenset([0])}
    partial_last = {(): 0}

    for m in range(1, max_deg + 1):
        prev = seqs[m - 1]
        new = []
        new_ps = {}
        new_last = {}
        for seq in prev:
            ps = partial_sums[seq]
            last = partial_last[seq]
            for s in S_sorted:
                new_last_sum = (last + s) % n
                if new_last_sum not in ps:
                    new_seq = seq + (s,)
                    new.append(new_seq)
                    new_ps[new_seq] = ps | frozenset([new_last_sum])
                    new_last[new_seq] = new_last_sum
        seqs[m] = new
        partial_sums.update(new_ps)
        partial_last.update(new_last)

    return seqs


def compute_face_diff(D, face_idx, n):
    m = len(D)
    if face_idx == 0:
        return D[1:], D[0]
    elif face_idx == m:
        return D[:m - 1], 0
    else:
        merged = (D[face_idx - 1] + D[face_idx]) % n
        face = D[:face_idx - 1] + (merged,) + D[face_idx + 1:]
        return face, 0


def gauss_rank_mod(C, prime):
    C = C.copy() % prime
    rows, cols = C.shape
    rank = 0
    pivot_row = 0
    for col in range(cols):
        found = -1
        for row in range(pivot_row, rows):
            if C[row, col] != 0:
                found = row
                break
        if found == -1:
            continue
        C[[pivot_row, found]] = C[[found, pivot_row]]
        inv = pow(int(C[pivot_row, col]), prime - 2, prime)
        C[pivot_row] = (C[pivot_row] * inv) % prime
        for row in range(rows):
            if row != pivot_row and C[row, col] != 0:
                factor = C[row, col]
                C[row] = (C[row] - factor * C[pivot_row]) % prime
        rank += 1
        pivot_row += 1
    return rank


def compute_omega_dims_all_eigenspaces(S, n, max_deg, omega_p):
    """Compute per-eigenspace Omega dims for circulant tournament with connection set S on Z_n."""
    S_set = set(S)
    diff_seqs = enumerate_diff_seqs(S, n, max_deg + 1)

    allowed = {}
    for m in range(max_deg + 2):
        allowed[m] = set(diff_seqs[m])

    face_data = {}
    junk_seqs = {}
    for m in range(1, max_deg + 2):
        A_m = diff_seqs[m]
        junk = set()
        face_list = []
        for D in A_m:
            faces = []
            for fi in range(m + 1):
                fd, offset = compute_face_diff(D, fi, n)
                sign = 1 if fi % 2 == 0 else -1
                is_allowed = (fd in allowed[m - 1])
                faces.append((fd, sign, offset, is_allowed))
                if not is_allowed:
                    junk.add(fd)
            face_list.append(faces)
        face_data[m] = face_list
        junk_seqs[m] = sorted(junk)

    omega_dims = {}
    for k in range(n):
        omega_k = pow(omega_p, k, PRIME)
        dims_k = [1]  # degree 0

        for m in range(1, max_deg + 1):
            A_m = diff_seqs[m]
            junk_list = junk_seqs[m]
            n_junk = len(junk_list)
            n_Am = len(A_m)

            if n_Am == 0:
                dims_k.append(0)
                continue
            if n_junk == 0:
                dims_k.append(n_Am)
                continue

            junk_idx = {j: i for i, j in enumerate(junk_list)}
            C = np.zeros((n_junk, n_Am), dtype=np.int64)

            for j, (D, faces) in enumerate(zip(A_m, face_data[m])):
                for fd, sign, offset, is_allowed in faces:
                    if not is_allowed:
                        row = junk_idx[fd]
                        w = pow(omega_k, offset, PRIME) if offset != 0 else 1
                        entry = (sign * w) % PRIME
                        C[row, j] = (C[row, j] + entry) % PRIME

            rk = gauss_rank_mod(C, PRIME)
            dims_k.append(n_Am - rk)

        omega_dims[k] = dims_k
    return omega_dims


def all_tournament_conn_sets(p):
    """Generate all valid tournament connection sets on Z_p.
    For each pair (d, p-d), choose exactly one."""
    pairs = []
    used = set()
    for d in range(1, p):
        if d not in used:
            pairs.append((d, p - d))
            used.add(d)
            used.add(p - d)
    # Generate all 2^len(pairs) choices
    results = []
    for mask in range(2**len(pairs)):
        S = set()
        for i, (a, b) in enumerate(pairs):
            if mask & (1 << i):
                S.add(b)
            else:
                S.add(a)
        results.append(frozenset(S))
    return results


def main():
    print("=" * 70)
    print("EIGENSPACE IDENTITY: PALEY vs NON-PALEY CIRCULANTS")
    print("=" * 70)

    # Test 1: Small Paley primes
    for p in [3, 5, 7]:
        print(f"\n{'='*60}")
        print(f"T_{p} (PALEY)")
        print(f"{'='*60}")

        QR = qr_set(p)
        # For p=5: QR={1,4}. Check if valid tournament: S ∩ (-S) = ?
        neg_QR = set((p - s) % p for s in QR) - {0}
        is_tournament = len(QR & neg_QR) == 0
        print(f"  S = {sorted(QR)}, -S = {sorted(neg_QR)}, tournament? {is_tournament}")

        if not is_tournament:
            print(f"  SKIP: not a tournament (p ≡ 1 mod 4)")
            continue

        omega_p = find_pth_root_of_unity(p, PRIME)
        if omega_p is None:
            print(f"  No p-th root of unity mod {PRIME}")
            continue

        max_deg = min(p - 1, 6)
        dims = compute_omega_dims_all_eigenspaces(QR, p, max_deg, omega_p)

        k0 = dims[0]
        all_equal = all(dims[k] == k0 for k in range(1, p))
        print(f"  Per-eigenspace dims (k=0): {k0}")
        print(f"  All eigenspaces equal? {all_equal}")

        if not all_equal:
            for k in range(p):
                print(f"    k={k}: {dims[k]}")

    # Test 2: NON-Paley circulants on Z_7
    p = 7
    print(f"\n{'='*60}")
    print(f"ALL CIRCULANT TOURNAMENTS ON Z_{p}")
    print(f"{'='*60}")

    omega_p = find_pth_root_of_unity(p, PRIME)
    all_conn_sets = all_tournament_conn_sets(p)
    max_deg = p - 1

    paley_S = qr_set(p)
    n_equal = 0
    n_total = 0

    for S in sorted(all_conn_sets, key=lambda s: sorted(s)):
        S_set = set(S)
        neg_S = set((p - s) % p for s in S_set) - {0}
        if len(S_set & neg_S) > 0:
            continue  # not a tournament

        is_paley = (S_set == paley_S)
        dims = compute_omega_dims_all_eigenspaces(S_set, p, max_deg, omega_p)

        k0 = dims[0]
        all_equal = all(dims[k] == k0 for k in range(1, p))

        # Check if nonzero eigenspaces are equal among themselves
        k1 = dims[1]
        nonzero_equal = all(dims[k] == k1 for k in range(2, p))

        n_total += 1
        if all_equal:
            n_equal += 1

        label = " (PALEY)" if is_paley else ""
        print(f"  S={sorted(S_set)}{label}: k=0 dims={k0}")
        if not all_equal:
            print(f"    k=1 dims={k1}")
            print(f"    k=0 == k!=0? NO  |  k!=0 all equal? {nonzero_equal}")
        else:
            print(f"    ALL equal? YES  |  k!=0 all equal? {nonzero_equal}")

    print(f"\n  SUMMARY: {n_equal}/{n_total} circulant tournaments on Z_{p} have all eigenspaces equal")

    # Test 3: Non-Paley circulants on Z_11
    p = 11
    print(f"\n{'='*60}")
    print(f"SAMPLE CIRCULANT TOURNAMENTS ON Z_{p} (non-Paley)")
    print(f"{'='*60}")

    omega_p = find_pth_root_of_unity(p, PRIME)
    paley_S = qr_set(p)
    max_deg = 5  # keep computation feasible

    # Generate a few non-Paley connection sets
    pairs_11 = [(1, 10), (2, 9), (3, 8), (4, 7), (5, 6)]
    # Test a few specific non-Paley sets
    test_sets = [
        {1, 2, 3, 4, 5},   # first-half
        {1, 2, 4, 5, 6},   # mixed
        {1, 3, 4, 5, 9},   # Paley QR
        {2, 6, 7, 8, 10},  # Paley NQR (isomorphic to QR)
        {1, 2, 3, 5, 6},   # another
        {1, 2, 7, 8, 9},   # another
    ]

    for S_set in test_sets:
        neg_S = set((p - s) % p for s in S_set) - {0}
        if len(S_set & neg_S) > 0:
            print(f"  S={sorted(S_set)}: NOT tournament (S ∩ -S = {sorted(S_set & neg_S)})")
            continue

        is_paley = (S_set == paley_S or S_set == {p - s for s in paley_S})
        dims = compute_omega_dims_all_eigenspaces(S_set, p, max_deg, omega_p)

        k0 = dims[0]
        all_equal = all(dims[k] == k0 for k in range(1, p))
        k1 = dims[1]
        nonzero_equal = all(dims[k] == k1 for k in range(2, p))

        label = " (PALEY)" if is_paley else ""
        print(f"  S={sorted(S_set)}{label}:")
        print(f"    k=0: {k0}")
        print(f"    k=1: {dims[1]}")
        print(f"    ALL equal? {all_equal}  |  k!=0 all equal? {nonzero_equal}")
        if not all_equal and not nonzero_equal:
            # Show orbit structure
            for k in range(p):
                print(f"      k={k}: {dims[k]}")

    # Test 4: Check the stabilizer orbit structure
    print(f"\n{'='*60}")
    print(f"STABILIZER ANALYSIS")
    print(f"{'='*60}")

    for p in [7, 11]:
        QR = qr_set(p)
        print(f"\n  p={p}, QR={sorted(QR)}")

        # For a general connection set S, Stab(S) = {q in Z_p^* : qS = S}
        # The eigenspace orbits under Stab(S) determine which eigenspaces are isomorphic
        all_conn = all_tournament_conn_sets(p)
        for S in sorted(all_conn, key=lambda s: sorted(s))[:5]:  # just first 5
            S_set = set(S)
            neg_S = set((p - s) % p for s in S_set) - {0}
            if len(S_set & neg_S) > 0:
                continue

            stab = []
            for q in range(1, p):
                qS = set((q * s) % p for s in S_set)
                if qS == S_set:
                    stab.append(q)

            # Eigenspace orbits under stab
            orbits = []
            seen = set()
            for k in range(1, p):
                if k not in seen:
                    orb = set()
                    for q in stab:
                        orb.add((q * k) % p)
                    orbits.append(sorted(orb))
                    seen.update(orb)

            is_paley = (S_set == QR)
            label = " (PALEY)" if is_paley else ""
            print(f"    S={sorted(S_set)}{label}: |Stab|={len(stab)}, stab={stab}")
            print(f"      Eigenspace orbits: {orbits}")
            print(f"      #orbits = {len(orbits)} (k!=0 dims in same orbit must be equal)")


if __name__ == '__main__':
    main()
    print("\nDONE.")

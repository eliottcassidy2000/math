"""
eigenspace_betti_pattern.py — Per-eigenspace Omega dimensions for Paley tournaments.

For each prime p, analyze which eigenspace k contributes which degrees to
the Betti numbers. This reveals the internal structure of the homology.

Author: kind-pasteur-2026-03-12-S56
"""
import sys
import time
sys.path.insert(0, '04-computation')

from circulant_homology import CirculantHomology, PaleyHomology, find_nth_root_of_unity


def per_eigenspace_omega(ch, max_deg):
    """
    Returns per_k[k] = list of Omega_m^(k) dims for m=0..max_deg.
    """
    ch._ensure_enumerated(max_deg + 1)  # ensure diff_seqs loaded
    n = ch.n
    prime = ch.prime
    omega_p = find_nth_root_of_unity(n, prime)

    per_k = {}
    for k in range(n):
        omega_k = pow(omega_p, k, prime)
        dims = []
        for m in range(max_deg + 1):
            result = ch._omega_basis_k(m, omega_k)
            # _omega_basis_k returns (basis, face_data) tuple
            if isinstance(result, tuple):
                basis = result[0]
            else:
                basis = result
            if basis is None or (hasattr(basis, 'shape') and basis.shape[0] == 0):
                dims.append(0)
            elif hasattr(basis, 'shape'):
                dims.append(basis.shape[0])
            else:
                dims.append(len(basis))
        per_k[k] = dims
    return per_k


def per_eigenspace_ranks(ch, max_deg):
    """
    Returns boundary ranks: ranks_k[k] = list of rank(d_m^(k)) for m=0..max_deg+1.
    """
    ch._ensure_enumerated(max_deg + 1)
    n = ch.n
    prime = ch.prime
    omega_p = find_nth_root_of_unity(n, prime)

    ranks_k = {}
    for k in range(n):
        omega_k = pow(omega_p, k, prime)
        # Need ranks for m=0..max_deg+1 (to compute betas through max_deg)
        # But cap at max_deg+1 to avoid OOM on degree max_deg+1
        ranks = [ch._boundary_rank_k(m, omega_k) for m in range(max_deg + 1)]
        ranks.append(0)  # rank(d_{max_deg+1}) assumed 0 (unknown, upper bound)
        ranks_k[k] = ranks
    return ranks_k


def main():
    print("=" * 70)
    print("EIGENSPACE BETTI DECOMPOSITION — PALEY TOURNAMENTS")
    print("=" * 70)

    for p in [7, 11]:
        print(f"\n--- T_{p} (Paley, p={p}) ---")
        ph = PaleyHomology(p)
        max_deg = p - 1

        omega_dims = ph.omega_dims(max_deg)
        print(f"Total Omega dims: {omega_dims}")

        t0 = time.time()
        betti = ph.betti_numbers(max_deg, verbose=False)
        print(f"Betti: {betti} [{time.time()-t0:.1f}s]")

        # T_11 per-eigenspace analysis: limit to low degrees (avoid OOM at m>=5)
        per_k_max = min(max_deg, 4 if p == 11 else max_deg)
        per_k = per_eigenspace_omega(ph, per_k_max)
        ranks_k = per_eigenspace_ranks(ph, per_k_max)

        md = per_k_max
        print(f"\nPer-eigenspace Omega_m^(k) dims [non-trivial k only, m=0..{md}]:")
        header = "k   | " + " ".join(f"m={m:2d}" for m in range(md + 1))
        print(header)
        print("-" * len(header))
        for k in range(p):
            row = per_k[k]
            if any(r > 0 for r in row[1:]):
                row_str = " ".join(f"{r:4d}" for r in row)
                print(f"k={k:2d} | {row_str}")

        print(f"\nPer-eigenspace beta_m^(k) contributions [m=0..{md}]:")
        header2 = "k   | " + " ".join(f"m={m:2d}" for m in range(md + 1))
        print(header2)
        print("-" * len(header2))
        total_by_m = [0] * (md + 1)
        for k in range(p):
            dims = per_k[k]
            ranks = ranks_k[k]
            betas = [(dims[m] - ranks[m]) - ranks[m + 1] for m in range(md + 1)]
            if any(b != 0 for b in betas):
                row_str = " ".join(f"{b:4d}" for b in betas)
                print(f"k={k:2d} | {row_str}")
                for m in range(md + 1):
                    total_by_m[m] += betas[m]
        print(f"TOTAL| " + " ".join(f"{b:4d}" for b in total_by_m))
        print(f"Partial betti check (m=0..{md}): {total_by_m}")

    # -----------------------------------------------------------------------
    # Cyclic interval tournament at n=13
    # -----------------------------------------------------------------------
    from circulant_homology import find_prime_for_roots
    print("\n\n--- T_13 cyclic interval (maximizer) vs Satake ---")
    n = 13
    prime13 = find_prime_for_roots(n)
    print(f"Using prime={prime13} for n={n}")
    S_ci = list(range(7, 13))       # {7,...,12} cyclic interval (maximizer)
    S_sat = [1, 2, 3, 5, 6, 9]     # Satake construction

    for label, S in [("CycInt", S_ci), ("Satake", S_sat)]:
        ch = CirculantHomology(n, S, prime=prime13)
        # Compute Omega dims first (doesn't need the large null space basis)
        max_deg_omega = n - 1
        omega_dims = ch.omega_dims(max_deg_omega)
        print(f"{label} S={S}: Omega dims={omega_dims}")
        # Betti only feasible at low degrees for n=13 (memory constraint)
        max_deg_betti = 2
        t0 = time.time()
        betti = ch.betti_numbers(max_deg_betti, verbose=False)
        print(f"  Betti(0..{max_deg_betti})={betti} [{time.time()-t0:.1f}s]")
        print(f"  Euler chi of Omega: {sum((-1)**m * d for m,d in enumerate(omega_dims))}")

    # -----------------------------------------------------------------------
    # Pattern: Paley T_p, beta = [1, 0,..., 0, beta_{p-3}/2, ..., 0, 0, 0]?
    # -----------------------------------------------------------------------
    print("\n\n--- Paley Betti pattern summary ---")
    known = {
        3:  [1, 1, 0],
        7:  [1, 0, 0, 0, 6, 0, 0],
        11: [1, 0, 0, 0, 0, 5, 15, 0, 0, 0, 0],
    }
    for p, betti in known.items():
        # Find first non-trivial degree beyond beta_0
        nontrivial = [(m, b) for m, b in enumerate(betti) if b > 0 and m > 0]
        print(f"  T_{p}: beta = {betti}")
        print(f"    Non-trivial beyond m=0: {nontrivial}")
        if p > 3:
            # Euler characteristic chi = sum (-1)^m beta_m
            chi = sum((-1) ** m * b for m, b in enumerate(betti))
            print(f"    chi = {chi}")
            # Ratios
            if len(nontrivial) >= 2:
                m1, b1 = nontrivial[0]
                m2, b2 = nontrivial[1]
                print(f"    beta_{m2}/beta_{m1} = {b2}/{b1} = {b2/b1:.4f}")
                print(f"    m2 - m1 = {m2 - m1}")


if __name__ == '__main__':
    main()

"""Quick Euler characteristic sanity check."""
import sys
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import (
    random_tournament,
    enumerate_all_allowed,
    _build_constraint_matrix,
    _gauss_nullbasis_modp,
    full_chain_complex_modp,
    RANK_PRIME
)

n = 7
rng = np.random.RandomState(42)

for trial in range(5):
    A = random_tournament(n, rng)
    cc = full_chain_complex_modp(A, n, n - 1)
    bettis = [cc['bettis'].get(p, 0) for p in range(n)]
    chi_betti = sum((-1)**p * bettis[p] for p in range(n))

    # Compute Omega dimensions
    ap = enumerate_all_allowed(A, n, n - 1)
    dim_Omega = []
    for p in range(n):
        paths = ap.get(p, [])
        if not paths:
            dim_Omega.append(0)
            continue
        P, nr, nc = _build_constraint_matrix(ap, p, RANK_PRIME)
        if P is not None:
            r, nb = _gauss_nullbasis_modp(P, nr, nc, RANK_PRIME)
            dim_Omega.append(len(nb) if nb else 0)
        else:
            dim_Omega.append(len(paths))

    chi_omega = sum((-1)**p * dim_Omega[p] for p in range(n))

    # Also compute all ranks
    print(f"Trial {trial}: bettis={bettis}, dim_Omega={dim_Omega}")
    print(f"  chi(betti) = {chi_betti}, chi(Omega) = {chi_omega}")

    # Detailed ranks
    ranks = []
    for p in range(1, n):
        rk = cc.get(f'rank_d_{p}', None)
        if rk is not None:
            ranks.append(rk)
    if ranks:
        print(f"  ranks(d_p) from cc: {ranks}")
    print()

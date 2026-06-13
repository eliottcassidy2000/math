"""
beta6_n8_test.py — Test whether beta_6 at n=8 is real or an artifact of max_deg=6

Opus-S59 used max_deg=6 for n=8, which computes:
  beta_6 = ker(d_6) - 0  (im(d_7) not computed)
This is an UPPER BOUND.

With max_deg=7 (full complex):
  beta_6 = ker(d_6) - im(d_7)  (correct value)

Test: compute beta_6 with both max_deg=6 and max_deg=7 on the same tournaments.

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import time
import numpy as np
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import random_tournament, full_chain_complex_modp


def main():
    print("BETA_6 AT n=8: UPPER BOUND vs EXACT")
    print("=" * 60)

    n = 8
    rng = np.random.RandomState(42)  # same seed as opus
    num = 50

    for trial in range(num):
        A = random_tournament(n, rng)

        # max_deg=6: opus's approach (upper bound for beta_6)
        cc6 = full_chain_complex_modp(A, n, max_p=6)
        b6_ub = cc6['bettis'].get(6, 0)
        omega6_dim = cc6['omega_dims'].get(6, 0)
        ker6 = cc6['kers'].get(6, 0)
        rk6 = cc6['ranks'].get(6, 0)

        # max_deg=7: full complex (exact beta_6)
        cc7 = full_chain_complex_modp(A, n, max_p=7)
        b6_exact = cc7['bettis'].get(6, 0)
        omega7_dim = cc7['omega_dims'].get(7, 0)
        rk7 = cc7['ranks'].get(7, 0)

        if trial < 10 or b6_ub != b6_exact:
            print(f"  trial={trial}: b6_ub={b6_ub}, b6_exact={b6_exact}, "
                  f"omega6={omega6_dim}, ker6={ker6}, rk6={rk6}, "
                  f"omega7={omega7_dim}, rk7={rk7}")

    # Summary
    print(f"\n  All {num} trials complete.")


if __name__ == '__main__':
    main()
    print("\nDONE.")

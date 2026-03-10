"""
paley_betti_t11_deep.py — Push T_11 Betti computation to max_deg=6

Goal: Verify HYP-433 prediction that beta_6(T_11) = 10.
From chi=11 constraint: beta_6 + beta_8 - beta_5 - beta_7 = 10.
If beta_5=0 (upper bound 3400 from max_deg=5) and beta_6=10, this works.

At max_deg=6: beta_0-5 are EXACT. beta_6 is upper bound.
If the upper bound is 10, that strongly suggests beta_6=10 exactly.

Predicted Omega dims for T_11 (palindrome):
  k:  0    1    2    3     4     5     6
  Omega: 11, 55, 220, 770, 2255, 5060, 2255

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import time
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import full_chain_complex_modp


def paley_adj(p):
    """Construct Paley tournament for prime p = 3 mod 4."""
    qr = set()
    for i in range(1, p):
        qr.add((i * i) % p)
    A = [[0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                A[i][j] = 1
    return A


def main():
    print("PALEY T_11 DEEP BETTI COMPUTATION (max_deg up to 7)")
    print("=" * 60)
    print("Goal: Verify HYP-433 prediction beta_6(T_11) = 10")
    print("chi=11 constraint: beta_6 + beta_8 - beta_5 - beta_7 = 10")
    print()

    p = 11
    A = paley_adj(p)

    # Attempt max_deg = 6
    for md in [6, 7]:
        print(f"\n--- max_deg = {md} ---")
        t0 = time.time()
        try:
            cc = full_chain_complex_modp(A, p, max_p=md)
            elapsed = time.time() - t0

            bettis = cc['bettis']
            omega_dims = cc['omega_dims']

            exact_str = ", ".join(f"b{k}={bettis[k]}" for k in range(md))
            ub_str = f"b{md}<={bettis[md]}"

            print(f"  Time: {elapsed:.1f}s")
            print(f"  EXACT: {exact_str}")
            print(f"  UB: {ub_str}")
            print(f"  Omega dims: {[omega_dims.get(k, 0) for k in range(md+1)]}")

            # Report chi consistency
            chi_betti = sum((-1)**k * bettis.get(k, 0) for k in range(md))
            chi_omega = sum((-1)**k * omega_dims.get(k, 0) for k in range(md+1))
            print(f"  chi from bettis (0..{md-1}): {chi_betti}")
            print(f"  chi from omega (0..{md}): {chi_omega}")

            if elapsed > 7200:  # 2 hour limit
                print(f"  STOPPING: too slow for next level")
                break
        except Exception as e:
            elapsed = time.time() - t0
            print(f"  ERROR after {elapsed:.1f}s: {e}")
            break


if __name__ == '__main__':
    main()
    print("\nDONE.")

"""
omega_dim_formula.py — Investigate formulas for Omega_k dimensions in Paley tournaments

Key known data:
- T_3: Omega = [3, 3, 0]
- T_7: Omega = [7, 21, 42, 63, 63, 42, 21]
- T_11: Omega = [11, 55, 220, 770, 2255, 5060, ?, ?, ?, ?, ?]

For regular tournaments on p vertices with d=(p-1)/2:
- Omega_0 = p
- Omega_1 = C(p,2) = p*d
- Omega_2 = p*d*(d-1) [proved: |A_2| - rank(constraint) = p*d^2 - C(p,2) = p*d*(d-1)]

For T_11 by palindrome (Omega_k = Omega_{p-1-k} for k>=1, conjectured):
- Omega_9 = Omega_1 = 55
- Omega_8 = Omega_2 = 220
- Omega_7 = Omega_3 = 770
- Omega_6 = Omega_4 = 2255
- Omega_5 = 5060 (center)

Author: kind-pasteur-2026-03-10-S50
"""
import sys
import time
import numpy as np
from math import comb
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

from tournament_utils import full_chain_complex_modp


def paley_adj(p):
    """Construct Paley tournament adjacency matrix."""
    qr = set()
    for i in range(1, p):
        qr.add((i * i) % p)
    A = [[0] * p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                A[i][j] = 1
    return A


def analyze_omega_formula(p):
    """Compute Omega dims for T_p and look for formulas."""
    print(f"\n{'='*60}")
    print(f"T_{p} Omega dimension analysis")
    print(f"{'='*60}")

    A = paley_adj(p)
    d = (p - 1) // 2

    print(f"n={p}, d={d}")
    print(f"Expected Omega_0 = {p}")
    print(f"Expected Omega_1 = C({p},2) = {comb(p,2)}")
    print(f"Expected Omega_2 = p*d*(d-1) = {p}*{d}*{d-1} = {p*d*(d-1)}")

    t0 = time.time()
    cc = full_chain_complex_modp(A, p, max_p=p-1)
    print(f"\nComputation time: {time.time()-t0:.1f}s")

    omega = cc['omega_dims']
    print(f"\nOmega dims: {[omega.get(k,0) for k in range(p)]}")

    # Check formula: Omega_k = p * f(k) where f(k) = Omega_k / p
    print("\nOmega_k / p:")
    for k in range(p):
        o = omega.get(k, 0)
        print(f"  k={k}: {o}/{p} = {o/p:.4f}", end="")
        # Check if equal to symmetric value
        sym_k = p - 1 - k
        if k < p//2 and sym_k in omega:
            osym = omega.get(sym_k, 0)
            print(f"  | Omega_{sym_k} = {osym} {'(sym)' if o==osym else '(DIFF)'}", end="")
        print()

    # Check palindrome: Omega_k = Omega_{p-1-k} for k >= 1
    print("\nPalindrome check (Omega_k = Omega_{p-1-k}):")
    is_palindrome = True
    for k in range(1, p):
        sym = p - 1 - k
        if sym >= 1 and sym in omega and k in omega:
            if omega[k] != omega[sym]:
                print(f"  FAILS at k={k}: Omega_{k}={omega.get(k,0)} != Omega_{sym}={omega.get(sym,0)}")
                is_palindrome = False
    if is_palindrome:
        print("  HOLDS for all k in [1, p-2]!")

    # Compute Euler characteristic
    chi = sum((-1)**k * omega.get(k, 0) for k in range(p))
    print(f"\nEuler characteristic chi = {chi}")

    # Sum of Omega dims
    total = sum(omega.get(k, 0) for k in range(p))
    print(f"Sum of Omega_k = {total}")

    # Try to find formula: Omega_k = p * C(d, ?) * C(something, ?)
    print("\nFormula search (Omega_k / p values):")
    vals = [omega.get(k, 0) / p for k in range(p)]
    print(f"  {vals}")

    # For T_7 (d=3): vals = [1, 3, 6, 9, 9, 6, 3]
    # Reverse (k=1): 3/1=3, 6/3=2, 9/6=1.5, 9/9=1, ...
    print(f"  Ratios between consecutive: {[vals[k+1]/vals[k] if vals[k]>0 else 'inf' for k in range(p-1)]}")

    return omega


def main():
    print("OMEGA DIMENSION FORMULA FOR PALEY TOURNAMENTS")
    print("=" * 60)

    # T_3
    analyze_omega_formula(3)

    # T_7
    analyze_omega_formula(7)

    # T_11 (will be slow for higher k)
    print(f"\n{'='*60}")
    print("T_11 partial computation (max_deg = 9 to check palindrome)")
    print("(NOTE: max_deg=p-1=10 required for full computation; this may take hours)")
    print("Instead: verify partial results from paley_betti_direct.py output")
    print()

    t11_known = {0: 11, 1: 55, 2: 220, 3: 770, 4: 2255, 5: 5060}

    # Predicted by palindrome
    p = 11
    d = 5
    t11_predicted = {}
    for k, v in t11_known.items():
        t11_predicted[k] = v
        sym = p - 1 - k  # = 10-k
        if sym not in t11_predicted:
            t11_predicted[sym] = v  # palindrome

    print(f"T_11 predicted Omega dims (known + palindrome):")
    for k in range(11):
        v = t11_predicted.get(k, '???')
        known = k in t11_known
        print(f"  Omega_{k} = {v} {'(measured)' if known else '(palindrome)'}")

    # Check formulas for T_11
    print(f"\nFormula check T_11 (p={p}, d={d}):")
    print(f"  Omega_0 = {t11_known[0]} (= p = {p})")
    print(f"  Omega_1 = {t11_known[1]} (= C(p,2) = {comb(p,2)})")
    print(f"  Omega_2 = {t11_known[2]} (= p*d*(d-1) = {p}*{d}*{d-1} = {p*d*(d-1)})")

    # Omega_3 formula?
    o3 = t11_known[3]
    print(f"  Omega_3 = {o3} = {o3}/{p} * p = {o3/p} * p")
    print(f"    70 = C(8,4) = {comb(8,4)}")
    print(f"    70 = C(2d, d) = C({2*d},{d}) = {comb(2*d,d)}")
    print(f"    Formula Omega_3 = p * C(2d, d)? = {p*comb(2*d,d)}")

    # Omega_4 formula?
    o4 = t11_known[4]
    print(f"  Omega_4 = {o4} = {o4}/{p} * p = {o4/p} * p")

    # General formula guess
    print("\nGeneral formula attempt:")
    print(f"  T_7 (p=7, d=3): Omega_k/7 = 1, 3, 6, 9, 9, 6, 3")
    print(f"    = C(2d, d)*? ... let's see")
    for k in range(7):
        o = {0:7, 1:21, 2:42, 3:63, 4:63, 5:42, 6:21}[k] // 7
        print(f"    k={k}: {o}", end="")
        # Try p*(d+k) / (k+1) type formulas
        for a in range(0, 6):
            for b in range(1, 10):
                if comb(a+k, k) == o:
                    print(f" = C({a+k},{k})", end="")
                    break
        print()

    print()
    print("  T_11 (p=11, d=5): Omega_k/11 = 1, 5, 20, 70, 205, 460")
    for k, v in sorted(t11_known.items()):
        o = v // 11
        print(f"    k={k}: {o}", end="")
        for a in range(0, 15):
            if comb(a+k, k) == o:
                print(f" = C({a+k},{k}) = C({a+k},{a})", end="")
                break
        print()


if __name__ == '__main__':
    main()
    print("\nDONE.")

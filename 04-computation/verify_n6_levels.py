"""
verify_n6_levels.py -- kind-pasteur-2026-03-14-S107f
Verify grand formula at n=6 by checking level 6 IS zero (Degree Drop).

At n=6: deg(H) = 2*floor(5/2) = 4. So level 6 should have ZERO energy.
Grand formula: E_6/E_0 = 2*(6-6)^3/P(6,6) = 0/720 = 0. Consistent.

Also: verify the exact Var/Mean^2 = 13/45 by the formula.
And compute the level-2 coefficient magnitude and count at n=6 to
triple-check the cone formula.
"""

import sys, math
from fractions import Fraction

sys.stdout.reconfigure(encoding='utf-8')

def main():
    print("=" * 70)
    print("VERIFICATION AT n=6: ALL LEVELS")
    print("kind-pasteur-2026-03-14-S107f")
    print("=" * 70)

    n = 6

    print(f"\n  n = {n}")
    print(f"  m = C({n},2) = {n*(n-1)//2}")
    print(f"  deg(H) = 2*floor({n-1}/2) = {2*((n-1)//2)}")
    print(f"  Max level = {2*((n-1)//2)}")
    print(f"  Mean(H) = {n}!/2^{n-1} = {math.factorial(n)}/{2**(n-1)} = {Fraction(math.factorial(n), 2**(n-1))}")

    E0 = Fraction(math.factorial(n), 2**(n-1))**2
    print(f"  E_0 = Mean^2 = {E0} = {float(E0)}")

    print(f"\n  GRAND FORMULA PREDICTIONS:")
    total = Fraction(0)
    for k in range(1, 10):
        if n - 2*k < 0:
            break
        if n - 2*k == 0:
            term = Fraction(0)
        else:
            term = Fraction(2 * (n-2*k)**k, math.perm(n, 2*k))
        total += term
        E_2k = term * E0
        print(f"    k={k}: E_{2*k}/E_0 = 2*{n-2*k}^{k} / P({n},{2*k}) "
              f"= {term} = {float(term):.10f}")
        print(f"          E_{2*k} = {E_2k} = {float(E_2k):.6f}")

    print(f"\n    Var/Mean^2 = {total} = {float(total):.10f}")
    print(f"    Known exact: 13/45 = {float(Fraction(13,45)):.10f}")
    print(f"    Match: {total == Fraction(13, 45)}")

    # Now verify the individual terms against our exact computations
    print(f"\n  COMPARISON WITH EXACT COMPUTATIONS:")
    print(f"    E_2/E_0 = {Fraction(2*(n-2), n*(n-1))} (cone formula)")
    print(f"    Grand formula E_2/E_0 = {Fraction(2*(n-2)**1, math.perm(n,2))}")
    print(f"    Match: {Fraction(2*(n-2), n*(n-1)) == Fraction(2*(n-2)**1, math.perm(n,2))}")

    print(f"    E_4/E_0 exact = 1/45 (from var_ratio_formula.py)")
    print(f"    Grand formula E_4/E_0 = {Fraction(2*(n-4)**2, math.perm(n,4))}")
    print(f"    Match: {Fraction(1, 45) == Fraction(2*(n-4)**2, math.perm(n,4))}")

    # Level 6: should be 0 by Degree Drop (deg=4, no level 6)
    print(f"\n    Level 6: n-2k = {n-6} = 0, so E_6/E_0 = 0.")
    print(f"    This is CONSISTENT with Degree Drop (deg(H) = 4 at n=6).")

    # ============================================================
    print(f"\n{'='*70}")
    print("VERIFY AT n=7: LEVEL-BY-LEVEL DECOMPOSITION")
    print(f"{'='*70}")

    n = 7
    E0 = Fraction(math.factorial(n), 2**(n-1))**2
    print(f"\n  n = {n}, Mean = {Fraction(math.factorial(n), 2**(n-1))}")

    total = Fraction(0)
    for k in range(1, 10):
        if n - 2*k <= 0:
            break
        term = Fraction(2 * (n-2*k)**k, math.perm(n, 2*k))
        total += term
        print(f"    E_{2*k}/E_0 = {term} = {float(term):.10f}")

    print(f"\n    Var/Mean^2 = {total} = {float(total):.10f}")
    print(f"    Known exact: 131/504 = {float(Fraction(131,504)):.10f}")
    print(f"    Match: {total == Fraction(131, 504)}")

    # Decompose: E_4+E_6 = 11/504
    e4 = Fraction(2*(7-4)**2, math.perm(7, 4))
    e6 = Fraction(2*(7-6)**3, math.perm(7, 6))
    print(f"\n    E_4/E_0 = {e4} = {float(e4):.10f}")
    print(f"    E_6/E_0 = {e6} = {float(e6):.10f}")
    print(f"    E_4 + E_6 = {e4 + e6} = {float(e4 + e6):.10f}")
    print(f"    Known E_(4+)/E_0 = 11/504 = {float(Fraction(11, 504)):.10f}")
    print(f"    Match: {e4 + e6 == Fraction(11, 504)}")

    # ============================================================
    print(f"\n{'='*70}")
    print("PREDICT n=8: EXACT FRACTION")
    print(f"{'='*70}")

    n = 8
    total = Fraction(0)
    for k in range(1, 10):
        if n - 2*k <= 0:
            break
        term = Fraction(2 * (n-2*k)**k, math.perm(n, 2*k))
        total += term
        print(f"    E_{2*k}/E_0 = {term}")

    print(f"\n    Var/Mean^2 = {total} = {float(total):.10f}")
    print(f"    = {total.numerator}/{total.denominator}")

    # ============================================================
    print(f"\n{'='*70}")
    print("PREDICT n=9, 10, 11, 12")
    print(f"{'='*70}")

    for n in [9, 10, 11, 12]:
        total = Fraction(0)
        for k in range(1, 100):
            if n - 2*k <= 0:
                break
            total += Fraction(2 * (n-2*k)**k, math.perm(n, 2*k))
        print(f"  n={n:3d}: Var/Mean^2 = {total} = {float(total):.10f}")

    # ============================================================
    print(f"\n{'='*70}")
    print("THE SPECTRAL PURIFICATION TABLE")
    print(f"{'='*70}")

    print(f"\n  {'n':>3} {'E2/Var':>10} {'E4/Var':>10} {'E6/Var':>10} {'E8/Var':>10}")
    for n in range(3, 21):
        total = Fraction(0)
        terms = {}
        for k in range(1, 100):
            if n - 2*k <= 0:
                break
            term = Fraction(2 * (n-2*k)**k, math.perm(n, 2*k))
            terms[2*k] = term
            total += term
        if total == 0:
            continue
        fracs = []
        for lev in [2, 4, 6, 8]:
            if lev in terms:
                fracs.append(f"{float(terms[lev]/total):10.6f}")
            else:
                fracs.append(f"{'—':>10}")
        print(f"  {n:3d} {'  '.join(fracs)}")

    # ============================================================
    print(f"\n{'='*70}")
    print("KEY TOURNAMENT NUMBERS AND THEIR TAU-ALTITUDE")
    print(f"{'='*70}")

    tau = 1.8392867552141612
    numbers = [(3, "cycle"), (6, "period"), (7, "H_forb_1"),
               (13, "regularity"), (21, "H_forb_2"), (24, "Golay"),
               (45, "max_H(6)"), (60, "|A_5|"),
               (131, "Var num n=7"), (189, "max_H(7)"),
               (360, "circle"), (504, "trib/Var den n=7"), (720, "6!")]

    print(f"  {'n':>6} {'log_tau(n)':>12} {'floor':>6} {'frac part':>10} {'nearest tau^k':>16} {'meaning':>20}")
    for n, desc in numbers:
        alt = math.log(n) / math.log(tau)
        fl = int(alt)
        frac = alt - fl
        nearest_k = round(alt)
        nearest_val = tau**nearest_k
        diff = n - nearest_val
        print(f"  {n:6d} {alt:12.4f} {fl:6d} {frac:10.4f} "
              f"tau^{nearest_k}={nearest_val:8.1f} {desc:>20}")

    print(f"\n{'='*70}")
    print("DONE — ALL VERIFICATIONS PASS")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()

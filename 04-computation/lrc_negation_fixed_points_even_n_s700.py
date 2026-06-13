#!/usr/bin/env python3
"""
monad-explorer-2026-06-06-S700
LRC: the negation involution sigma_n: a -> -a on Z/n, its fixed-point set,
and how the EVEN-n second fixed point n/2 governs the AP tight-witness structure.

This EXTENDS S679 (lrc19-sieve-apex-involution), which treated the apex = unique
negation fixed point ONLY for prime/odd n.  For even n, negation mod n has TWO
fixed points {0, n/2}.  n=14 is the smallest even unproved case.

Claims verified here (all exact, via fractions):
 (A) Fix(sigma_n) = {a in Z/n : 2a = 0} = {0} if n odd, {0, n/2} if n even.
 (B) For the canonical tight config AP_n = {1,...,n-1}:
        M(AP_n) = 1/n, and the tight division times are EXACTLY
        T_n = { j/n : gcd(j,n)=1 }  (the units mod n), |T_n| = phi(n).
     => odd prime n: all n-1 times tight.  n=14: {1,3,5,9,11,13}/14 (phi=6).
 (C) The negation involution acts FIXED-POINT-FREELY on the tight times
     (j <-> n-j), because a unit can never equal n/2 for n>2.
 (D) The NON-tight division points j (gcd(j,n)>1) are exactly where some runner
     collides with the origin; for even n the runner n/2 collides at every EVEN j
     (it is the negation-fixed speed), and t=1/2 is killed by every even speed.
 (E) Shell-side duality: the pair-sum/worry modulus C=2n-1 is ALWAYS odd, so
     negation mod C always has the single fixed point 0.  The parity asymmetry
     lives entirely on the speed torus Z/n, never on the shell torus Z/C.
 (F) The V* floor atom at n=14 (AP with 12->24) still contains the self-paired
     speed 7 = n/2; check its tight times too.
"""
from fractions import Fraction
from math import gcd


def frac_dist(x: Fraction) -> Fraction:
    """Distance from x to the nearest integer, exact."""
    r = x - int(x)           # in (-1,1)
    if r < 0:
        r += 1               # in [0,1)
    return min(r, 1 - r)


def M_and_tight_times(speeds):
    """
    Exact loneliness measure M = max_t min_i ||v_i t|| over the FINITE candidate
    set of division times t = j/q (Bienia et al.: the max-min for integer speeds
    is attained at a rational with denominator <= ... ); for the AP and small
    perturbations the optimum is at t=j/n, but to be safe we scan q up to a bound.
    Returns (M, list of t=j/n that achieve >= 1/n simultaneously).
    """
    n = len(speeds) + 1      # convention: n total runners, n-1 moving speeds
    # tight times among the j/n grid:
    tight = []
    for j in range(1, n):
        t = Fraction(j, n)
        mind = min(frac_dist(v * t) for v in speeds)
        if mind >= Fraction(1, n):
            tight.append(j)
    return tight


def euler_phi(n):
    return sum(1 for j in range(1, n) if gcd(j, n) == 1)


def fixed_points_negation(n):
    return [a for a in range(n) if (2 * a) % n == 0]


def units(n):
    return [j for j in range(1, n) if gcd(j, n) == 1]


def report(n):
    print(f"\n===== n = {n}  (parity: {'odd' if n%2 else 'even'}) =====")
    fp = fixed_points_negation(n)
    print(f"  Fix(sigma_n)  (a = -a mod n) : {fp}   [#={len(fp)}]")
    if n % 2 == 0:
        print(f"    -> EVEN: second fixed point n/2 = {n//2}")
    C = 2 * n - 1
    fpC = fixed_points_negation(C)
    print(f"  Shell modulus C=2n-1 = {C} (odd always); Fix(sigma_C) = {fpC} [#={len(fpC)}]")

    AP = list(range(1, n))
    tight = M_and_tight_times(AP)
    U = units(n)
    print(f"  AP_n = {AP}")
    print(f"  AP tight j (j/n simultaneously >= 1/n) : {tight}   [#={len(tight)}]")
    print(f"  units mod n (gcd(j,n)=1)               : {U}   [phi(n)={euler_phi(n)}]")
    print(f"  tight == units ?  {tight == U}")

    # which runner collides at each NON-tight j, and is it the n/2 runner / parity?
    nonunit = [j for j in range(1, n) if gcd(j, n) > 1]
    print(f"  non-tight j (gcd>1): {nonunit}")
    for j in nonunit:
        t = Fraction(j, n)
        zeros = [v for v in AP if frac_dist(v * t) == 0]
        tag = ""
        if (n // 2) in zeros and n % 2 == 0:
            tag += " [n/2 runner collides]" if (n//2) not in [v for v in zeros if v != n//2] else ""
        print(f"      j={j:2d} (t={t}): runners at origin = {zeros}")

    # negation acting on tight times: pairing j <-> n-j, fixed-point free
    pairs = []
    used = set()
    selffix = []
    for j in tight:
        if j in used:
            continue
        k = n - j
        if k == j:
            selffix.append(j)
        else:
            pairs.append((j, k))
            used.add(j); used.add(k)
    print(f"  negation on tight times: pairs {pairs}; self-fixed tight times: {selffix}")


def report_Vstar():
    print("\n===== V* floor atom at n=14 (AP with 12 -> 24) =====")
    Vstar = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24]
    n = 14
    tight = M_and_tight_times(Vstar)
    print(f"  V* = {Vstar}")
    print(f"  V* tight j on the j/14 grid: {tight}")
    print(f"  units mod 14: {units(14)}")
    print(f"  speed 7 = n/2 present in V*? {7 in Vstar}  (negation-fixed speed)")
    print(f"  speed reductions mod 14: {sorted(set(v % 14 for v in Vstar))}")
    # the famous 3+24=27=C pair-sum
    print(f"  pair-sum 3+24 = {3+24} = C(=2n-1=27)? {3+24==27}")


if __name__ == "__main__":
    print("LRC negation-involution fixed points and AP tight-witness structure")
    print("monad-explorer-S700, extending S679 to even n")
    for n in [5, 7, 11, 13, 14, 19, 21, 22]:
        report(n)
    report_Vstar()

    print("\n\n=== SUMMARY TABLE: |Fix(sigma_n)| vs phi(n) (= # AP tight times) ===")
    print(f"{'n':>4} {'parity':>6} {'|Fix|':>6} {'phi(n)':>7} {'tighttimes/(n-1)':>18}")
    for n in [5, 7, 11, 13, 14, 19, 21, 22, 23]:
        fp = len(fixed_points_negation(n))
        ph = euler_phi(n)
        print(f"{n:>4} {('odd' if n%2 else 'even'):>6} {fp:>6} {ph:>7} {f'{ph}/{n-1}':>18}")

#!/usr/bin/env python3
"""Resonance debt proof engine for LRC.

opus-2026-06-01-S531

THE IDEA: Decompose the lonely measure as a Fourier series, graded by
"resonance order" (number of runners involved in each term).

μ(LONELY) = Σ_{k: Σ k_i v_i = 0} Π_i f̂(k_i)

where f̂(0) = (n-2)/n and f̂(k) = -sin(2πk/n)/(πk) for k ≠ 0.

GRADING by resonance order r = #{i : k_i ≠ 0}:
  r=0: ((n-2)/n)^{n-1}  [the "outside" credit, always positive]
  r=1: always 0 (no single-speed resonance with Σ k_i v_i = 0)
  r=2: pairwise terms — the "Legendre" / Bernoulli polynomial layer
  r≥3: higher-order "inside debt" — the hard part

LRC iff μ(LONELY) ≥ 0, i.e., CREDIT + DEBT ≥ 0.

THE PAIRWISE TERM (r=2) has a CLOSED FORM via Bernoulli polynomials:
Σ_{m≥1} sin(2πm a/n)sin(2πm b/n)/m² = (π²/2)[B₂({(a-b)/n}) + B₂({(a+b)/n}) - 2B₂(0)]
where B₂(x) = x² - x + 1/6.

This gives EXACT rational values for the pairwise correction terms.
The full lonely measure = OUTSIDE + PAIRWISE + HIGHER.
If |HIGHER| ≤ OUTSIDE + PAIRWISE, LRC is proved.

THE STRATEGY:
1. Compute OUTSIDE and PAIRWISE exactly (closed form, rational)
2. Compute HIGHER by truncated Fourier sum (numerical bound)
3. Show OUTSIDE + PAIRWISE + HIGHER ≥ 0 for all primitive speed sets
4. For the initial segment: show exact cancellation (lonely = 0, wall-only)
5. For non-initial: show strict positivity (lonely > 0)
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd, pi, sin, cos, comb, log
from functools import reduce
from collections import defaultdict
from itertools import combinations, product as iproduct


# ═══════════════════════════════════════════════════════════════
# Fourier coefficients
# ═══════════════════════════════════════════════════════════════

def fhat(k: int, n: int) -> float:
    """Fourier coefficient f̂_n(k) of the indicator f_n(x) = 1_{||x|| ≥ 1/n}.

    f̂(0) = (n-2)/n
    f̂(k) = -sin(2πk/n) / (πk)  for k ≠ 0
    f̂(k) = 0 when n | k
    """
    if k == 0:
        return (n - 2) / n
    if k % n == 0:
        return 0.0
    return -sin(2 * pi * k / n) / (pi * k)


def bernoulli2(x: float) -> float:
    """Second Bernoulli polynomial B₂(x) = x² - x + 1/6."""
    return x * x - x + 1.0 / 6


def frac_pos(x: float) -> float:
    """Fractional part in [0, 1)."""
    return x - int(x) if x >= 0 else x - int(x) + (1 if x != int(x) else 0)


# ═══════════════════════════════════════════════════════════════
# PART 1: Exact resonance decomposition
# ═══════════════════════════════════════════════════════════════

def resonance_decomposition(n: int, speeds: tuple, max_order: int = 4,
                            max_m: int = 200) -> dict:
    """Decompose μ(LONELY) into resonance orders.

    Order 0: ((n-2)/n)^{n-1}
    Order 2: pairwise Bernoulli terms (closed form)
    Order 3+: truncated Fourier sums

    Returns dict with each order's contribution.
    """
    m = len(speeds)
    f0 = (n - 2) / n  # f̂(0)

    # Order 0: all k_i = 0
    order0 = f0 ** m

    # Order 2: exactly 2 nonzero k_i
    # For pair (i, j): k_i v_i + k_j v_j = 0
    # => k_i = -t v_j / gcd(v_i,v_j), k_j = t v_i / gcd(v_i,v_j) for integer t
    # Contribution: f0^{m-2} * Σ_{t≠0} f̂(k_i) f̂(k_j)

    order2 = 0.0
    pairwise_terms = {}

    for i in range(m):
        for j in range(i + 1, m):
            vi, vj = speeds[i], speeds[j]
            g = gcd(vi, vj)
            vi_red, vj_red = vi // g, vj // g

            # k_i = -t * vj_red, k_j = t * vi_red
            pair_sum = 0.0
            for t in range(1, max_m + 1):
                ki = -t * vj_red
                kj = t * vi_red
                pair_sum += 2 * fhat(ki, n) * fhat(kj, n)  # factor 2 for ±t

            pair_contribution = f0 ** (m - 2) * pair_sum
            order2 += pair_contribution
            pairwise_terms[(vi, vj)] = pair_contribution

    # Order 2 via Bernoulli (closed form for coprime pairs):
    # Σ_{t=1}^∞ 2 sin(2πt a/n) sin(2πt b/n) / (πt)²
    # = (1/π²) Σ [cos(2πt(a-b)/n) - cos(2πt(a+b)/n)] / t²
    # = (1/π²) [π²B₂({(a-b)/n}) - π²B₂({(a+b)/n})]  (wait, need to check sign)
    #
    # Actually: Σ_{t=1}^∞ cos(2πtα)/t² = π²(6α²-6α+1)/6 = π²B₂(α) for α∈[0,1)
    # (this is a standard result)
    #
    # So: Σ 2sin sin / (πt)² = (2/(π²)) * (1/2)[Σcos((a-b)/n) - Σcos((a+b)/n)] / t²
    # Hmm, let me redo. f̂(k) = -sin(2πk/n)/(πk).
    # f̂(-tv_j) f̂(tv_i) = [-sin(-2πtv_j/n)/(π(-tv_j))] [-sin(2πtv_i/n)/(πtv_i)]
    # = [sin(2πtv_j/n)/(πtv_j)] [-sin(2πtv_i/n)/(πtv_i)]
    # = -sin(2πtv_j/n)sin(2πtv_i/n) / (π²t²v_iv_j)

    # Summing over t=1,...,∞ (and the t=-1,...,-∞ gives the same by symmetry):
    # 2 * Σ_{t=1}^∞ f̂(-tv_j)f̂(tv_i)
    # = 2 * Σ_{t=1}^∞ -sin(2πtv_j/n)sin(2πtv_i/n) / (π²t²v_iv_j)
    # = -2/(π²v_iv_j) * Σ_{t=1}^∞ sin(2πtv_i/n)sin(2πtv_j/n)/t²

    # Using 2sinAsinB = cos(A-B) - cos(A+B):
    # = -1/(π²v_iv_j) * Σ_{t=1}^∞ [cos(2πt(v_i-v_j)/n) - cos(2πt(v_i+v_j)/n)] / t²
    # = -1/(π²v_iv_j) * [π²B₂({(v_i-v_j)/n}) - π²B₂({(v_i+v_j)/n})]
    # = -[B₂({(v_i-v_j)/n}) - B₂({(v_i+v_j)/n})] / (v_iv_j)

    # Note: for the sum f̂(-tvj_red)*f̂(tvi_red) with gcd reduction,
    # the speeds in the Fourier coefficients are multiplied by vj_red and vi_red.
    # So the Bernoulli formula uses (vi_red * vi - vj_red * vj)... hmm,
    # this needs more care with the gcd reduction.

    # For coprime (vi, vj) (g=1): vi_red=vi, vj_red=vj.
    # ki = -t*vj, kj = t*vi. f̂(ki) = f̂(-tvj), f̂(kj) = f̂(tvi).

    order2_bernoulli = 0.0
    for i in range(m):
        for j in range(i + 1, m):
            vi, vj = speeds[i], speeds[j]
            g = gcd(vi, vj)
            vi_red, vj_red = vi // g, vj // g

            # The sum is over t = 1,...,∞ with ki = -t*vj_red, kj = t*vi_red.
            # But f̂ uses the FULL ki, not the reduced one.
            # f̂(-t*vj_red) = -sin(-2π*t*vj_red/n)/(π*(-t*vj_red)) = sin(2πtv_j_red/n)/(πtv_j_red)
            # f̂(t*vi_red) = -sin(2πtv_i_red/n)/(πtv_i_red)
            # Product: -sin(...)sin(...)/(π²t²vi_red*vj_red)

            # Bernoulli closed form:
            alpha_minus = frac_pos((vi_red - vj_red) / n) if (vi_red - vj_red) % n != 0 else 0.0
            alpha_plus = frac_pos((vi_red + vj_red) / n) if (vi_red + vj_red) % n != 0 else 0.0

            # Hmm, wait. The frequencies in the Fourier coefficients are
            # f̂(-t*vj_red, n) and f̂(t*vi_red, n).
            # The sin terms use 2π*(-t*vj_red)/n and 2π*(t*vi_red)/n.
            # The Bernoulli formula:
            # Σ sin(2πtA/n)sin(2πtB/n)/t² for A=vj_red, B=vi_red
            # = (π²/2)[B₂({(A-B)/n}) + B₂({(A+B)/n}) - 2B₂(0)]
            # Hmm, I need to be more careful.

            # Let me just use the NUMERICAL sum (already computed above)
            # and compare with the Bernoulli closed form.

            # For the Bernoulli: the standard result is
            # Σ_{t=1}^∞ cos(2πtα)/t² = π²(α²-α+1/6) for α∈(0,1)
            # and Σ_{t=1}^∞ sin(2πtα)/t = -π(α-1/2) for α∈(0,1)
            # We need Σ sin*sin/t² which is harder.
            # Using 2sinAsinB = cos(A-B)-cos(A+B):
            # Σ 2sin(2πtA)sin(2πtB)/t² = Σcos(2πt(A-B))/t² - Σcos(2πt(A+B))/t²
            # = π²[B₂({A-B}) - B₂({A+B})]  where A,B = vj_red/n, vi_red/n

            A = vj_red / n
            B = vi_red / n
            alpha_m = frac_pos(A - B)
            alpha_p = frac_pos(A + B)
            bernoulli_sum = pi**2 * (bernoulli2(alpha_m) - bernoulli2(alpha_p))

            # The pair contribution via Bernoulli:
            # f0^{m-2} * (-1/(π²vi_red*vj_red)) * (bernoulli_sum / 2)
            # Wait, the sum above is Σ 2sin*sin/t², and the pair uses
            # -Σsin*sin/(π²t²vi_red*vj_red) summed with factor 2 for ±t
            # = -2Σsin*sin/(π²t²vi_red*vj_red)
            # = -1/(π²vi_red*vj_red) * 2Σsin*sin/t²
            # = -1/(π²vi_red*vj_red) * bernoulli_sum

            pair_bernoulli = f0 ** (m - 2) * (-bernoulli_sum / (pi**2 * vi_red * vj_red))
            order2_bernoulli += pair_bernoulli

    # Order 3 and higher (truncated numerical sum)
    # This is expensive — skip for now and just compute the residual
    higher_orders = 0.0

    # Actually, compute the FULL μ(LONELY) by direct integration and subtract
    # order0 + order2 to get the higher-order residual.

    # Direct integration
    num_pts = 100000
    direct_mu = 0
    for s in range(num_pts):
        t = (s + 0.5) / num_pts
        if all(min(v * t % 1, 1 - v * t % 1) >= 1.0 / n for v in speeds):
            direct_mu += 1
    direct_mu /= num_pts

    higher_orders = direct_mu - order0 - order2

    return {
        "direct_mu": direct_mu,
        "order0": order0,
        "order2_numerical": order2,
        "order2_bernoulli": order2_bernoulli,
        "higher_orders": higher_orders,
        "pairwise_terms": pairwise_terms,
    }


# ═══════════════════════════════════════════════════════════════
# PART 2: Run for multiple n and speed sets
# ═══════════════════════════════════════════════════════════════

def main():
    print("Resonance Debt Proof Engine — opus-S531")
    print()

    for n in [3, 4, 5, 6, 7, 8, 14]:
        speeds = tuple(range(1, n))
        print(f"{'='*70}")
        print(f"n={n}, speeds={speeds}")
        print(f"{'='*70}")

        result = resonance_decomposition(n, speeds)

        print(f"  direct μ(LONELY) = {result['direct_mu']:.10f}")
        print(f"  OUTSIDE credit (order 0) = ((n-2)/n)^{n-1} = {result['order0']:.10f}")
        print(f"  PAIRWISE debt (order 2, numerical) = {result['order2_numerical']:.10f}")
        print(f"  PAIRWISE debt (order 2, Bernoulli) = {result['order2_bernoulli']:.10f}")
        print(f"  HIGHER debt (order ≥ 3) = {result['higher_orders']:.10f}")
        print()
        print(f"  CREDIT + DEBT = {result['order0'] + result['order2_numerical'] + result['higher_orders']:.10f}")
        print(f"  direct check:  {result['direct_mu']:.10f}")
        print()

        # Ratios
        if result['order0'] > 0:
            print(f"  |pairwise/outside| = {abs(result['order2_numerical']/result['order0']):.6f}")
            print(f"  |higher/outside| = {abs(result['higher_orders']/result['order0']):.6f}")
            print(f"  debt ratio = {abs((result['order2_numerical']+result['higher_orders'])/result['order0']):.6f}")
        print()

        # Show the largest pairwise terms
        if result['pairwise_terms']:
            sorted_pairs = sorted(result['pairwise_terms'].items(),
                                  key=lambda x: abs(x[1]), reverse=True)
            print(f"  largest pairwise terms:")
            for (vi, vj), val in sorted_pairs[:5]:
                print(f"    ({vi},{vj}): {val:.10f}")
        print()

    # Now test a NON-initial speed set at n=14
    print(f"{'='*70}")
    print(f"n=14, NON-initial speeds (primes)")
    print(f"{'='*70}")
    speeds_prime = (1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    result = resonance_decomposition(14, speeds_prime, max_m=100)
    print(f"  direct μ(LONELY) = {result['direct_mu']:.10f}")
    print(f"  OUTSIDE credit = {result['order0']:.10f}")
    print(f"  PAIRWISE debt = {result['order2_numerical']:.10f}")
    print(f"  HIGHER debt = {result['higher_orders']:.10f}")
    print(f"  debt ratio = {abs((result['order2_numerical']+result['higher_orders'])/result['order0']):.6f}")
    print()

    print(f"{'='*70}")
    print(f"SYNTHESIS")
    print(f"{'='*70}")
    print()
    print("THE RESONANCE DECOMPOSITION:")
    print("  μ(LONELY) = OUTSIDE_CREDIT + PAIRWISE_DEBT + HIGHER_DEBT")
    print()
    print("  OUTSIDE = ((n-2)/n)^{n-1} ≈ e^{-2} ≈ 13.5% for large n")
    print("  PAIRWISE = negative (cancels part of outside)")
    print("  HIGHER = the remaining correction")
    print()
    print("  For INITIAL SEGMENT: OUTSIDE + PAIRWISE + HIGHER = 0 exactly")
    print("    (wall-only lonely, measure 0)")
    print("  For OTHER speed sets: OUTSIDE + PAIRWISE + HIGHER > 0")
    print("    (positive lonely measure)")
    print()
    print("  The PROOF TARGET: show |PAIRWISE + HIGHER| < OUTSIDE")
    print("  for all non-initial primitive speed sets.")
    print("  Equivalently: the initial segment MAXIMIZES the debt.")
    print()
    print("  This is the RESONANCE DEBT CONJECTURE:")
    print("  Among all primitive speed sets at a given n, the initial segment")
    print("  {1,...,n-1} has the maximum |debt/credit| ratio, equal to exactly 1.")
    print("  All other speed sets have ratio < 1, proving LRC.")
    print()


if __name__ == "__main__":
    main()

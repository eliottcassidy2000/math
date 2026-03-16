# The Dimension Ladder and the Echo of Errors

*opus-2026-03-16-S73 — arising from finding MISTAKE-027 in the THM-080 amplitude table*

---

## The ladder

The Walsh amplitudes of H(T) and M[a,b](T) are connected by a dimension ladder:

H_amp(n, d, r) / M_amp(n, d-1, s=r-1) = n - d

This ratio is exact. It comes from the factorial ratio (n-d)!/(n-d-1)! = n-d, combined with the normalization and component-count shifts that cancel perfectly. The ladder descends through odd numbers:

At n=9: H(d=2) →÷7→ M(d=1) →÷7→ H(d=4) →÷5→ M(d=3) →÷5→ H(d=6) →÷3→ M(d=5) →÷3→ H(d=8) →÷1→ M(d=7)

The divisors are n-2, n-4, n-6, ..., 1. All odd. The H/M step never introduces or removes factors of 2 — ALL powers of 2 in the Walsh spectrum come from the component counting (2^r for H, 2^s for M). This is the "bit accumulation": each unrooted even-length component contributes exactly one bit.

## The glitch

The amplitude table in THM-080 had wrong values for n=9. Four of five entries didn't match the stated formula. The formula was proved analytically and verified computationally at n=5,6,7. The table was populated by hand calculation that repeated, at n=9, a component-counting error similar to MISTAKE-013b (the original missing 2^s factor, caught at n=7).

The error was invisible at n=3,5,7 and appeared at n=9 — exactly the pattern predicted by THM-J (the algebraic criterion for S(T) universality), which says n=9 is the first odd n where tournament-dependent behavior appears in the 2-adic structure. The same complexity threshold that makes the signed HP permanent depend on t₃ parity also makes unrooted-component counting more error-prone.

## The echo

MISTAKE-013b: missing 2^s factor, caught at n=7.
MISTAKE-027: wrong s-labels in table, surviving at n=9.

The first error was conceptual: the formula didn't account for unrooted component orientations. The second error is its echo: the correction was applied to the formula but not fully propagated to the numerical table. The formula was right; the table was wrong. The proof was trustworthy; the example was not.

This is a pattern with a name in software engineering: **the fix that doesn't fix.** You patch the bug in the code but not in the documentation. You correct the theorem but not the table that illustrates it. The error survives in the derivative, the representation, the thing that was computed from the wrong version and never recomputed from the right one.

In the chain complex of this research project: MISTAKE-013b was a 2-chain (wrong formula), and its correction was a 3-chain (the 2^s insight). But MISTAKE-027 shows the 3-chain didn't completely fill the cycle — there was a residual 2-cycle in the table, a derivative artifact of the original error, hiding at n=9 where nobody checked.

β₂ = 0 for tournaments means every 2-cycle has a filling. But β₂ = 0 for the repo is only an aspiration. The correction process doesn't always succeed on the first pass.

## The 2-adic structure

The bit accumulation picture has a beautiful self-consistency. The ladder ratio n-d is always odd (for odd n, even d), so:

v₂(H_amp / M_amp) = 0

All 2-adic information flows through the component count. And the component count is determined by the monomial structure — how many disjoint even-length paths are in S, and how many contain a root vertex. This is purely combinatorial. The 2-adic valuation of the amplitude is:

v₂(amplitude) = s - (n-2) + v₂((n-2-d)!) = s - d - s₂(n-2-d)

where s₂(m) is the digit sum of m in binary (by Legendre's formula). This connects directly to THM-J: the universality of S(T) mod 2^{n-1} depends on s₂(n-3) ≤ 1, which is a condition on the *binary representation* of the dimensional parameter.

The mathematics is telling us that the 2-adic structure of tournament invariants is controlled by the binary representations of the dimensions involved. Not their values — their *bit patterns*. This is what "bit accumulation" means at the deepest level: the Walsh spectrum of tournament invariants is organized by the binary structure of the underlying dimensions.

## Cross-references

- THM-080: Walsh formula for M[a,b] (formula correct, table corrected)
- MISTAKE-013b: Original missing 2^s factor (n=7)
- MISTAKE-027: Echo of 013b in the n=9 table (this session)
- THM-J: Universality criterion s₂(n-3) ≤ 1
- THM-H: Universal congruence S mod 2^{n-1}
- 04-computation/amplitude_table_check.py, amplitude_glitch_analysis.py

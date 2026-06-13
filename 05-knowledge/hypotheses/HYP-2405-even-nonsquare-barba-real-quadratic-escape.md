# HYP-2405 — The even non-square Barba branch escapes into real quadratic fields

**Status:** OPEN (claimed claudebox-2026-06-11-S1); evidence at n = 10, 18.
**Companions:** THM-476 (the square law that creates this branch), THM-475 (odd branch),
HYP-2389 (n ≡ 1 mod 4), OPEN-Q-058.

## Statement

For n ≡ 2 (mod 4) with 2n−3 NOT a perfect square (so tournaments cannot attain the
Ehlich–Wojtas bound, THM-476), the tournament maxdet B_t(n) = max det(I+S) is attained at
SSᵀ spectra that are irrational — conjugate pairs in a real quadratic field — rather than at
the nearest integer two-level spectrum.

Concretely at n = 10: B_t(10) = 64000 = 2⁹·5³, attained with S-char-poly
(x²−18x+61)⁴(x−9)² ... i.e. SSᵀ eigenvalues 9 ± 2√5 (multiplicity 4 each) and 9 (mult 2):
the maximizer lives in ℚ(√5) (note 9 ± 2√5 = 3 + (1±√5)², a golden-ratio shift). Best of
20 annealing restarts × 150k steps, never exceeded; conjectured exact.

At n = 18 (2n−3 = 33): best found 131 533 373 440 = 2²⁹·5·7² ≈ 0.901·EW(18) (8 restarts);
spectrum not yet identified; the 2013 dihedral-cocyclic literature records ~0.853·EW at n=18,
so this already beats the published record fraction if it survives verification — flagged for
follow-up, NOT yet a claim (different normalizations possible; verify against
Đoković/Kotsireas-era tables before claiming).

## Why interesting

The disproof-of-grid-conjecture motif (extremal objects escape the rational/lattice structure
into a number field) recurs INSIDE the maximal determinant problem: when the integer two-level
spectrum is arithmetically barred (no k² = 2n−3), the optimum splits the excited level into a
Galois-conjugate pair. The repo's flat-vs-peaked autocorrelation dichotomy (THM-441) acquires
a third regime: flat (DRT/skew-Hadamard) > split-conjugate (this branch) > peaked.

## Tests

1. Exhaust n = 10 per iso class (9 733 056 classes — feasible with the mac-mini numpy
   pipeline; claudebox lacks numpy) to settle B_t(10) = 64000.
2. Identify the n = 18 maximizer spectrum; check whether the ℚ(√d) pattern persists and which
   d is selected (d = 5 at n = 10; predict d | disc-related to 2n−3?).
3. n = 22 (2n−3 = 41): the literature records 0.877·EW from 2013 cocyclic searches; anneal +
   spectral analysis.

# The maxdet ladder: flags, square laws, and the golden escape (claudebox-2026-06-11-S1)

User directive: tournaments × Hadamard matrices, leverage the repo's machinery toward famous
open questions; consider Tao's teorth/optimizationproblems; hunt hidden gems in cross-area
connections.

## What happened

The session landed on the **tournament maximal-determinant ladder** — the skew/tournament
shadow of the Hadamard maximal determinant problem (Tao constants C23a/b/c are the external
anchors; C23a = smallest open Hadamard order 668). mac-mini-2026-06-10-S2 had just built the
floor/ceiling/average canon (THM-468/472/473) and left the two non-attainable residues as
HYP-2389/OPEN-Q-058. This session closed the construction half and found the law of the
even residue:

1. **THM-475 (the flag).** At n ≡ 1 mod 4 the integrality obstruction (no DRT) materializes
   as ONE excited eigenpair, and the minimal-defect spectrum is REALIZED: take a DRT on n−2
   and add two stacked apexes (u, v beat everything, u → v). A 3×3 char-poly computation
   (eigenvalues {2n−3, 2n−3, 0} on the flag-coupled space) gives det(I+S) = 2(n−1)^((n−1)/2)
   exactly — the conjectured maximum, with the conjectured spectrum. The two flag vertices
   are literally apexes: the Barba maximizer is "a DRT watched by two stacked observers."
   The construction works through every DRT species we have: Paley primes, the THM-448
   doubling tower (DRT(15)), and GF(27) (where circulants provably fail — the ramified
   3³ shell needs the field group, not the cyclic group; the LRC ramification motif again).

2. **THM-476 (the square law).** At n ≡ 2 mod 4 a tournament attains the Ehlich–Wojtas
   bound (the TRUE maxdet over ALL ±1 matrices) only if 2n−3 is a perfect square. Proof in
   five lines: EW Gram rigidity pins the top eigenspace to two signed indicators a, b;
   skewness forces Sa = μb with μ ∈ ℤ and μ² = 2n−3. Equivalently: the classical
   two-squares datum 2n−2 = a₁² + b₁² degenerates to a₁ = 1 — the +1/observer again, now
   as an arithmetic obstruction. Attainable orders = (k²+3)/2 = 6, 14, 26, 42, 62, 86, …
   — exactly the known/open frontier of the skew E-W literature (Greaves–Suda; first open
   86). We constructed explicit two-circulant witnesses (exact integer PAF conditions) at
   n = 14, 26, 62 — every candidate ≤ 62 — and showed n = 86 RESISTS all
   multiplier-symmetric two-circulant ansätze (subgroup orders 3, 7, 21 on ℤ/43).

3. **HYP-2405 (the golden escape).** At the non-square orders the optimum leaves ℚ: the
   n = 10 best (64000, 20 restarts stable) has SSᵀ eigenvalues 9 ± 2√5 — a Galois-conjugate
   excited pair in ℚ(√5) replacing the barred integer level. The grid-disproof motif
   (extremal objects escape the lattice into a number field) shows up INSIDE maxdet.

4. **Doubling corollary.** det(I + S_{D(T)}) = 2ⁿ·det(I+S_T)² (one line from THM-447's
   spectral law): the skew-Sylvester double preserves the Hadamard ceiling exactly and is
   strictly suboptimal at every other residue — the ladder is doubling-closed only at the
   Hadamard rung.

## The ladder (the session's one picture)

| n mod 4 | ceiling                       | attainment law            | status |
|---------|-------------------------------|---------------------------|--------|
| 0       | n^(n/2) (Hadamard)            | skew-Hadamard conjecture  | open frontier 4k = 167·4 = 668 (Tao C23a) |
| 1       | 2(n−1)^((n−1)/2) (conj.)      | DRT(n−2) + flag (THM-475) | lower PROVED; upper = OPEN-Q-058 |
| 2       | 2(n−1)(n−2)^((n−2)/2) (= EW)  | iff 2n−3 = k² (THM-476)   | necessity PROVED; witnesses ≤ 62; open 86 |
| 3       | (n+1)^((n−1)/2)               | DRT ⟺ skew-Hadamard (THM-472) | = the skew-Hadamard conjecture |

Two of the four rungs are governed by the observer/+1 (the flag apexes; the forced 1² in the
two-squares datum); the other two ARE the skew-Hadamard conjecture. The whole ladder is the
repo's flat-vs-peaked autocorrelation dichotomy made arithmetic: flat (DRT/conference) where
allowed, minimal-defect (one excited pair) where barred, Galois-split where even that is barred.

## Honest scope

- The n ≡ 1 mod 4 UPPER bound remains open (OPEN-Q-058); n=9 is settled only by mac-mini's
  exhaustion, n=13 by 1M+ restarts of evidence.
- THM-476 necessity may be implicit in the Greaves–Suda char-poly characterization — our
  proof is elementary and self-contained, and the witness constructions + the n=86 negative
  are new data either way; exact-source check emitted as part of t-0114.
- The n=18 anneal value 131 533 373 440 ≈ 0.901·EW(18) would beat the 2013 cocyclic record
  fraction (0.853) if normalizations match — flagged in HYP-2405, NOT claimed.
- Annealing lower bounds are not maxima; all witness determinants are exact (Bareiss over ℤ).

## Files

THM-475, THM-476, HYP-2405 (+ INDEX), OPEN-Q-058 update,
04-computation/tournament_barba_flag_cbx1.py (+ .out),
04-computation/skew_ew_two_circulant_cbx1.py (+ .out), tasks t-0114/0115/0116.

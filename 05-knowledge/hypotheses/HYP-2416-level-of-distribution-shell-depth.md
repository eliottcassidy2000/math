# HYP-2416 — The level-of-distribution / shell-depth dictionary, and additive representation systems as transfer operators

**Status:** OPEN (synthesis frame, claimed claudebox-2026-06-11-S5). A DICTIONARY plus provable
sub-pieces (THM-485/486); NOT a resolution of the famous problems. Honest scope stated throughout.
**Companions:** THM-485/486, THM-406/S561 (LRC sieve ρ), the S625 window lemma, GoldbachLemoine.lean
(S630), goldbach_polygonal_zeckendorf_bridge_s501. **External:** Elliott–Halberstam / Bombieri–
Vinogradov (level of distribution θ); von Mangoldt Λ = μ * log; Helfgott (ternary Goldbach); Cauchy
(polygonal); Euler (pentagonal number theorem).

## A. The Möbius/sieve face: von Mangoldt ↔ the LRC covering-depth sieve

von Mangoldt Λ = μ * log: Λ(n) = Σ_{d|n} μ(d) log(n/d) — a Möbius convolution that extracts the
PRIME-POWER atoms from the smooth function log. The LRC covering-depth / Erdős-2 covering density
(THM-406, S561, S724) is ρ = Σ_{T} (−1)^|T| / lcm(T) — a Möbius/inclusion–exclusion sum over the
divisor (lcm) lattice that extracts the IRREDUCIBLE-RESONANCE atoms (the tight/collapse cores).
**Both are Möbius transforms over a divisibility lattice that isolate atoms from a smooth/uniform
baseline.** Λ ↔ ρ is a structural parallel (the prime-distribution face of the repo's sieve).

## B. The level-of-distribution face: Elliott–Halberstam θ ↔ the LRC shell-window depth

EH: the primes have level of distribution θ if Σ_{q ≤ x^θ} max_a |ψ(x;q,a) − x/φ(q)| ≪ x/log^A x.
θ = how DEEP into the modulus range one controls prime-equidistribution: Bombieri–Vinogradov gives
θ = 1/2 unconditionally (the √-barrier), EH conjectures θ → 1 (false at θ = 1). The LRC window lemma
(S625) is the SAME shape: among the φ(m) multipliers of a shell m = (ℤ/m)*, the width-2 danger band
blocks ≤ 2 per runner, so a good multiplier survives once the shell is deep enough (m ≳ 2n). The
"level" = the shell DEPTH one can use; the critical shell is 2n−1. The √-barrier θ = 1/2 mirrors the
gap between the easy constant-factor bound M > 1/(2n) (shell ~ 2n, S625) and the OPTIMAL 2/(2n−1)
(the critical shell itself, THM-415) — closing to the critical shell is the LRC analogue of pushing
θ past 1/2 to 1. The parity/perspective obstruction that floors prime gaps at 6 (GEH) mirrors the
single-corrector → multi-sieve ceiling (the perspective key). **HONEST: heuristic analogy, not a
proven equivalence; the value of θ for the multiplier orbits is the content of OPEN-Q-062.**

## C. Additive representation systems as a redundancy ladder of transfer operators

Every "represent n as a constrained sum" system is counted by a transfer operator / partition
function; they form a ladder by redundancy:
- **Zeckendorf (unique):** non-consecutive Fibonacci = golden-mean SFT = the path independence
  operator at x=1, growth φ (THM-485). Zero redundancy — a bijective numeration.
- **Fermat polygonal (bounded depth, Cauchy):** every n = sum of k k-gonal numbers. Triangular
  T_k = k(k+1)/2 = C(k+1,2) = the EDGE COUNT of the tournament on k+1 vertices (and the repo tiling
  dimension m = C(n−1,2) is itself triangular); Gauss's "every n = Δ+Δ+Δ" = every n is a sum of 3
  tournament-edge-counts. PENTAGONAL numbers g_j = j(3j−1)/2 give Euler's pentagonal number theorem
  ∏(1−x^k) = Σ(−1)^j x^{g_j} = the recurrence for the PARTITION FUNCTION p(n) — the repo's
  "partition functions everywhere" (S626). (Contrast: p(n) is sub-exponential / NOT C-finite, the
  opposite regime from the exponential, C-finite tournament families — S716.)
- **Goldbach / Lemoine (conjectural):** even = p+q (σ-symmetric / additive), odd = p+2q (σ-broken /
  multiplicative) — the S630 pair, with 6 = 2·3 the unique σ-fixed commutation point
  (GoldbachLemoine.lean). Maximal redundancy (Hardy–Littlewood r(n) → ∞), unproven existence.

## D. The disorder axis (Viswanath) over the whole ladder

THM-485's second temperature (quenched sign-disorder) applies to every transfer operator above:
randomizing signs replaces the growth eigenvalue by a Lyapunov exponent (Viswanath at the Zeckendorf
node) and introduces a phase transition at the Embree–Trefethen activity β*. So the ladder has TWO
coordinates: redundancy (which system) and disorder (ordered eigenvalue ↔ quenched Lyapunov).

## Tests / open sub-questions

1. OPEN-Q-062: make the LRC "level of distribution" precise (a Bombieri–Vinogradov-type average over
   shells) and locate its θ = 1/2 analogue.
2. Is the LRC sieve ρ literally a (twisted) von Mangoldt average over the shell tower? (compare the
   Euler products: shell tower ∏(1 − …) vs −ζ'/ζ = ∏ local factors.)
3. Compute Lemoine/Goldbach representation counts r(n) through the repo's σ-pair lens; confirm the
   Hardy–Littlewood growth (not a proof, a sanity/structure check).
4. Does the polygonal → theta-function → self-dual-code chain (triangular/square sums = θ-powers;
   even unimodular lattice θ-series at dim 8,16,24) connect Fermat polygonal directly to THM-481's
   Golay/Leech gauge codes? (the 24 of THM-486.)

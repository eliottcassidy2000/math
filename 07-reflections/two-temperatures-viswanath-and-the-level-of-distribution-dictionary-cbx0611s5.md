# Two temperatures, Viswanath, and the level-of-distribution dictionary (claudebox-2026-06-11-S5)

User dispatch: continue the Gleason/involution-modulus line; integrate von Mangoldt, the
Elliott–Halberstam level-of-distribution exponent, Viswanath's constant, and the repo's
Zeckendorf / Fermat-polygonal / Goldbach-Lemoine work; make concrete progress on famous problems.

## The organizing idea

Every object on the list is a TRANSFER OPERATOR read at some temperature. The repo already had the
deterministic side (S720/721: the temperature ladder, the tournament H = I(Ω,2), φ/tribonacci/3.383
growth constants, the Fibonacci=independence-on-paths identity). What was missing, and what these
new ideas supply, is a SECOND temperature axis and the arithmetic that ties the two sides together.

## What is genuinely new and solid (THM-485/486)

1. **Two temperatures of one operator (THM-485).** The path/hard-core operator
   I(P_n,x)=I(P_{n-1},x)+x·I(P_{n-2},x) has a deterministic FUGACITY axis — x=1 is Zeckendorf
   (golden-mean SFT, growth φ), x=2 is the repo's H (Jacobsthal) — and a QUENCHED SIGN-DISORDER
   axis. Randomizing the recurrence signs replaces the growth eigenvalue by a Lyapunov exponent,
   and at x=1 that is exactly **Viswanath's constant 1.13198824…** (reproduced to 4 digits).
   φ and Viswanath are the ordered and disordered endpoints of ONE operator — the same operator
   whose x=2 reading is the tournament partition function. NEW DATA: the disordered constants of
   the repo's families (tribonacci 1.839→1.223, base-path 3.383→2.979), and a DISORDER-INDUCED
   PHASE TRANSITION at the Embree–Trefethen activity β*≈0.70258 — the quenched hard-core gas decays
   below β* and grows above, a transition the deterministic operator (always growing) never has.
   This is the transfer-operator face of the S637 glass transition: disorder = the spin-glass /
   Lyapunov regime; the eigenvalue ≥ Lyapunov gap is the glass cost.

2. **24 is the Pisano modulus (THM-486).** The same fact 24 | p²−1 that makes (ℤ/24)* all-involutions
   (THM-484, the Gleason/Golay orders) also governs the Fibonacci period: π(p) | p²−1, and π(24)=24
   with π(n)=n exactly on {1}∪{24·5ᵏ} — 24 is the base Pisano-fixed point. α(24)=12, F₁₂=144=12²
   (the unique nontrivial Fibonacci square). This is the bridge from the code/involution side
   (THM-481/484) to the Zeckendorf/Fibonacci side (THM-485): ℚ(√5) controls the Pisano 5-tower,
   ℚ(√−3) the code side, and 24 = lcm of the 8-part (involutions) and the 3-part. Caught and logged
   a false recalled-folklore claim ("24 = largest n with π(n)≤n") rather than canonize it.

## The dictionary (HYP-2416, honest framing)

von Mangoldt Λ=μ*log and the LRC covering sieve ρ=Σ(−1)^|T|/lcm are both Möbius transforms extracting
atoms (primes / irreducible resonances). The Elliott–Halberstam level-of-distribution θ — how deep
into the modulus range primes equidistribute, BV's θ=1/2 vs EH's θ→1 — mirrors the LRC shell-window
DEPTH (S625: a good multiplier survives once the shell is deep enough; the √-barrier ↔ the gap
between the easy M>1/(2n) and the optimal 2/(2n−1)). Additive representation systems form a
redundancy ladder: Zeckendorf (unique, φ) → Fermat polygonal (bounded; triangular = tournament edges
C(k+1,2); pentagonal = Euler's partition function) → Goldbach/Lemoine (conjectural, the σ-pair).
OPEN-Q-062 makes the LRC level-of-distribution a precise Bombieri–Vinogradov-type question.

## Honest scope (the part that matters)

The famous problems (binary Goldbach, Lemoine, Elliott–Halberstam, bounded gaps) get a DICTIONARY
and precise reframings, NOT proofs — and I say so. The concrete NEW mathematics is: the disordered
transfer constants (genuinely new objects extending S720/721), the disorder-induced phase transition
at β*, and the 24=Pisano-modulus bridge. The deterministic-vs-quenched two-temperature structure is
the real conceptual gain — it places Viswanath's constant inside the repo's transfer-spectrum thread
and gives the glass-transition theme (S637) an exact small-operator instance. The polygonal→theta→
self-dual-code chain (test 4 of HYP-2416) is the most promising lead back to the Golay/Leech 24.

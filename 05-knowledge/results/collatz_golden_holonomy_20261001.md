# Collatz XXI: Collatz parity read in base φ — one polynomial identity behind `7·Φ_T(n) ∈ Z`, the golden β-map, the pentagon cosines, and a holonomy mirror of Apéry

**opus S15, twenty-first note, 2026-10-01.** Worktree `collatz-functional-uniqueness-20261001`.

- Script: `04-computation/experiments/collatz_golden_holonomy_20261001.py` (+ `.out`, ALL CHECKS PASSED).
- Canon: THM-4528.

**Status.**
- PROVED (elementary): Theorems G1–G3.
- FINITE-EXACT: the lattice census in G2, and the checks for `n ≤ 20000`.
- **NO PROOF of Collatz.** These are exact restatements (§6).
- Independent audit: DONE (2026-10-01, blind re-derivation with exact arithmetic). G1–G3 are SOUND; the census and
  every identity were reproduced. Wording corrections are applied (semi-conjugacy, the exceptional point, the
  measure; MISTAKE-556), and the record is in §8.

## 0. The owner's prompt

> Collatz holds iff 7·Φ_T(n) is an integer for every n.
>
> φ⁸ + 1 = 7φ⁴, φ¹⁰ = 1 + 11φ⁵. In base φ, "11" equals "100" because φ + 1 = φ². Consecutive 1s always rewrite
> that way. Think of 11 as bit shift with memory.

The first line is Proposition 3 of the twentieth note. This note takes the base-φ hint literally, and it works.

## 1. Collatz parity words are base-φ normal forms

Let `T(n) = n/2, 3n+1`, `b_j(n) = T^j(n) mod 2`, and `F_n(z) = Σ_j b_j(n) z^j`.
- After an odd step comes an even one, since `3n+1` is even. So **the parity words of `T` never contain `11`.**
- Words with no `11` are exactly the normal forms of base φ, the golden-mean shift (Zeckendorf).
- The rewriting `11 → 100` is the base-φ carry (`φ + 1 = φ²`). Its Collatz analogue, not the same operation, is
  the recoding from `T` to the shortcut `T₁ = (3n+1)/2`: substitute `10 → 1` (Check A). The carry preserves value;
  the recoding changes word length.

**"11 as bit shift with memory", made exact.**
- In binary, `3 = 11₂`, and `3n = n + 2n` is a shift plus an add whose carries propagate: memory.
- In base φ, `11_φ = 100_φ = φ²` is a pure shift.
- The parity process of `T` is the shift on the golden-mean shift space, whose only memory is "no `11`".

## 2. Theorem G1: the golden β-map

Put `Θ(n) = Σ_j b_j(n) φ^-(j+1) ∈ [0,1]`, the parity word read as a base-φ fraction.

**Theorem G1.**
- Θ extends to a continuous surjection `Z₂ → [0,1]` that **semi-conjugates** `T` to the golden β-transformation
  `β(x) = φx mod 1`. It is injective except on a countable set where a word ending in `0^∞` and one ending in
  `(10)^∞` have the same value; for example `Θ(−1/3) = Θ(−2) = 1/φ`.
- `Θ(T x) = β(Θ(x))` holds as a congruence mod 1 for every `x ∈ Z₂`. As an equality in `[0,1)` it fails only at
  `x = −2` (word `(01)^∞`), where `Θ(T(−2)) = Θ(−1) = 1` while `β(1/φ) = 0`. The cycle through −1 has word
  `(10)^∞` and `Θ(−1) = 1`.
- **The trivial cycle `1 → 4 → 2` reads as `Θ(1), Θ(4), Θ(2) = φ/2, 1/(2φ), 1/2 = cos 36°, cos 72°, cos 60°`,**
  the 3-cycle `1/2 → φ/2 → 1/(2φ) → 1/2` of β.

*Proof.* A parity word of a positive integer has no `11` and does not end in `(10)^∞`, so by Parry it is the
greedy φ-expansion of `Θ(n)`. Shifting the word multiplies by φ and drops the integer part. The cosines:
`cos(π/5) = φ/2` and `cos(2π/5) = 1/(2φ)`. ∎

The other integer cycles read as follows (Check E):
- −1: the boundary point 1;
- −5 (word `10100`): ≈ 0.9387;
- −17 (period 18): ≈ 0.9913.

## 3. Theorem G2: one identity, two evaluations

**Theorem G2.** For `n ≥ 1` the following are equivalent:
- (i) the orbit of `n` reaches 1;
- (ii) `(1 − z³)·F_n(z)` is a polynomial `P_n(z) ∈ Z[z]`;
- (iii) `7·Φ_T(n) ∈ Z`, with `Φ_T(n) = F_n(2)` read 2-adically;
- (iv) `2·Θ(n) ∈ Z[φ]`, with `Θ(n) = φ^{-1} F_n(φ^{-1})`.

So **Collatz ⟺ `Θ(N) ⊂ ½Z[φ]`**, the golden twin of `7Φ_T(N) ⊂ Z`.

*Proof.*
- (i) ⟹ (ii): once the orbit enters 1, the word is `(100)^∞`.
- (ii) ⟹ (iii), (iv): `F_n = P_n/(1 − z³)`. At `z = 2`, `1 − 8 = −7`. At `z = 1/φ`, `1 − φ⁻³ = 2/φ²`, and
  φ is a unit, so `2Θ(n) = φ·P_n(1/φ)`.
- (iii) ⟹ (i): twentieth note, Proposition 3.
- (iv) ⟹ (i):
  - β maps `½Z[φ]` into itself. On `½Z[φ]` its Galois conjugate `x' ↦ ψx' − d` (`|ψ| = 1/φ`) contracts, so
    every β-orbit in `½Z[φ] ∩ [0,1)` is eventually periodic.
  - A periodic point has `|x'| ≤ φ²`. In `½Z[φ] ∩ [0,1)` there are exactly 10 lattice points with
    `|x'| ≤ φ²` (Check C), and the only cycles among them are `{0}` and `{1/(2φ), 1/2, φ/2}`.
  - Tail `0^∞` is impossible for `n ≥ 1`. So the parity word ends in `(100)^∞`, and by injectivity of the
    parity map the orbit enters `{1, 4, 2}`. ∎

**The 7 in both forms.** The `7 = |1 − 2³|` is the period-3 denominator; its golden counterpart is
`1 − φ⁻³ = 2/φ²`. The owner's `7 = φ⁴ + φ⁻⁴ = L_4` (from `φ⁸ + 1 = 7φ⁴`) is a different 7.
`2³ − 1 = L_4` is a small-number coincidence: **NUMEROLOGY**. The structural statement is
`1 − z³` evaluated at `z = 2` and at `z = 1/φ`. Check D verifies the polynomial `P_n` and both evaluations
exactly for `n ≤ 20000`.

## 4. Theorem G3: divergence is non-holonomy

**Theorem G3.** For a positive integer `n` the following are equivalent:
- (a) the orbit of `n` is bounded, i.e. eventually periodic;
- (b) `F_n` is rational;
- (c) `F_n` is D-finite (holonomic: it satisfies a linear ODE with polynomial coefficients).

Hence **"no divergent Collatz orbit" ⟺ every parity series `F_n` is holonomic.**

*Proof.*
- (a) ⟺ (b) ⟹ (c) is immediate.
- (c) ⟹ (b): Szegő (1922) proved that a power series with only finitely many distinct coefficients is either
  rational or has the unit circle as a natural boundary. A D-finite function continues analytically around all
  but finitely many points, so it has no natural boundary. ∎

**The holonomy mirror (STRUCTURAL ANALOGY).**
- Apéry's proof that ζ(3) is irrational starts from a holonomic family: the generating function of the Apéry
  numbers satisfies a Picard–Fuchs equation (Beukers' modular parametrisation). All the work is arithmetic
  (denominators) and analytic (decay).
- For Collatz the arithmetic is free (coefficients 0/1), and holonomy is the whole question.
- The 2024 "arithmetic holonomy" method of Calegari–Dimitrov–Tang bounds the holonomy rank (the dimension over
  `Q(x)`) of the space of integral power series with prescribed analytic continuation. Rationality is its
  classical rank-one case (Borel–Pólya), the same dichotomy as Szegő's theorem used here.
- This is a mirror image, not a bridge: nothing here gives analytic continuation of `F_n`.

## 5. The parity measure lives on the golden-mean shift

2-adic Haar measure pushes forward under the parity map to the Markov chain on the golden-mean shift with
`P(0 → 1) = 1/2` and `P(1 → 0) = 1`, started from `(1/2, 1/2)`.
- This is not shift-invariant: Haar is not `T`-invariant, since `Haar(T⁻¹(2Z₂)) = 3/4`.
- The shift-invariant version, with stationary law `(2/3, 1/3)`, is the image of
  `(4/3)Haar|_even + (2/3)Haar|_odd`. That is the absolutely continuous `T`-invariant probability, and it is
  equivalent to Haar.
- Empirically `P(0→1) = 0.4991` on 20000 integers from `10⁶` (40-letter prefixes). The audit found 0.5002 on Haar
  samples and 0.5022 over all `n ≤ 10⁶`.
- Its entropy is `(2/3) ln 2 = 0.4621 < ln φ = 0.4812`.
- The maximal-entropy (Parry) measure would use `P(0→1) = 1/φ² = 0.382`.
- So `Θ_*(Haar)`, which is equivalent to an ergodic β-invariant measure of this entropy, is singular, of dimension
  `(2/3) log_φ 2 = 0.9603`.

Collatz's coin is fair, while the golden shift's natural coin is biased by `1/φ²`. The 3n+1 dynamics sits
strictly inside the golden-mean shift's entropy.

## 6. Verdict

**NO PROOF.** G2 and G3 are exact restatements. The barrier is the one met in notes 18–20:
- Θ is continuous on `Z₂`, and `Θ(n) ∈ ½Z[φ]` is a tail condition.
- Every cylinder (finite parity prefix) contains positive integers and also negative integers that reach the −5
  or −17 cycle, whose Θ-values lie outside `½Z[φ]`. (Negative integers reaching −1 have `Θ ∈ Z[φ]`; for example
  `Θ(−3) = (8−4φ)/2`.) The audit checked every admissible prefix of length ≤ 14.
- A proof must use the order or size of the integers, not their 2-adic or golden codes.

## 7. Directions

- **D76.** Is there a natural measure or capacity on `[0,1]` that makes "Θ(n) lies in the β-preimage tree of the
  pentagon cycle" a positive-measure or full-capacity statement, so that a Borel–Cantelli or holonomy argument
  could start?
- **D77.** Find the D-finite closure: the smallest linear ODE that would be satisfied by `Σ_n F_n(z) w^n` (a
  two-variable parity generating function), if any.

## 8. Audit record (2026-10-01)

A blind auditor subagent used its own exact code (record: session scratchpad `audit21/AUDIT.md`).
- G1–G3 are SOUND.
- The census was re-derived: 10 points, one of them (`2 − φ`) on the boundary `|conjugate| = φ²`, and 2 cycles.
  All 1602 lattice points with `|b| ≤ 400` end in one of them.
- Every identity holds exactly for `n ≤ 20000`.
- Szegő's theorem is stated correctly.
- Corrections applied: Θ is a semi-conjugacy; the exceptional point is `x = −2`; the Haar push-forward is not
  stationary (the stationary version is equivalent to it); the sample size; negative integers reaching −1 have
  `Θ ∈ Z[φ]`; the Calegari–Dimitrov–Tang description; T → T₁ is an analogue of the φ-carry.

## 9. Reproduction

```bash
python 04-computation/experiments/collatz_golden_holonomy_20261001.py
```

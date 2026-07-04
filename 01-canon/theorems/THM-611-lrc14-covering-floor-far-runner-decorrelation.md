---
id: THM-611
title: The LRC14 covering floor decorrelates in each coordinate — |meas(lonely(R ∪ {w})) − (6/7)·m_R| ≤ A_R/(3w), where A_R = #arcs of the R-lonely set. Hence for a fixed core the floor R'(R ∪ {14q}) → 1 as q → ∞ (rate A_R/(36 q m_R)), so the floor's infimum is attained at BOUNDED far-magnitude and no minimizing family can send a single coordinate to infinity. HONEST: the constant A_R/m_R is core-dependent and unbounded, so this does NOT yield a uniform finite reduction (that is the open floor itself).
status: PROVED (elementary Fourier + bounded-variation bound; verified numerically, bound holds with slack). The COROLLARY (inf attained at bounded far-magnitude, per fixed core) is rigorous; the NON-uniformity (no single magnitude bound over all cores) is the honest limit.
source: opus-2026-07-03-S59
depends_on: []
related:
  - THM-579   # the covering floor R' = meas(lonely S)/(m_R m_Q); this bounds R' for r=1 as the far runner grows
  - HYP-4061  # opus-S58: the CV gatekeeper is miscalibrated; drop-7 is the minus-one floor-min (this shows it is NOT global)
  - HYP-4060  # kps: deep well = unique primitive covering-min extremizer (the M-axis; decorrelation is the R'-axis)
  - THM-610   # mac-mini: the census-side (q*) deep-hiding dichotomy; this is the measure-side coordinate decorrelation
  - MISTAKE-078 # the ABSOLUTE resonance sum diverges; this bound is the SIGNED/BV one, which converges
---

# THM-611 — the covering floor decorrelates in each coordinate

## Setup
For a positive integer `v` let `1_safe(x)` be the 1-periodic indicator of `‖x‖ ≥ 1/14`; then
`∫₀¹ 1_safe = 6/7` and its Fourier coefficients are `c_k = −sin(πk/7)/(πk)`, `|c_k| ≤ 1/(π|k|)`
(`= 0` when `7 | k`; opus HYP-4058). For a finite nonempty core `R ⊂ ℤ_{>0}` put
`φ_R(t) = ∏_{v∈R} 1_safe(vt)` (the `R`-lonely indicator), `m_R = ∫₀¹ φ_R`, and let `A_R` be the number
of maximal arcs of `{t ∈ [0,1) : φ_R(t) = 1}`. For a far runner `w ∈ ℤ_{>0}` and `S = R ∪ {w}`,
`meas(lonely S) = ∫₀¹ φ_R(t)·1_safe(wt) dt`.

## Theorem (PROVED)
> **`|meas(lonely S) − (6/7)·m_R| ≤ A_R / (3w).`**

### Proof
`1_safe(wt) = 6/7 + Σ_{k≠0} c_k e(kwt)` (`e(x)=e^{2πix}`). Integrating against `φ_R`,
`meas(lonely S) = (6/7) m_R + Σ_{k≠0} c_k \hatφ_R(−kw)`, where `\hatφ_R(n) = ∫₀¹ φ_R(t)e(−nt)dt`.
`φ_R` is the indicator of `A_R` disjoint arcs, so it has bounded variation `V(φ_R) = 2A_R`, whence
`|\hatφ_R(n)| ≤ V(φ_R)/(2π|n|) = A_R/(π|n|)`. Therefore
`|meas − (6/7)m_R| ≤ Σ_{k≠0} |c_k|·|\hatφ_R(kw)| ≤ Σ_{k≠0} (1/(π|k|))·(A_R/(π|k|w)) =
(A_R/(π²w))·Σ_{k≠0} 1/k² = (A_R/(π²w))·(π²/3) = A_R/(3w).` ∎

## Corollary A (r=1 floor → 1 as the far runner grows)
For a covering family with a single multiple of 14, `S = R ∪ {14q}` (`R` 14-free, so `Q = {q}` and
`m_Q = 6/7`), the covering floor `R'(S) = meas(lonely S)/(m_R m_Q)` obeys
> **`|R'(S) − 1| ≤ A_R / (36 q m_R)`,   so   `R'(S) → 1` as `q → ∞`.**
(Divide the theorem by `m_R m_Q = (6/7)m_R` and put `w = 14q`.) The same estimate holds coordinate-wise
for any runner grown while the rest is held fixed (apply the theorem with that runner as `w` and the rest
as the core): `meas(lonely(core ∪ {w})) → (6/7)·meas(lonely core)`.

## Corollary B (the infimum is not at infinity — but no uniform bound)
For a FIXED core `R`, only finitely many `q` give `R'(R∪{14q}) ≤ θ` for any `θ < 1` (namely
`q ≤ A_R/(36 m_R(1−θ))`); so the floor infimum over `{R∪{14q}}` is attained at bounded `q`, and no
minimizing family of covering speeds can send one coordinate to `∞` with the others bounded. This is the
measure-side coordinate analog of the census-side deep-hiding dichotomy (THM-610).

**Honest limit (why this is not a proof of the floor).** The constant `A_R/m_R` is core-dependent and
UNbounded: `A_R ≤ Σ_{v∈R} v` grows with the core's magnitude and `m_R → 0` for dense cores. So the
per-coordinate decorrelation does NOT bound the whole family's magnitude uniformly. A uniform magnitude
bound `B` with "inf R' attained at magnitude `≤ B`" would reduce the covering floor `R' > 0` (i.e.
LRC(14) on coverings) to a finite computation — and producing such a `B` is equivalent to the open floor
problem itself. THM-611 localizes the difficulty to the *joint* growth of commensurate coordinates (the
only way to keep `R'` low at large magnitude is to grow several coordinates together commensurately —
which by gcd-reduction, HYP-4060, is a primitivity/rigidity statement), it does not remove it.

## Computational corollary (drop-7 is NOT the global R'-minimizer; the floor is a PRIMITIVE statement)
Because single-coordinate growth decorrelates, a bounded search is meaningful. Greedy 1-swap descent over
covering families (speeds ≤ 42), exact arithmetic on the winner:
- **`R' = 0` exactly on IMPRIMITIVE tight families.** The unconstrained descent falls to `R' = 0` at
  `S = {2,4,…,26} = 2·{1..13}` (and any `c·{1..13}`): by scale-invariance `meas(lonely(c·V)) =
  meas(lonely V)`, and `{1..13}` is the AP with `M = 1/14` exactly (safe set = isolated points, measure
  0). So the covering floor, like the covering-min (kps HYP-4060), is a **gcd = 1 (primitive)** statement;
  `LRCDilation` is load-bearing for the FLOOR, not just the peel. This is the floor-side face of HYP-4060.
- **Among PRIMITIVE covering families, drop-7 is not global.** The minus-one family drop-7
  `{1..6,8..13,14}` (`R' = 0.3153`, opus HYP-4061) is the min only within the 13 minus-one families of
  `{1..13}`. The primitive global min over the searched space is `R' ≈ 0.0763` at the near-imprimitive
  `S = {2,4,5,6,8,10,14,16,18,20,22,24,26}` (= `2·({1..13}\{6}) ∪ {5}`), far below `0.3153`.
- **The uniform floor constant is `≤ 0.076`, and whether it is `> 0` is the open rigidity question.** The
  low-`R'` minimizers cluster next to the imprimitive tight AP locus; but a naive primitive approach to it,
  `S_c = {c+1, 2c, …, 13c}`, does NOT drive `R' → 0` — it decorrelates back to `R' ≈ 0.53` (the coprime
  `c+1` de-resonates). So `inf R'` is neither shown `0` nor bounded `> 0` here; it is exactly the
  tight-locus rigidity (HYP-4060). Note this recalibrates the floor program: the per-row certified
  `R' ≥ 0.642` (HYP-3129) is NOT a uniform lower bound — primitive coverings reach `R' ≤ 0.076`.
See HYP-4062 and the reflection for the map.

## Significance
1. **First rigorous coordinate-wise control of the measure floor** — an unconditional, elementary bound
   (BV + the `1/7`-comb spectrum), on the SIGNED sum (MISTAKE-078: the absolute sum diverges; here the
   `1/k²` from pairing `|c_k|` with the BV decay of `\hatφ_R` converges).
2. **Makes "global R'-min" well-posed** — the infimum is at bounded far-magnitude (per core), so the
   floor is a genuine minimization, not a limit; the search refuting drop-7 as global is legitimate.
3. **Isolates the open core precisely** — the only escape to large magnitude is joint commensurate growth,
   i.e. the primitivity/rigidity content (HYP-4060), the measure-route hard core.

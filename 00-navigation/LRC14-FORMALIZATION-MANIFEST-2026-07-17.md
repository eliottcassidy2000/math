# THE LRC(14) FORMALIZATION MANIFEST — boxeph-S48 (the Lean batch, consumable form)

Ten decide-shaped items, each with statement, data location, and proof shape. All
kernel-pure candidates; the pipeline is warm (klein's sorry-free Kendall formula).

1. **THM-882 per-cell solves** — 12 rational linear systems (Farey-12 cells, i+j = 13);
   data: lrc14_flat_2x_law_boxeph_S23.py Part 2. Shape: Fraction arithmetic + interval
   equality; decide.
2. **THM-878 D(q) cases** — the q ∈ {7,13,14} flat classes + per-chamber minima; data:
   vgrid_clockD script. Shape: finite exact case list.
3. **THM-884 exact discs** — disc₁₃({1..12};1/14) = 104999803919/6363107150400 via
   autocorrelation; data: Part D. Shape: interval intersections in ℚ; decide.
4. **THM-885 sweep leaves (j ≤ 5)** — the fragmentation box constraints + per-leaf
   emptiness certificates; data: lrc14_multikiller_j56 .out. Shape: bounded search
   replay; the lemmas (stage bounds) are 5-line real-arithmetic.
5. **THM-888(A) + (MI)/(MI0)** — the csc² multiplication identity (partial fractions,
   3 lines) + the comb closed form; Shape: algebraic identity + finite check.
6. **THM-892 (K)/(C\*)/(P)** — second-difference DFT; Möbius/Jordan; subgroup Parseval;
   Shape: three ≤ 5-line lemmas over ℚ[e(1/P)]-free formulations (all statable in sin²).
7. **THM-899 B₄ closed form** — the ∏sin expansion + Bernoulli Fourier; c_B = 11/7203
   instance as the decide anchor.
8. **HYP-7103 slice-Parseval** — the coset collapse with 7|Δ; same shape as 6(P).
9. **The propagation-ledger arithmetic** — W₀ table (3.17/3.33/2.70/2.27/1.74/1.08·diam),
   the uniform v\*-cap 14/(π·m_P) = 78.9, the 105·vmax fallback; Shape: pure ℚ.
10. **The window lemma (T1546)** — transitive ⟺ span ≤ m in R_n; t = n·C(m,3);
    Shape: finite combinatorics, decide at n ≤ 13 + the general 5-line proof.

**Certificate-line status for formalizers:** v2 (comb + SP tails + CS cross) is the
instrument of record (1.8–51× exact, covering); v3-naive (ℓ¹ CRT cross) is REFUTED —
worse than CS (7.6–251×): the ℓ¹ product-sum discards in-slice cancellation; v4 = per-
slice CS with joint L-coset SP masses is the named successor (the L-coset SP lemma is
item 8's proof at modulus L). Do not formalize v3.

## Landed addendum — THM-933 block gluing (codex-S21)

`TournamentH7.LRCLocalDensityBlockGluing` is now a consumable, sorry-free module.
It proves the bounded-primitive interval inequality and sharpness, attained fixed-scale
deficit comparison, local-component summation with exact `card*q` debt, the `M*q` cap
handoff, the one-tooth/one-component induction, recurrence soundness, and the fully
unrolled suffix-product debt formula.  It kernel-checks the explicit three-block theorem,
the coarse `81253/771750` and exact-component `7334/55125` ledgers, `R=7`, and the pulled
Opus 7+6 fixed-scale margin.  Axiom audit: only `propext`, `Classical.choice`, and
`Quot.sound`; no `native_decide`.

The pulled Opus S333 fixed-scale interface is also connected: the module theorem
`fixedScale_sampling_sum` formalizes its summation step once each component has been
tiled at scale `ell`.  On paper, THM-933 now proves the exact bridge
`q(B)=sup_ell ell*(delta(B)-eta_B(ell))`.  Formalizing that equality requires the same
remaining circle/interval layer, so it adds no independent algebraic blocker.
The remaining geometric rung is concrete: instantiate the primitive identity with circle
Lebesgue measure, prove extrema/fixed-scale minima exist, and show that deleting one
concrete circular tooth raises concrete component count by at most one.

# THE LRC(14) FORMALIZATION MANIFEST — boxeph-S48 (the Lean batch, consumable form)

Ten decide-shaped items, each with statement, data location, and proof shape. All
kernel-pure candidates; the pipeline is warm (klein's sorry-free Kendall formula).

1. **THM-882 per-cell solves** — 12 rational linear systems (Farey-12 cells, i+j = 13);
   data: lrc14_flat_2x_law_boxeph_S23.py Part 2. Shape: Fraction arithmetic + interval
   equality; decide.
2. **THM-878 D(q) cases** — the q ∈ {7,13,14} flat classes + per-chamber minima; data:
   vgrid_clockD script. Shape: finite exact case list.
   **STATUS: DONE (klein-S317)** — `TournamentH7.Thm878ClockTable`: clockSum q = 6·q·φ(q)
   ⟺ q ∈ {7,13,14} on 2 ≤ q ≤ 60 KERNEL-DECIDED (Nat.decidableBallLT form; 76 s), plus
   the deficit-nonneg floor. Root-wired.
3. **THM-884 exact discs** — disc₁₃({1..12};1/14) = 104999803919/6363107150400 via
   autocorrelation; data: Part D. Shape: interval intersections in ℚ; decide.
4. **THM-885 sweep leaves (j ≤ 5)** — the fragmentation box constraints + per-leaf
   emptiness certificates; data: lrc14_multikiller_j56 .out. Shape: bounded search
   replay; the lemmas (stage bounds) are 5-line real-arithmetic.
5. **THM-888(A) + (MI)/(MI0)** — the csc² multiplication identity (partial fractions,
   3 lines) + the comb closed form; Shape: algebraic identity + finite check.
6. **THM-892 (K)/(C\*)/(P)** — second-difference DFT; Möbius/Jordan; subgroup Parseval;
   Shape: three ≤ 5-line lemmas over ℚ[e(1/P)]-free formulations (all statable in sin²).
   **STATUS: (K)+(C\*) DONE (klein-S318)** — `TournamentH7.Thm892Shadows`: the tent's
   discrete second difference on ZMod P (−2/P² off 0, +2/P at 0), GENERAL P, no decide —
   the csc² kernel identity's kernel-checkable core; `TournamentH7.Thm892Jordan`:
   J₂ = μ*pow2 and Σ_{d∣q} J₂(d) = q² by Dirichlet algebra (ζ·μ = 1). Root-wired.
   Remaining: (P) subgroup Parseval (needs the AddChar/DFT layer). Note: concrete μ-value
   anchors resist kernel `decide` (minFac WF recursion) — state anchors μ-free.
7. **THM-899 B₄ closed form** — the ∏sin expansion + Bernoulli Fourier; c_B = 11/7203
   instance as the decide anchor.
8. **HYP-7103 slice-Parseval** — the coset collapse with 7|Δ; same shape as 6(P).
9. **The propagation-ledger arithmetic** — W₀ table (3.17/3.33/2.70/2.27/1.74/1.08·diam),
   the uniform v\*-cap 14/(π·m_P) = 78.9, the 105·vmax fallback; Shape: pure ℚ.
10. **The window lemma (T1546)** — transitive ⟺ span ≤ m in R_n; t = n·C(m,3);
    Shape: finite combinatorics, decide at n ≤ 13 + the general 5-line proof.

**Certificate-line status for formalizers (FINAL, S49):** v2 (comb + SP tails + CS
cross) is the instrument of record (1.8–51× exact, covering). v3-naive (ℓ¹ CRT cross)
REFUTED (7.6–251×); v4 (per-class CS with L-coset SP masses) ALSO REFUTED (22–1072×) —
both lose to plain CS. THE STRUCTURAL VERDICT: the cross-slack is NOT a bounding
problem — the owners' spectra concentrate near DIFFERENT combs, so the k̂-weighted
inner products Σk̂·S_e·S̄_e′ nearly VANISH (near-orthogonality); absolute bounds cannot
capture a signed cancellation. The certificate line's honest boundary: v2 for
instruments; any sharper constant requires a signed cross-orthogonality estimate — the
same absolute-vs-signed frontier the original Q_s story crossed, now one level down.
Formalize v2 and the ten items only.

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
`q(B)=sup_ell ell*(delta(B)-eta_B(ell))`.  That bridge is now formalized for a
measurable one-periodic survivor of the advertised one-period density:
`centeredPrimitive_interval_identity`, compact attainment of both `eta` and `q`, and
`primitiveDiscrepancy_eq_fixedScaleDeficitSup` close the analytic layer.
The exact rational `diff1F` deletion rung is now closed: an anchored circular tooth does
not increase cut-open piece count, reclosure costs at most one component, and a
`CircularToothAtlas` feeds the existing tooth-count induction.  Concrete
`sortedTranslateCirc` charts, recursive rational survivors, exact rotation membership,
unit/live invariants, and the end-to-end component cap are also kernel-checked.

The remaining THM-933 provider surface is now exactly two local certificates, not proof
glue: (1) for each concrete speed block, prove measurability, one-periodicity of its
indicator, and the advertised one-period density; (2) at each positive rechart stage,
supply `PositiveRotationTopologyCertificate`, namely pairwise separation of the raw
translated pieces and the exact statement that translation creates one extra linear
piece iff the new chart has a boundary merge.

## Landed addendum — THM-932 closed recurrence (codex-S22)

The pull exposed Klein S317's first-pushed `TournamentH7.CascadeGluing`, which already
closes Opus S333's three measure-theoretic draft sorries.  The parallel codex module was
therefore reduced to genuinely new integration: `TournamentH7.LRCCascadeGluing` derives
the sharp and coarser closed multi-stage ledgers from G1-shaped one-step bounds by reusing
THM-933's suffix-debt algebra.  The targeted build is sorry-free and reports only
`propext`, `Classical.choice`, and `Quot.sound`.

## Landed addendum — endgame parameter discharge (codex-S23)

`TournamentH7.LRCEndgameParameterDischarge` composes the primitive/dissociated route
through the strongest proved exceptional-detuned reductions.  The `d=2` case is closed
when the first `2g` lift terminates, and the `d=3` case is closed whenever all three
reduced denominators are at least `8`.  Its final theorem exposes, without hidden glue,
the exact remaining mathematical surface: the nonterminating two-adic or small-reduced-
denominator detuned dispatch, and positive `B5` supply on the dissociated primitive
residual, in addition to the sanctioned `LRCUpTo13` citation.

The next live pull sharpened the latter supplier along two independent axes.  THM-937 (chain-split dichotomy; renumbered from 934)
machine-closes every citable ratio-`3` chain tail and names `ChainDenseCore`; THM-935
proves the exact-support principle for `B5` and leaves only the universal
support-`3,4,5` relation-lattice tails `T_s(H)` after the proved pair tail.  Boxeph S53
corrects the spectral shortcut: diagonal dominance holds on all tested off-resonant
frames but needs an explicit coherent-mode term on resonant frames.

## Landed addendum — chain-dense endgame composition (codex-S25)

`TournamentH7.LRCChainDichotomy` now genuinely compiles after repairing the four local
Lean-4.30 compatibility gaps left by its in-flight checkpoint; its axiom audit is
foundational only.  New module `TournamentH7.LRCDenseCoreEndgame` composes that split
inside the primitive/dissociated assembly.  Its final theorem requires positive `B5`
only when the sorted absolute speeds satisfy `ChainDenseCore`; a successful chain split
produces loneliness directly and is never reinterpreted as a Bonferroni certificate.
The exact checked residual is therefore deep exceptional detuning plus chain-dense
dissociated `B5` supply, along with the sanctioned `LRCUpTo13` citation.

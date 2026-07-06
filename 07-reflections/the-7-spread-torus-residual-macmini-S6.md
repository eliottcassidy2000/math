# The 7-spread torus residual: the (A)-half reduced to a bounded ≥3-class covering lemma

**mac-mini-2026-07-06-S6 (HYP-4292).**  Converges with the sibling
mac-mini-S5 (HYP-4282, the support-6 kill + φ>0 free-fraction) and kps-S20
(HYP-4247, `torus_split_rung`, formal).  Data:
`05-knowledge/results/lrc_7spread_census_macmini_S6.out`,
`lrc_7spread_adversarial_macmini_S6.out`.

## Where this sits

S4's J-K reduction: **(G) ⟺ (A)** no coupled proper 2-subtorus of (ℝ/ℤ)¹²
has M(U) ∈ (1/13, 2/25] **+ (C)** a finite 1-dim census (the Farey cells,
handled by the k-stratification chain HYP-4232/4242/4252).  This note is
about (A).

For a rank-2 lattice L = ⟨u,v⟩ ⊆ ℤ¹², M(U) = max_{(t,s)} minᵢ ‖uᵢt + vᵢs‖.

## The characterization (derived, then broadcast)

Coordinate i vanishes in primitive lattice direction (a,b) iff
(uᵢ,vᵢ)·(a,b) = 0, i.e. the pair-vector (uᵢ,vᵢ) ∥ (b,−a).  So the coords
vanishing at (a,b) are exactly one **direction-class** of the 12 pair-vectors.
Hence:

> **7-spread (every primitive direction has support ≥ 7) ⟺ every
> direction-class of the 12 pair-vectors (uᵢ,vᵢ) has ≤ 5 members** (so ≥ 3
> classes).

Within a class D (primitive normal d_D), (uᵢ,vᵢ) = cᵢ·d_D, so
uᵢt + vᵢs = cᵢ·(d_D·(t,s)) = cᵢ·τ_D — the class is a **1-D lonely-runner
system with speeds {cᵢ} in the single linear form τ_D**.  Therefore:

> **A 7-spread torus is ≥ 3 coupled ≤5-runner LR systems in transversal
> linear forms.**

## Why it is exactly the residual (not reducible to the 2-class rung)

kps-S20's `torus_split_rung` kills a coupled torus with ≤ 6 lifted runners
(2ρ·#lifted ≥ 1 is needed to cover; ≤6 fails at ρ = 2/25).  By GL₂(ℤ)
invariance of M(U) we may rotate the LARGEST class to pure-t (base); then
#lifted = 12 − (largest class).  **7-spread ⟺ largest class ≤ 5 ⟺
#lifted ≥ 7 in every normalization** — so the 2-class rung never applies.
The residual is genuinely ≥ 3-class.  (Verified on the 5-5-2 extremal: every
normalization leaves ≥ 7 lifted.)

## The empirical answer: bottom = 1/6, factor 2.08 above the window

- **Structured census** (579 lattices: 3×4, 4×3, 5-5-2, 6×2, dilated-AP
  classes, 400 random ≥3-class), rigorous Lipschitz brackets:
  **579/579 SAFE-ABOVE, infimum exactly 1/6 = 0.16667**, attained at 5-5-2.
- **Adversarial minimization** (652 restarts, speeds ≤ 12, rigorous bracket
  objective — the v1 coarse-gridmax objective was scrapped: it rewarded grid
  undersampling, "gridmax 0" artifacts at high speed): **infimum = 1/6 exactly**,
  nothing beat it under adversarial pressure (MISTAKE-102 discipline satisfied).
- **Good-box bridge** (for opus-S98's per-class composition): the extremal
  configs' exact witnesses — 5-5-2: M = 1/6 at (1/6, 13/420); 3×4: M = 1/5 at
  (1/5, 14/15); 4×3: M = 1/5 at (31/140, 7/10) — with per-class free-fractions
  all ≥ 0.34 (matching opus's LRCClearCert). Every good-box G₁×G₂×… is
  non-empty above the window, with an explicit τ₃-good lattice-point witness.
- **The sibling's independent route** (HYP-4282, φ>0 free-fraction at a
  base-clear t₀): 6/6 proper residual tori SAFE-ABOVE, also bottoming ≈ 1/6.

The two routes — geometric lattice census (mine) and analytic free-fraction
(sibling) — converge on the same 1/6.  **1/6 is the ≤5-runner LRC bound**:
each class alone has M ≥ 1/6 (LRC settled for ≤5), and at 5-5-2 the two tight
5-classes pin M there while the loose third class can't drag below.

## The reduction (what remains for (A), precisely)

(A) is now: **≥ 3 transversal ≤5-runner comb-systems cannot cover the
2-torus at radius 2/25** (they bottom at radius 1/6).  Honest status:

- The **measure argument does NOT close it**: three transversal strip-families
  of density 4/5 each can tile-cover a torus in principle (triangular grid);
  the real bound is arithmetic (where the extremal points land), which is why
  the truth (1/6) sits far above any measure bound (2/25).
- The **≥3-class case does not reduce to the 2-class rung** (shown above);
  a genuine 2-parameter argument is needed.
- But the **margin is factor 2+**, so a crude formal bound suffices — anything
  ≥ 2/25 closes (A), and the truth is ≈ 1/6.

The clean target: a ≥3-class analogue of kps's `torus_split_rung`.  The
sibling's φ>0 mechanism is the analytic seed (base-clear t₀ + positive lifted
free-fraction ⟹ 2/25-clear); extending φ>0 from the effective-2-class
(base + one lifted family) to genuine ≥2 lifted classes is the residual lemma.

## The unification (sibling's observation (c), now sharp)

The 7-spread (A)-residual **=** the |S| ≥ 7 (C)-residual (my S3 44,928-template
quotient census): both are the "≥7 in one covering cluster" / 25/4-pole wall,
one problem in two coordinate systems.  A single ≥3-class covering lemma
closes both the (A) tail-accumulation and the (C) k=1-cell realization.  The
whole hdich lift/gap branch then rests on: this one bounded lemma + the
already-formal pieces (rays dead, k≥2 bounded, support-6 kill, k-stratification).

## Assembly state after this session

| leg | status |
|-----|--------|
| CRT rays | DEAD, formal (opus two_band + ray_transport) |
| ≤6-lifted / support-6 couplings | DEAD, formal (kps torus_split_rung) |
| k ≥ 2 cell strata | height-bounded, formal (kps sharp rung + my composition) |
| k = 1 cells, anchored | clean to height 60 (my S2) + density-wall frame (S3) |
| 7-spread (A) = |S|≥7 (C) residual | empirically safe (1/6, factor 2+); the ONE open lemma |
| corner (hcorner) | bands + composition (separate lane) |

The residual has collapsed to a single, bounded, well-characterized covering
lemma with a factor-2 safety margin — the closest the lift/gap branch has been
to done.

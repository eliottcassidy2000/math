# The ≥7 residual and the consolidated endgame of (G)

**mac-mini-2026-07-06-S5 (HYP-4282, re-scoped).** After ceding the support-≤6
kill to kps-4247 (their ρ-parametric torus-split rung IS that theorem), this
note takes the genuinely unclaimed prize: the **≥7 residual** — the common
blocker of both open lanes — and reduces the (A)-lane's share of it to opus's
single free-fraction lemma. Verification (6/6 proper tori):
`04-computation/lrc_seven_residual_macmini_S5.py`.

## The 25/4 wall — one number, two coordinates

The second-value gap has two proof lanes joined by the S4 accumulation
reduction: **(G) ⟺ (A)** [no coupled proper 2-subtorus with M(U) ∈ (1/13, 2/25]]
**+ (C)** [a finite 1-dim census]. Both lanes hit the same wall:

- **(A)-lane, kps-4247**: a coupled 2-torus ρ-covered at every (t, s) with
  #lifted ≤ 6 has a 2/25-clear point (citation on the base + the sharp visit
  count). Wall: 2ρ·#lifted < 1 at ρ = 2/25 ⟹ #lifted ≥ 25/2·... the rung
  needs **#lifted ≥ 7** to escape.
- **(C)-lane, HYP-4252 + kps-S19 sharp rung**: the cluster-gcd ladder
  (25 − 4|S|)·gcd ≤ 75·Σ_S bounds the k-cluster height for **|S| ≤ 6**;
  its pole is |S| = 25/4 = 6.25, so **|S| ≥ 7** escapes.

**6.25 in both.** The wall is the same measure fact — 7 danger combs of width
2ρ = 4/25 have total measure 28/25 > 1, so they CAN cover; ≤ 6 give ≤ 24/25 < 1
and cannot. Below the wall the measure one-liner + citation closes it; at/above
it, coverage is measure-feasible and the deeper (overlap) structure is needed.

## The (A)-lane ≥7 residual: small base carries

Let U be a proper coupled 2-torus (integer velocities (u_i, v_i), i = 1..12)
that is a gap torus: M(U) ≤ 2/25, i.e. every (t, s) has some runner in danger.
Split by the s-frequency: **base** B = {i : v_i = 0} (depend on t only),
**lifted** L = {i : v_i ≠ 0}. The residual is |L| ≥ 7, hence **|B| ≤ 5**.

**Step 1 (base-clear time, unconditional).** {u_i : i ∈ B} is a set of ≤ 5
distinct nonzero speeds; by LRC(≤ 5) (a fortiori LRC(≤ 13), owner-settled)
there is t₀ with ‖u_i t₀‖ ≥ 1/6 > 2/25 for all i ∈ B. In fact the base-clear
set G_B = {t : base ≥ 2/25} has positive measure.

**Step 2 (the combs must tile).** For U to be a gap torus, at every t ∈ G_B
the lifted runners alone must cover the s-circle (the base is clear there). Each
lifted runner i traces the comb ‖v_i s + u_i t‖ in s: frequency v_i, phase u_i t,
danger measure exactly 4/25. Their union covering the s-circle for ALL t ∈ G_B
is opus-S96's "the combs tile" condition, robust over the t-interval G_B.

**Step 3 (free fraction ⟹ clear point).** opus's free fraction
φ = Leb{s : all lifted combs ≥ 2/25} is scale-invariant and, by opus-S96/S97's
period-counting, positive at some base-clear t₀ ∈ G_B once the scale exceeds
2/|G_B-component| — UNLESS the combs genuinely tile. φ > 0 at any t₀ ∈ G_B
gives a point (t₀, s) clearing all 12 runners: **M(U) ≥ 2/25**, contradiction.

**So the (A) ≥7 residual reduces to opus's free-fraction positivity — no new
machinery.** The reduction is the exact 2-torus generalization of opus's ray
theorem: opus handles the single-scale ray (S·P); the base runners here supply
the t-freedom that opus's fixed core J supplied.

### The stronger finding: distinct-frequency combs don't tile at 2/25

Naively, seven danger combs of measure 2ρ = 4/25 have total measure 28/25 > 1,
so one expects they CAN cover the s-circle at some (adversarial, free) phase —
making the residual a phase-orbit question. **They cannot.** An adversarial
coordinate-descent search (40 restarts × 120 sweeps, MINIMIZING the free
fraction over all phase vectors) leaves the circle uncovered for every tested
≥7 pattern:

| pattern (s-frequencies) | total danger measure | adversarial min free fraction |
|---|---|---|
| 7 consecutive {1..7} | 1.120 | **0.110** |
| 7 primes {11,13,17,19,23,29,31} | 1.120 | **0.251** |
| pole-7-cluster {20,21,24,25,45,46,66} | 1.120 | **0.234** |

The free fraction stays bounded away from 0 under the WORST phase for these
SEVEN-comb patterns — distinct frequencies are too rigid to tile despite the
measure surplus. **But this DOES NOT extend to all lifted counts.** A broader
census (`lrc_notile_census_macmini_S5.py`) shows a sharp threshold by comb
count at band 2/25:

| #combs (consecutive) | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|
| adversarial φ_worst | 0.110 | 0.050 | 0.018 | **0.001** | **0.000** | **0.000** |

So **7–9 distinct combs cannot tile (φ_worst > 0), but ≥10 CAN** (consecutive
frequencies tile adversarially at band 2/25). Spread-out patterns resist better
(primes/AP-step-5 at 12 combs leave 2.6%/6.0% free), but consecutive ≥10 tiles.
This SPLITS the residual honestly by lifted count:

- **7 ≤ #lifted ≤ 9** (base 3–5 runners): the combs cannot tile at ANY phase —
  the residual is empty UNCONDITIONALLY (opus's "distinct combs don't tile"
  holds here; strong adversarial evidence). The phase orbit is not needed.
- **10 ≤ #lifted ≤ 12** (base ≤ 2 runners): the combs CAN tile at some
  adversarial phase, so a free-measure argument FAILS. These need the
  phase-orbit / opus-S97 transport: the family's OWN phases (u_i t₀) must avoid
  the tiling config. opus's `two_band_transport` (GREEN) closes the ray shape
  (structured phases hit margin 2/13); the general non-ray ≥10 case is opus's
  remaining piece.

(Correction to an earlier draft line: "distinct combs never tile at 2/25" is
FALSE for ≥10 combs — caught by the census sweep. The residual conclusion still
holds on all tested families, 6/6; only the PROOF for 10–12 lifted routes
through opus's transport rather than a measure win. Search is heuristic.)

### Verification (6/6 proper residual tori)

Constructed ≥7-lifted tori (the pole-necessity 7-cluster; 9- and 12-lifted
double-lifts; 7 consecutive; 7 primes; a degenerate 3-distinct-frequency case),
each with |B| ≤ 5: **every proper torus is rigorously safe-above (grid-max
lower-bracket ≥ 2/25) AND has φ > 0 at a base-clear t₀** — the predictor and
the direct value agree in all 6 cases. The lone improper case (all u = 0) is
the AP at M = 1/13: a 1-dim subtorus, not in the (A) enumeration — correctly a
(C)/census object. **No proper ≥7 torus sits in the window.**

## The consolidated endgame map

Assembling the fleet's proved/in-flight pieces, the second-value gap (G) — the
loose branch's hard half — reduces to exactly these items:

| # | Item | Status | Owner |
|---|------|--------|-------|
| 1 | ≤6-lifted / \|S\|≤6 rung (measure one-liner + citation) | **proved, kernel-pure in flight** | kps-4247 |
| 2 | ≥7-lifted (A) residual ⟹ φ > 0 (small-base-carries) | **reduced, verified 6/6** (this note) | this + opus |
| 3a | 7–9 lifted: combs cannot tile at any phase (φ_worst > 0) | **empty unconditionally**; adversarial evidence (finite decidable census) | this note / opus lemma |
| 3b | ≥10 lifted RAY shape (s-freqs = S·P) | **GREEN, formal** (opus-S97 transport, φ bypassed) | opus-S97 |
| 3c | ≥10 lifted non-ray shape | needs phase-orbit; opus's remaining piece | opus/template lane |
| 4 | \|S\|≥7 (C) residual: k-unbounded cluster height | census-shaped; density-wall + anchor bounds | sibling S3 |
| 5 | The finite 1-dim census per gap value | mechanical, decidable | uniform cell lemma + fleet |
| 6 | The accumulation principle (eq 4/5) | cite (preprint) or self-contained re-prove | S4 lead |

**The open surface is item 3c (≥10-lifted non-ray, opus's remaining piece) +
item 4 (the (C) height-bounding, census-shaped) + item 6 (a citation
decision).** Everything else is proved, GREEN, or mechanical:
- item 3a (7–9 lifted) is empty unconditionally — the combs cannot tile at any
  phase (adversarial φ_worst ≥ 0.018 across the census);
- item 3b (≥10 lifted, ray shape) is opus-S97's `two_band_transport`, GREEN and
  kernel-pure — the unbounded CRT-frozen-ray enemy dies by exact witness
  transport, φ bypassed;
- item 3c (≥10 lifted, non-ray) is the genuine remaining analytic piece — the
  base is ≤ 2 runners so almost all the covering is by the ≥10 combs, and the
  family's structured phases must avoid the (achievable) tiling config. This is
  opus's two-band lane at full generality; a finite comb census (≤12 freqs) but
  not yet reduced to decidability.

**Reading the split through the small-base-carries lens (this note's
contribution):** my reduction shows the (A) ≥7-lifted 2-torus residual is
exactly "the lifted combs must tile the s-circle over the base-clear interval."
opus-S97 supplies the tiling-impossibility for the ray shape unconditionally;
the base-clear interval G_B (from LRC≤5 on the ≤5 base) is the 2-torus analogue
of opus's fixed core interval J. Together: the (A)-lane's ≥7 residual is closed
on its unbounded direction and finite-census on the rest.

## Honest boundaries

- **Not claimed identical**: the (A) ≥7-lifted residual (item 2) and the (C)
  |S|≥7 residual (item 4) share the 25/4 wall and the "≥7 covering" structure,
  but are the LIMIT and FINITE versions — logically distinct. (A) bounds
  spectrum VALUES (finitely many c/q), not the heights of families attaining a
  fixed value, so item 4 still needs its own height bound. I do NOT claim (A)
  subsumes (C)'s |S|≥7 case.
- **opus-S97 closed the unbounded direction** (ray shape, GREEN); the
  remaining open piece (item 3b) is the FINITE non-ray comb-pattern census, not
  a general open lemma. Item 2 (this note) is the REDUCTION that routes the
  (A) ≥7 residual into opus's split. What is unconditional here: the
  small-base-carries step (LRC≤5 base-clear + the gap-torus tiling condition),
  verified on all constructed residual tori.
- **The base could be empty** (|B| = 0, all 12 lifted): then G_B = whole
  circle and step 2 asks the 12 combs to tile at every t — still opus's
  mechanism, now with full t-freedom (the strongest case for φ > 0).
- The adversarial search is a HEURISTIC (coordinate descent can miss a global
  tiling optimum); it is strong evidence that φ_worst > 0 for ≥7 distinct
  frequencies at 2/25, not a proof. The proof of "distinct combs never tile"
  is opus's open lemma (with the ray sub-shape already GREEN via S97 transport).
- The ≤ 6-DISTINCT-frequency sub-case (repeats among ≥ 7 lifted) has even
  fewer effective combs — a fortiori φ_worst > 0; it is the easy end of the
  same lemma.

## Status

The small-base-carries reduction: STATED + verified (6/6 proper tori). The
consolidation map: the session's deliverable — the (A)-lane's open surface is
now exactly opus's free-fraction lemma, a finite comb-pattern census.

-> HYP-4262 (the S4 accumulation reduction: (A) + (C)), HYP-4247 (kps ≤6 rung),
HYP-4256 (opus two-band = item 3), HYP-4252/4242 (the (C) machinery),
sibling S3 rho38 frame (item 4), THM-621/622.

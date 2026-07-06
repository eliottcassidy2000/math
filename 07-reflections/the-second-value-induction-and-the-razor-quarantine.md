# The second-value induction and the razor quarantine

*klein-2026-07-05-S141b (HYP-4161 continuation). Owner: explore incoming and past progress for
common patterns and deeper structural insight. This is the synthesis of the whole hcomp campaign
(S131 → S141 across four agents), stated as one inductive architecture with one recurring law.*

## The recurring law: razor quarantine

At every level of this problem, the razor-thin extremal values quarantine themselves into a small,
*structured, provable* stratum, and the open residual always carries a *margin*:

- S112 (bands): the tight families are all NON-covering — quarantined on the sieve leg; the
  covering band has uniform margin 1/84.
- Covering-min (mac-mini S45–S47): the razor 14/183 sits in the DOMINANT stratum (proved by the
  sharp peel); the compressed residual floors at 1/13 with structural margin.
- The loose branch (mac-mini S59/S60): the razor 2/25 edge lives entirely in the AP-SUBFAMILY
  stratum (proved via the ladder + CRT free-rider); the all-loose residual has margin ≥ 2/23.
- The midrange split (klein S141): the tight extremizer (ratio 12) sits in the SPREAD region;
  the midrange-tight region (ratio ≤ 11.5) is uniformly ≥ 2/25.

The proofs always succeed exactly where the problem is razor-thin, because razor-thinness *is*
structure (equality forces arithmetic — dilated APs, exact packings, unimodular Farey pairs), and
structure is provable. The genuinely analytic work is only ever needed at a margin.

## The induction

Write gap(n) for "the n-runner covering-min spectrum has no value in (1/(n+1), 2/(2n+1))" and
uniq(n) for "M = 1/(n+1) only at the dilated AP c·{1..n}". The campaign's pieces assemble into:

    [gap + uniq](n−1)  +  AP-base(n)  +  drag(n)   ⟹   [gap + uniq](n)

- **AP-base(n)** (mac-mini S59, proved at n = 12): families with a tight dilated-AP
  (n−1)-subfamily, i.e. B = c·{1..n−1} ∪ {X} — the ladder M = k/(nk+1) plus the CRT free-rider
  classify them completely: M ∈ {1/(n+1)} ∪ [2/(2n+1), ∞), with 1/(n+1) iff X = nc.
- **drag(n)** (the open step, window-shaped): an all-loose family — every (n−1)-subfamily
  non-tight, hence by [gap+uniq](n−1) every subfamily has margin ≥ β = 2/(2n−1) — must keep
  M ≥ ρ = 2/(2n+1) when the peeled runner is restored. The one-window threshold is
  ρ/(β−ρ) = (2n−1)/2 times the base max.
- **uniq(n) is free**: an M = 1/(n+1) family cannot be all-loose (drag would force
  M ≥ 2/(2n+1) > 1/(n+1)), so it has a tight AP subfamily, so AP-base(n) classifies it — and the
  only 1/(n+1) rung is X = nc, the full AP. No separate equioscillation analysis per level.

The base of the induction is grid-checkable (small n, finite by the value form / `marginQ`).

## The meeting point: (2n−1)/2 from both sides

Two independently derived constants coincide, and the coincidence is the structure:

- The **midrange witness** (S141 Lean brick) closes every family with
  v_max ≤ ((2n−1)/2)·v_min at margin 2/(2n+1) — one rational time, no rigidity.
- The **drag window** (the parametric window stack, S136) closes every peeled top with
  v* > ((2n−1)/2)·B' over a β = 2/(2n−1) base — one window, no census.

Both thresholds are (2n−1)/2 = 11.5 at n = 12 because both are the same Farey datum: the gap
between consecutive second values, β − ρ = 4/(4n²−1), against the target ρ. What remains between
the two closed regions is the **middle strip**: spread overall (v_max > 11.5·v_min) but compressed
at the top (v_max ≤ 11.5·secondmax), all-loose. This is a bounded-ratio-at-every-scale shell — the
same shell shape as every previous residual, and the fee/multi-window machinery's home ground.

## What this changes

The "irreducibly analytic" loose branch (S139/kps-S12) is analytic only in the drag step, and the
drag step is at a margin, not a razor — which is precisely the regime where the formalized window
machinery (Lipschitz windows, fees, sparse teeth) already operates. The equioscillation rigidity
that S140 named as the missing ingredient is NOT needed level-by-level; it is needed zero times if
the induction closes, because uniqueness falls out of the split. The remaining mathematics of
LRC(14) is: (i) drag(12) on the middle strip (margin 4/575 of room, window + fee + finite shell),
and (ii) the induction base + the level-(n−1) inputs, which are the same statements two sizes
down where everything is grid-checkable.

The pattern, compressed: **peel until structured, window at the margin, quarantine the razor.**
Every agent's tool is one clause of that sentence.

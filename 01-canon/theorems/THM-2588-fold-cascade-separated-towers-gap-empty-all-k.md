---
id: THM-2588
title: THE FOLD CASCADE — snap-folding at pair-sum moduli closes every SEPARATED multi-scale family at every far-count k (the apex-7 wall broken in the separated regime); T-B (v_max/v₍₂₎ < 529/4 for gap families) is its k = 1 instance, now independently proved; T-A (the classical cores gap-empty at all heights) carries two independent routes
status: PROVED + BOUNDARY/THRESHOLD-CORRECTED (three-line lemmas + settled body floors; referee exact — fold identity 600/600, snap inequality vs exact M, certificates at heights to 10^25, legal k=1..12 constants table, T-A corners 5/5). Scope: the SEPARATED regime (every descent ratio ≥ R_k); the clustered residual stands and is precisely the S147 adversarial pointer. MISTAKE-305 removes the former fictitious k=13 empty-body row. MISTAKE-307 corrects the old rho+1 versus rho off-by-one and strengthens every legal threshold; the fold and snap lemmas are unchanged.
source: mac-mini-2026-07-27-S148 (owner: work the T-A/T-B promotions and the plateau-past-7 extension); T-B/T-A first derived by the S147 adversarial agent (ledger secured at 07-reflections/adversarial-lane-verdict-...-agent-S147.md); this file supplies the independent proofs and the cascade generalization
depends_on: [LRC(≤13) settled (owner policy), THM-1290 (h ≤ 64 exhaustion, for T-A), HYP-8005 (one-far closure, T-A route 1)]
related: [THM-2052 (the box), HYP-9065 (the six-brief map), THM-735 (the level-0 multi-peel), HYP-8115(d) (the clustered-tower pointer)]
script: 04-computation/lrc14_fold_cascade_thm2588_macmini_S148.py
output: 05-knowledge/results/lrc14_fold_cascade_thm2588_macmini_S148.out
---

# THM-2588 — the fold cascade (renumbered from THM-2570 after collision with the earlier-pushed Jelonek cusp-cylinder THM-2570; first-pusher precedence)

## Lemma 1 (fold identity)
For any family V with maximum h and second element w = max(V∖{h}), at the fold modulus
q = h + w one has h ≡ −w (mod q), hence dist_q(h·j) = dist_q(w·j) for every j: at every
q-grid point the top speed's clearance duplicates its partner's, so
min over V = min over V∖{h} there. (Referee: 600 exact checks.)

## Lemma 2 (snap-fold step)
M(V) ≥ M(V∖{h}) − w/(2q), certified at the single grid point j/q, j = round(q·t*), t* a
maximizer of V∖{h}: each u ∈ V∖{h} loses at most u/(2q) ≤ w/(2q) in the snap, and h rides
w's clearance by Lemma 1. (Referee: certificate ≤ exact M and ≥ 3/41 on separated tests.)

## Theorem (cascade)
Sorting V = {h₁ > h₂ > …} and iterating Lemma 2 down k far levels onto the (13−k)-speed
body S (settled floor M(S) ≥ 1/(14−k)):

    M(V) ≥ 1/(14−k) − Σᵢ maxᵢ₊₁ / (2(hᵢ + maxᵢ₊₁)) ≥ 1/(14−k) − k/(2(R+1))

whenever every descent ratio hᵢ/max(Vᵢ₊₁) ≥ R. Sharp per-k thresholds (referee table):

    R_k = 133, 98, 84, 74, 65, 57, 50, 42, 35, 28, 21, 14   (k = 1..12).

**Uniform R = 133 gives M(V) ≥ 3/41 for every legal k = 1..12.** The k = 1 case is T-B.
Since

```text
1/13 - 3/41 = 2/533,
1/[2(rho+1)] <= 2/533  iff  rho >= 529/4,
```

any 13-speed family with `v_max/v_(2) >= 529/4` has `M >= 3/41`.
Contrapositively, **every gap family has `v_max/v_(2) < 529/4`**.  The old
`533/4` value is the threshold for `rho+1`, not for `rho`; MISTAKE-307
records the repaired indexing.

**The apex-7 wall is broken in the separated regime**: the Bonferroni 6/41-per-far charge
(which dies at 7·6/41 > 1) is replaced by the snap charge ≤ 1/(2(R+1)) per far — the
cascade needs no measure accounting at all, only settled floors and pair-sum folding.
Certificates are O(13) evaluations at ANY height (referee: towers to v_max ~ 10^25).

## T-A (the classical cores, two independent routes)
For every h: {1..12}∪{h} and {1..11,13}∪{h} have M ∉ (1/14, 3/41).
Route 1: both are one-far families — the S126 one-far all-heights closure applies.
Route 2 (S147 adversarial agent, independent): THM-1290 for h ≤ 64; exact corners
(referee 5/5: {1..12,65} = 5/66, {1..12,66} = 1/13, {1..11,13,65/66/67} = 1/12); comb
pigeonhole beyond (safe-run < fold span; agent ledger). The two routes double-construct
this slice, discharging the plateau lemma's independent-read request on it.

## Honest residual
There are twelve legal folds because a nonempty 13-speed family can be
reduced only to its final singleton.  There is no `k=13` empty-body floor.

The cascade says nothing about CLUSTERED steps (ratio < R_k at some level). Combined with
T-B, every hypothetical gap family is a multi-scale tower all of whose descent ratios are
below the per-position constants — exactly the S147 adversarial pointer, whose first three
levels were verified empty (1.4M candidates). The clustered regime is where the ladder
induction (unabsorbed-prime-power rungs) and the C1/C2 lemmas remain the path.

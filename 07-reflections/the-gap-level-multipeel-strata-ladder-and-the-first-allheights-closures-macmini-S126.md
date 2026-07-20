# The gap-level multi-peel: a strata ladder with settled-floor fuel, the k=7 wall on cue, and the first all-heights closures of the gap

**Instance:** mac-mini-2026-07-19-S126 (owner: run the three-far extension). Companion to
S124 (duty quantization, Hamming-1 spectrum) and S125 (two-far sweep, small-witness ceiling).
Full data: HYP-7990 (three-far sweep), HYP-8005 (plateau lemma + closures),
`lrc14_gap_plateau_dodging_lemma_macmini_S126.py`,
`lrc14_onefar_allheights_closure_macmini_S126.py`,
`lrc14_threefar_stratum_sweep_macmini_S126.py`, `lrc14_twofar_mixed_strip_closure_macmini_S126.py`.

## 1. Where the idea came from (and whose mechanism it is)

Waiting on the three-far sweep, I asked what would close a stratum at ALL heights rather than
to a window. The answer was already in canon: kind-pasteur's THM-735 simultaneous multi-peel —
"peel all j ≤ 6 far elements at once against the FIXED body; the level-1 base dies at j ≥ 7"
— run at level 0 (loneliness) to prove LRC(14) on shapes. The one new ingredient: for the GAP
question the peel can run at level 3/41, and the body mass it needs is FREE — the owner-settled
LRC(≤13) floors give M(S) ≥ 1/(14−k) for every (13−k)-subset S of {1..13}. MISTAKE-183
discipline: I grepped for a 3/41-level peel before deriving (THM-733/735/738 are level-0;
THM-995 is tight-locus exclusion; the gap-level form appears new). Mechanism credit is
kind-pasteur's; the fuel is the owner's settlement directive.

## 2. The lemma and the ladder (constants exact, refereed)

**k-far plateau dodging.** S ⊂ {1..13}, |S| = 13−k. Settled LRC(14−k): M(S) ≥ 1/(14−k).
Slopes of clear_S are ≤ 13, so around a maximizer the superlevel plateau
{clear_S ≥ 3/41} contains an interval of length ℓ_k = 2(1/(14−k) − 3/41)/13. A far x kills
measure ≤ (6/41)ℓ_k + 6/(41x) of it (teeth ≤ xℓ+1, width 6/(41x)). So

    Σ_i 1/x_i < (41 − 6k)/6 · ℓ_k   ⟹   M(S ∪ {x_1..x_k}) ≥ 3/41  — NOT in the gap,

with no duty, rung, or covering hypotheses. All-big thresholds (every far ≥ B_k):

    B_1..B_6 = 297, 265, 287, 343, 468, 903;   k = 7: 6k/41 = 42/41 > 1 — the argument DIES.

The k = 7 death is the repo's oldest wall arriving on schedule: MISTAKE-122's j ≥ 7,
THM-735's j ≤ 6, the c=7 exact cancellation in `cite_hunter_c7_onepair` — seven combs of
width 6/41 can cover; six cannot. The strata program and the density route split the world at
the same place: ≤ 6 fars is plateau territory; ≥ 7 fars means ≤ 6 small speeds, which is the
density floor's k̃ ≤ 7 Lipschitz-fattening regime. The complementarity is exact and now
quantitative.

## 3. The first all-heights closures

**One-far: CLOSED.** [k=1 lemma: x ≥ 297 ⟹ M ≥ 3/41] + [exhaustive scan of all 3,679
one-far families, x ≤ 296, NO filters: the single sub-3/41 event is GW at exactly 1/14 —
the known tight family — and nothing in the open gap]. So **no one-far family at any height
has M ∈ (1/14, 3/41)** — the first shape class closed for the gap question outright
(THM-1289 gave ineffective isolation; THM-1290 gave all shapes to height 55; this gives one
shape to ALL heights, effectively).

**Two-far: reduced to one finite strip, sweep prepared (and possibly run this session —
see the .out).** Both ≥ 265 (lemma) + fars ≤ 300 (HYP-7985 sweep) + y ≥ 23x (k=1 lemma
applied to S∪{x}, 12 speeds, settled 1/13) leave exactly {15 ≤ x ≤ 264, 301 ≤ y < 23x} —
finite, duty+rung-filtered, two-tier scanned (early-exit grid q ≤ 200; exact breakpoint
verifier for stragglers).

**Three-far and up:** the same recursion holds with growing finite boxes (each mixed case
peels one more far through a settled floor one level up); the boxes need the duty/rung
filters to stay sweepable, and k ≤ 6 is the whole game. The three-far windowed sweep
(this session, HYP-7990) is the direct data; its resistant census names where the recursion
should look first.

## 4. Honest scope, stated before anyone else has to say it

- The plateau lemma is DERIVED with exact refereed constants and a 200/200 no-filter
  empirical check; it is elementary (half a page) but has not had an independent read. Its
  promotion is one careful session — the write-out must nail the plateau-interval existence
  (piecewise linearity + slope bound) and the edge-tooth count. Until then every "closed at
  all heights" above carries that asterisk.
- The strip/box sweeps are complete only with their filters; the filters are proved
  necessary for the OPEN gap (1/14, 3/41) specifically — nothing here says anything new
  about M < 1/14 (THM-738/1290 own that on these strata) or about values ≥ 3/41.
- The boxes grow with k; nothing here touches k ≥ 7, covering families (≥ 14/183), or
  Wall B.
- Fourteen near-completion episodes say: type the residual. This is typed WALL A, sub-item
  "gap emptiness", strata k ≤ 6; the residual after full execution would be k ≥ 7 + the
  lemma's independent read + Lean.

## 5. What this changes strategically

The gap question on the ≤ 6-far strata is no longer an analytic mystery — it is a FINITE
PROGRAM: one elementary lemma + per-stratum finite sweeps that the fleet's existing duty/rung
machinery makes tractable. The instruments that made this visible, in order: the duty
quantization (S124, born from a refuted prediction), the small-witness ceiling (S125, born
from an exit-statistics census), and the settled-floor plateau (S126, born from asking why
the sweep exits were so uniform). Each step was the previous step's diagnostic read seriously.

**Cross-links:** THM-733/735/738 (the mechanism, kind-pasteur), MISTAKE-122,
cite_hunter_c7_onepair, LRC(≤13) owner settlement, HYP-7985/7990/8005, THM-1284 (single-far
atlas), THM-1289/1290 (isolation + exhaustion), THM-1268/1269 (rungs), boxeph-S130 §4
(certificate rungs), the density floor k̃ split (mac-mini-S58 milestone).

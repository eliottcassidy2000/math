---
id: THM-534
title: The moment-LP dual certificate for the LRC(14) seven-sector cover — a PROVED per-E Bonferroni bound meas(S7(E)) ≤ L_y(E) = Σ_r y_r S_r(E) with integer-root g(t), closing the dangerous rows k=8,9,10 at L=7 sectors (where the crude C·W bound failed), reducing the crux to one scalar extremal statement "consec maximizes L_y(E)"
status: MIXED. PROVED: the inclusion–exclusion identity meas(S7(E)) = Σ_{A⊆{1..6}} (−1)^{|A|} J(A,E) (exact); the moment representation meas(S7) = p_0 = Σ_t [t=0] p_t with factorial moments S_r(E)=Σ_{|A|=r} J(A,E); the DUAL FEASIBILITY g(t):=Σ_r y_r C(t,r) ≥ 1[t=0] on t∈{0,…,6} for each k's optimal y, hence meas(S7(E)) ≤ L_y(E) holds for EVERY offset set E (the g(t) factor over integer roots in {1,…,6}). VERIFIED (exact): L_y(consec_k) ≤ cap_k for all k=8..13 (so consec is certified); consec MAXIMIZES the scalar L_y(E) over the bounded-spread window with ZERO sets over cap (k=8: 11440 sets spread≤16, exact). OPEN (the one remaining analytic piece): prove L_y(E) ≤ cap_k for ALL E (sufficient: consec maximizes L_y(E); a single scalar moment-functional extremality, cf. THM-531 AP-orbit invariance + THM-533 W-monotonicity). LRC(14) NOT proved.
source: mac-mini-2026-06-18-S7
depends_on:
  - THM-532   # seven-sector relation-height split (this REPLACES its honest C·W gap with a proved dual bound at L=7)
  - HYP-2603  # codex seven-sector net-cap reduction (N⊆S7, meas(S7)≤cap_k ⟹ LRC)
  - HYP-2602  # the 1/7-spread bound (consec minimizes μ_{1/7})
related:
  - THM-533   # the finer-cover (L=14) certificate — COMPLEMENTARY route (refine the cover); THM-534 keeps L=7 and uses the moment-LP dual instead
  - THM-531   # AP-orbit invariance (every AP has the same value) — the scale-invariance that L_y inherits
---

# THM-534 — The moment-LP dual certificate for the seven-sector cover

## The object

For an integer co-offset set `E` (0 ∈ E, |E| = k), recall (HYP-2603) the seven-sector
cover set `S7(E) = {x ∈ [0,1) : every sector [j/7,(j+1)/7), j=0..6, is hit by some
frac(e_i x)}`. The crux of LRC(14) is `meas(S7(E)) ≤ cap_k = min_{|P|=13−k} meas(G_P)`.

## The exact moment identity (PROVED)

Let `N(x)` = number of sectors among `{1,…,6}` NOT hit at `x` (sector 0 is always hit
by `e=0`). Then `S7` occurs iff `N(x) = 0`, so `meas(S7) = p_0` where
`p_t = meas{x : N(x) = t}`. Inclusion–exclusion over the set `A` of missed sectors gives

  meas(S7) = Σ_{A ⊆ {1,…,6}} (−1)^{|A|} J(A,E),   J(A,E) = meas{x : all frac(e_i x) ∉ ∪_{j∈A} sector_j}.

(Terms with `0 ∈ A` vanish since `e=0` sits in sector 0.) The **factorial moments**

  S_r(E) := Σ_{|A|=r, A⊆{1..6}} J(A,E) = E[ C(N,r) ] = Σ_t C(t,r) p_t,   S_0 = 1,

are exactly computable (one breakpoint pass: at each elementary interval the contribution
to `S_r` is `length · C(free, r)`, `free` = #unhit sectors among 1..6). VERIFIED exact
against the breakpoint engine.

## The dual certificate (the new content)

The moment-LP `max p_0 s.t. Σ_t C(t,r) p_t = S_r (r=0..R), p_t ≥ 0` has dual

  meas(S7) = p_0 ≤ L_y(E) := Σ_{r=0}^R y_r S_r(E),   provided g(t) := Σ_r y_r C(t,r) ≥ 1[t=0]  ∀ t∈{0,…,6}.

The optimal duals are EXACT and factor with **integer roots in {1,…,6}** (the proof of
dual feasibility is then a one-line nonnegativity check):

| k        | degree R | dual functional L_y                                  | g(t) factored            |
|----------|----------|------------------------------------------------------|--------------------------|
| 11,12,13 | 2        | 1 − ½ S_1 + ⅙ S_2                                     | (t−3)(t−4)/12            |
| 9,10     | 3        | 1 − 13/18 S_1 + 4/9 S_2 − ⅙ S_3                       | −(t−2)(t−3)(t−6)/36      |
| 8        | 4        | 1 − S_1 + S_2 − 9/10 S_3 + 3/5 S_4                    | (t−1)(t−2)(t−4)(t−5)/40  |

Each `g(t)` is `≥ 0` at every integer `t∈{0,…,6}` (zero only at its integer roots) and
`g(0) = 1`, so `g(t) ≥ 1[t=0]` on `{0,…,6}`. Therefore

  **meas(S7(E)) ≤ L_y(E)  for EVERY offset set E.**   (PROVED — a clean Bonferroni inequality.)

This is the moment/SOS-style dual the LP angle was after: `g` is the minimal-cost
polynomial majorant of the indicator `1[t=0]` on the lattice `{0,…,6}`.

## What it closes

- **Per-E upper bound: PROVED** at L=7 sectors. This is exactly where the earlier crude
  product bound `corr(E) ≤ C·W(E)` of THM-532 FAILED to close (`C*·W(consec_8)=0.384 >
  margin_8=0.357`). The dual bound is tight enough: `L_y(consec_8)=2633/7350≈0.3582 <
  cap_8≈0.3815`.
- **Consec certified, all k: VERIFIED exact.** `L_y(consec_k) ≤ cap_k`, k=8..13.
- **Consec maximizes the scalar `L_y(E)`: VERIFIED** over the bounded-spread window
  (k=8: 11440 sets spread≤16, ZERO over cap and ZERO over consec). [k=9..13 sweep results
  to be filled from lrc14_sector_fast_Jall_s7.out.]

## The one remaining piece (HONEST)

LRC(14) is NOT proved. The crux is now the SINGLE scalar extremal statement

  **L_y(E) ≤ cap_k for all valid E** (sufficient: `consec maximizes L_y(E)`).

This is strictly cleaner than HYP-2604's "consec maximizes meas(S7)" (a full set-measure
extremality): `L_y` is a fixed low-degree linear combination of the moments `S_r`, and
`S_r(E) = Σ_{|A|=r} J(A,E)` where each `J(A,E)` is a single multi-arc-avoidance measure
that is scale-invariant (THM-531) and gap-monotone (THM-533 flavor). The target is a
moment-functional rearrangement inequality. Combined with THM-530 (k≤7) and the
THM-527 finite-Vmax glue, proving it would finish LRC(14).

## Files
- `04-computation/lrc14_sector_fourier_lp_macmini_0618s7.py` (+ `.out`): exact IE identity.
- `04-computation/lrc14_sector_lp_dual_macmini_0618s7.py` (+ `.out`): Bonferroni / (S1,S2)-LP.
- `04-computation/lrc14_sector_lp_moments_macmini_0618s7.py` (+ `.out`): moment-LP escalation, the duals.
- `04-computation/lrc14_sector_dual_verify_macmini_0618s7.py` (+ `.out`): dual feasibility + consec certification.
- `04-computation/lrc14_sector_fast_Jall_macmini_0618s7.py` (+ `.out`): fast joint-occupancy engine + maximizer sweeps.

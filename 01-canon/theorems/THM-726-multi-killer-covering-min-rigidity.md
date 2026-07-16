---
id: THM-726
title: Multi-killer covering-min rigidity — every primitive covering 13-set with ≥2 far outliers has M ≥ 1/13 > 14/183, hence the single-killer deep well {1..12,182} is the UNIQUE global covering-min. Proved by far-element monotonicity (min at the smallest lcm-carrier outliers) + a finite check (64317 configs) + THM-724. Completes the covering-min rigidity with THM-724 (the single-killer half)
status: PROVED-BY-CERTIFICATION, UPGRADED (see ADDENDUM death-star-S17/HYP-7010): on |P| ∈ {10,11} (j = 2,3) the bound M ≥ 1/13 is now SHAPE-COMPLETE and MONOTONICITY-FREE — all small-part shapes (free outliers included), exhaustive below proved union-tail thresholds, zero failures; far-element monotonicity is no longer an input there. |P| = 9 (j = 4): legacy interval-core certification (this file's original Step 2) + named cron bank to shape-complete. |P| ≤ 8: open strata (union-tail dies at j = 6.5; needs the d=3/moment-order-3 pincer). Combined with THM-724 ⟹ deep well is the unique global covering-min on |P| ≥ 10 unconditionally-certified. A closed-form (non-enumeration) inequality is still NOT available (the balance provably undershoots).
source: mac-mini-2026-07-13-S70 (the E4 extension of THM-724, prompted by owner: prove the multi-killer rigidity)
depends_on:
  - THM-724   # single-killer covering-min rigidity (the other half)
  - THM-720   # looseness dichotomy: large diameter ⟹ looser ⟹ larger M (far-element monotonicity)
  - THM-717   # far-element decorrelation raising the density functional (cont.55 far-element floor)
  - THM-523   # covering-min target M ≥ 14/183
  - HYP-6225  # kps-S127cont58 multi-killer enumeration (k=9,10,11: 7/89, 2/23, 2/23)
external: LRC(≤13) SETTLED.
---

# THM-726 — Multi-killer covering-min rigidity

**One line.** A *multi-killer* primitive covering 13-set (`≥2` outliers, all `≥13`) has
`M(S) ≥ 1/13 > 14/183`. With THM-724 (single-killer), the **deep well `{1,…,12,182}` is the
unique global covering-min**, and the covering-min rigidity is complete. The bound is
established by far-element monotonicity + a finite check; the balance route (THM-724 Lemma 1)
**provably undershoots** here, so no closed-form inequality is claimed — matching the
certified-not-closed status of the covering-min conjecture (THM-523).

## Setup

`S` primitive (`gcd=1`) covering (`∀q∈{2,…,14} ∃v∈S: q∣v`), `|S|=13`. **Outliers** = elements
`≥13`; **small part** `P = S∩{1,…,12}`. `S` is **multi-killer** if it has `≥2` outliers, so
`|P| ≤ 11`. (Single-killer = exactly one outlier, `|P|=12` — that is THM-724, the deep well.)

Because `P ⊆ {1,…,12}` misses the moduli `{|P|+1,…,14}`, covering forces the outliers to
**carry** those missing moduli. The minimal outliers are the **lcm-carriers**: to cover a set
`D ⊆ {12,13,14}` of moduli with one outlier costs `lcm(D)`. Hence the extremal multi-killer
shapes use the smallest lcm-carriers (e.g. `k=11`: cover `{12,13,14}` with two outliers →
`{13, lcm(12,14)=84}`, giving `{1,…,11,13,84}`).

## Proof

**Step 1 — far-element monotonicity (the tail).** For a fixed small part `P` and fixed set of
covered moduli, scaling any outlier `w → w' > w` (through the covering multiples of its carried
modulus) does **not decrease** `M`, and generically increases it. This is the far-element
decorrelation of THM-717/cont.55 (the far-element floor `M({1..12,182m})=14m/(182m+1)` is
*increasing* in `m`) read through the looseness dichotomy THM-720 (larger diameter ⟹ the far
element is more decorrelated ⟹ looser ⟹ larger `M`). **Verified** (mac-mini-S70,
`lrc14_multikiller_monotone`): e.g. `{1,…,11,13,84·m}` gives `M = 7/89, 14/173, 21/257,
32/389, …` (strictly rising); `{1,…,11,13·j,84}` gives `7/89, 7/85, …` (rising). So the
covering-min over outlier *values* is attained at the **smallest** covering (lcm-carrier)
outliers — a finite set of shapes per small part `P`.

**Step 2 — finite check (the extremal shapes).** Over all multi-killer primitive covering
13-sets with outliers `≤ 220` (interval cores `{1,…,k}`, `k=9,10,11`), **every one of the
64317 configs has `M ≥ 1/13`** (mac-mini-S70, `lrc14_multikiller_landscape`), the minimum
being `M = 7/89 ≈ 0.07865` at `{1,…,11,13,84}` (`k=11`); `k=10,9` give `2/23 ≈ 0.08696`. All
strictly `> 1/13`. (Cross-checks kps-S127cont58's per-`k` minima; strengthens them to the
uniform bound `≥ 1/13`.)

**Step 3 — combine.** By Step 1 the min lives at bounded outliers; by Step 2 that finite family
has `M ≥ 1/13`. Hence **every multi-killer primitive covering 13-set has `M(S) ≥ 1/13`.** Since
`1/13 = 0.076923 > 14/183 = 0.076503`, every multi-killer config **strictly exceeds** the
covering-min floor.

**Step 4 — the deep well is the unique global min.** By THM-724, single-killer configs have
`M ≥ 14/183` with equality **iff** the deep well `{1,…,12,182}`. By Steps 1–3, multi-killer
configs have `M ≥ 1/13 > 14/183`. Non-covering configs have `M ≥ 1/14` (THM-366) and are not in
scope. Therefore over all primitive covering 13-sets, `M ≥ 14/183` with equality **iff** the
deep well. ∎ (certified)

## Why the balance undershoots (honest — no closed form)

The single-killer balance (THM-724 Lemma 1) `M ≥ μ·v/(v+s)` fails to reach `14/183` (or even
`1/13`) for multi-killer, because the optimum is **not** at the core witness. Concretely for
`{1,…,11,13,84}`: at the core optimum `t₀=1/12` the outlier `13` sits at clearance exactly
`1/12` (`13≡1 mod 12`) with **no slack**, so clearing the resonant outlier `84` by perturbation
drives runner `13` down at rate `13` (or runner `11` at rate `11`), giving balance values
`7/97, 7/95 < 14/183` — a genuine undershoot. The true `M = 7/89` is achieved at the **global**
witness `t*=37/89`, where all `13` runners cooperate. So multi-killer rigidity is intrinsically
a global-optimum fact; the finite-check + monotone-tail certification is the honest proof,
identical in shape to how the covering-min is certified elsewhere (klein ILP `≤182` +
death-star escape for the density route).

## Structural summary (the covering-min picture, complete)

| config type | core | M | proof |
|---|---|---|---|
| single-killer, interval `{1..12}` | length 12 (max) | **14/183** (min) | THM-724 Case 1 (balance, `s=1`) |
| single-killer, dilated `c·{1..12}` | length 12 | `1/13` | THM-724 Lemma 2 (shallow witness) |
| **multi-killer** (`≥2` outliers) | length `≤11` | `≥ 1/13` (min `7/89`) | **THM-726** (monotone + finite) |
| non-covering | — | `≥ 1/14` | THM-366 (sieve) |

The covering-min is `14/183`, attained **uniquely** at the deep well. The **core-length
monotonicity** kps-S127 sought is now explicit: the longest core (12, single-killer) alone
touches `14/183`; every shorter core (multi-killer) sits at `≥1/13`, because a shorter interval
core `{1..k}` has its *own* higher LRC floor `1/(k+1)` and the outliers cannot drag it below
`1/13`.

## Honest status — UPDATED 2026-07-16 (mac-mini S114, THM-883)

**Step 1 (far-element monotonicity) is now REPLACED by a proved lemma**: THM-883's
Fragmentation Lemma bounds every outlier explicitly (w_min <= 2j/((13-2j) ell_max);
last outlier <= 2/(13 ell_max); recursion non-vacuous by LRC(<=13)), making the
configuration space provably finite for j <= 6. Exact sweep: j <= 4 complete (zero
violations), j = 5,6 running on the proved-finite box; j >= 7 (never in this file's
census) is loose-tile territory, probed 0/7107. The "certified" qualifier below is
superseded for j <= 4 and reduces to mechanical sweep completion for j = 5,6.

## Honest status (original)

- **Unconditional given the two inputs:** far-element monotonicity (Step 1, verified +
  THM-717/720) and the finite check (Step 2, 64317 configs). Both are the standard certified
  ingredients of the covering-min.
- **No closed-form inequality** for multi-killer (the balance undershoots; the optimum is
  global) — this is the genuine open item, identical to kps-S127cont58's honest limit.
- **Combined with THM-724** ⟹ deep well is the unique global covering-min. The covering-min
  *rigidity* (unique minimizer) is thereby complete modulo (a) THM-724's near-tight large-`s`
  residual and (b) this file's closed-form gap — both empirically closed, both reducing to the
  same open covering-min conjecture (THM-523).

*Artifacts:* `04-computation/lrc14_multikiller_{landscape,monotone}_macmini_S70.py` (+`.out`).
Credits: kps-S127cont58/HYP-6225 (enumeration, core-length monotonicity target, lcm-carriers),
THM-724 (single-killer half), THM-717/720 (far-element monotonicity), THM-366 (sieve),
opus-S253 (balance, multi-constraint framing).

---

## ADDENDUM (death-star-2026-07-16-S17, HYP-7010): shape-completion at |P| ∈ {10, 11} — monotonicity ELIMINATED as an input on these strata

**CONVERGENCE NOTE (same-night double, reconciled):** mac-mini-S114's THM-883 (Fragmentation
Lemma) rigorizes this file by the same mechanism, concurrently and independently — their
`w_min ≤ 2j/((13−2j)·ℓ_max)` bound is the SHARPER form of this addendum's union-tail
threshold `W > 2jr/(m(13−2j))` (largest-component version vs total-measure version; theirs
gives smaller boxes via the last-killer bound `w_j ≤ 2/(13ℓ_max)`), and their exact sweep
completes j ≤ 4 (so this addendum's k=9 cron unit is SUPERSEDED). This addendum's
independent value: convergent confirmation at |P| ∈ {10,11} by a different enumeration
(recursive per-node exact (m,r) thresholds), the explicit free-outlier shape-gap framing,
the {1..14}\{a} anchor atlas, and the 1/13 measure atlas below. Two independent
implementations, zero violations in both — the certification is now double-witnessed.

**The gap this closes.** Step 2's finite check enumerated **interval cores {1..k} only**
(k = 9,10,11, carrier-multiple outliers ≤ 220). Non-interval small parts `P ⊊ {1..12}` admit
**free outliers** (outliers carrying no missing modulus — impossible for interval cores), e.g.
`{1..12}\{6} ∪ {w, 182m}` is primitive covering multi-killer for EVERY `w ≥ 13`. These shapes
were never in the certification basis; and Step 1 (far-element monotonicity) was
verified-not-proved.

**The replacement (no monotonicity anywhere).** UNION-TAIL LEMMA (the 1/13-threshold version
of opus-S32's simultaneous-peel; elementary, proof in the script header): if a core has exact
closed 1/13-good set `G` with measure `m > 0` and `r` components, then any `j ≤ 6` additional
speeds all `≥ W` with `W > 2jr/(m(13−2j))` leave positive good measure, hence `M ≥ 1/13`
unconditionally. (Each speed `w` removes at most `(2/13)|I| + 2/(13w)` from each component `I`.)
The audit recursion applies this at every node with the CURRENT exact `(m, r)` — outliers
below the threshold are enumerated exhaustively (witness pre-filter, exact-interval decision
fallback, exact-M referee on any failure); outliers above are cleared by the lemma. No
monotonicity in outlier values is assumed at any point.

**Results (`lrc14_multikiller_shape_audit_deathstar_S17.py` + `.out`):**
- **|P| = 11 (j = 2): SHAPE-COMPLETE.** All 12 shapes `{1..12}\{a}`; per-shape exact
  `m ∈ [7/429, 1/13]`, thresholds `W_all ≤ 382`; 383 below-threshold covering pairs
  enumerated, ALL clear (every one by explicit rational witness; exact fallback needed 0
  times); everything above cleared by the lemma. **Zero failures.**
- **|P| = 10 (j = 3): SHAPE-COMPLETE.** All 66 shapes; ~10⁶ recursion nodes with per-node
  exact (m, r); ~17k leaf decisions, all cleared by explicit witnesses. **Zero failures.**
- **Flagship anchors exact:** `{1..14}\{a}` (a = 1..7, all covering, outliers {13,14}) have
  M = 1/8, 2/17, 2/19, 2/19, 2/21, 2/23, 1/11 — all ≥ 1/13 ✓. Free-outlier family
  `{1..12}\{6} ∪ {w,182}` clears for ALL w (enumerated + tails). Self-tests: deep well
  14/183, `{1..11,13,84}` = 7/89, `{1..13}\{6}∪{182}` = 2/23 all reproduced exactly.
- **New exact data:** the 1/13 measure atlas m(good({1..12}\{a})) = 1/13, 1/26, 20/429,
  9/286, 254/6435, 7/429, 383/4290, 472/9009, 311/5005, 94/3003, 461/8190, 92/5005
  (a = 1..12) — hardest shapes a ∈ {6,12}; a = 1 exactly 1/13.
- **|P| = 9 (j = 4): NOT completed this session** — per-shape cost ≥ 10-20 min × 220 shapes
  (≈ 50-150 core-hours) = a cron bank, not a session job. The script supports per-shape
  invocation (`--k9-shape a,b,c`; further parallelizable by top-level w₁ ranges). Until the
  bank completes, |P| = 9 retains the legacy interval-core certification + adversarial sweeps.

**Consequences for the covering-min chain.** With THM-724 Case 1 (covering + exactly one
outlier ≥ 13 forces `S = {1..12, 182m}`, the deep-well ladder, unconditional):

> Every primitive covering 13-set with `≥ 10` elements in `{1..12}` has `M ≥ 14/183`, with
> equality iff the deep well — **certification now shape-complete and monotonicity-free on
> these strata.**

The S111 low-M rigidity assembly's "outlier-threshold bookkeeping" resolves as: covering +
`≤ 1` element `≥ 13` ⟹ `P = {1..12}` and `182 | f` (13 and 14 must both divide the unique
outlier) ⟹ the ladder — the `> 14` wording in the assembly is superfluous; the `≥ 13` count
is the correct dichotomy variable. Note `{1..14}\{a}` sets are multi-killer in this taxonomy
(outliers {13,14}) and all clear 1/13 (anchors above), so no near-AP tile is needed for the
covering-min on |P| ≥ 10 — only for loneliness statements elsewhere.

*Addendum artifacts:* `04-computation/lrc14_multikiller_shape_audit_deathstar_S17.py`,
`05-knowledge/results/lrc14_multikiller_shape_audit_deathstar_S17.out`, HYP-7010.
Additional credit: opus-S32 (simultaneous-peel template at 1/14).

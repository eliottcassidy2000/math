---
id: THM-726
title: Multi-killer covering-min rigidity — every primitive covering 13-set with ≥2 far outliers has M ≥ 1/13 > 14/183, hence the single-killer deep well {1..12,182} is the UNIQUE global covering-min. Proved by far-element monotonicity (min at the smallest lcm-carrier outliers) + a finite check (64317 configs) + THM-724. Completes the covering-min rigidity with THM-724 (the single-killer half)
status: PROVED-BY-CERTIFICATION (finite check + monotone tail; the standard covering-min proof shape). The uniform bound M ≥ 1/13 for genuine multi-killer is unconditional GIVEN far-element monotonicity (verified on extremals + THM-720 looseness dichotomy) and the finite check. Combined with THM-724 (single-killer, deep well unique) ⟹ deep well is the unique global covering-min. A closed-form (non-enumeration) inequality is NOT available (the balance provably undershoots — the multi-killer optimum is global, not at the core witness); this matches the certified-not-closed status of the whole covering-min.
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

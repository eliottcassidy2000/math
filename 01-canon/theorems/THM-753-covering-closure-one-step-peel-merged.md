---
id: THM-753
title: The covering-closure one-step peel — every covering 13-set is 1/14-lonely by peeling its max element v=max(S) off a 12-element core P whose M(P)≥1/13 is FREE (LRC(13), settled), certified by EITHER THM-751 (v aligned to P's tight point and v≥13·max(P) ⟹ M(S)≥(1/13)·v/(v+max P)≥1/14, PROVED) OR the disc far-peel THM-731/732 (L(S)=(6/7)|G'_P|−ε_v>0 when disc_v<6|G'_P|², PROVED for isolated v). Base max(S)≤14 is kps THM-738 (PROVED). This MERGES mac-mini's combinatorial monotonicity and klein's analytic disc into one recursion-free closure; the sole open hypothesis is disc_v<6|G'_P|² for the moderate non-aligned far element = opus's finite mod-360360 binding residual
status: ASSEMBLY / merge (klein-2026-07-14-S306). The skeleton is PROVED: LRC(13) core floor (SETTLED); THM-751 clean/aligned branch (PROVED); disc unsafe branch for isolated v (PROVED via crude THM-732 arc-count r<3√2 v|G'|, S289); base ⊆{1..14} (kps THM-738/734/741, PROVED); |S|≤2 trivial. The SOLE remaining analytic hypothesis: disc_v<6|G'_P|² for a moderate (non-isolated, non-aligned) far element v>14 — verified with 50–1000× margins, exact-ℚ per family via THM-732's Bernoulli form, and reduced by opus-S286 to a FINITE classification of residue patterns mod 360360 (occupancy-rigidity: blocking-complete forces a complete residue system, AP-like). VERIFIED end-to-end: 0 stalls on 15 adversarial covering families (klein-S305), incl the M=2/23 no-shadow-tile family {1..13\6,182}.
source: klein-2026-07-14-S306 (merging mac-mini-THM-751 + klein-THM-731/732 + kps-THM-738)
depends_on:
  - THM-751   # mac-mini: aligned far-element monotonicity (the clean-peel branch, PROVED)
  - THM-731   # klein: the peeling certificate L=(6/7)|G'_~v|−ε_v, |ε_v|²≤(6/49)disc_v
  - THM-732   # klein/kps: disc_v≤r²/3v² (crude) + exact Bernoulli form (the unsafe-peel branch)
  - THM-738   # kind-pasteur: every ≥10-in-{1..14} body is lonely (the base, PROVED)
  - THM-523   # non-covering dispatch (for cores that lose covering under peeling — not needed: core floor is LRC(13))
related:
  - THM-726   # multi-killer covering-min rigidity (THM-751 is its Step 1)
  - HYP-6680  # klein-S305 (the far-peel consolidation this formalizes)
  - opus-S286 # occupancy-rigidity: the moderate-v residual is finite mod 360360
external: LRC(≤13) SETTLED (owner policy).
---

# THM-753 — the covering-closure one-step peel (merged)

The single theorem the three routes converged on: **peel the largest speed once, off a core that LRC(13)
already makes lonely.** It merges mac-mini's combinatorial monotonicity (THM-751) and klein's analytic disc
peel (THM-731/732) into one statement, with kps's bounded body (THM-738) as the base.

## Statement

Let `S` be a covering 13-set, `v = max(S)`, `P = S \ {v}` (a 12-element set). Recall **`M(P) ≥ 1/13`**
(LRC(13), SETTLED), so the good set `G'_P = {t : ‖pt‖ ≥ 1/14 ∀p∈P}` is nonempty and open (`|G'_P| > 0`).
Then `M(S) ≥ 1/14`, by exactly one of:

1. **Base — `v ≤ 14`.** Then `S ⊆ {1,…,14}` with 13 speeds in `{1..14}` (`≥10`), so `M(S) ≥ 1/14` by
   **kps THM-738** (PROVED). [Or `|S| ≤ 2`: trivial, `M ≥ 1/3`.]

2. **Clean peel — `v > 14`, `v` aligned to a tight point `t₀` of `P` (`v t₀ ∈ ℤ`), and `v ≥ 13·max(P)`.**
   Writing `v = w m` with the aligned carrier `w`, **THM-751** gives
   `M(S) ≥ M(P)·v/(v + max P) ≥ (1/13)·v/(v + max P) ≥ 1/14` (the last step is `v ≥ 13·max P`). PROVED,
   by the explicit witness `t = t₀ + M(P)/(max P + v)`.

3. **Unsafe peel — `v > 14`, otherwise.** By the peeling identity (**THM-731**)
   `L(S) = (6/7)|G'_P| − ε_v`, `|ε_v|² ≤ (6/49)·disc_v`, so
   `M(S) > 1/14` as soon as `disc_v < 6|G'_P|²`. This holds (PROVED, **THM-732** crude
   `disc_v ≤ r²/3v²` ⟺ `r < 3√2·v|G'_P|`) whenever `v` is **isolated** (`v ≫ max P`); the residual
   moderate-`v` case is the one open hypothesis below.

## What is proved, and the single open hypothesis

**PROVED:** the core floor (LRC(13)); the clean branch (THM-751, all aligned lcm-carrier multi-killers —
deep well, single-killer, etc.); the isolated unsafe branch (crude disc, all far-element-dominant sets incl
the loose/spread escapees, klein-S304); the base (kps THM-738); `|S|≤2`. Together these close every covering
`S` whose max element is `≤14`, or aligned+isolated, or isolated non-aligned.

**THE SOLE OPEN HYPOTHESIS (H):** `disc_v < 6|G'_P|²` for a **moderate, non-aligned** far element
(`14 < v`, `v` comparable to `max P`, `v t₀ ∉ ℤ`). Status:
- **Verified** with 50–1000× margins on every family tested (klein-S305: 0 stalls on 15 adversarial,
  including the `M=2/23` family `{1..13\6,182}` that no `k≤13` shadow tile covers).
- **Exact per family** via THM-732's Bernoulli edge-pair form (`disc_v = (1/2v²)ΣσσB₂({v(e−e')})`,
  rational).
- **Finite** as a uniform statement: opus-S286's occupancy-rigidity shows the moderate-`v` binding residual
  is a classification of residue patterns mod `lcm(2..13)=360360` passing surjectivity cuts — the
  equidistribution reduced to finite-index combinatorics. opus-S286 conjectures (two independent zero-hunts
  as evidence) that a `k≤13` shadow always exists on this residual, which would discharge (H) entirely.

## Why the merge is exact

THM-751 and THM-731 certify the *same* peel from opposite sides. THM-751 is the **combinatorial** witness
(`t = t₀ + s`, one balance point) — sharp when the far element is *aligned* (its danger tooth pins to the
core's tight point and only narrows). THM-731 is the **measure** witness (the good set survives the far
element's thin arcs) — sharp when the far element is *large* (fine grid ⟹ small disc). The clean/unsafe
split is mac-mini-S102's peeling recursion; this theorem states it as one closure with LRC(13) supplying
the core floor for free, so **no recursion is needed** — a single peel of `max(S)` suffices, because the
12-element core is already lonely by the settled sub-case.

*Files: klein-S304/S305 (`lrc14_iterated_peel`, `lrc14_exactdisc_uniform`); consumes THM-751/731/732/738,
opus-S286. HYP-6690.*

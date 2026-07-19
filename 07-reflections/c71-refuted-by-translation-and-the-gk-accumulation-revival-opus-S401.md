# Two S400-residue angles: Conjecture 7.1 refuted by translation, and the Giri–Kravitz accumulation revival

**Instance:** opus-2026-07-19-S401 (owner: "keep working creative angles like HYP-7920").
**Claims:** THM-1288 (the refutation lemma — self-contained, proved), HYP-7930 (the G-K
revival — derived, citation-pending). Scripts + frozen outs:
`lrc_c71_refutation_divisor_alignment_opus_S401.py`.

Both angles are HYP-7920-shaped: external statements mined with repo lore. The first is a
complete, self-contained refutation; the second revives the *legitimate* half of a citation
the fleet dropped entirely after burning itself on the illegitimate half.

---

## Angle B — THM-1288: S-T Conjecture 7.1 is FALSE for every k ≥ 2

**Their conjecture (arXiv:2604.23906 §7, verbatim):** "Let k+1 be a positive integer. There
exists a constant D such that, for any integer d ≥ D, every non-tight speed tuple
v ∈ ℤ_{>0}^k with gcd(v) = 1 has a witness time in (1/d)ℤ." (D uniform over v; their
motivating observation: non-tight tuples have an interval of witness times.)

**Refutation (three steps, fully self-contained — no dependence on their paper beyond the
conjecture's statement):** fix any d with (k+1) ∤ d and let V_d = {d+1, d+2, …, d+k}.

1. *V_d is coprime and non-tight.* Consecutive integers are coprime; at t* = 1/(2d+k+1) the
   positions (d+i)t* fill a band symmetric about 1/2 with min distance
   (d+1)/(2d+k+1) = 1/(k+1) + d(k−1)/((k+1)(2d+k+1)) — strictly above threshold (margin
   → (k−1)/(2(k+1)) ≈ 1/2), so strict witnesses exist and V_d is non-tight.
2. *At the d-grid, V_d is indistinguishable from the AP.* For t = j/d:
   ‖(d+i)·j/d‖ = ‖i·j/d‖ — the translation by d vanishes mod 1 on the grid.
3. *No j/d is an AP-witness when (k+1) ∤ d.* If min_{1≤i≤k}‖i·j/d‖ ≥ 1/(k+1), then among
   the k+1 points {0, t, 2t, …, kt} two lie within 1/(k+1) (the k+1 circular gaps sum
   to 1), so ‖m·j/d‖ ≤ 1/(k+1) for some 1 ≤ m ≤ k; with the hypothesis this is equality,
   1/(k+1) = ‖m j/d‖ ∈ (1/d′)ℤ with d′ | d, forcing (k+1) | d — contradiction.

So V_d has no witness in (1/d)ℤ; since d can be any (k+1)-nonmultiple ≥ any D, no uniform D
exists. ∎ (Exact verification k = 12, 13 at d up to 100,003, with the (k+1)|d controls
correctly producing witnesses — frozen output.)

**Why the repo saw this instantly:** the mechanism is the repo's translation lore, fifth
appearance. The (1/d)ℤ-restricted problem is blind to translations by multiples of d; the
cluster is the AP translated by d — M jumps from 1/(k+1) to ≈ 1/2, but the grid can't see
it; and the AP's witness set is the measure-zero {s/(k+1)} (their own Remark 3.2), which
the grid misses entirely unless (k+1) | d. Loose-but-grid-blind families are exactly what
the repo's detection-floor discipline trains one to construct.

**What survives of their observation (the repaired conjecture):** per-family, D(v) < ∞
always (witness intervals have positive length). The two blowup mechanisms for D(v) are now
cataloged: (i) translation-alignment (this refutation — D(v) ≥ largest (k+1)-free divisor
structure of the family ~ v_min); (ii) witness-width collapse near the floor (the S400
probe — D(v) ~ v_max/(M − 1/(k+1))). A version of C7.1 with D allowed to depend on v_max
(or on min-speed) is the honest candidate; and by Angle C below, mechanism (ii) cannot
occur along value-sequences (no accumulation at the floor), which sharpens what a repaired
conjecture can hope for. This refines — partially supersedes — the S400 "entanglement"
framing: the uniform-D version dies for a cheap translation reason before the gap program
even enters; the gap program governs the v_max-normalized version.

## Angle C — HYP-7930: the Giri–Kravitz accumulation theorem, revived on its legitimate side

**History:** MISTAKE-117 (opus-S130, 07-06) correctly retracted Route 2's use of
Giri–Kravitz to bound the *sup* — their paper governs *accumulation points*, not extrema.
The fleet then dropped the paper entirely (grep: zero uses since). That was an
overcorrection: the accumulation content is real, citable, and bears directly on the live
(1/14, 3/41) question — whose second horn (opus-S396/THM-1268: "…or 1/14 is an
accumulation point from above and there is no gap at all") is exactly accumulation-type.

**Their statement (arXiv:2304.01462, v4; abstract + Theorem 1.4 chain as extracted):**
S(n) := {D(T) : T a 1-dim subtorus of (ℝ/ℤ)ⁿ not in a coordinate hyperplane}, D(T) = inf
L∞-distance (torus metric) from T to (1/2,…,1/2); main result "the set of accumulation
points of S(n) is precisely S(n−1)", refined in-paper (Thm 1.4, k=1 chain) to
S(n−1) ⊆ acc(S(n)) ⊆ S₂(n) ⊆ S₂*(n) = S*(n−1), with "only upper accumulation points".

**Translation (exact):** a nondegenerate 1-dim subtorus is {t·v} for a primitive integer
v with all coordinates nonzero (n = number of SPEEDS), and circle distance to 1/2 is
1/2 − dist-to-0, so D(T) = 1/2 − M(v): S(n) is the M-spectrum in the variable 1/2 − M.
Accumulation structure is preserved under x ↦ 1/2 − x. The upper inclusion (the only
direction used below) gives: **acc(M-spectrum of k speeds) ⊆ {M-values of starred (k−1)-
directions}** — where the star at worst allows zero coordinates (fewer effective runners)
and repeats (M unchanged), so every such value is ≥ 1/13 for k−1 = 12 by the settled
LRC(≤13). Hence **every accumulation point of the 13-speed M-spectrum is ≥ 1/13**, and
likewise ≥ 1/12 at 12 speeds.

**Corollaries (derived; citation-pending on the exact Thm 1.4 wording + the S* definition):**
- **(C1) The n=14 floor is isolated from above:** ∃δ > 0 with no primitive 13-speed family
  having M ∈ (1/14, 1/14 + δ). The accumulation horn of THM-1268's dichotomy is DEAD.
- **(C2) (1/14, 3/41] contains only FINITELY many attained values.** (3/41 = 1/13 − 2/533,
  so [1/14, 1/13 − ε] applies with ε < 2/533; Bolzano–Weierstrass on a compact interval
  avoiding the accumulation set.) In particular the slack-1 ladder D/(14D−1) has only
  finitely many realized rungs — retro-explaining why the repo's structured searches keep
  finding nothing between 3/41 and 1/14, and consistent with death-star's cross-N gate
  tower (rung realizability as a sparse arithmetic condition).
- **(C3) At 12 speeds: (1/13, 2/25) ∩ spectrum is a finite list.** HYP-7310's gap
  statement weakens, by citation, from "an interval is empty" to "**a finite list is
  empty**" — a type change for Wall A's shadow question. (The rigidity core at the attained
  floor is untouched.)
- **(C4) Composition with HYP-7920:** the cage empties the micro-gap (1/13, 113/1466)
  effectively but only below height 258,276; G-K empties all-but-finitely-many values at
  all heights but ineffectively. The sandwich isolates a NEW well-posed target:
  **effectivize the bottom of the G-K argument** (their degeneration mechanism — box
  geometry/dimension counting — applied only near the floor) to bound the heights of the
  finite list. Even a crude effective bound composes with the cage and the repo's
  enumeration machinery.

**Honest scope:** (i) the G-K constants are ineffective as extracted — no δ, no list bound;
(ii) the extraction of Theorem 1.4's chain and the S* definition needs one PDF-pinning
session (same discipline as HYP-7920's F1–F3; the failure mode of MISTAKE-117 was sloppy
citation, and the fix is precise citation, not no citation); (iii) versions: v2 (2023) vs
v4 (2025-03) — pin against v4; (iv) the lower inclusion S(n−1) ⊆ acc(S(n)) additionally
predicts infinitely many 13-speed values just above 1/13 — no conflict with anything in
canon, and its *constructive* content (how they build approaching families) is a machine
the rung-realizability program should read.

## What changed for the repo, in one paragraph

The (1/14, 3/41) question was a dichotomy with an infinite horn; the horn is now dead
modulo one citation-pinning session, and the question becomes "determine a finite list" —
with an effective-height companion target named (C4). Their Conjecture 7.1 is settled in
the negative with a three-line translation argument, which is simultaneously a
literature-facing note (the owner may wish to communicate it to the authors alongside the
HYP-7920 cage — both are contributions TO their program) and a fresh confirmation that the
repo's invariance triage (S399 §3, axis 1) has predictive power outside the repo. Wall A's
rigidity core is untouched by both angles; what both do is shrink and retype everything
around it.

## Cross-links

HYP-7920 + S400 reflection (the cage; the entanglement framing refined here) · THM-1288
(this refutation) · HYP-7930 (this revival) · MISTAKE-117 (the sup-misuse; the boundary of
legitimacy) · THM-1268/1269 (the dichotomy whose horn dies) · THM-1230/1235 (the 3/41
witness/ladder) · HYP-7310 (Wall A; gap → finite list) · death-star S59 gate tower
(rung-realizability; the finite list's natural describer) · S399 synthesis §3 axis 1
(translation blindness, fifth strike) · scripts + frozen outs (S401).

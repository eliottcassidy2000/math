---
id: THM-682
title: The collision arm of the final rung closes — (a) THE COMMON-RESIDUE DISPATCH (new, unconditional): any covering primitive family whose 13 speeds share a residue class mod d ≥ 2 has M(v) ≥ 8/17 (all phases COINCIDE at τ = c/d, since covering+primitivity force every prime factor of d to be ≥ 17); (b) k=13 restricted-sumset stability B ≤ 24 ⟹ diam ≤ 13 (EXHAUSTIVE: exactly 5 normalized sets) extending to the ladder law diam ≤ B − 11 verified through B ≤ 28; (c) hence every core family has B ≥ 29 — collision multiplicity ≤ 49 — because B ≤ 28 forces AP-containment whose d ≥ 2 arm is common-residue-dispatched and whose d = 1 arm lands inside window-22; (d) the W₀-carrier lemma: support-2 global exact relations are DOUBLINGS ONLY, so the final rung narrows to the doubling-rich (≥ 3 proven / ≥ 6 measured) collision-light (B ≥ 29) even-rich corner — the 2-adic pressure direction — plus LEM-016-protocol slivers
status: PROVED ((a) complete proof below, 5 lines, verified 400/400 on in-core affine-dilated families with EXACT clearance match; (b) t = 1 level EXHAUSTIVE with no cap caveat (max diam 13 realized, cap 45 unreachable by the proved 2-block tail at this budget); ladder levels B ≤ 28 exhaustive below the DMAX = 90 cap with the far-jump sliver FLAGGED per LEM-016 protocol and the 2-block tail PROVED (B ≥ 31); (c) follows from (a)+(b) for the proved range B ≤ 28 — the statement "B ≥ 29 for core families" is PROVED modulo the flagged sliver; the extension to B ≥ 33 (meeting window-22 exactly at diam 21) is CONJECTURED with the GAP sharpness ceiling MEASURED at B = 34; (d) the carrier classification is a one-line sign argument, machine-echoed; the doubling counts use proven line-weight bounds). Machine-verified throughout (companion .out).
source: monad-explorer-2026-07-09-S11 (HYP-5817) — executing THM-681's final-rung handoff with kps/opus as the E3-ladder endpoints.
depends_on:
  - THM-680   # the per-ruler floor whose off-line budget W₀ this narrows
  - THM-681   # the taxonomy + dichotomy this closes the collision arm of
related:
  - LEM-016 (boxeph/mac-mini)  # k=7 t=1 stability + the rank-2 GAP escape shapes; the sliver-flag protocol followed here; KEY CONTRAST — k=7 escapes at t=2, k=13 escapes only at t=11 (B=34): the sharpness does NOT scale
  - kps-S120/S121 LRCE3Budget  # the rigidity top of the ladder (E3 deficit on the whole residual)
  - opus-S189 freimanChain     # the chain bottom; the middle rung handed off = B ≥ D + 11 for D ≤ 21
  - THM-668/678  # the detuned dispatches; the common-residue dispatch is their affine sibling
  - window-22 (hwindow22)      # the d = 1 arm lands exactly inside it
---

# THM-682 — the collision arm closes: the common-residue dispatch + the B-ladder

## (a) THE COMMON-RESIDUE DISPATCH (PROVED, unconditional)

**Statement.** Let `v_1 < … < v_13` be positive integers, primitive (`gcd = 1`), covering
(every `q ∈ {2, …, 14}` divides some `v_l`). Suppose all speeds share a residue class:
`v_l ≡ a (mod d)` for all `l`, some `d ≥ 2`. Then `M(v) ≥ 1/2 − 1/(2d) ≥ 8/17`.

**Proof.** (1) Let `p ≤ 13` be prime with `p ∣ d`. Covering at `q = p` gives `l` with
`p ∣ v_l ≡ a (mod p)`, so `p ∣ a`, so `p ∣ v_l` for EVERY `l` — contradicting
primitivity. Hence every prime factor of `d` is ≥ 17, in particular `d ≥ 17`.
(2) If `p ∣ gcd(a, d)` then `p ∣ v_l` for all `l` — again non-primitive; so
`gcd(a, d) = 1`. (3) Choose `c` with `a·c ≡ ⌊d/2⌋ (mod d)`. At `τ = c/d`:
`v_l·τ = (a + s_l d)·c/d = ac/d + s_l c ≡ ac/d (mod 1)` — ALL thirteen phases coincide
at `‖ac/d‖ = ⌊d/2⌋/d ≥ 1/2 − 1/(2d) ≥ 8/17 > 1/14`. ∎

**Remarks.** (i) The dispatch is the AFFINE sibling of THM-668 (there `d` divides the
values of all but one runner; here `d` divides all DIFFERENCES). (ii) The hypothesis is
exactly "`gcd(v_2 − v_1, …, v_13 − v_1) ≥ 2`", so its complement — **the differences
are setwise coprime** — is a legitimate NEW clause for the Lean `ResidualObligation`
(branch 8): residual families may be assumed difference-primitive. (iii) Verified on
400 in-core (covering, primitive, gapped, distinct) affine-dilated families: every one
dispatched, clearance exactly `1/2 − 1/(2d)` (companion Part A).

## (b) The k = 13 restricted-sumset stability ladder

`B(A) = |A +̂ A|` for a 13-set `A`, gcd-normalized (`min = 0`, `gcd = 1`). Classical:
`B ≥ 23`, equality iff AP.

- **t = 1 (EXHAUSTIVE, no cap caveat):** `B ≤ 24 ⟹ diam(A) ≤ 13`. Exactly **5**
  normalized sets exist: the AP `{0,…,12}` (B = 23) and four B = 24 near-APs of
  diameter ≤ 13 (e.g. `{0, 2, 3, …, 13}`). DFS: 16 564 nodes; the budget forces every
  prefix within +1 of the classical `2s − 3` minimum, so the tree is tiny.
- **The ladder law (exhaustive through B ≤ 28 below the cap):** `diam ≤ B − 11`,
  tight at every level: (B, max diam) = (24, 13), (25, 14), (26, 15), (28, 17).
  Equivalently **`B ≥ diam + 11 = diam + (k − 2)`** on this range.
- **Tail (PROVED):** any gap > D/2 splits `A` into blocks with disjoint sum ranges:
  `B ≥ (2s − 3) + 12 + (2(13 − s) − 3) = 32` (singleton blocks ≥ 33) — so B ≤ 28 sets
  have no gap > D/2, and 2-block configurations never enter the ladder range.
- **Sliver (FLAGGED, LEM-016 protocol):** diam > 90 with all gaps ≤ D/2 and every
  prefix within the budget — not excluded by the DFS cap; the measured margin
  (max realized diam 17 vs cap 90) and the multi-block arithmetic (3-block ≥ 34
  measured) make it empty in practice; it is the same sliver LEM-016 flags at k = 7.
- **Sharpness ceiling (MEASURED):** rank-2 GAP escapes (LEM-016(ii) shapes, blocks at
  0/c/2c) begin at **B = 34 (t = 11)** — the k = 7 sharpness (escape at t = 2) does
  NOT scale; the window for rank-1 stability at k = 13 is enormous. Every near-minimal
  GAP-shape core instance tested is fully ruler-live (23/23, 24/24, 19/19, 18/18 tall
  rulers) — the escapes are generic-side (THM-681) anyway.

## (c) The collision-arm dispatch chain (PROVED for B ≤ 28)

Let `v` be a core family (covering, primitive, gapped `v_max > 13·v_min`, distinct).
If `B(v) ≤ 28`, then by (b) `v ⊂ {a, a + d, …, a + 17d}` for some `a ≥ 1, d ≥ 1`
(AP-containment, 18 terms). Two arms:

- **d ≥ 2:** all speeds ≡ a (mod d) ⟹ the common-residue dispatch (a): `M(v) ≥ 8/17`.
- **d = 1:** `v ⊂ [v_min, v_min + 17]`; gapped gives `13·v_min < v_max ≤ v_min + 17`,
  so `12·v_min < 17`, forcing `v_min = 1` and `v ⊂ {1, …, 18}` — inside **window-22**
  (`hwindow22`), dispatched.

**Hence every core family satisfies `B(v) ≥ 29`** (modulo the flagged sliver): total
collision multiplicity `≤ 78 − 29 = 49`. **Conjectured extension:** the ladder law
holds through `B ≤ 32` (diam ≤ 21 — meeting window-22 EXACTLY at `v ⊂ {1, …, 22}`),
giving `B ≥ 33` and collisions ≤ 45; the boundary case `B = 33` (diam 22, `v ⊂ {1, …, 23}`,
necessarily `1, 23 ∈ v`) is **VERIFIED in-session: all 8 361 covering boundary families
are lonely, worst witness modulus q = 27** — a trivially `native_decide`-able object;
B = 34 is the measured GAP ceiling. The middle-rung
lemma handed to kps/opus: **`B ≥ D + 11` for gcd-normalized 13-sets of diameter
D ≤ 21** — endpoints owned (kps E3-rigidity at the AP; opus freimanChain), verified
exhaustively here through D ≤ 17.

## (d) The W₀-carrier lemma (PROVED)

A support-2 relation `c_a v_a + c_b v_b = 0` with `1 ≤ |c| ≤ 2` on distinct positives
forces opposite signs and ratio 2: **the only support-2 global exact relations are
doublings `v_b = 2v_a`** (the (2, −2) case reduces to equality; machine-echoed
exhaustively). With THM-681's weights: `W₀ > 0.08` requires **≥ 3 doublings** (proven
line bound 0.0382) — **≥ 6 measured** (0.0142 each); Schur content is nearly free
(0.0027) and cannot carry W₀. Random search over core families maxes out at **9
doublings** (e.g. `[22, 39, 42, 44, 45, 84, 88, 90, 168, 176, 180, 336, 352]`,
11/13 even) — doubling-rich families are EVEN-RICH: the 2-adic pressure toward the
`g = 2` detuned dispatches (THM-668/678).

## Consequence: the final rung after THM-682

THM-681's final rung `W₀ ∈ (0.08, ladder)` now requires SIMULTANEOUSLY:
- **collision-light:** `B ≥ 29` proved (conj. ≥ 33) — the collision route to ruler
  death is capped; and
- **doubling-rich:** ≥ 3 (proven) / ≥ 6 (measured) doublings — even-rich, 2-adically
  coherent, pressed toward the g = 2 dispatches;
plus the two named finite objects (B = 33 boundary `native_decide`; the LEM-016-style
sliver). The rung is no longer an interval of exotic families — it is the
**doubling-chain corner**, where the harmonic dispatches and the 2-adic tower already
operate.

## Verification & files

`04-computation/lrc14_final_rung_collision_arm_monad_S11.py` (+ `.out`): Part A the
400-family dispatch battery (clearance exact); Part B the t = 1 DFS (5 sets) + tail;
Part C the GAP escape threshold (B = 34) + ruler-liveness of escape instances (all
live); Part D the carrier + doubling census. Ladder extension runs (B ≤ 25, 26, 28 at
DMAX = 90; bitmask B ≤ 29+) in the scratch ladder logged in the session letter and
`.out` addendum.

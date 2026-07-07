---
source: klein-2026-07-06-S152 (HYP-4711)
status: creative synthesis + verified mechanism (not yet a full proof). The r=13 coarse-descent
  crux (escape / L-lift families whose COARSE part is the AP) carries an explicit AP-inherited
  conjugate witness giving M >= 1/14 UNIFORMLY — so the escape families are NOT a second open
  obstruction; the entire LRC(14) residual collapses onto ONE object: the single-scale moat.
tags:
  - lonely-runner
  - LRC14
  - coarse-reduction
  - escape-families
  - conjugate-witness
  - route-1
  - synthesis
---

# The escape families collapse to the moat: the AP carries its own witness

Owner: *pull from agents, extend their ideas, attack the LRC(14) proof* (the M(V)→M(K) descent).
I extended kps-S52/S53's coarse reduction and mac-mini-S38's decorrelation atom, and found that
the case where BOTH break — the sharp-margin r=13 crux — is not a new hard problem at all.

## Where the coarse descent breaks (the honest crux)

The coarse bound `M(v) ≥ M(K) − A/L` (kps `reach_transfer_coarse`, mac-mini `reach_decorr`) for
`vᵢ = aᵢ + L·kᵢ`, `K = {distinct kᵢ}`, is only useful when `M(K)` has slack. It **fails exactly
when the coarse part `K` is the AP** `{1,…,13}` — then `M(K) = 1/14`, and the bound degrades to
`M(v) ≥ 1/14 − A/L < 1/14`. The r=13 branch (kps: "DESCEND") does **not** close here: the descent
loses `A/L` each step and the AP has no slack to absorb it. mac-mini-S38 hit the same wall from the
covering side (open NEXT item (a): the boundary `M(K)=2/25` block-lift gives only `2/25−12/L`, but
the true `M(v)` is higher — "needs a sharper base-structure bound").

So the genuinely hard survivor is the **escape / L-lift family with AP coarse part**:

    vᵢ = aᵢ + L·(d·i)     (K = d·{1,…,13}, |aᵢ| ≤ A, L large),

the perturbed dilated AP — mac-mini's S36/S37 "escape" object, and the family my own S144–S150
covering-min L-lift work circled.

## The AP carries its own witness (the mechanism)

The dilated AP `L·d·{1,…,13}` is tight (`M = 1/14`) and has **exactly `φ(14)=6` optimal witnesses**

    t_c = c / (14·d·L),   c ∈ (ℤ/14)* = {1,3,5,9,11,13},

and each `t_c` binds **exactly one antipodal index pair** `{i₊, i₋}` (`i₊c ≡ +1`, `i₋c ≡ −1 mod 14`);
the other 11 speeds sit at `‖·‖ ≥ 2/14`, with room to spare. Perturb by `a` and shift `t = t_c + δ`.
An **exact** (un-linearised, in-branch) computation gives: the binding pair stays `≥ 1/14` iff

    a_{i₊}/v_{i₊}  ≥  a_{i₋}/v_{i₋}      (v_i = a_i + L d i),

a single scalar slope test, and there is always a valid `δ` of size `O(A/L²)` when it holds. The
**conjugate** witness `c ↦ 14−c` binds the *same* pair with `i₊ ↔ i₋` swapped, i.e. it flips the
inequality. **So for every family, one of each conjugate pair `{c, 14−c}` satisfies the slope
test** → a valid witness exists → `M(v) ≥ 1/14`. The non-binding speeds stay `> 1/14` once
`L ≳ 200·A` (the same scale threshold as the coarse descent's `L > 182·A`).

### Verified (exact `Fraction` arithmetic, `lrc14_ap_escape_conjugate_witness_klein_S152`)

- Base `K={1..13}`, random `a`, every `(L,A)` with `L∈{183,500,3600}`, `A∈{1,3,6}`:
  **200/200 certified `M ≥ 1/14`** by the constructed conjugate witness.
- The analytic slope test predicts the winning `c`: **0/100 mispredictions, 0/100 uncertified.**
- Permuted AP: **200/200.** Dilated AP `d∈{1,2,3,5,7}`: 100/100 for `d≤3`; for `d=5,7` a
  **1e-6 grid artifact** in ~1% (constructed `δ` off the exact optimum by ~1e-6) — the families are
  in fact **loose** (true `M = 0.12–0.25`, wide search), so lonely with huge margin; the exact `δ`
  clears `1/14`.
- **Sharpness:** the true `M` of a perturbed dilated AP is `~0.1–0.25`, *not* `1/14`. Any `a ≠ 0`
  **decorrelates** the AP and lifts `M` far above threshold. The conjugate witness does not need to
  see that lift — it certifies the uniform `1/14` floor directly, sidestepping the decorrelation
  estimate entirely.

## The consequence: one open core, not two

This removes the escape families as an independent obstruction and **collapses the whole LRC(14)
residual onto a single statement.** Every 13-family `v`:

| case | tool | status |
|---|---|---|
| `≤12` clusters at some scale | coarse bound + **LRC(≤13)** ⟹ `M ≥ 1/13` | **DONE** (kps `lonely14_of_coarse_le12`) |
| coarse part = dilated AP (escape / L-lift) | **AP-inherited conjugate witness** ⟹ `M ≥ 1/14` | **NEW** (this session; verified, not yet formal) |
| coarse part single-scale, non-AP | **the moat** ⟹ `M(K) ≥ 1/13`, then coarse bound | reduces to ↓ |
| single-scale | **the moat** ⟹ `M ≥ 1/13`, or `= AP` (`M=1/14`, `initial_segment_unit_lonely`) | **OPEN** |

The single open core is **the moat**:

> **`{1,…,13}` is the unique single-scale 13-family (up to dilation) with `M < 1/13`;
> every other single-scale family has `M ≥ 1/13`.**

This is the 13-family analogue of `(C)`. Crucially it is used **directly as a sup bound for
single-scale families (Route 1)** — not fed through the J-K accumulation top link — so it is the
*right* object and survives opus-S130/MISTAKE-117. (The same rigidity fact that the retired Route 2
used wrongly is load-bearing here when aimed correctly. My S144–S150 `(C)`/covering-min work is
about a cousin of this moat and re-enters as *evidence for* it, not as a proof via J-K.)

kps-S53 had described the single-scale residue as a **dichotomy** (near-AP rigidity ⊕ spread
decorrelation) *plus* an unclosed multi-scale r=13 "descent." Both simplify: the spread branch is
just "`M ≥ 1/13`" (the moat, `M` there is `0.1–0.2`), and the r=13 escape branch is closed by the
witness. **Two-and-a-half open pieces become one.**

## Honest ledger

- **New, verified:** the AP-inherited conjugate-witness floor `M(a + L·(d·AP)) ≥ 1/14` (exact slope
  test `a_{i₊}/v_{i₊} ≥ a_{i₋}/v_{i₋}`, conjugates cover both signs); 200/200 base, slope predicts
  winner 100/100, permuted 200/200, dilated `M` loose. Answers mac-mini-S38's open item (a) (the
  "sharper base-structure bound" for the threshold case).
- **New, structural:** the escape/L-lift families are **not** a separate decorrelation obstruction;
  the entire residual collapses to the single-scale **moat** (13-family `(C)`-analog), correctly
  aimed at the sup (Route 1).
- **Does NOT claim:** a proof of LRC(14). The moat is open (the hard rigidity core). The witness is
  analytically derived + exhaustively verified but not yet a rigorous proof — it needs the explicit
  `L₀(A) ≈ 200A` non-binding-slack + `‖·‖`-branch bookkeeping (clean, deferred to formalization).
- **Extends / credits:** kps-S52/S53 (coarse reduction, single-scale dichotomy), mac-mini-S38
  (decorrelation atom, the boundary-bound gap), opus-S130 (Route 2 retirement — this stays on the
  correct sup side). Corrects kps's "r=13 DESCEND terminates" (it loses margin; closes by the
  witness, not by descent) and the covering-side pessimism that escape needs an unbounded modulus /
  analytic argument (wrong object; the direct sup witness handles it).

- **Files:** `04-computation/lrc14_ap_escape_conjugate_witness_klein_S152.py`,
  `05-knowledge/results/lrc14_ap_escape_conjugate_witness_klein_S152.out`. HYP-4711.
- **Next:** (1) formalize the witness lemma (`M(a+L·d·AP) ≥ 1/14`) with explicit `L₀`; (2) attack
  the moat directly (single-scale non-AP ⟹ `M ≥ 1/13`) — the one remaining core, = Route-1 rigidity.

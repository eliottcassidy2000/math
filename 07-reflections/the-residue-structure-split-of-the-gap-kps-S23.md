# The residue-structure split of the gap: the multiplicative side has a *wider* gap, and the razor-thin edge lives entirely on the additive side

*kind-pasteur-2026-07-06-S23 — synthesizing the many lenses (additive/multiplicative,
odd/even, ±, rational/irrational) into a single decomposition of (G) by residue
structure mod the prime, and connecting the LRC tight locus to the tournament
project's `Z_2` complement involution.*

## The finding, in one line

Split every `n = 13` covering family by its residue set mod 13:

- **Full residue system** (residues `= {1,…,12}`, "residue-pinned"):
  `M = 1/13` (the AP orbit, unique) **or** `M ≥ 1/12`. **Never in `(1/13, 1/12)`.**
- **Residue collision** (some class repeated, another missing): this is where the
  entire near-gap ladder `m/(12m+1)` lives, in `(2/25, 1/12)`.

Verified: 25 000 residue-pinned lifts (bounded height) → **zero** in
`(1/13, 1/12)` (the general statement is the multiplicative M-minimizer, still
open — see below); every ladder family `{1,…,11, 12m}` (`m ≥ 2`) is a collision
(`12m ≡ −m mod 13`, colliding with runner `−m`). The AP (`m = 1`, residue `12`)
is the sole full-system member.

The ladder value itself is mac-mini's **sum-lever** (HYP-4422) in the flesh: in
`{1,…,11, 12m}` the pair `(1, 12m)` sums to exactly `12m + 1`, so the witness
denominator is `q = 12m + 1` and `M = m/(12m+1)`. The additive lever *produces*
the ladder; the collision residue `−m` is its mod-13 shadow.

So the two sides of opus's sum-product coincidence (HYP-4396) are **cleanly
separated by residue structure**, and the crux gap `(1/13, 2/25)` sits exactly at
their junction.

## Why this is the additive/multiplicative dichotomy made discrete

opus (HYP-4396) and I (HYP-4407) framed the AP as the coincidence of the
**additive** interval `[1,12]` and the **multiplicative** group `(ℤ/13)*`. The
residue split makes the two sides literal and disjoint:

**Multiplicative side = full residue system.** Residue-pinned means the runners
hit every nonzero class once — at `t = 1/13` they sit at all the nonzero 13th
**roots of unity** (margin exactly `1/13`). The rigidity here is `(ℤ/p)*`-flavored:
perturb within the full system and `M` jumps past a **wide** gap, all the way to
`1/12` — not the razor-thin `2/25`. The multiplicative case of (G) is the
**generous** case. (`1/12 − 2/25 = 1/300`, but the *gap* it skips,
`(1/13, 1/12)`, strictly contains the crux gap.)

**Additive side = residue collision.** A collision means two runners share a
class and one class is vacant — the set is no longer a full residue system but an
**arithmetic interval with a defect**. This is where the theta/Farey/covering
mechanics live. The ladder is the signature: `{1,…,11, 12m}` has big-runner
residue `−m mod 13`, and as `m` climbs `2,3,4,…` the **collision class walks down**
`11, 10, 9, …` while `M` walks **up** the ladder `2/25, 3/37, 4/49, …`. The
razor-thin `2/25` edge (the `m = 2` family `{1,…,11,24}`, mac-mini's second-best
covering) is a **collision** family — the crux's hardest point is purely additive.

**The junction.** The crux gap's two edges are the two sides meeting:
`1/13` is the full-system AP (multiplicative), `2/25` is the first collision
(`m = 2`, additive). The gap `(1/13, 2/25)` is empty because you cannot be *between*
"full system, roots of unity" and "first arithmetic defect" — the residue set is
either complete (jump to `1/12`) or has a defect (start the ladder at `2/25`).

## What it buys the proof (the leverage)

This decomposition **splits (G) into two independent, differently-flavored
obligations**, and one of them is essentially done:

1. **Residue-pinned ⟹ `M ∉ (1/13, 1/12)`** — the roots-of-unity / `(ℤ/p)*`
   rigidity (this is mac-mini's M-minimizer HYP-4392, *but with a wider target*:
   the pinned families avoid the whole `(1/13, 1/12)`, so a **cruder** bound
   closes them than the crux `2/25` demands). The multiplicative case is the
   easier half.
2. **Residue-collision ⟹ `M ∉ (1/13, 2/25)`** — the collision families start
   their ladder at `2/25` and climb. The `{1,…,11, v}` collisions are **already
   closed** (`slice11_loose`, HYP-4357 + opus HYP-4366): every `v ≥ 13` is loose.
   The remaining collisions are the non-`{1..11}`-shaped defects — the genuine
   additive residual, where the density floor / three-gap quantization
   (mac-mini HYP-4412) does its work.

So the razor-thin edge is **isolated to the additive/collision side**, and the
multiplicative/pinned side is the generous `1/12` case. That is a real narrowing:
prove the pinned case with room to spare, and the hard analysis is confined to the
collision families, whose simplest and sharpest members I have already formalized.

## The `Z_2` involution: the LRC tight locus is "self-complementary"

The dichotomy has an odd/even (order-2) heart that ties the LRC crux to the
**tournament half of this very project**. The tight locus is invariant under the
involution `σ : i ↦ −i ≡ 13 − i` — the unique order-2 element of `(ℤ/13)*`
(verified: every dilated AP is `σ`-invariant; the AP's witness points `{i/13}` are
reflection-symmetric through `1/2`). Under `σ` the resonance sum pairs `k ↔ −k`
into conjugate Fourier terms — **real, sign-balanced** (the ± dichotomy: the
covering cancellation is `σ`-symmetric).

This is the **same `Z_2`** that governs the tournament side: the complement
involution `T ↔ T^op`, the self-complementary tournaments (THM-024: every SC
tournament has a `−1` anti-automorphism), the merged metagraph `G_n / Z_2`. The
AP is the LRC's **self-complementary fixed point** — the roots-of-unity
configuration invariant under `−1`, exactly as an SC tournament is a complement
fixed point invariant under arc-reversal. One order-2 element of the governing
group pins the extremal configuration on **both** halves of the project. The LRC
tight locus and the SC tournaments are two faces of the same `Z_2`.

## Convergence with mac-mini's witness-denominator lever (HYP-4422)

The owner dispatched the four-dichotomy prompt to several machines at once;
mac-mini's HYP-4422 is the **additive** companion to this reflection, and the two
weld exactly. Their proven lever: if `M(S) = c/q` in lowest terms then
`q | (v_i ± v_j)` or `q | 2v_i` — so `q ≤ 2·max v_i`, and **bounding height bounds
the witness denominator** (the finite-check that answers MISTAKE-110's
unboundedness, and the additive realization of my S22 "bound the off-13-grid
denominators" guidance).

The residue split *is* their lever read at `q = 13`:

- **Collision** `v_i ≡ v_j (mod 13)` ⇔ `13 | (v_i − v_j)` — the **difference**
  lever at `q = 13`. The near-gap ladder families are exactly those with a
  difference-collision.
- **Residue-pinned** (all classes distinct) has **no** difference divisible by 13
  — its `q = 13` witness comes instead from the **sum** lever: antipodal pairs
  `v_i + v_j ≡ i + (13−i) = 13 ≡ 0 (mod 13)`, which is mac-mini's POS/NEG reading
  and the `σ`/`T^op` reflection below. To beat `1/13`, a pinned family must move
  to a *different* denominator `q ≠ 13`, and the nearest one lands it at `≥ 1/12`
  — the wider gap. So "full system ⇒ sum-lever `q=13` only ⇒ jump to `1/12`" is
  the denominator-lever explanation of the residue split.

mac-mini also reached the antipodal/`T^op` involution independently — good
cross-confirmation of the `Z_2` bridge. This reflection's distinct contribution is
the **residue-structure decomposition itself**: full system (wider `1/12` gap,
multiplicative) vs collision (the ladder, additive), the collision-class *walk*
down as `M` climbs, and the observation that the razor-thin `2/25` edge is
confined to the collision side — where `slice11_loose` already lives.

## Pointers

- `lrc_involution_signbalance_kps_S23.py`, `lrc_residue_split_kps_S23b.out` — the
  residue split (pinned avoids `(1/13,1/12)`; ladder = collision) and the
  `σ`-invariance / sign-balance of the tight locus.
- opus HYP-4396 (sum-product), kps HYP-4407 (`(ℤ/p)*`-orbit), mac-mini HYP-4412
  (three-gap / CF quantization), HYP-4392 (M-minimizer, the pinned obligation).
- kps HYP-4357 (`m/(12m+1)` ladder), `slice11_loose` / opus HYP-4366 (the
  `{1..11,v}` collision families, closed).
- Tournament side: THM-024 (SC ⟹ `−1` anti-aut), the merged metagraph `G_n/Z_2`,
  `07-reflections/merged-metagraph-invariants.md`.

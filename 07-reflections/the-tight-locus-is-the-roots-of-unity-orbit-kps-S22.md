# The tight locus is the `(ℤ/p)*`-orbit of the roots of unity — and that is why primality gives rigidity

*kind-pasteur-2026-07-06-S22 — a synthesis reflection on the open remainder of
LRC(14), pulling together mac-mini's prime/composite dichotomy (HYP-4382), opus's
discrepancy inversion (HYP-4074/4013), the Farey ladder (HYP-4357), and the
strict lift-rigidity (HYP-4392) under a single multiplicative-group picture.*

## The one picture

Peel LRC(14) to its `n = 13` base. The **tight locus** — the covering families
with `M = 1/13`, the gap's lower edge — is exactly one object:

> the orbit of the nonzero 13th roots of unity under the multiplicative group
> `(ℤ/13)*`.

At `t = 1/13` the AP `{1,…,12}` places its runners at `e^{2πik/13}`, `k = 1…12`
— every nonzero 13th root of unity, once each. That is the **perfectly
equidistributed** configuration: minimal star-discrepancy (opus's inversion,
HYP-4074), minimal covering threshold, minimal `M`. Dilating the speeds by
`d ∈ (ℤ/13)*` sends `k ↦ dk mod 13` — a **permutation of the same roots**
(the multiplicative action of `(ℤ/13)*` on `μ₁₃ ∖ {1}`), so `M(d·AP) = 1/13`
still. Verified: `d = 1,…,12` all tight.

So the tight locus is a **single `(ℤ/p)*`-orbit**, and the rigidity question
"`M = 1/13 ⇒ dilated AP`" (mac-mini HYP-4392) is:

> the roots-of-unity configuration is an **isolated** minimum of the covering
> threshold — no lifted (non-orbit) configuration reaches it.

Verified this session (lrc_nearap_rigidity_kps_S22): over 15 000 residue-fixed
`13`-lifts of the AP, **zero** tight non-dilations and **zero** gap members; the
only `M = 1/13` families in the `{0,13}`-lift cube are the AP and its `d = 2`
dilation. The orbit is isolated.

## Why primality is the whole story

`residue_pinning_13` (formal, mac-mini) says `M = 1/13` forces the residues mod
13 to be a **full residue system**. Here is the multiplicative punchline:

- **`p` prime** ⇒ `(ℤ/p)*` is cyclic of order `p−1`, and a full residue system
  of nonzero classes is `{1,…,p−1}` — forced to be the AP (up to the `(ℤ/p)*`
  dilation). One orbit, one shape. **Rigidity is clean.**
- **`n` composite** ⇒ `(ℤ/n)*` has proper subgroups and zero-divisors; a
  "full-ish" residue set can avoid a coset and still cover. mac-mini's `n = 6`
  witness `{1,3,4,5,9}` is exactly a **non-AP tight set** living off the AP
  orbit — the orbit is no longer unique, and the AP-tower induction that would
  climb through composite `n = 12` is **unsound** (HYP-4382).

So "the AP is the unique tight locus **because 13 is prime**" (mac-mini
OPEN-Q-108 R3) is the statement that `(ℤ/p)*` acts **transitively with a single
orbit-shape** on full residue systems. Composite `n` breaks transitivity; prime
`p` restores it. The endgame's insistence on peeling to the **prime** `n = 13`
(not the composite `12`) is not a convenience — it is the arithmetic that makes
the tight locus a single rigid orbit.

## Three classical shadows this casts

The same object appears in three classical theories, and each lends the crux a
tool:

1. **Cyclotomy / equidistribution.** The AP = `μ_p ∖ {1}`. Discrepancy is
   minimized exactly at roots of unity (Erdős–Turán: discrepancy is controlled
   by the exponential sums `Σ e(k v_i t)`, whose resonances `Σ k_i v_i = 0` are
   the AP's rich additive relations). The **density floor** ("AP = min
   star-discrepancy", opus HYP-4013) is the quantitative isolation of this
   minimum — a discrepancy-lower-bound problem, classical in shape.

2. **Covering systems.** `M(S) < β ⇔` the combs `A_i = {t : ‖v_i t‖ < β}`
   **cover** `[0,1]`. The AP-covering at `β → 1/13⁺` is the continuous analog of
   an **exact covering system** with the tightest moduli, and
   Davenport–Mirsky–Newman–Rado rigidity (an exact cover forces a repeated
   largest modulus) is the discrete shadow of "distinct-frequency combs can't
   tile" — the covering-impossibility I formalized (`torus_split_rung`,
   `CircleClearFloor`). The razor-thin cancellation mac-mini found (no low-order
   Bonferroni certificate at `2/25`, HYP-4402) is exactly the DMN phenomenon:
   the covering is a **global** cancellation, invisible to any finite-order
   local expansion.

3. **The Farey ladder.** The near-gap spectrum is the Stern–Brocot chain
   `m/(12m+1)` (HYP-4357): `1/13` (the AP orbit, `m=1`), `2/25` (`m=2`), `3/37`,
   … → `1/12`. The AP is the **deepest rung**, and the gap `(1/13, 2/25)` is the
   Farey interval between the two deepest — empty because the ladder is the
   quantization of the covering threshold near the roots-of-unity minimum.

## The convergence with opus's sum-product reframe (HYP-4396)

Concurrently, opus reframed the whole crux as a **sum-product rigidity**: at the
prime, `{1,…,12}` is *simultaneously* the additive interval `[1,p−1]` (maximal
relation lattice `L(S) = {m : Σ mᵢvᵢ = 0}`, theta-extremal — the **additive**
side: Farey, mediant, grid-attainment, my ladder) **and** the multiplicative
group `(ℤ/13)*` (character-rigid, residue-pinned — the **multiplicative** side).
`(G)` is the rigidity of that coincidence: a set can't be maximally additive
*and* multiplicative unless it is (essentially) a subfield, and `[1,p−1]` is the
extreme case. Perturb and you break additive richness (theta stops cancelling ⇒
`safe > 0` ⇒ loose) **or** multiplicative pinning (residues leave the units),
never neither.

This reflection is the **multiplicative face** of opus's coincidence, made
geometric: `(ℤ/13)*` acting character-rigidly *is* the dilation orbit of the
roots of unity `μ₁₃ ∖ {1}`, and "residue-pinned" *is* "sits on the orbit." What
I add to the governing frame: (i) the roots-of-unity **geometry** and the
discrepancy-minimum **isolation** language; (ii) the **covering-system**
(Davenport–Mirsky–Newman–Rado) shadow of the additive side, which is where the
razor-thin no-low-order-Bonferroni cancellation comes from; (iii) the concrete
**orbit-isolation confirmation** (strict rigidity verified); (iv) a bridge
between the two sides, below.

## The character-sum bridge: where the two sides meet

opus's proof-strategy hypothesis is that the hard **additive** density floor is
the *shadow* of the clean **multiplicative** pinning — "push the character-sum
link on `(ℤ/13)*`, not either side alone." Here is exactly that link, concrete:

> **A `13`-lift of the AP is invisible at every `t = a/13`.** For
> `vᵢ = i + 13cᵢ`, `vᵢ·(a/13) = i·a/13 + cᵢa ≡ i·a/13 (mod 1)`, so
> `‖vᵢ·a/13‖ = ‖i·a/13‖` — lift-independent.

At the `13`-rational points (the multiplicative points, where the roots-of-unity
/ character structure lives) **the lift does nothing**: margin is exactly `1/13`,
pinned. So a lift can only raise `M` at **other denominators** — the additive
points. This is why the multiplicative pinning is "done" (it is a finite
`(ℤ/p)*` fact, opus's `residue_pinning_13` and `margin_of_residue_witness`) while
the density floor is "open" (it is an additive/theta fact at the non-`13`
denominators). The bridge — and the guidance for mac-mini's height-bound
reduction (HYP-4392) — is: **the rigidity witness lives off the `13`-grid; bound
those denominators (Lipschitz + the `1/300` seam) and it is a finite check.**

## What this buys the open remainder

The three open pieces are one isolation statement seen three ways:

- **Strict lift-rigidity** (HYP-4392) = the `(ℤ/p)*`-orbit is isolated among
  lifts. The clean route is mac-mini's **height-bound reduction**: the roots-of-
  unity behavior is *identical* for every lift at the rational points `t = a/13`
  (there `v_i a/13 ≡ i a/13 mod 1`, lift-independent), so a lift can only help at
  **other** denominators — bound those denominators (Lipschitz + the `1/300`
  seam) and the rigidity is a finite check, exactly the certificate shape
  (`rational_point_margin`) my ladder used. The orbit picture says *where* to
  look: away from the `13`-power-of denominators.

- **Density floor** (HYP-4013) = the discrepancy minimum at `μ_p ∖ {1}` is
  strict, with a quantitative gap to the next configuration. Erdős–Turán with
  the AP's resonance set is the classical hammer.

- **Gap-emptiness (G)** = the covering threshold **quantizes** onto the Farey
  ladder near the minimum (skips `(1/13, 2/25)`), which is the DMN rigidity of
  the near-exact cover. mac-mini's two-scale decorrelation already isolates this
  to the single-scale near-AP kernel; the orbit isolation is its heart.

The endgame is not three problems. It is one: **the roots-of-unity configuration
is a rigid, isolated optimum, and it is unique because 13 is prime.** Cyclotomy,
discrepancy, and covering systems are three classical theories all describing
that one isolation — and the fleet's three open pieces are their three faces.

## Pointers

- `lrc_nearap_rigidity_kps_S22.py` / `.out` — the tight-locus orbit + strict-
  rigidity confirmation (dilations tight; no tight non-dilation; gap empty).
- mac-mini HYP-4382 (prime/composite tight locus), HYP-4392 (strict lift-
  rigidity), HYP-4402 (two-scale descent to the single-scale kernel).
- opus HYP-4074/4013 (discrepancy inversion, AP = roots of unity = min
  discrepancy), HYP-4376 (`M ≥ 1/13` pinned floor).
- kps HYP-4357 (the Farey ladder), `CircleClearFloor` / `torus_split_rung` (the
  covering-impossibility), `slice11_loose` (the `{1..11,v}` face).
- Classical: Erdős–Turán discrepancy; Davenport–Mirsky–Newman–Rado covering
  systems; the three-distance theorem.

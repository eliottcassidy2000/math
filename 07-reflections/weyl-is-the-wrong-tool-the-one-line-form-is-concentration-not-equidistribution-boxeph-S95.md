# Weyl is the wrong tool: the one-line rigidity is concentration, not equidistribution

*boxeph-2026-07-18-S95. Attacking the S94 one-line form (`M<1/13 ⟹` twelve of thirteen residues
`v_i·a mod(13val+1)` are multiples of `val`) with an extremal Weyl bound, as directed. Outcome: a clean
exponential-sum reformulation, a precise proof that **Weyl bounds cannot reach it** (wrong direction),
and the redirect — Weyl belongs to the *density* route, not the rigidity. LRC(14) not closed. Verified
S95 exponential-sum computation.*

## The one-line form as an exponential sum

Let `T(ℓ) = Σ_i e(ℓ·r_i/val)`, `r_i = v_i·a mod q`, `q = 13val+1`. The offset-zero count is
`N₀ = #{i : val ∣ r_i} = (1/val)·Σ_{ℓ=0}^{val−1} T(ℓ)` (orthogonality). Since `T(0)=13`,

> **One-line form ⟺ `N₀ ≥ 12` ⟺ `Σ_{ℓ≠0} T(ℓ) ≥ 12·val − 13` ⟺ `T(ℓ)` is LARGE for all `ℓ`.**

Verified: for the deep well and the ladder `{1..12,182m}`, `N₀ = 12` and
`min_{ℓ≠0}|T(ℓ)| = 11`, `avg ≈ 12` — the residues are **maximally concentrated** on `val·ℤ`. A generic
covering family (`N₀=4`) has `|T| ≈ 2.7`, near the equidistribution scale `√13 ≈ 3.6`.

## Why an extremal Weyl bound cannot reach it (the direction is wrong)

A Weyl / Erdős–Turán bound proves **equidistribution**: it delivers an **upper** bound `|T(ℓ)| ≤ B`
with power-saving. Feeding that into the identity gives
`N₀ ≤ (13 + (val−1)·B)/val` — an **upper bound on `N₀`.** But the target is a **lower** bound
(`N₀ ≥ 12`). Equidistribution machinery is not just too weak here; it points the wrong way:

- The one-line form asserts the residues are **maximally non-equidistributed** (`|T| ≈ val-scale
  maximum`, all thirteen almost on one lattice). Weyl exists precisely to *rule out* such concentration
  for generic inputs — it can never *force* it.
- This is the sharp/soft mismatch klein-S279 already flagged for the covering side ("needs the sharp
  constant; naive Erdős–Turán is ~700× too weak"): a soft upper bound cannot certify a sharp lower
  bound / exact cancellation.

## What the extremality actually buys (2 of the 13, not the core)

The maximizer's structure (death-star THM-999: opposite-slope active pair) pins exactly **two**
residues: `v₊` at `r = val` (offset 0) and `v₋` at `r = q−val = val·12+1` (offset 1). So extremality
contributes `1 + e(ℓ/val)` to `T(ℓ)` — the active pair, `2/13` of the concentration. The remaining
**eleven** core residues being on `val·ℤ` is *not* a consequence of this maximizer's local (extremal)
geometry — those runners sit at distance strictly `> M` (slack), so the extremality does not touch
them. Their concentration is the **global** content (the family is covering with no better `t`
anywhere), which is exactly the crux (S94: equivalent to LRC(14)). An "extremal Weyl bound" splits as
[extremal part → active pair, gives 2] ⊕ [Weyl part → equidistribution, wrong direction for the other
11]; neither half reaches the core.

## The redirect: Weyl belongs to the density route

Weyl is the *right* tool on the **other** face of LRC(14) — the **density** route (Route A), which is
the **soft** side: `μ₀ > 0 ⟺ M > 1/14`, provable by a soft oscillatory Weyl bound on the arc midpoints
where **"any power-saving suffices"** (klein finish-map). There, equidistribution (an upper bound on an
exponential sum) is exactly what's needed and Weyl delivers. The rigidity (Route B / the one-line form)
is the **sharp** side and needs the opposite — a **concentration / lower-bound / inverse** tool
(sharp Freiman `3k−4` / PFR), not Weyl. Pairing Weyl with the rigidity is the wrong marriage; pairing
it with density is the right one.

## Net (honest)

- **New:** the one-line form is `min_{ℓ≠0}|T(ℓ)| ≳ val` (maximal concentration of the maximizer
  residues on `val·ℤ`), verified `≈ 12` on the deep-well tower.
- **Proved:** an extremal Weyl bound cannot attack it — it bounds `N₀` from **above**; the residues are
  maximally non-equidistributed, so equidistribution is vacuous. The extremal part reaches only the
  active pair (`2/13`).
- **Redirect:** Weyl is the density route's tool (soft, any power-saving); the rigidity needs a sharp
  concentration/inverse theorem. This is consistent with S94 (the extremality content is global, not
  local) and klein-S279 (covering needs the sharp constant).

LRC(14) is not closed. The Weyl attack on the rigidity is ruled out with its reason; the concentration
must come from the additive inverse theorem, and the soft Weyl tool should be aimed at the density
route instead.

Cross-links: [[the-crux-reduced-to-bedrock-j1-is-zero-and-the-offset-vanishing-IS-LRC14-boxeph-S94]],
[[both-lrc14-routes-bottom-on-the-same-multilinear-cancellation-not-one-mollification-klein-S279]],
[[the-abandoned-attempts-decoded-the-crux-is-additive-dimension-two-not-any-scalar-boxeph-S92]],
[[THM-1017-ap-core-bridge-reduction]].

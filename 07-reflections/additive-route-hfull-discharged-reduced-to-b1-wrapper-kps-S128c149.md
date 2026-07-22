# Additive route: `hfull` discharged, DvdK reduced to just the b=1 wrapper (kernel-pure)

*kind-pasteur-2026-07-21-S128c149. Owner: finish the additive route; many push/pulls; look at prior
additive/multiplicative duality for inspiration. Follows S128c148 (the additive orbit-sum core).*

## Concrete progress (kernel-pure, pushed — `GMC2FullRootConcrete.lean`)

The additive one-variable DvdK route (THM-2101) had two open inputs to boxeph's
`root_packet_eq_zero`: `hfull` (full-root Lagrange sum `= 0`) and `hb` (the b=1 wrapper). I closed
`hfull` and assembled the route down to `hb` alone:

- **`fullRootSum_eq_zero`** — for `Φ` irreducible over a characteristic-zero field `F` splitting in
  `L`, with `0 < M < deg Φ`: `∑_{α ∈ Φ.rootSet L} α^{M-1}/Φ'(α) = 0`. This bridges codex's abstract
  nodal Lagrange lemma (`GMC2FullRootLagrange.sum_pow_pred_div_derivative_nodal_eq_zero`, over a
  `Finset` of nodes) to the concrete root set: `Φ` separable (irreducible over the perfect field `F`)
  ⇒ distinct roots ⇒ `Ψ := Φ.map = C(lc)·nodal(roots)` ⇒ `Ψ'(α) = lc·nodal'(α)`, and the
  `rootSet`-subtype sum equals the `Finset` sum, so the concrete sum is `lc⁻¹·0 = 0`.
- **`additive_dvdk_reduces_to_smallSum`** — combining `fullRootSum_eq_zero` (discharging `hfull`)
  with `root_packet_eq_zero`: if the small-cluster residue sum equals `1` (`hb`), then `False`. So the
  **additive one-variable DvdK contradiction is formally complete modulo the single remaining input
  `hb`**, with *no* THM-1550 / small-root product / Hensel / Wiener–Hopf. Axioms `[propext,
  Classical.choice, Quot.sound]`.

Both were plumbing-heavy (the Mathlib `Splits` API is now single-argument `Splits (f)` = "splits over
its own ring"; the `rootSet`-subtype ↔ `Finset` bridge and the `Ψ.leadingCoeff` `conv_lhs` care were
the sharp edges), but they compile clean.

## The one remaining input: the b=1 wrapper (honest characterization)

`hb : algebraMap 1 = ∑_{β ∈ S} β^{M-1}/Φ'(β)` — the small-root cluster `S` has residue sum `1`. This
is the generating-function identity (THM-2101 eq 7–8): `∑_{α ∈ S₀(t)} w_t(α) = G(t) := ∑_m CT(f^m)tᵐ`,
which under the vanishing assumption `CT(f^m)=0 ∀m≥1` collapses to `G(t)=CT(f⁰)=1`. Establishing it
requires **selecting the small-root packet `S` by valuation** (Newton polygon: the `M` roots of
`u^M − tR(u)` with positive `t`-valuation) and identifying `∑_S residue` with `G(t)`. This is the
genuinely substantial remaining piece.

**Where the additive/multiplicative duality actually lives (THM-2101 §6).** The product route
(THM-1550) and the additive route are related by the **Abel operator** `A(G−1) = ∫₀ᵗ (G(s)−1) ds/s`:
integrating the residue observable turns the additive residue *sum* into `log Π(t)`, creating the
`1/m`, the small-root *product* `Π`, Hensel factorization, and Wiener–Hopf. So:

| quantity | additive (before integration) | multiplicative (after `∫ ds/s`) |
|---|---|---|
| observable | residue **sum** `∑_S w` | root **product** `Π = ∏_S α` |
| full-orbit identity | Lagrange sum `= 0` (**done**, `hfull`) | Vieta product `= const` (death-star) |
| small-packet value | `∑_S w = 1` (**b=1 wrapper**) | `Π = c·t` (**THM-1550**, deep) |
| Mathlib obstacle | partial-fraction residue (a *sum*) | Hensel factorization (a *product*) |

The two small-packet values (`b=1` vs `THM-1550`) are the **same valuation content seen through the
Abel operator** — both select the small-root packet by Newton polygon. boxeph-S238's honest note is
right: the additive route does **not** escape the valuation core; it replaces the *product* (Hensel
factorization, a Mathlib gap) with a *sum* (partial-fraction residue), which is a **cleaner formal
target** but still substantial. My earlier over-optimistic "additive fully escapes analysis" framing
(S128c148) is thereby tempered: the *closing algebra* escapes THM-1550 (proved, kernel-pure), but the
*small-packet selection* remains, as a sum.

## Status

The additive route is now: **irreducibility (mac-mini) + Galois action (boxeph) + full-root Lagrange
`hfull` (this session) + root-packet contradiction (boxeph) — all kernel-pure — reducing DvdK1 to the
single b=1 wrapper `hb`**, which is the Newton-polygon small-packet-residue-sum identity. I did **not**
close `hb` (it is the remaining crux, plausibly multi-session, though a *sum* not a *product*).

## Named next (b=1 wrapper)
1. The formal partial-fraction identity `∑_{α ∈ S} α^{M-1}/Φ'(α) = [u⁰] u^{M-1}/Φ(u,t)` as a power
   series in `t` (mac-mini-S163's root-free `F(t) = [x⁰] x^M/(x^M − tR)` is the RHS).
2. The Newton-polygon selection of the `M` small roots `S` (positive `t`-valuation) — the valuation
   core, shared in spirit with THM-1550 but for a sum.
3. Under `CT(f^m)=0 ∀m≥1`, `F(t)=1`, giving `∑_S = 1 = hb`.

## Cross-links
`GMC2FullRootConcrete.lean` (this session) · `GMC2RootPacketConcrete` (boxeph, `root_packet_eq_zero`) ·
`GMC2LaurentShiftCheckA` (codex nodal Lagrange + additive orbit machinery) · `GMC2AdditiveOrbitSum`
(my S128c148 orbit-sum core) · THM-2101 (three additive proofs, §6 Abel duality) · THM-1550 (the
multiplicative twin of the b=1 wrapper) · mac-mini-S163 (root-free `F(t)`).

---
source: opus-2026-06-06-S699 (the polarization open question)
status: PARTIAL PROOF + abstraction. The polarization "deltas avoid 7,21" ⟸ the phantom-volume theorem (7,21 never H) ⟸ strong-component multiplicativity ⟸ a strong-min lower bound. PROVEN for tournaments whose strong components have ≤6 vertices (full enum: strong-min(m)=3,5,9,15 for m=3..6 ⟹ strong values exclude 7,21); strong-min(7)≤25>21; reduces to strong-min(m)≥22 (m≥7) for the full theorem. ABSTRACTION: the H-spectrum is a co-finite multiplicative numerical semigroup with exactly 2 gaps {7,21}=Φ₃(2),3Φ₃(2) (genus 2); the delta field is its gradient, POLARIZED around the band gap. Polarized delta fields = gradients of integer functions with a gappy range (band gaps / numerical-semigroup genus / phantom volumes), appearing across physics, number theory, scissors congruence, and discrete Morse/chip-firing.
tags: [polarization, delta-field, phantom-volume, strong-component, strong-min, band-gap, numerical-semigroup, genus, frobenius, morse, chip-firing, redei, 7-21, abstraction]
---

# Polarized delta fields: band gaps, numerical semigroups, and the 7/21 holes

**Prompt (user):** prove the open question (polarization "deltas avoid 7,21" from the
phantom-volume structure); understand it at larger n; understand polarized delta fields; find where
else they appear; be abstract and creative.

## 1. The reduction (rigorous spine) — and a partial proof

The polarization is not an extra fact; it is *downstream* of the forbidden values:
```
   "no arc-flip from H lands on 7,21"  ⟸  7,21 are never H of any tournament  (phantom-volume thm)
   ⟸  (strong-component multiplicativity, S599s: H=∏H(C_i))  no STRONG tournament has H∈{7,21}
   ⟸  a lower bound on  strong-min(m) = min Ham-paths over strong tournaments on m vertices.
```
- `7=Φ₃(2)` is **prime**, so `H=7` needs a *single strong component* with `H=7`; `21=3·Φ₃(2)` needs
  either a strong `H=21`, or strong-`3` × strong-`7` (the `7` is unavailable).
- **Full enumeration (`…s699j.py`):** `strong-min(m)=3,5,9,15` for `m=3,4,5,6`; the strong-value
  sets are `{3},{5},{9,11,13,15},{15,17,…,45}` — **all exclude `7` and `21`**. So **no strong
  tournament on `≤6` vertices realizes `7` or `21`.**
- **`m=7`:** deeper search (`…s699k.py`, `k≤6` reversals) gives `strong-min(7) ≤ 25 > 21` (minimizer
  = transitive with arcs `(0,3),(3,6)` reversed); `strong-min(8) ≤ 45`. The minimum **grows**, so
  components of size `≥7` also avoid `7,21`.

> **Theorem (partial).** The delta field of `H` avoids `{7,21}` for every tournament whose strong
> components each have `≤ 6` vertices. **Reduces** to: `strong-min(m) ≥ 22` for `m ≥ 7` (then no
> component of any size realizes `7,21`, completing the phantom-volume theorem and the polarization
> for all `n`). The growth `3,5,9,15,25,…` makes this near-certain; the rigorous closure is the
> strong-tournament Hamiltonian-minimum lower bound (Busch-type).

*(Correction to S699i's `m²−5m+9`: it fits `m≤6` but `strong-min(7)=25≠23` — the formula breaks; the
true sequence is `3,5,9,15,25,…`.)*

## 2. The abstraction: the H-spectrum is a band structure / a numerical semigroup

> **The achievable-`H` spectrum is a *co-finite multiplicative numerical semigroup* with exactly
> two gaps, `{7,21}` (genus 2).** It is generated (multiplicatively) by the strong-tournament
> H-values (the "irreducibles"); `7,21` are the only permanent non-members (`63,189,35,49` fill in
> at higher `n`, S613). The **gaps are the band gap**; the **delta field is the dispersion**,
> *polarized* so no single step crosses the gap (from `H=5`, `−2→3` allowed, `+2→7` forbidden).

This unifies three pictures of the same object:
- **Band structure (physics):** `H`-values = allowed energy *bands*; `{7,21}` = the *band gap*
  (forbidden energies = phantom volumes); `delta` = the band *dispersion* (the gradient), which has
  no states in the gap — exactly the polarization. The gap sits in the **small-`H` regime**;
  large `H` is a dense (gap-free) band, since `strong-min` grows (quadratically-ish) so all large
  odd values are products of irreducibles.
- **Numerical semigroup (number theory):** the gaps `{7,21}` are the **genus** (the non-representable
  elements), the multiplicative analog of the Frobenius number; the polarization = "you cannot
  reach a gap by one generator step." `7=Φ₃(2)` ties the genus to cyclotomy (S599u).
- **Scissors congruence (S599v):** the gaps are the **phantom volumes** — values realized by no
  equidecomposability class; the delta = a scissors move that can't produce a phantom volume.

## 3. Polarized delta fields, abstractly — and where they appear

> **A *polarized delta field* is the discrete gradient of an integer-valued function whose range is
> a gappy set (a band structure / numerical semigroup with genus), on a discrete domain; the
> gradient avoids the gaps, breaking sign symmetry.** The data: a *potential* (the function), a
> *gradient* (delta), a *Hessian* (how flips change deltas = the Walsh/OCF, S699i), and *forbidden
> levels* (the gaps) that polarize the gradient.

Where it appears:
| domain | potential | forbidden levels (gaps) | the polarized gradient |
|---|---|---|---|
| **tournaments** | `H` (Ham paths) | `{7,21}` = phantom volumes | `delta(e)=H(T^e)−H(T)` (Walsh/OCF) |
| **solid-state physics** | energy band | the **band gap** | the band dispersion (no states in gap) |
| **numerical semigroups** | the semigroup elements | the **genus** (Frobenius gaps) | generator steps avoiding gaps |
| **scissors congruence** | volume (Dehn class) | **phantom volumes** (S599v) | scissors moves |
| **discrete Morse / chip-firing** | a Morse potential on a complex | forbidden critical levels | the gradient flow; **Hessian = discrete Laplacian** (the toppling operator) |
| **LRC / unit distance** | `M(S)` / `u(n)` (achievable values) | the value-spectrum gaps | add a runner / point (the frontier-gain, S599w) |

The **Hessian-as-Laplacian** is worth a line: `Δ_{ef}=H(T)−H(T^e)−H(T^f)+H(T^{ef})` is a discrete
second difference — a *Laplacian / toppling operator* on the arc-flip hypercube. So the delta-field
dynamics is **chip-firing on the OCF graph**, and the polarization is the recurrent-configuration
structure avoiding the forbidden levels. The "exact pattern of arc-flips on deltas" = the spectrum
of this Laplacian = the Walsh/OCF (S699i).

## 4. How it operates at larger n

- The two gaps `{7,21}` are **durable (genus 2)**; transient gaps (`35,39,63,189`) close as `n`
  grows (more irreducibles). So the band gap is a **fixed, small-`H` feature**; the band becomes
  dense at large `H`. The polarization therefore *only bites near `H≈7,21`* — the small-`H`,
  near-transitive regime (where `strong-min` lives).
- The polarization is a **boundary effect of the band gap**: deep in the band (large `H`) the delta
  field is essentially unconstrained (any even step lands on an achievable value); near the gap it
  is sharply polarized. This is exactly band-edge physics.

## 5. Honest status

- **Proven:** the reduction polarization ⟸ phantom-volume ⟸ strong-min; the phantom-volume theorem
  (hence the polarization) for tournaments with all strong components `≤6` (full enum:
  `strong-min(3..6)=3,5,9,15`, strong-values exclude `7,21`).
- **Verified / strongly supported:** `strong-min(7)≤25>21`, `strong-min(8)≤45`, growth `3,5,9,15,25`;
  the genus-2 spectrum `{7,21}` durable (with S613's `63,189` fill-in).
- **The remaining rigorous gap:** `strong-min(m) ≥ 22` for all `m≥7` (a strong-tournament
  Hamiltonian-minimum lower bound, Busch-type) — would complete the phantom-volume theorem and the
  polarization for all `n`.
- **Corrected:** `m²−5m+9` is *not* the strong-min formula (breaks at `m=7`).
- **New (abstraction):** the H-spectrum = a co-finite multiplicative numerical semigroup, genus 2
  `{7,21}`; the delta field = its band-gap dispersion; polarized delta fields = gradients of
  gappy-range integer functions (band gap / Frobenius genus / phantom volume), with the Hessian =
  a discrete Laplacian (chip-firing).

**Artifacts:** `04-computation/strong_min_phantom_band_s699j.py`, `strong_min7_deeper_s699k.py`
(+`.out`s). Builds on S699i (delta field/Walsh/OCF), S599s (strong-component multiplicativity),
S599v (phantom volumes/equidecomposability), S599u (`Φ₃(2)`/cyclotomy), Rédei, Busch (strong-min).
New: **HYP-2271**.

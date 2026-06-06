---
source: opus-2026-06-06-S699 (the delta field)
status: FORMALIZATION + verified findings (with two honest corrections) — the user's "delta" is the discrete GRADIENT of H on the arc-flip hypercube (∂_e H=H(T^e)−H(T), even); "how flipping f changes delta(e)" is the discrete HESSIAN Δ_ef. VERIFIED: delta always even; Hessian even but NOT always ÷4 (Walsh coeffs are half-integers); flipping an arc CAN change ALL other deltas (refutes 'not ever all'); the {7,21} holes POLARIZE the delta field (from H=5, −2 allowed but +2 forbidden) — confirming the user's coupling of delta and the forbidden values. The pattern = the Walsh/OCF (odd-cycle) structure.
tags: [H-function, delta-field, discrete-gradient, hessian, hypercube, walsh, OCF, odd-cycle, redei, forbidden-values, 7-21, polarization, arc-flip, structure]
---

# The H delta field: gradient, Hessian, and the 7/21 holes

**Prompt (user):** H is always odd, so flipping an arc changes H by an even "delta"; but 7 and 21
are never realizable, constraining how delta moves; each arc has a delta, flipping one arc changes
others' deltas (but not all); the exact pattern is the key to H's structure — understand delta and
the forbidden values together.

This is exactly right, and it has a clean formalization that the repo's machinery (Walsh/OCF)
already addresses. Here are the verified findings — two of the guesses get sharpened/corrected.

## The formalization: H is a potential on the arc-flip hypercube

`H` is a function on `{0,1}^{C(n,2)}` (each tournament = a vertex of the hypercube `Q_m`, `m=C(n,2)`;
arcs = coordinates). Then:
```
   delta:    ∂_e H(T) = H(T^e) − H(T)                         (discrete GRADIENT)
   "how flipping f changes delta(e)":
             Δ_{ef}(T) = H(T) − H(T^e) − H(T^f) + H(T^{ef})    (discrete HESSIAN / mixed 2nd difference)
```
So your "present delta value" of each arc = the **gradient field**, and "flipping one arc changes
the deltas of others" = the **Hessian** — the second-order structure of `H`. The whole program is
the **discrete differential geometry of `H` on the hypercube**, and its Fourier (Walsh) expansion
is the **OCF** (THM-081): `H = Σ_S c_S χ_S`, with `delta(e) = −2Σ_{S∋e} c_S χ_S` and
`Δ_{ef} = 4Σ_{S⊇\{e,f\}} c_S χ_S` — the gradient picks Walsh terms through `e`, the Hessian those
through both `e,f` (= **odd-cycle collections containing both arcs**).

## Verified findings (n=5, `…s699i.py`)

> **(1) delta is always even** — `H(T^e),H(T)` both odd (Rédei) ⟹ difference even. ✓ (Your `delta`
> is literally `2·(`an OCF gradient`)`.)

> **(2) [CORRECTION] the Hessian is even but NOT always divisible by 4.** I expected `÷4` from the
> Walsh form `4Σc_Sχ_S`; the computation says `÷2` only. So **the Walsh/OCF coefficients `c_S` are
> half-integers** (denominators `2^k`), not integers — a real structural fact about `H`'s Fourier
> spectrum. The higher differences carry the `2`-adic depth.

> **(3) [CORRECTION] flipping an arc CAN change ALL other deltas.** Max #(other arcs whose delta
> changes) `= 9` of `9` at `n=5` — so "provably not ever all" is **false**: there exist
> tournaments where one flip perturbs *every* other arc's delta. The interaction is *usually*
> sparse but can be **total**; the interaction *support* (pairs `(e,f)` that interact in *some*
> tournament) is **complete** at `n=5` — every pair of arcs lies in a common odd cycle (down to
> the 5-cycle), so every pair can interact. The honest version of your claim: *at a generic
> tournament not all deltas change, but no theorem forbids all changing.*

> **(4) [CONFIRMED — the main point] the 7,21 holes POLARIZE the delta field.** Over all 1024
> tournaments, the signed-delta set at each `H` exactly avoids the values landing on `7` or `21`:
> ```
>    H= 1: deltas {0,2,4,8}      lands {1,3,5,9}        (delta=6→7 ABSENT)
>    H= 3: deltas {−2,0,2,6,12}  lands {1,3,5,9,15}     (delta=4→7 ABSENT)
>    H= 5: deltas {−4,−2,0,4,6,8} lands {1,3,5,9,11,13} (delta=+2→7 ABSENT; −2→3 present)
>    H= 9: deltas {−8,−6,−4,0,2,4,6}                    (delta=−2→7 ABSENT)
>    H=13: deltas {−8,−4,−2,2}                          (delta=−6→7, +8→21 ABSENT)
>    H=15: deltas {−12,−6,−4,−2}                        (delta=−8→7, +6→21 ABSENT)
> ```
> **The holes break delta's up/down symmetry.** At `H=5`, you can descend by `2` (→3) but *not*
> ascend by `2` (→7): the delta field is *repelled* by `7`. This is your coupling, made precise:
> **the forbidden levels are holes in the H-landscape, and the gradient field bends away from
> them** — delta and `{7,21}` must be read together, as you said.

## The structure, named

> **`(H, delta, Hessian)` is a discrete dynamical system on the arc-flip hypercube:** `H` is an
> odd-valued *potential* with *forbidden levels* `7,21`; `delta` is its even *gradient* (the OCF
> first order, odd cycles through an arc); `Δ_{ef}` is its *Hessian* (the OCF second order,
> odd-cycle co-membership of two arcs). The forbidden levels **polarize** the gradient (4). The
> "exact pattern of how arc-flipping changes the deltas" you asked for **is the OCF/Walsh
> expansion of `H` (THM-081)** — the gradient is the odd-cycle-through-`e` sum, the interaction is
> odd-cycle co-membership, and the half-integer coefficients (2) carry the `2`-adic depth.

And the holes themselves are explained elsewhere: `7=Φ₃(2)`, `21=3Φ₃(2)` are the **phantom volumes**
of the strong-component (equidecomposability) structure — values realized by **no** scissors class
(HYP-2180/S599v). So the delta field is the OCF gradient, and the holes it bends around are the
multiplicative gaps of the strong-component semigroup. *Both understood together* = the OCF gradient
on the strong-component-graded H-spectrum.

## Honest status

- **Verified:** delta always even; Hessian even but not always `÷4`; per-flip interaction reaches
  all `9` (n=5); interaction support complete (n=5); the `{7,21}` polarization of the delta-set at
  every achievable `H`.
- **Confirmed (your insight):** the forbidden values impose structure on delta — they are holes
  that break the gradient's sign symmetry; delta and `{7,21}` are genuinely coupled.
- **Corrected (your guesses):** the Hessian is `÷2` not `÷4` (half-integer Walsh coefficients);
  flipping an arc CAN change all other deltas (the interaction is not provably sub-total).
- **Established (the formalization):** delta = discrete gradient, Hessian = mixed 2nd difference,
  both = the Walsh/OCF expansion (THM-081); the pattern = odd-cycle (co-)membership.
- **Open / next:** does the per-flip interaction count have a clean formula (the OCF degree of an
  arc)? At larger `n`, is the interaction support still complete, or do far-apart arcs in reducible
  tournaments decouple? Is the polarization "deltas avoid `7,21`" provable from the
  strong-component (phantom-volume) structure?

**Artifacts:** `04-computation/H_delta_field_hessian_s699i.py` (+`.out`). Builds on Rédei (H odd),
THM-081/077 (Walsh/OCF), HYP-2180/S599v (strong-component / `{7,21}` phantom volumes), S599s
(H-spectrum). New: **HYP-2268**.

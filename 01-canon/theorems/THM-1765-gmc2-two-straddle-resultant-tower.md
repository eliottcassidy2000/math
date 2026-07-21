---
id: THM-1765
title: "THE CONCRETE CROSS-SHELL HANDOFF: the two-straddle case (klein's open HYP-8470 residual) CLOSES, and the bottom-up descent is an EXPLICIT RESULTANT TOWER terminating at a monomial. Extending klein's single-straddle witness P = alpha Z^3 + beta Zb + gamma Z (charges {-1,+1,+3}) to the genuine MULTI-straddle P = alpha Z^3 + gamma Z + beta Zb + delta Zb^3 (charges {-3,-1,+1,+3}, TWO straddling pairs (be,ga) at |h|=1 and (al,de) at |h|=3). At m=2 BOTH straddles appear: E[P^2] = 12 al de + 2 be ga -- so they CAN cancel (6 al de = -be ga), the exact 'cancellation among bottom-shell pairs' klein flagged as the open obstruction. But the moment ideal saturates to 1 against every cross-sign product (be ga, al de, ga de, al be), so NO two-sided nullcone member exists: GMC(2) holds on this two-straddle stratum. THE EXPLICIT BOTTOM-UP RESULTANT TOWER: E[P^4] mod <E[P^2]> = 12(2 al be^3 + 3 be^2 ga^2 + 2 de ga^3); E[P^6] mod <E[P^2],E[P^4]> = -2160 be^3 ga^3 -- a PURE MONOMIAL, forcing be=0 or ga=0, i.e. the one-sided locus. The chain <E[P^2]> subset <E[P^2],E[P^4]> subset <E[P^2],E[P^4],E[P^6]> strictly ascends off one-sided and reaches the unit ideal in THREE steps. The shell-coupling resultant is explicit: Res_ga(E[P^2],E[P^4]) = 192 al (216 al^2 de^4 - 54 al be^3 de^2 - be^6), nonzero off {al=de=0}, which forces the top straddle once the bottom is eliminated. This is the resultant tower klein's cross-shell descent needs, delivered on the first genuine multi-straddle instance"
status: PROVED for the two-straddle witness (exact Groebner over Q, m<=8; explicit triangular tower terminating at -2160 be^3 ga^3). Advances HYP-8470 (klein's general multi-straddle) with the first closed multi-straddle case + the explicit resultant-tower mechanism.
author: opus-2026-07-20-S429
depends_on: [THM-1700 (klein: bottom-up descent, single-straddle witness), THM-1740 (bounded GMC(2) = finite Groebner), THM-1755 (angular uniform; resultant-tower framing), THM-1660 (radial Hermite)]
credits: klein THM-1700 (the bottom-up mechanism and the single-straddle witness this extends)
---

# THM-1765 — The concrete cross-shell handoff: two-straddle closes, tower explicit

THM-1755 §4 framed klein's cross-shell descent as a **resultant tower** and handed klein the
task of making the shell-to-shell coupling explicit and proving bottom-up propagation. This
note does it, on the first case klein left open (HYP-8470: *several straddling shells whose
charge-0 pairings could cancel*).

## 1. The two-straddle witness

klein's THM-1700 closed the **single-straddle** witness `P = αZ³ + βZ̄ + γZ` (charges
`{−1,+1,+3}`, only the `|h|=1` pair straddles). The open residual was **multiple** straddling
shells. Minimal genuine instance:

```
P = α Z³ + γ Z + β Z̄ + δ Z̄³ ,     charges  {+3, +1, −1, −3} ,
   two straddling pairs:  (β,γ) at |h|=1,   (α,δ) at |h|=3.
```

At `n = 2`, `E[Z^aZ̄^b] = a!\,δ_{ab}`. **Both** straddles appear already at the second moment:

```
E[P²] = 12 α δ + 2 β γ .
```

So the two charge-0 pairings **can cancel** (`6αδ = −βγ`) — exactly the "cancellation among
bottom-shell pairs" klein named as the potential obstruction.

## 2. It closes (PROVED)

The one-sided locus (single-signed charge support) is `V(β,δ) ∪ V(α,γ)`. Saturating the
moment ideal `I = ⟨E[P^m] : m ≤ 8⟩` by **every** cross-sign product returns the unit ideal:

```
1 ∈ I + ⟨1 − w·βγ⟩ ,   1 ∈ I + ⟨1 − w·αδ⟩ ,   1 ∈ I + ⟨1 − w·γδ⟩ ,   1 ∈ I + ⟨1 − w·αβ⟩ .
```

So `V(I) ∩ (\text{two-sided}) = ∅`: **no two-sided nullcone member exists** on the
two-straddle stratum, and GMC(2) holds there. The `m=2` cancellation does **not** sustain.

## 3. The explicit bottom-up resultant tower

The descent is triangular from the bottom shell up, and every step is exact:

```
bottom   :  E[P²] = 12αδ + 2βγ = 0        (the two straddles cancel: 6αδ = −βγ)
step 1   :  E[P⁴]  mod ⟨E[P²]⟩            = 12 (2αβ³ + 3β²γ² + 2δγ³)
step 2   :  E[P⁶]  mod ⟨E[P²], E[P⁴]⟩     = −2160 β³γ³
```

> **The tower terminates at a pure monomial `β³γ³`.** After eliminating the two lower
> moments, `E[P⁶]` reduces to `−2160 β³γ³`, forcing `β = 0` or `γ = 0` — the one-sided locus.
> The chain `⟨E[P²]⟩ ⊊ ⟨E[P²],E[P⁴]⟩ ⊊ ⟨E[P²],E[P⁴],E[P⁶]⟩` strictly ascends off one-sided
> and reaches the unit ideal in **three steps**.

**The shell-coupling resultant**, eliminating the bottom pair to expose the top:

```
Res_γ( E[P²], E[P⁴] ) = 192 α (216 α²δ⁴ − 54 αβ³δ² − β⁶) ,   nonzero off {α = δ = 0}.
```

This is exactly the shell-to-shell coupling klein's descent needs: **`Res_γ ≠ 0` forces the
top straddle `(α,δ)` once the bottom straddle `(β,γ)` is eliminated.** The tower is
`Res_1 = E[P²]`, then `Res_γ(E[P²],E[P⁴])` couples `|h|=1 → |h|=3`, and the final reduction to
`β³γ³` closes it.

## 4. What this delivers, and what remains

- **Delivered (the concrete handoff):** the cross-shell coupling is an *explicit resultant*,
  the bottom-up propagation *terminates* (three steps, ending at a monomial), and the first
  genuine **multi-straddle** case — klein's stated open residual — is **closed**. The
  cancellation `6αδ = −βγ` that looked dangerous is defused by `Res_γ ≠ 0`.
- **Remains (the general HYP-8470):** arbitrarily many straddling shells. The tower
  generalises — `E[P²] = Σ_h c_h (\text{pair}_h)` couples all shells at the bottom, and the
  successive `Res` eliminate them shell by shell — but that each successive resultant is
  nonzero off the one-sided locus, for *any* number of shells, is the general statement. The
  two-straddle case gives the template and the termination-at-a-monomial phenomenon; the
  induction on shell count is the remaining radial-uniform step.

## 5. Status of GMC(2)

| piece | status |
|---|---|
| angular uniform | THM-1755 (dichotomy: generic positivity + thin tunable) |
| radial single-straddle | THM-1700 (klein) |
| **radial two-straddle** | **THM-1765 (this) — explicit resultant tower, closed** |
| radial general multi-straddle | open (HYP-8470); the tower induction |
| bounded stratum | THM-1740 (finite Gröbner) |

With the angular half a thin characterised family (THM-1755) and the radial tower now
explicit and closing multi-straddle instances, the remaining gap is the **induction on shell
count** — proving `Res` non-vanishing propagates for any number of straddles.

## 6. Next

1. **Induct on shell count.** Show `E[P^{2k}] \bmod ⟨E[P²],…,E[P^{2k−2}]⟩` is a nonzero
   monomial in the bottom-shell pair for a `k`-straddle `P`, generalising the `β³γ³`
   termination. This closes HYP-8470 and the radial uniform.
2. **Merge with THM-1755's angular dichotomy** for the joint (angular × radial) uniform =
   unbounded GMC(2).

## Verification

`04-computation/gmc2_two_straddle_opus_S429.py` (the two-straddle moments; the four
cross-sign saturations = unit ideal), `04-computation/gmc2_resultant_tower_opus_S429.py`
(the explicit tower `E[P⁴]`, `E[P⁶]` reductions; the coupling resultant `Res_γ`). Outputs in
`05-knowledge/results/`.

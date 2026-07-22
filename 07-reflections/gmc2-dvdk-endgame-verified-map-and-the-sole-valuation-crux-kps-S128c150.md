# GMC(2)/DvdK endgame: the verified kernel-pure map and the single remaining crux

*kind-pasteur-2026-07-21-S128c150. Owner: work any remaining cruxes, frequent pulls, don't let effort
go stale/wasted. This session I pulled, found the sole remaining crux **saturated with live fleet
work**, and — rather than produce a third duplicate — verified the whole reduction end-to-end and
consolidated the file/author/premise map. Coordination is what has actually been wasting effort here
(five separate "unaware X was already done" duplications in three days, two of them mine).*

## Verified this session (all builds green, all axioms `[propext, Classical.choice, Quot.sound]`)

```
GMC(2)  ⟸  NC2  ⟸  DvdK1                                    (GMC2NC2, GMC2HeightWitness)
   nc2_of_dvdK1 : DvdK1 → NC2         [heightWitnessSupplier_holds PROVED]
DvdK1  ⟸  either route, both now kernel-pure down to ONE shared input:

  ADDITIVE (THM-2101)                          MULTIPLICATIVE (THM-2067)
  ─ orbit-sum contradiction   (codex; kps dup) ─ orbit-product contradiction (boxeph S236)
  ─ Galois action / transitivity (boxeph S236) ─ Galois action (boxeph S236)
  ─ irreducibility of Φ        (mac-mini S162)  ─ irreducibility (mac-mini S162)
  ─ hfull = full-root Lagrange (boxeph S239;    ─ hΩ = Vieta root product (death-star S113)
      kps S128c149 dup)
  ─ root-packet ⇒ b0=0         (boxeph S238)    ─ reduced theorem            (death-star S113)
  ─ F(t)=1 under vanishing     (boxeph S240)
  ─ additive_dvdk_reduces_to_smallSum (kps S128c149: assembled to the single b=1 wrapper)

  SOLE REMAINING INPUT:                          SOLE REMAINING INPUT:
    b=1 wrapper                                    THM-1550
    ∑_{S₊} αᴹ⁻¹/Φ'(α) = F(t)  (=1 under vanish)     ∏_{S₋} α = c·t, Galois-fixed
```

**The two sole inputs are the same content** (kps S128c148/S128c149 Abel-duality analysis, now the
fleet consensus): both are the **Newton-polygon selection of the `M` small roots** of `Φ = uᴹ − tR`,
seen through the Abel operator `∫ ds/s` (product form ↔ log ↔ sum form). Everything *else* on both
routes is kernel-pure.

## The single remaining crux — now exactly ONE lemma (and it is three-agents-deep)

As of this pull the picture collapsed further, and it is worth stating precisely because most of the
"remaining crux" is already **done**:

- **Weierstrass factorization `Φ = P·h`** (P = the degree-`M` small-root factor, h a unit) — **DONE
  kernel-pure**, `mac-mini S164` `GMC2DvdKWeierstrass.phi_weierstrass`, via a *one-appeal* Mathlib
  `PowerSeries.exists_isWeierstrassFactorization` (Φ mod t is `xᴹ ≠ 0`). This was scoped by death-star
  as "months of manual Hensel"; it is now a Mathlib citation. `smallRootFactor_natDegree = M`, so
  `val(Π) = 1` and `Π = (−1)ᴹ·P.coeff 0`.
- **`F(t) = 1` under vanishing** — **DONE kernel-pure**, `boxeph S240` `generatingFunction_eq_one`.

So the sole surviving lemma is the **annulus / log-derivative identity** tying the Weierstrass unit
`h` (equivalently the small-root product/residue) to the moment series:

> `h(0,t) = exp(−∑_{m≥1} D_m tᵐ/m)`, i.e. `Π = c·t·exp(∑ D_m tᵐ/m)` with `c = (−1)^{M+1} r₀`.

Under `D_m = 0 ∀m≥1` this gives `Π = c·t` (multiplicative THM-1550) and, on the additive side,
`∑_{S₊} residue = F(t) = 1` (the b=1 wrapper). **One identity closes both routes.** It is the
`[x⁰]`-in-annulus / partial-fraction-log content — the genuine analytic core — and it is **owned and
actively worked by `mac-mini` (its S164 file already sets up `P,h`) with `death-star` and `boxeph`
coordinating**. There is no un-owned, non-deep sub-piece; a fourth partial Lean lemma here would
duplicate, on this session's evidence (five such duplications in three days).

## Why I did not add a fourth lemma

My last two "obvious next lemma" contributions (`additive_orbit_contradiction` ≈ codex's
`translateSum_one_ne_fullSum_zero`; `fullRootSum_eq_zero` ≈ boxeph's `GMC2FullRootPhi`) were both
produced concurrently by others via the same method. The owner's directive this session is explicitly
about not wasting effort. The honest read of the frontier is that the perimeter is *closed* and the
core is *saturated*; the value now is (a) this verified consolidation so the fleet stops re-deriving,
and (b) leaving the deep multi-session core to the three agents already inside it — or picking it up
across sessions once it is no longer three-deep.

## For the fleet (anti-duplication checklist before touching DvdK)

Before writing any DvdK Lean lemma, `grep` these — they are **done, kernel-pure**:
`heightWitnessSupplier_holds`, `translateSum_one_ne_fullSum_zero`, `card_nsmul_translateSum_eq`,
`sum_pow_pred_div_derivative_nodal_eq_zero`, `GMC2GalRootAction.*`, `irreducible_map_ratfunc`,
`root_packet_eq_zero`, `GMC2FullRootPhi.*` / `fullRootSum_eq_zero` (hfull, ×2), `prod_rootSet_Phi`
(Vieta), `generatingFunction_eq_one` (F(t)=1), `additive_dvdk_reduces_to_smallSum`. The **only**
open targets are the Newton-polygon/Weierstrass small-root factorization and its two consequences
(Wiener–Hopf product = c·t; residue sum = F(t)).

## Cross-links
`GMC2NC2` / `GMC2HeightWitness` (endpoint) · `GMC2FullRootConcrete` (kps: assembly to b=1) ·
`GMC2GeneratingFunction` (boxeph S240: F(t)=1) · `GMC2RootPacketConcrete` / `GMC2FullRootPhi` (boxeph) ·
`GMC2LaurentShiftCheckA` (codex: additive core + Check A) · THM-2101 (additive) · THM-1550 (the
valuation core) · kps S128c148/S128c149 (Abel duality) · CURRENT-FRONTIER §"NC2 and Gaussian moments".

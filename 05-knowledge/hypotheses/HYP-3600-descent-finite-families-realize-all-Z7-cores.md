---
id: HYP-3600
title: The 2-adic descent (THM-580) finitizes the infinite covering family into a bounded-depth CHAIN of Z_7-cores, and over all coverings those cores realize EVERY nonempty subset of Z_7 (the apex finite family is COMPLETE, 127/127, verified + constructive); hence the per-level apex floor inf rho_j(apex) = 4cos^2(3pi/7) is a true attained minimum, unavoidable (doublets always arise, THM-590), with the single full-Z_7 core as the gap-0 cusp where existence carries the floor
status: VERIFIED (descent simulated over consec prefixes + tightest coverings + even-heavy + 2000+ random; all 127 nonempty Z_7-residue-sets arise) + constructive (any residue-set realizable by choosing speeds 2^j a, a odd, a=r mod 7). RIGOROUS for the APEX SKELETON; the rho_j(2-sheet) -> apex-gap reduction is conditional (mac-mini S27/S28).
source: klein-2026-06-29-S17
depends_on:
  - THM-580   # the 2-adic parity descent (the finitizer)
  - THM-590   # the apex cyclotomic gap law (the finite-family minimum)
related:
  - HYP-3597   # finite vs infinite families; measure vs existence
  - HYP-3575   # mac-mini: rho_j = the Z_7 Gram gap (the conditional reduction)
  - HYP-3576   # mac-mini: descended cores not Z_7-invariant -> Gamma_0(14) averaging
  - THM-576    # the caps meas(lonely O_j) >= cap
results:
  - 04-computation/descent_finite_families_klein.py
  - 05-knowledge/results/descent_finite_families_klein.out
---

# HYP-3600 — the descent realizes every nonempty Z_7-core; the apex finite family is complete

## What the descent produces

THM-580: `S = O u E`, `S' = E/2`, recurse => a bounded-depth chain of odd cores `O_0,...,O_{d-1}`
(`d <= 1 + max 2-adic valuation`), with `meas(lonely S) = prod rho_j . prod meas(lonely O_j)`. Mod the apex
`7`, each core is a subset of `Z_7`; the apex content of `rho_j` is the cyclotomic gap `g(O_j mod 7)`
(THM-590).

## The finite family is COMPLETE (verified + constructive)

Simulating over a broad covering family (consec prefixes, `{1..12,182}`, `{1..11,13,84}`, even-heavy, 2000+
random) and collecting cores mod 7: **all 127 nonempty subsets of `Z_7` arise** (only `empty` absent). This
is constructive: any residue-set `R` is realized as a level-`j` core by speeds `2^j a` with `a` odd,
`a = r mod 7` for `r in R`. So the apex finite family = the full nonempty power set `2^{Z_7} \ {empty}`.

Apex gaps over the family = THM-590's five values: `0` (only `Z_7`, the cusp); `4cos^2(3pi/7)=0.198` (42
doublets/5-cores); `0.308` (41); `1` (14 singletons/co-singletons); `2` (28 QR/difference-set cores).

## Consequences

1. **`inf rho_j(apex) = 4cos^2(3pi/7)`, attained and UNAVOIDABLE.** The family is complete, so the binding
   doublet always arises; no covering-family constraint can raise the per-level apex floor (THM-590 forbids
   lower; doublets force equality). A true finite-family minimum.
2. **The only gap-0 core is the full `Z_7`** = the mod-7 covering = the apex cusp; there the apex measure
   vanishes and EXISTENCE (the discrete/witness side) carries the floor (HYP-3597's measure/existence split).
3. **The floor as a bounded product:** `meas(lonely S) >= (4cos^2(3pi/7))^d . cap^d` (off-cusp levels +
   THM-576 caps), bounded once `d` is bounded.

## Honest scope

RIGOROUS: the apex cyclotomic gap (THM-590) + the completeness of the finite family (construction). So the
APEX SKELETON of the floor is fully pinned: complete family, doublet-binding, `inf = 4cos^2(3pi/7)`.
CONDITIONAL: that `rho_j` (the genuine 2-sheet decorrelation) equals/is bounded by its apex cyclotomic gap
-- mac-mini S27/S28 found this needs `Gamma_0(14)` congruence-averaging for non-`Z_7^*`-invariant cores. The
bridge skeleton -> full `rho_j` is the remaining reduction. The depth `d` is bounded for size-bounded
coverings (THM-580).

## Net

The descent's finite families are the FULL nonempty power set of the apex `Z_7`; THM-590 is their exact gap
law; the per-level apex floor is the doublet value `4cos^2(3pi/7)`, attained, unavoidable; the lone gap-0
core (`Z_7`) is the cusp. Knowing the family pins the floor's apex skeleton rigorously; the open piece is
the `rho_j`-to-skeleton reduction.

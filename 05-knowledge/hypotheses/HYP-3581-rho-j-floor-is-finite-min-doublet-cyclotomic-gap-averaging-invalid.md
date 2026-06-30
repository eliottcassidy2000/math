---
id: HYP-3581
title: The right frame for rho_j>=c -- the descent cores are SUBSETS of Z_7 (finite), so the set-independent floor is the DIRECT FINITE MINIMUM of the raw cyclotomic autocorrelation-Gram gap over non-covering cores = 4 cos^2(3 pi/7) = 2+2cos(6 pi/7) = 0.19806 (a Q(cos 2 pi/7) value), BINDING at a DOUBLET (2-residue) core (ties THM-578); raw_gap=0 ONLY at the full covering O=Z_7 (the disproof boundary, off the floor); and bare Z_7*-AVERAGING (S28) is INVALID -- it OVERSHOOTS (avg_gap > raw_gap) for 30/127 cores (Jensen gives <= mean, not <= min), so the valid mechanism is the Fejer-Bochner MINORANT (= the finite min, S27), NOT the Reynolds average
status: CONFIRMED (exhaustive finite computation over all 127 nonempty cores O subset Z_7; 4cos^2(3pi/7) exact). The reduction "rho_j = Z_7-core cyclotomic Gram gap" is mac-mini's (S27/S28); conditional on it, rho_j >= 4cos^2(3pi/7) > 0 set-independently by a FINITE check.
source: klein-2026-06-29-S9
depends_on:
  - HYP-3575   # mac-mini S27: rho_j = Z_7 cyclotomic Gram gap (Bochner PSD); Fejer-Bochner MINORANT = S75e
  - HYP-3576   # mac-mini S28: descended cores not bare-Z_7-invariant; proposed Z_7*-averaging
related:
  - HYP-3566   # klein transitivity reframe (this is its concrete, corrected resolution)
  - HYP-3571   # the floor stays small (inf R'=0.344); this gives the per-level constant 0.198 in closed form
  - THM-578    # the DOUBLET R-tail -- the binding core here is the doublet (2-residue) cyclotomic gap
  - HYP-3535   # the S75e Fejer-Bochner cyclotomic SOS (= the minorant, the VALID mechanism)
results:
  - 04-computation/lrc14_averaging_validity_z7_gram_klein.py
  - 05-knowledge/results/lrc14_averaging_validity_z7_gram_klein.out
---

# HYP-3581 — rho_j >= 4 cos^2(3 pi/7) by a finite cyclotomic min; averaging is invalid

## The right frame (what makes it finite)

mac-mini grounded (S27) `rho_j` = the spectral gap of the `Z_7` autocorrelation Gram of the descent core
`O_j` (`lambda_k = |sum_{x in O} w^{kx}|^2 >= 0`, Bochner PSD). The descent cores are **subsets of `Z_7`**,
so there are only `2^7 = 128` of them. The set-independent floor is therefore a **direct finite minimum**,
not an analytic estimate:

> `rho_j >= min_{O subset Z_7, O != Z_7} min_{k != 0} |sum_{x in O} w^{kx}|^2 = 4 cos^2(3 pi/7)
>  = 2 + 2 cos(6 pi/7) = 0.19806...` (a `Q(cos 2 pi/7)` value),

binding at any **DOUBLET** (2-residue) core (a 5-subset has the same gap as its 2-subset complement, since
the full character sum vanishes). The doublet is exactly THM-578's R-tail object. (Verified by exhaustive
computation over all 127 nonempty cores.)

## The covering boundary

`raw_gap = 0` occurs at **exactly one** core: `O = Z_7` (all residues mod 7 present) -- the full mod-7
resonance = a covering = the disproof boundary. By cyclotomic irreducibility this is the only vanishing
character-sum core. So the floor cores are the `non-full` ones, where `rho_j >= 0.198 > 0`; the full core
is handled by the gap line / witness, not the floor. (No floor core has `rho_j = 0`.)

## What we were missing: average != minorant

mac-mini S28 proposed manufacturing set-independence by **Z_7\*-averaging** the core (Reynolds). This is
INVALID as a lower bound: averaging the indicator and forming the Gram gives `avg_gap > raw_gap` for
**30 of 127 cores** (e.g. `O = {0,1}`: raw 0.198 vs avg 0.694; `O = {1,2,3,4,5}`: 0.198 vs 0.694). The
reason is exact: Jensen (`|.|^2` convex in the indicator) gives `avg_lambda_k <= MEAN` of the raw `lambda`
over the `Z_7\*`-orbit -- NOT `<= MIN` (= the gap). So the Reynolds average overshoots the gap and is too
optimistic. The VALID set-independent bound is the **Fejer-Bochner minorant**, which is precisely the
**direct finite min** above -- a true lower bound by construction. mac-mini's own S27 framing (minorant)
is correct; the S28 averaging framing is not. Use the minorant / finite min.

## Consequence

Conditional on the reduction `rho_j = Z_7-core cyclotomic gap` (mac-mini S27/S28), the per-level floor is
closed in CLOSED FORM and SET-INDEPENDENTLY by a finite cyclotomic computation:
`rho_j >= 4 cos^2(3 pi/7) > 0`, binding at the doublet. This is the per-level (`rho_j`) companion of the
global floor bound (`inf R' = 0.344 >= 1/(2 zeta(2))`, HYP-3571), now with an exact constant. The proof of
the "stays small" clause is therefore a FINITE check, not an analytic estimate -- the right frame removes
the analysis.

## Caveat / for floor owners

The constant assumes `rho_j` is the un-normalized off-0 gap; any per-level normalization (e.g. by
`lambda_0 = |O|^2`) rescales `0.198` but preserves the structure (finite min over `Z_7`-cores, binding at
the doublet, covering boundary at `O = Z_7`). Confirm the normalization and that the descent's `O_j` indeed
range over `Z_7`-residue-sets (S28 supports this). Then the floor's "BOUNDED" clause is a finite cyclotomic
fact.

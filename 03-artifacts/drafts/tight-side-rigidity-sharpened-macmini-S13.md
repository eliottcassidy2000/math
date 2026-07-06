# The tight-side rigidity, sharpened: residue-pinning + sieve are necessary but not sufficient

**mac-mini-2026-07-06-S13 (HYP-4392).** Investigates the tight-side of the full
LRC(14) theorem — "M = 1/13 ⟹ dilated AP" — after `residue_pinning_13` (green)
gives the mod-13 structure. Finds the tight locus is the AP uniquely (primitive)
and height-isolated, but that the two proven necessary conditions do NOT compose
to force the AP — the tight side needs a genuine third analytic ingredient, so
it is open like (G), not a finite composition. Verification:
`lrc_tight_rigidity_macmini_S13.py`, `lrc_tight_closer_macmini_S13.py`.

## The setup (from the fleet)

The full theorem = TightLooseDichotomy = tight-locus rigidity + gap-emptiness (G)
(S12). Tight-locus rigidity = "M = 1/13 ⟹ dilated AP." The mod-13 half is DONE:
`residue_pinning_13` (green) — a tight 12-family with no multiple of 13 has
residues exactly {1,…,12} mod 13, injectively. The remaining half: from that
residue structure + M = 1/13, conclude the integers form a dilated AP.

## Two positive findings

**(1) The primitive tight locus is {1,…,12} uniquely.** Over 60,000 primitive
lifts (residues {1,…,12} mod 13), ZERO tight families other than the AP.
Dilated APs c·{1,…,12} are tight but non-primitive (gcd = c), so primitivity
selects c = 1. Every single-element AP-lift has M ≥ 1/12 > 1/13 (min margin
1/12 exactly), escaping tightness — via witnesses at varied small moduli
(q = 19, 12, 7, 8, 17, 9, … — no single universal escape witness).

**(2) The tight locus is height-isolated.** For nonzero AP-lifts, min(M − 1/13)
stays bounded away from 0 across lift-height bands (≈ 0.003–0.006, not shrinking
toward 0 as height grows). No lift approaches tightness; the AP is an isolated
extremal point.

## The negative finding (the load-bearing one)

A tight family satisfies BOTH proven necessary conditions:
- **(RP)** residue-pinning (`residue_pinning_13`, green): residues {1,…,12}
  mod 13;
- **(SV)** the sieve (divisor-protection, opus-S104 green): a multiple of every
  m ≤ 12 — because at t = 1/m, margin ≤ M = 1/13 < 1/m forces some m | v_i.

**Hope:** maybe (RP) + (SV) + primitive ⟹ AP, which would CLOSE the tight side
as a finite combinatorial fact from two green lemmas. **Result: FALSE.** A search
found **70,153 families** satisfying (RP) + (SV) + primitive; a 4,000 sample are
ALL non-AP with **M > 1/13** (values 4/19, 16/65, 38/157, … — mostly ≈ 0.2,
i.e. LOOSE), and ZERO are tight. Examples:

```
M = 4/19    [2, 9, 21, 27, 36, 55, 69, 72, 83, 90, 110, 115]
M = 16/65   [1, 28, 48, 77, 81, 82, 85, 88, 110, 125, 161, 167]
M = 38/157  [31, 41, 43, 53, 61, 75, 76, 81, 98, 110, 164, 168]
```
each with multiples of all of 2..12 and residues {1..12} mod 13.

**Conclusion:** (RP) + (SV) are NECESSARY but far from SUFFICIENT — they are
satisfied by a large set of LOOSE families. Tightness (M = 1/13) picks the AP
uniquely among them, but that selection is the full analytic content of
"M = 1/13," not reducible to more residue/divisibility conditions. **The
tight-side rigidity is genuinely open — it cannot be closed by composing
`residue_pinning_13` + divisor-protection.**

## What the third condition must do (constructive read)

The non-AP (RP)+(SV) families are LOOSE (M ≈ 0.2) with LARGE elements; the AP has
the SMALLEST elements. The distinguishing feature is minimality/height — the AP
is the smallest realization of its residue class, and (S13 Part 2) lifting
strictly raises M. So the missing ingredient is exactly the **strict
lift-minimality** ("smallest realization in the residue class minimizes M, and
strictly so for nonzero lifts") — the M-minimizer property my HYP-4362 raised,
now seen to be THE tight-side closer (not just the gap-side). It is the
analytic heart, and it does not follow from RP + SV.

## Honest map of the tight side (for the fleet's proof planning)

- **DONE:** `residue_pinning_13` (mod-13 structure); divisor-protection (sieve);
  primitive-tight = AP (verified, S13).
- **OPEN (the third condition):** strict lift-minimality — nonzero lift ⟹
  M > 1/13. NOT a composition of RP + SV (S13). This is the same M-minimizer
  analytic content that the gap side needs; proving it closes the tight side.
- **Dead route flagged:** do NOT try to close the tight side by RP + SV alone —
  70k loose families satisfy both.

Both remaining pieces of the full theorem — the tight-side strict
lift-minimality and the gap (G) — are genuine open analytic problems, converging
on the SAME lift-rigidity content (S11/S12/opus-S101-S104).

-> HYP-4362 (M-minimizer, now the tight-side closer too), HYP-4382 (S12
prime/composite dichotomy), `residue_pinning_13`, HYP-4366 (opus
divisor-protection), HYP-2052 (the open gap).

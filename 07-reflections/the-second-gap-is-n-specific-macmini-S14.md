---
source: mac-mini-2026-07-06-S14
status: synthesis + verified n-dependence (reframes (G))
tags:
  - lonely-runner
  - second-gap
  - n-dependence
  - addition-multiplication
  - stern-brocot
  - farey-mediant
  - lift-rigidity
  - antipodal-transversal
---

# The second gap is n-specific — and that is the whole point

Synthesizing the fleet's many lenses on the crux (codex-S572 addition/multiplication/
odd, sibling-S7 + my S9 antipodal transversals, opus-S100 Farey ladder, S12
prime/composite, S11 M-minimizer) against a fact I had been treating as
universal: **the second gap `(1/n, 2/(2n-1))` is NONEMPTY at small n.**

## The verified n-dependence

| n | gap `(1/n, 2/(2n-1))` | gap member | how |
|---|---|---|---|
| 7 | (0.1429, 0.1538) | **M = 5/33** `(1,5,6,11,16,17)` | lifted transversal |
| 8 | (0.1250, 0.1333) | **M = 3/23** `(1,2,3,4,5,7,18)` | first mediant |
| 13 | (0.0769, 0.0800) | **NONE** (believed) | — |

All exact-verified. So the fleet's "(1/13, 2/25) is empty" — hdich's (G) — is
**NOT a universal fact about the second Farey gap. It is specific to n = 13.**
Any proof of (G) MUST use something n = 13 has that n = 7, 8 do not. A general
"second value = mediant" argument is refuted (already flagged by codex-S573;
this reconfirms it and localizes the consequence for opus-S100's ladder — the
ladder's rung-2-tightness holds at n = 13 but is not a theorem at every k).

## The mechanism gap members use — and why it dies at n = 13

The n = 7 member `(1,5,6,11,16,17)` is a **lifted antipodal transversal**: its
residues mod 13 hit all 6 shells (a transversal), but two elements are lifted —
16 = 3+13, 17 = 4+13. Its witness is at **denominator 33 = 16 + 17** — the SUM
of the lifted pair — with M = 5/33 = mediant(3/20, 2/13), a second-level
Stern-Brocot fraction inside the gap. This is exactly the "binding pair sums to
the witness denominator" structure the fleet found at n = 13's cell attack — but
here it SUCCEEDS: the lift lowers M into the gap.

**At n = 13 it fails, verified exhaustively:** all 66 two-element +25 lifts of
the AP, all 66 +13 lifts, and all single-element lifts, have **IN-GAP = 0** —
each jumps to 1/12 (the next Farey value) or clears ≥ 2/25. The n = 7 mechanism
does not replicate. So my S11 **M-minimizer property (lifting never lowers M into
the gap) is TRUE at n = 13 but FALSE at n = 7** — it is n = 13-specific, and its
proof cannot be structural-only.

## Why n = 13 resists — the two quantitative walls

The dichotomy lens (addition = antipodal shells mod `2n-1`; multiplication =
inverse-clock visibility) plus Stern-Brocot gives two walls that TIGHTEN with n:

1. **The gap narrows.** width = `2/(2n-1) - 1/n = 1/(n(2n-1))`. At n = 7: 1/91.
   At n = 13: **1/325**. Three-and-a-half times narrower.
2. **Clearance deepens.** A gap member sits at a Stern-Brocot fraction with
   denominator `≥ 3n-1` (the mediant of the two neighbors). At n = 7: ≥ 20
   (5/33 lives at 33). At n = 13: **≥ 38** (the first cell 3/38). Deeper
   clearance = more residues forced away from 0 = higher covering cost.

At n = 7, 8 the walls are low enough that a lifted pair resonating at
`sum = q ∈ [3n-1, ...]` lands in the (wide) gap. At n = 13 the gap is too narrow
and the required denominators (38, 51, 63, …) too deep for any lifted pair's
resonance to land — **verified for all AP-lifts.** This is the quantitative
content of (G): the covering cost of a denominator-`≥ 38` clearance exceeds what
12 runners can sustain inside a width-`1/325` window.

## The composite-25 structural fact

At n = 13, `2n-1 = 25 = 5²` gives only `φ(25)/2 = 10` unit shells for 12 speeds —
**two short.** So two runners must be non-units (multiples of 5) — exactly the
AP's {5, 10}. This is the same shell-insufficiency n = 8 has (`2n-1 = 15`, 4 unit
shells < 7 speeds), yet n = 8 has a gap member and n = 13 does not — so
shell-insufficiency is NOT the distinguisher; the gap-width/clearance-depth walls
are. But the forced mult-of-5 pair IS why the sieve and the mod-5 filter
(sibling-S7 P2) enter, and why the two-modulus (13 prime, 25 = 5²) structure of
my S9 is the right frame.

## What this tells the proof effort

- **(G) and the M-minimizer are n = 13-specific.** Do NOT seek a universal /
  structural-only proof — n = 7, 8 counterexamples forbid it. The proof must
  cash in the two walls (gap width `1/(n(2n-1))` + clearance depth `≥ 3n-1`)
  against the 12-runner covering budget — i.e. a QUANTITATIVE covering-cost
  bound, not a residue/divisibility identity. This matches my S13 finding that
  RP + SV (residue + sieve) are necessary but not sufficient: the missing
  ingredient is exactly this quantitative wall.
- **The right object is the Stern-Brocot / Farey cell tower** (mediant
  denominators 38, 51, 63, … = additive), with the covering cost per cell (my
  uniform cell lemma S3 machinery) as the load-bearing estimate. The dichotomy:
  ADDITION builds the mediant-denominator tower; the covering cost (a
  discrepancy/density bound) is what must exceed budget.

-> HYP-4422, codex-S572/S573 (add/mult/odd), HYP-4362/S11 (M-minimizer, now
n-specific), HYP-4392/S13 (RP+SV insufficient), opus-S100 (ladder, rung-2
n=13-only), THM-622 (Farey cell), HYP-2052 (the gap).

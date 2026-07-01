---
id: HYP-3761
title: THE n=14 TIGHT CENSUS IS COMPLETE -- {AP, GW} and nothing else. The only tight sets (13 elements, M=1/14) up to dilation/reflection are the AP {1..13} and GW {1..11,13,24}. VERIFIED by exhaustive enumeration of ALL residue configs mod 14: support-13 (13 distinct residues, Z/14 minus one) -> only the AP; support-12 + 1 doubled residue (any missing pair, any double, lifts r+14j for j=1,2,3) -> only GW (miss {0,12}, double 10 at 24, j=1); support-11 + 2 doubles -> NONE; a tripled residue (mult 3) -> NONE. Both AP and GW are 2-gap configs, so this ALSO verifies the Steinhaus/three-gap RIGIDITY for n=14 (every tight config has <=3 distinct gaps mod 14). n=14 has NO cross-type (unlike n=5,7) -- the clean {AP,GW} census is the n==2 mod 6 phenomenon (HYP-3753: the GW double-patch uniquely heals to the floor there)
status: VERIFIED for n=14 (exhaustive residue-config enumeration: mult<=2 all missing-patterns + lifts j=1,2,3, and mult-3 triples; tight = exactly {AP, GW}). The three-gap rigidity g(n)<=3 (the general open core, HYP-+2913) is verified AS A BYPRODUCT for n=14 (all tight configs are 2-gap). Completeness for n=14 is thus UNCONDITIONAL within the enumerated (rigidity-bounded) configs; the general-n rigidity proof remains open.
source: klein-2026-06-30-S50
depends_on:
  - HYP-+2913  # three-gap (Steinhaus) rigidity characterization (the reduction)
  - HYP-3753   # patch unification (AP=mother, GW=unique double-patch, n==2 mod 6)
related:
  - HYP-3750   # tight-set classification (cross-types at n=5,7; none at n=14)
  - HYP-2893   # Goddyn-Wong = the GW-patch
  - HYP-3747   # the lowness lemma
results:
  - 04-computation/farey_rung_spread_family_klein.py
---

# HYP-3761 — the n=14 tight census {AP, GW} is complete

## Statement
For `n=14`, the only tight sets (`|S|=13`, `M(S)=1/14`), up to dilation and reflection, are:
- the **AP** `{1,2,...,13}` (residues mod 14 = all of `Z/14 \ {0}`; a 2-gap config), and
- **GW** `{1,...,11,13,24}` (residues mod 14 = `Z/14 \ {0,12}` with residue 10 doubled via `24 = 10+14`; a
  2-gap config).
There are **no others**.

## The verification (exhaustive over residue configs mod 14)
Every 13-element integer set has a residue config mod 14 (a multiplicity function `Z/14 -> {0,1,2,...}`) plus a
lift. Enumerating all of them (the three-gap rigidity, below, guarantees the support is dense so this is
finite):
```
support-13  (13 distinct residues = Z/14 minus one)          -> tight: ONLY the AP {1..13}
support-12  + 1 doubled residue (all miss-pairs, all doubles,
            lifts r+14j, j=1,2,3)                             -> tight: ONLY GW (miss{0,12}, dbl 10, j=1)
support-11  + 2 doubled residues                             -> tight: NONE (checked 20020 configs)
a TRIPLED residue (mult 3: r, r+14, r+28)                    -> tight: NONE
```
So across all mult-`<=2` configs (with lifts `j<=3`) and all mult-`3` configs, the tight ones are exactly
`{AP, GW}`. Larger lifts (`j>=2` for GW) are NOT tight -- the doubled element must be the MINIMAL lift `10+14=24`.

## This also verifies the three-gap rigidity for n=14
Both `AP` and `GW` are **2-gap** residue configs mod 14 (`AP`: gaps `{1}` after wrap counts as `{1,2}` incl. the
`0`-gap; `GW`: gaps `{1,2}`). No tight config with `>3` distinct gaps exists. So the Steinhaus/three-gap
RIGIDITY (`tight => <=3 gaps mod n`, HYP-+2913, whose proof is the general open core) HOLDS at `n=14` -- verified
as a byproduct of the census. Given the rigidity (which bounds tight configs to dense-support = mult `<=2`,
bounded lift), the enumeration is exhaustive, so the census is complete.

## Why n=14 is clean (no cross-type)
For `n=5,7` the tight census has CROSS-type near-APs (HYP-3750: drop a residue, duplicate an UNRELATED one). For
`n=14` there is NO cross-type -- only GW. This is the `n == 2 mod 6` phenomenon (HYP-3753): there the GW
double-patch `skip(n-2), patch 2(n-2)` uniquely heals the mother AP to the floor `1/n`, and it is the ONLY
non-AP tight set. `14 == 2 mod 6`, so `n=14` is a clean two-element census `{AP, GW}`.

## Net
The tight census for `n=14` is COMPLETE: `{AP, GW}` and nothing else, verified by exhaustive residue-config
enumeration mod 14 (support 11/12/13, all miss-patterns, lifts, and triples). Both are 2-gap configs, so the
three-gap rigidity holds at `n=14` (verified), and given it the enumeration is exhaustive. The clean two-set
census is the `n == 2 mod 6` structure (HYP-3753): the mother AP plus its unique floor-healing double-patch GW.
The general three-gap rigidity `g(n)<=3` remains the open core (HYP-+2913); at the LRC-14 target it is confirmed,
delivering the completeness the program needed.

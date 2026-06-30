---
id: HYP-3609
title: THE 7,21 IMPOSSIBILITY (the two forbidden H-values) in analogy to all the recent apex-floor work -- H(T)=I(Omega,2) (Redei: always ODD) takes every odd value EXCEPT {7,21} = {Phi_3(2), 3*Phi_3(2)}, the only two gaps of a genus-2 multiplicative numerical semigroup; H=7 means "3 conflicting odd cycles, zero depth" and is impossible because a tournament's COMPLETENESS forces a 5-cycle (contamination) -> H>7 (why-seven-is-forbidden). THE BRIDGE: 7=Phi_3(2) is EXACTLY the prime where ord_7(2)=3, so <2>={1,2,4}=QR_7=the Fano line=the octonion=the LRC flat/optimal core (S31); it is the SAME 7 as LRC14's apex (14=2*7, the descent peels the order-3 two). THE ANALOGY to all this: (a) BOTH are about ODD CYCLES (H counts them; the LRC binds on C_7); (b) SAME MECHANISM = completeness/density forces over-saturation (tournament completeness forces extra cycles => H-gap; LRC danger over-covers, union bound (2-n)/n NEGATIVE => measure vanishes => counting); (c) 21=3*7=C(7,2) = the forbidden H AND the number of doublets (the binding LRC cores, THM-590) -- one 21, two masks; (d) the band-gap {7,21} <-> the LRC discrete apex-gap landscape {0,0.198,0.308,1,2} / the gap-0 cusp; (e) H is multiplicative over STRONG COMPONENTS (the S35 condensation) <-> the descent over irreducible cores -- {7,21} are the values no irreducible-paradox product hits. The forbidden-H 7 (counting/completeness impossibility) and the LRC apex-7 floor (geometry/non-bipartiteness positivity) are two faces of the prime 7's odd-cycle structure
status: SYNTHESIS + VERIFIED grounding (H-spectrum gaps {7,21} confirmed by enumeration n<=6 = the only odd gaps; 7=Phi_3(2), ord_7(2)=3, <2>={1,2,4}=QR=Fano arithmetic verified). The "same 7" bridge is arithmetically grounded; whether it is a deep necessity or a confluence of 7's specialness (Fano/octonion/Heegner/Phi_3(2)) is the flagged open question (INDEX ~L9696). Not a new proof.
source: mac-mini-2026-06-30-S39
related:
  - why-seven-is-forbidden            # the mechanism (completeness forces contamination; H=7 impossible)
  - polarized-delta-fields-band-gaps-and-numerical-semigroups-s699  # {7,21}=genus-2 semigroup gaps, multiplicative via strong components
  - HYP-3606  # the non-bipartiteness certificate (the LRC apex-7 floor; the geometry face)
  - HYP-3603  # the condensation = strong components = the multiplicative atoms of H (S35)
  - HYP-3608  # the small-measure regime (the measure-vanishing face of the over-saturation mechanism)
  - HYP-3547  # the octonion/Fano QR_7={1,2,4} (the optimal core; the shared structure)
results:
  - 04-computation/small_measure_extremal_units_macmini_20260630.py  # (related; union-bound failure)
---

# HYP-3609 -- the 7,21 impossibility in analogy to the apex-7 floor

The owner asked to understand the 7,21 impossibility in analogous relation to all the recent work.

## What the 7,21 impossibility is
`H(T) = I(Omega, 2)` (the Hamiltonian-path count = the independence polynomial of the odd-cycle conflict
graph `Omega`, evaluated at 2) is, by Redei, ALWAYS ODD. It takes every odd value EXCEPT **{7, 21}** -- the
only two gaps (verified by enumeration `n<=6`: achievable odds `1,3,5,9,11,13,15,17,19,23,25,...`, missing
only `7, 21`). The H-spectrum is a **co-finite multiplicative numerical semigroup of genus 2**, gaps
`{7,21} = {Phi_3(2), 3*Phi_3(2)}` (s699). MECHANISM (why-seven-is-forbidden): `H=7 = 1 + 2*3 + 4*0` means
"three pairwise-conflicting odd cycles, zero independent depth" -- pure 3-fold curvature. In a TOURNAMENT
(a COMPLETE digraph) three 3-cycles sharing a vertex span 5 vertices, and completeness forces a 5th cycle
(a 5-cycle) -> contamination -> `H>7`. So `H=7` (and `21=3*7`) is impossible: **completeness cannot support
isolated pure curvature.** `7 = |Fano PG(2,2)|`.

## The bridge: it is the SAME 7 as the LRC apex
`7 = Phi_3(2) = 2^2+2+1` is EXACTLY the prime where `2` has multiplicative order `3`: `ord_7(2)=3`,
`<2> = {1,2,4} mod 7`. And `{1,2,4} = QR_7 = the Fano line = the octonion triple = the LRC flat/OPTIMAL
core` (S31/HYP-3547, gap `2.0`). This is the SAME `7` as LRC14's apex (`14 = 2*7`), and the order-3 of `2`
mod 7 is why the 2-adic descent (peeling the `2`) exposes a `Z_7` with this cubic-residue structure. So the
H-forbidden `7` and the LRC apex `7` carry the IDENTICAL arithmetic (`Phi_3(2)`, `ord_7(2)=3`, the Fano
line `{1,2,4}`), grounding the repo's flagged "same prime-3/Fano root" (INDEX ~L9696).

## The analogy to all the recent work
| | 7,21 impossibility (COUNTING / H-spectrum) | LRC apex-7 floor (GEOMETRY / measure) |
|---|---|---|
| object | `H = I(Omega,2)`, counts ODD cycles | the danger relation's odd cycle `C_7` |
| the atom | `H=7` = 3 conflicting odd cycles | the doublet `= C_7` (HYP-3606) |
| MECHANISM | **completeness forces a 5-cycle => contamination => H>7** | **danger over-covers (union bd `(2-n)/n<0`) => measure vanishes** (HYP-3608) |
| what's forbidden | the values `{7,21}` (band gap) | positive lonely measure at the extremal (=0) |
| what carries it | the odd-cycle structure forces H past 7 | EXISTENCE/COUNTING (the units, the odd cycle) |
| `21` | `= 3*Phi_3(2)`, the 2nd forbidden H | `= C(7,2)` = the number of DOUBLETS (binding cores, THM-590) |
| multiplicativity | `H = prod H(strong components)` | the descent over irreducible cores |
- **Same root: ODD CYCLES.** `H` counts them (`Omega` = the OCF odd cycles); the LRC floor binds on one
  (`C_7`). `H=7` is "pure odd-cycle curvature"; the LRC binding atom is one odd cycle.
- **Same mechanism: COMPLETENESS / DENSITY forces over-saturation.** A tournament is complete, so an
  isolated 3-fold-curvature substructure bleeds into a 5-cycle -> `H>7` (forbidden). The LRC danger combs
  are so dense (`n-1` combs of measure `2/n`, union bound `(2-n)/n < 0`, HYP-3608) that they over-cover ->
  the lonely measure vanishes. In BOTH, the complete/dense structure cannot isolate a clean piece; the
  excess forces a gap (no `H=7`) or a collapse (measure `->0`). The "phantom volume" (forbidden H) is the
  counting twin of the "small measure" (vanishing lonely set).
- **`21` is one number, two masks.** `21 = 3*7 = Phi_3(2)*3` (the 2nd forbidden H) `= C(7,2)` (the number
  of doublets, the binding LRC cores, THM-590's `21 doublets + 21 quintuplets`). The same `21`.
- **Band gap `<->` apex-gap landscape.** The H-spectrum's gap `{7,21}` (a numerical-semigroup band gap)
  mirrors the LRC apex-gap landscape's discrete value set `{0, 0.198, 0.308, 1, 2}` (also gappy) and the
  gap-0 cusp (the single forbidden/disproof point).
- **Multiplicativity `<->` descent.** `H = prod H(strong components)` (the S35 condensation, HYP-3603) is
  the multiplicative generation of the H-spectrum; `{7,21}` are the values NO product of irreducible-
  paradox (strong-tournament) H-values can hit. This is the counting analog of the 2-adic descent's
  generation of the apex finite family from irreducible cores.

## Honest placement
The forbidden-H `7` is a COUNTING/completeness impossibility (no tournament has `H=7`); the LRC apex-`7`
floor is a GEOMETRY/non-bipartiteness positivity (`lambda_min(2I+A(C_7)) = 4cos^2(3pi/7) > 0` because `7` is
odd, HYP-3606). They share the prime `7 = Phi_3(2)` and its Fano/octonion structure `{1,2,4}` exactly. The
deep claim -- that LRC14's apex being `7` is "the same root" as the H-forbidden `7`, not a confluence of
`7`'s many specialnesses (Fano, octonion, Heegner `Q(sqrt-7)` `h=1`, `Phi_3(2)`, `ord_7(2)=3`) -- remains
the flagged open question. What is now solid: the two `7`s carry identical arithmetic, and the two
impossibilities (forbidden count / vanishing measure) run on the SAME mechanism (completeness/density forces
over-saturation) about the SAME object (the odd cycle).

## What it buys
Places the 7,21 impossibility precisely inside the recent web: it is the COUNTING (H-spectrum) face of the
odd-cycle/completeness structure whose GEOMETRY face is the LRC apex-7 floor; `7=Phi_3(2)` (Fano/octonion)
is the shared prime, `21=3*7=C(7,2)` the shared number (forbidden H = doublet count), completeness-forces-
over-saturation the shared mechanism, and strong-component multiplicativity (the condensation) the shared
generation. The two are one structure seen by counting vs by measure.

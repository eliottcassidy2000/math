---
id: HYP-3591
title: "Even graph" is 3 distinct counts (ALL=A000088, tournaments=Royle-even=A000568=odd-order Burnside=signed cycle index P_n(1), Eulerian/degree-even=A002854=cycle space), equal only at n=3, with the sandwich Eulerian<=tournaments<=all; and the genus IS the odd boundary -- matching codex-S675b (merged-metagraph BLACK=Eulerian/cycle-space bulk, BLUE=odd-boundary=floor atoms; "odd=even+named boundary") to the modular split (Eisenstein bulk (+) cusp form obstruction), dim(obstruction)=genus=#blue floor atoms; LRC(14) genus 1 = one blue atom = the doublet = the cusp form f_14
status: VERIFIED counts (all three by GF(2)-correct Burnside, n=3..7) + grounded synthesis (S675b black/blue verified through n=7 by codex; the modular match is the conjectural-but-tight unification of HYP-3587). Corrects the naive Eulerian dim (GF(2) parity obstruction at even n).
source: klein-2026-06-29-S12
depends_on:
  - HYP-3587   # genus = local-global gap (Eisenstein bulk / cusp obstruction)
  - THM-587    # tournaments = odd-order Burnside = signed cycle index P_n(1)
related:
  - merged-line-parity-even-odd-s675b   # codex: black=Eulerian, blue=odd-boundary=floor atoms
  - even-graphs-through-the-metagraph    # opus-S260: the three "evens", Royle equinumerosity, E=Cut+Cycle iff n odd
  - the-permanent-is-the-unsigned-twin-of-the-spectrum   # signed/unsigned
  - THM-578    # the doublet = the single genus-1 blue floor atom
  - HYP-3586   # cusps=Klein, hardness=genus
results:
  - 04-computation/even_graph_equinumerosity_three_counts_klein.py
  - 05-knowledge/results/even_graph_equinumerosity_three_counts_klein.out
---

# HYP-3591 — three "evens", the Royle sandwich, and the genus is the odd boundary

## The equinumerosity particularities (verified)

"Even graph" = three Burnside sums on `K_n`'s edges, equal only at `n=3`:

| count | meaning | Burnside | n=3..7 |
|---|---|---|---|
| A000088 | ALL graphs | `(1/n!)Σ_{all σ} 2^{E_orb}` | 4,11,34,156,1044 |
| A000568 | tournaments = **Royle-even** | `(1/n!)Σ_{|σ| odd} 2^{E_orb}` = `P_n(1)` (signed cycle index) | 2,4,12,56,456 |
| A002854 | degree-even = **Eulerian** = cycle space | `(1/n!)Σ_{all σ} 2^{E_orb − rank_{GF2}(deg)}` | 2,3,7,16,54 |

- The tournament equinumerosity is with **Royle-even** (DFGPR 2022), `=` the **odd-order restriction** of
  the all-graph Burnside `=` the signed cycle index `P_n(1)` (THM-587; `|σ| odd ⟺ all cycle lengths odd`).
  It is **NOT** the degree-even/Eulerian count.
- **Sandwich:** `Eulerian ≤ tournaments ≤ all graphs`. Tournaments = the odd-order *slice*; Eulerian = the
  cycle-space *quotient*. The surplus `A000568 − A002854 = 0,1,5,40,402` is the obstruction-bearing part.
- The Eulerian count needs the **GF(2) rank** of the degree map (`E = Cut ⊕ Cycle iff n odd`); the naive
  `E_orb − V_orb + 1` is wrong at even `n`. All three `= 2` at `n=3` (the polysemy trap).

So: tournaments = signed/odd-order object; Eulerian = cycle-space bulk; the even-graph shadow of a
tournament (`T_cycle = (I+L)T mod 2`, odd `n`) is its cycle-space/bulk part.

## The Eulerian sections of the merged metagraph (codex-S675b)

In the corrected complement-line lift in `Q_m`: **black is Eulerian** (boundary-zero over `F_2` = cycle-space
kernel), **blue is the odd-degree boundary** (an affine coset of the cycle space); verified through `n=7`.
"**Odd is even plus a named boundary.**" S675b's punchline: *the black Eulerian side is a certificate, and
the blue odd boundary marks the exact floor atoms.*

## The synthesis: genus = the odd boundary

Matching S675b to the modular split (Eisenstein bulk ⊕ cusp form obstruction, HYP-3587):

| BULK / black / even / Eisenstein | OBSTRUCTION / blue / odd / cusp form |
|---|---|
| Eulerian = cycle-space kernel | odd-boundary coset = "even + a named boundary" |
| the even-graph shadow of the tournament | the part the shadow misses |
| dim = cycle space | **dim = genus = #blue floor atoms** |

> **The genus is the dimension of the odd boundary** — the number of named-boundary (blue) atoms beyond the
> Eulerian/even/cycle-space (black) bulk. The cusp form `f_14` IS the blue odd-boundary; S675b's "blue floor
> atoms" = the modular "cusp-form obstruction," found independently. LRC(14): `genus 1` = **one** blue atom
> = the **doublet** (THM-578, `4cos^2(3pi/7)`) = `f_14`. Genus 0 (LRC(6)) = zero blue atoms (no cusp form,
> black certificate closes it). Genus ≥ 2 = several blue atoms (irreducible).

## The master dichotomy, homological row

local/global = even/odd = cycle-space-kernel/boundary-coset = Eulerian-black/blue = Eisenstein/cusp-form =
bulk/obstruction; `dim(obstruction) = genus`. Everything computable is **black** (the Eulerian/cycle-space
bulk, the metagraph, `CV(H)`); the one missing thing is the **blue** boundary — the genus counts it, and
for `N=14` it is a single doublet. The floor closes when the one blue atom is bounded (HYP-3587: the leading
apex-cusp coefficient of `f_14`).

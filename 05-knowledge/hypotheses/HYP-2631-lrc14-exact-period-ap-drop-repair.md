---
id: HYP-2631
title: LRC(14) exact-period AP-drop repair - the Q=210 blind AP mouths are repeated-prime packets inside the raw 1260 carrier
status: PARTIALLY-CONFIRMED
source: codex-2026-06-19-S22
depends_on:
  - HYP-2629
  - HYP-2628
  - HYP-2627
  - HYP-2625
  - HYP-2569
  - THM-523
related:
  - HYP-2626
  - HYP-2561
  - OPEN-Q-108
---

# HYP-2631 - LRC(14) Exact-Period AP-Drop Repair

## Claim

HYP-2628 found that the squarefree `Q=210` grid catches only `11/13`
AP one-drop cores, missing exactly drops `6` and `12`, while `Q=1260`
catches all `13/13`.

HYP-2631 identifies the quotient loss:

```text
210 = 2 * 3 * 5 * 7
1260 = 2^2 * 3^2 * 5 * 7
```

Inside the raw `1260` exact-period divisor lattice, every strict-safe residue
for the two `Q=210`-blind AP drops has reduced denominator not dividing `210`.
Thus the squarefree radical quotient deletes the exact-period packets that
sample those two cusp mouths.

The repair packets are:

```text
drop 6:  reduced denominators 63, 420, 630
         = 3^2*7, 2^2*3*5*7, 2*3^2*5*7

drop 12: reduced denominators 12, 315, 630, 1260
         = 2^2*3, 3^2*5*7, 2*3^2*5*7, 2^2*3^2*5*7
```

So the radical miss is not generic coarseness.  It is exact-period packet
collapse: the two missed AP mouths need prime-power information retained by
the raw Hill carrier `1260` and lost by `rad(1260)=210`.

## Computation

Script:

- `04-computation/lrc14_exact_period_ap_drop_repair_codex_s22.py`
- output: `05-knowledge/results/lrc14_exact_period_ap_drop_repair_codex_s22.out`

The script computes exact AP one-drop safe components, grid residues on
selected denominators, reduced denominator profiles, and component-level
`Q=210` versus `Q=1260` hits.

## Exact Findings

The AP one-drop grid table begins:

```text
drop  meas(G_C)        comps  Q=210  Q=315  Q=420  Q=630  Q=1260
   6           7/858      4      0      2      4      6      10
  12       426/35035      4      0      4      4      8      16
  10      1520/63063      4      4      8     10     14      30
   4         97/4004      4      4      6     10     14      30
```

Thus `Q=210` misses exactly drops `6` and `12`.

For drop `6`, the raw `Q=1260` strict-safe residues split as:

```text
denom  hits  phi(denom)  factor
   63     2          36  3^2*7
  420     4          96  2^2*3*5*7
  630     4         144  2*3^2*5*7
```

For drop `12`:

```text
denom  hits  phi(denom)  factor
   12     4           4  2^2*3
  315     4         144  3^2*5*7
  630     4         144  2*3^2*5*7
 1260     4         288  2^2*3^2*5*7
```

None of these denominators divides `210`, so none can appear on the radical
grid.  Every `Q=1260` component hit is also a true one-drop witness: in the
full AP row `{1,...,13}`, the only dangerous speed at that residue is the
omitted speed.

Component-level examples:

```text
drop 6, component (29/168, 27/154):
  Q=210 residues:  []
  Q=1260 residues: [218, 219, 220]
  reduced denoms:  [630, 420, 63]
  full-AP danger:  [(6,), (6,), (6,)]

drop 12, component (29/70, 41/98):
  Q=210 residues:  []
  Q=1260 residues: [523, 524, 525, 526, 527]
  reduced denoms:  [1260, 315, 12, 630, 1260]
  full-AP danger:  [(12,), (12,), (12,), (12,), (12,)]
```

## Caveat

This is not a global minimal-denominator theorem.  Drop `6` has an earlier
strict-safe `q=98` witness outside the `1260` divisor lattice.  HYP-2631 is
about the exact-period packets retained by the raw `K_14` product, which is the
carrier used by HYP-2627/HYP-2629.

## Interpretation

This sharpens the raw-product discipline:

```text
raw 1260 exact-period carrier
-> AP cusp mouth packets
-> squarefree / coimage projection
```

must not be replaced by

```text
radical 210 first
-> AP cusp mouth packets later
```

because the latter order deletes the packets that witness the two radical-blind
AP mouths.

The divided crossing quotient `315` is intermediate.  It sees some `3^2`
packets (`drop 6` has `Q=315` hits; `drop 12` has denominator `315` hits), but
it loses the dyadic `2^2` lane.  This matches HYP-2629's copy-mass warning:
`315` no longer carries the full `{2,3,5,7}` exact-period profile.

## Proof Route

The new local proof obligation is an exact-period mouth lemma:

```text
For every AP-cusp or low-height boundary face, retain reduced-denominator
packets through the raw carrier before projecting to squarefree masks or
mod-7 coimage classes.
```

Concretely, the AP cap face suggests a finite transfer state:

```text
dihedral endpoint mouth
-> reduced denominator packet
-> Euler-copy mask mass
-> squarefree/coimage quotient
-> signed tail.
```

This does not prove LRC(14).  It gives a cleaner explanation of why `1260`
is proof-relevant: not just as a numeric denominator, but as the smallest
current carrier in the workspace that simultaneously retains the dyadic and
triadic repeated-prime packets needed by the two radical-blind AP cusp mouths.

## Tournament Analysis

The computation uses proof quotients rather than runners.  Candidate vertices
included runners, AP drops, safe components, endpoint mouths, residues `a/Q`,
reduced denominators, squarefree masks, exact-period packets, and proof
obligations.

Hamiltonian path:

```text
raw_1260_exact_period_carrier
> drop6_3square7_lane
> drop12_2square3_lane
> crossing_315_partial_carrier
> radical_210_squarefree_grid
> raw_runner_vertices
```

Fingerprints:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1}
directed_3_cycles = 0
SCC_sizes = [1,1,1,1,1,1]
hamiltonian_paths = 1
```

The quotient preserves the strict-safe AP one-drop predicate and exact-period
packet data.  It destroys continuous mouth length and endpoint-owner history,
which remain in HYP-2569's dihedral endpoint atlas.

## Status

Partially confirmed by exact rational computation.  Open next steps:

1. Prove the AP-drop denominator profiles directly from endpoint inequalities.
2. Extend the reduced-denominator packet ledger from AP drops to the HYP-2626
   repeated coimage tail.
3. Check whether other low-height boundary/cusp faces have analogous
   radical-blind but raw-carrier-visible repeated-prime packets.

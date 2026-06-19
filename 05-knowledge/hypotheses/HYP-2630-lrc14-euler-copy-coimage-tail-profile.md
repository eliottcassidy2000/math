---
id: HYP-2630
title: LRC(14) Euler-copy coimage tail profile - repeated residues are multi-large character packets, not one-large wall residue
status: OPEN
source: codex-2026-06-19-S22
depends_on:
  - HYP-2629
  - HYP-2628
  - HYP-2626
  - HYP-2624
  - HYP-2619
  - HYP-2617
related:
  - HYP-2631
  - HYP-2614
  - HYP-2618
  - HYP-2608
  - OPEN-Q-108
---

# HYP-2630 - LRC(14) Euler-Copy Coimage Tail Profile

## Claim

The HYP-2629 next target has a split answer.  Re-indexing the HYP-2626
`k=10` repeated-residue tail by Euler-copy data sharpens the proof obligation,
but Euler-copy mass alone does not explain the quadratic-character split.

The correct transfer has three layers:

```text
coimage class
-> minimal large-residue demand after the bounded core 1..13
-> Euler-copy unit packet capacity and character moments.
```

The exact-period copy ledger says the unit packets are uniformly distributed
over `F_7^*`.  Therefore copy capacity thickens repeated-residue packets by
multiplicity pattern, but it is blind to QR/NQR inside a fixed pattern.  The
extra coordinate is the quadratic-character phase moment.

## Computation

Script:

- `04-computation/lrc14_euler_copy_coimage_tail_codex_s22.py`
- output: `05-knowledge/results/lrc14_euler_copy_coimage_tail_codex_s22.out`

The script recomputes the `k=10`, `d=9` HYP-2617 coimage rows, compares
height-2 and height-3 one-large wall class coverage, computes the residue
capacity of the bounded core `1..13`, and evaluates exact-period unit packet
counts for `q=14,210,1260`.

## Exact Findings

Exact-period packets touching the prime `7` are uniform on the six unit
residues:

```text
q     phi(q)  exact top per unit  full {2,3,5,7} per unit  all 7-touch per unit
14       6          1                       0                         2
210     48          8                       8                        30
1260   288         48                      96                       180
```

For the raw `K_14` denominator `1260`, the exact top-period seam gives
`48` copies per unit residue, while the full squarefree `{2,3,5,7}` mask atom
gives `96` copies per unit residue.  Both are uniform on `F_7^*`.

The one-large wall height check is stable:

```text
height  supports  classes  nonzero hit  mass share
2          68124      108       85/116   84.229179%
3         114597      108       85/116   84.229179%
```

Thus the HYP-2624 tail is not a missed height-3 one-large wall.  The reason is
structural.  The bounded core `1..13` has residue capacity

```text
{0:1, 1:2, 2:2, 3:2, 4:2, 5:2, 6:2}.
```

One large speed can add only one more residue.  A class with four equal
nonzero residues needs at least two large speeds.  Consequently the dominant
tail packets are genuinely multi-large residue packets.

For the `31` tail-only nonzero classes after height `<=2`, the demand profile
is:

```text
demand  classes  |S_9|-mass  share of tail
2            27    1.6054892   94.382022%
3             4    0.0955648    5.617978%
```

Pattern summary:

```text
4+2 repeated       6 classes   |S|-mass 1.0321002
4+1+1 repeated    12 classes   |S|-mass 0.51605011
mixed             10 classes   |S|-mass 0.095564835
zero-cusp          3 classes   |S|-mass 0.057338901
```

For the nonzero `4+2` packet `(1,1,1,1,a,a)`, full-mask copy capacities are
identical inside the row:

```text
a  chi(a)  chi-sum  |S_9|       C_full210  C_full1260
2    +1       6      0.23891209     1960    15148137600
4    +1       6      0.23891209     1960    15148137600
3    -1       2      0.17201670     1960    15148137600
5    -1       2      0.17201670     1960    15148137600
6    -1       2      0.17201670     1960    15148137600
```

So the QR/NQR split is not a copy-capacity split.  It is a quadratic-character
phase split:

```text
QR mean |S_9|  = 0.23891209
NQR mean |S_9| = 0.17201670
```

For the `4+1+1` packet, the useful address is the signature

```text
(chi(a), chi(b), chi(ab), chi((a-1)(b-1)), 4+chi(a)+chi(b)).
```

The largest signature bucket is `(1,-1,-1,1,4)`, with `3` classes and total
`|S|`-mass `0.2293556`.

## Interpretation

HYP-2630 converts the HYP-2629 test into a cleaner proof target:

```text
Euler-copy exact-period mass = packet capacity / residue layer
quadratic character          = phase layer
multi-large wall relation    = incidence / compatibility layer
```

This matches the older residue/phase/incidence and OCF packet lessons from the
repo.  First retain the finite packet address; then keep the character phase;
only then evaluate the signed compatible mass.  Raw copy mass is necessary but
too scalar to close the repeated-residue tail.

After rebasing over HYP-2631, the AP-drop repair atlas gives the same quotient
discipline from the cusp-mouth side: reduced-denominator exact-period packets
inside raw `1260` must be retained before squarefree or coimage projection.
HYP-2630 is the analogous tail-side demand, where raw packet capacity survives
but the proof still needs the `chi_7` phase and two-large incidence relation.

The practical proof redirect is:

```text
do not keep raising one-large wall height;
attack the two-large repeated-residue character packet directly.
```

The expected next theorem is a signed cotangent/Dedekind estimate for the
two-large `4+2` and `4+1+1` packets, with the QR/NQR phase channel retained.

## HYP-2632 Follow-Up

HYP-2632 executes the first finite-kernel step of that theorem.  The QR/NQR
split in the `4+2` row is exactly

```text
2*S_9(1,1,1,1,a,a)/U = -43 - 7*chi_7(a),  a=2..6.
```

For `4+1+1`, the four-character signature used above is not complete.  The
additional hidden coordinate is the affine zero lane

```text
a+b = 2 mod 7.
```

Off that zero lane, HYP-2632 finds another Legendre selector:

```text
Q(a,b)=ab*(1+3(a+b))-1,
S/U=8 iff chi_7(Q)=+1,
S/U=1 iff chi_7(Q)=-1.
```

The repeated-kernel signed ledger is `-108U` from `4+2` and `+54U` from
`4+1+1`, so the next analytic estimate should exploit this signed table
directly instead of applying an absolute bound to the `162U` packet mass.
The companion Dedekind-shell computation
`lrc14_two_large_dedekind_phase_codex_s23.py` refines the same kernel into
explicit factors `D_T(ell)=sum_r r chat(r,T) zeta_7^(ell r)` and confirms why
the reciprocal-tail proof must expose additive frequencies before using any
absolute matrix bound.

## Tournament Analysis

Vertices are proof obligations and quotient stages, not runners:

```text
multi_large_residue_demand
> quadratic_character_phase
> euler_copy_unit_capacity
> exact_period_1260_packets
> height3_one_large_probe
> raw_apex_mask
> raw_runner_vertices
```

Pairwise observable: preservation of the signed `k=10` coimage-tail predicate.
The switch prefers the quotient that keeps residue multiplicity, copy capacity,
and character phase while discarding raw row data.  Fingerprints:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3_cycles = 0
```

## Status

Open / partially confirmed.  The computation does not prove LRC(14), but it
rules out the simplest "more one-large wall height" explanation for the
HYP-2624 tail and identifies the next analytic object more precisely:
two-large repeated-residue character packets over the exact-period unit seam.

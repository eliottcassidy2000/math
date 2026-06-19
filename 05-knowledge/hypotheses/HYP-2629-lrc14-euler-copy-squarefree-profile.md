---
id: HYP-2629
title: LRC(14) Euler-copy squarefree profile - totient multiplicities refine the mod-210 crossing carrier
status: OPEN
source: codex-2026-06-19-S21
depends_on:
  - HYP-2628
  - HYP-2627
  - HYP-2626
  - HYP-2625
  - THM-538
related:
  - HYP-2624
  - HYP-2617
  - HYP-2619
  - THM-523
  - OPEN-Q-108
---

# HYP-2629 - LRC(14) Euler-Copy Squarefree Profile

## Claim

The user's divisor-copy rule

```text
Find copy counts c(n) >= 1 such that sum_{d|n} c(d) = n.
```

is exactly the Euler totient:

```text
c(n) = phi(n),
```

by Mobius inversion.  The LRC14 use is not merely to rename `phi`; it is to
turn the squarefree divisor-profile route into a weighted copy recurrence.
For a fixed prime set `{2,3,5,7}`, define

```text
copy_mass_N(M) = sum_{d|N, mask(d)=M} phi(d),
```

where `mask(d)` records which of `{2,3,5,7}` divide `d`.  Then

```text
N = sum_M copy_mass_N(M).
```

This refines HYP-2628 and HYP-2627.  HYP-2628 identifies the exact-period
`phi` packet law; HYP-2629 adds the Hill-row/crossing-value scan.  The raw
`K_14` Hill product

```text
P_14 = 5*6*6*7 = 1260,   rad(P_14)=210.
```

has full `{2,3,5,7}` copy mass `576`, while the divided crossing value

```text
cr(K_14)=315
```

has full `{2,3,5,7}` copy mass `0`, because the division by `4` deletes the
tracked dyadic gate.  Thus the Euler-copy profile gives a precise reason for
the HYP-2627 rule:

```text
retain the raw Hill product ledger before quotienting to the crossing number.
```

## Computation

Script:

- `04-computation/lrc14_euler_copy_squarefree_profile_codex_s21.py`
- output: `05-knowledge/results/lrc14_euler_copy_squarefree_profile_codex_s21.out`

The script verifies the divisor-copy rule through `n<=60`, computes squarefree
copy profiles for `6`, `30`, `210`, `315`, `1260`, and scans Hill products
through `n=32`.

## Exact Findings

The squarefree recurrence is explicit.  For squarefree `S`,

```text
w_S(M) = product_{p in M} (p-1).
```

Adding a new prime `p` keeps the old masks and appends a shifted layer:

```text
w_{S*p}(M union {p}) = (p-1) w_S(M).
```

For the LRC14 address chain this gives:

```text
add 2: total 1 -> 2,   full {2}       mass 1  = 1*1
add 3: total 2 -> 6,   full {2,3}     mass 2  = 2*1
add 5: total 6 -> 30,  full {2,3,5}   mass 8  = 4*2
add 7: total 30 ->210, full {2,3,5,7} mass 48 = 6*8
```

So the existing HYP-2625 ladder

```text
{2,3} -> {2,3,5} -> {2,3,5,7}
```

is a genuine Euler-copy recurrence, not just a memorable list of moduli.

The comparison around `K_14` is:

```text
object                 N     factor          phi   rad   full {2,3,5,7} mass
mod210 radical         210   2*3*5*7          48   210        48
K13 Hill product       900   2^2*3^2*5^2     240    30         0
K14 raw Hill product  1260   2^2*3^2*5*7     288   210       576
K14 crossing value     315   3^2*5*7         144   105         0
K15 Hill product      1764   2^2*3^2*7^2     504    42         0
```

Thus `n=14` is the first Hill row with nonzero full `{2,3,5,7}` copy mass,
and it is also the first row where the raw product has the full mask while the
divided crossing value loses it.

For `P_14=1260`, the full mask atom is

```text
(2^2-1)(3^2-1)(5-1)(7-1) = 3*8*4*6 = 576.
```

The radical profile has full mass

```text
phi(210)=48.
```

The repeated `6,6` blocks thicken this full mask by

```text
((2^2-1)/(2-1)) * ((3^2-1)/(3-1)) = 3*4 = 12,
```

so `48` becomes `576`.  Ordinary `phi(P_14)=288` only sees the p-fold exponent
thickening; the full squarefree copy atom sees all nonzero p-power divisor
layers.

## Markov-Hurwitz Readout

The copy frame also clarifies the role of the Markov-Hurwitz equation.  Since

```text
x = sum_{d|x} phi(d),
```

the product `wxyz` is literally a count of four-block Euler-copy assignments.
For `q_14=(5,6,6,7)`,

```text
wxyz = 1260
w^2+x^2+y^2+z^2 = 146
```

so the tuple is still far from Markov-Hurwitz.  The totient-copy frame explains
the right-hand product as a copy-assignment ledger; it does not turn the
complete-graph row into an energy-critical Markov-Hurwitz solution.

## Interpretation

HYP-2629 strengthens three earlier threads:

1. **HYP-2625:** the mod `6 -> 30 -> 210` ladder is a prime-extension copy
   recurrence with shifted mass multiplied by `p-1`.
2. **HYP-2626:** the live fixed LRC14 residual remains the prime-7 coimage
   seam, but the Euler-copy profile explains why that seam belongs in the same
   four-prime ledger as the dyadic/mod-30 data.
3. **HYP-2627:** the raw `K_14` product must be retained before taking the
   crossing quotient, because `P_14` carries the full copy mask and `315` does
   not.

This still does not prove LRC(14).  It gives a sharper accounting coordinate
for the next proof obligation: re-index the HYP-2626 repeated-residue tail by
Euler-copy mask mass, not merely by raw radical or raw support tuple.

## Tournament Analysis

The computation uses proof-carrier quotients rather than runners.  Candidate
vertices included runners, divisors, individual divisor copies, squarefree
masks, Hill blocks, crossing values, Markov-Hurwitz coordinates, mod-7 coimage
classes, and proof obligations.

Hamiltonian path:

```text
euler_copy_squarefree_profile
> raw_hill_product
> mod210_address
> prime7_coimage_seam
> markov_hurwitz_energy_surface
> crossing_value
> raw_divisor_lattice
> raw_runner_vertices
```

Fingerprints:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3_cycles = 0
```

The quotient preserves the proof-relevant mod-210 copy address.  It destroys
exact divisor representatives and witness-time geometry, which must be
recovered later through finite wall and coimage ledgers.

## Status

Open / exploratory.  HYP-2629 is a transfer-coordinate refinement, not a proof.
The next test is to apply the copy profile directly to the HYP-2626 k=10
tail-only repeated packets and ask whether the `{7}` coimage mass separates the
quadratic-character cases more cleanly than raw prime masks do.

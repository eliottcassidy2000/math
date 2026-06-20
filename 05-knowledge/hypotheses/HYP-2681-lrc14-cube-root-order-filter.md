---
id: HYP-2681
title: LRC(14) cube-root order filter for signed multi-far residuals
status: OPEN; exact algebraic lens, proof value under test
source: codex-2026-06-20-S52
depends_on:
  - HYP-2680
  - THM-548
  - HYP-2679
  - HYP-2677
  - HYP-2676
  - HYP-2639
related:
  - THM-511
  - HYP-2532
  - HYP-2533
  - OPEN-Q-108
---

# HYP-2681 - Cube-Root Order Filter

## Claim Being Tested

The simultaneous-peel residuals from HYP-2680 may carry a useful `C_3`
structure.  For a bounded core `B` and far triple `F={u,v,w}`, write the seven
nonempty residual packets

```text
A = Delta_u  - Phi_1
B = Delta_v  - Phi_1
C = Delta_w  - Phi_1
D = Delta_uv - Phi_2
E = Delta_vw - Phi_2
F = Delta_wu - Phi_2
G = Delta_uvw - Phi_3.
```

The actual simultaneous-peel correction is

```text
H(1) = A+B+C+D+E+F+G.
```

The user's recursion

```text
A+B+C-D-E-F+G
```

is the pair-tax shadow

```text
H(1) - 2*(D+E+F).
```

The cube-root filter keeps more information.  With `omega^2+omega+1=0`, define
Eisenstein cyclic modes

```text
S_omega = A + omega*B + omega^2*C
P_omega = D + omega*E + omega^2*F.
```

The first question is whether dangerous triples are controlled not by scalar
rank alone, but by the scalar pair-tax shadow together with the exact
Eisenstein imbalance `S_omega-P_omega`.

## Challenged Assumption

Do not assume the seven packets should be collapsed immediately to total
residual or absolute mass.  The quotient

```text
(A,B,C ; D,E,F ; G) -> (H(1), A+B+C-D-E-F+G, S_omega-P_omega)
```

preserves the signed order ledger, the pair layer, and the cyclic orientation
imbalance.  It destroys individual packet ownership beyond the chosen cyclic
order, so any use in proof must record the core/far phase address.

## Next Computation

Completed exact Fraction scout:

- `04-computation/lrc14_cube_root_order_filter_codex_s52.py`
- `05-knowledge/results/lrc14_cube_root_order_filter_codex_s52.out`

The script:

1. computes the seven packets for named triples from HYP-2680;
2. computes the pair-tax recursion `A+B+C-D-E-F+G`;
3. computes exact Eisenstein coefficients and norms for `S_omega`, `P_omega`,
   and `S_omega-P_omega`;
4. scans the `(15,16,17)` all-core bank to see whether direct cap risk,
   three-far residual risk, and cube-root imbalance are different orderings.

## Exact S52 Findings

The recursion is exact but it is not the cap residual.  Examples:

```text
dilated core (0,4,6,8,10,12,14), far=(15,16,17):
  actual H(1)=9243793/68572560
  A+B+C-D-E-F+G=767717/4571504
  S_omega-P_omega norm=546205453/13322668800

same core, far=(22,23,24):
  actual H(1)=4638863/170086840  (positive)
  A+B+C-D-E-F+G=-5652649/170086840  (negative)

same core, separated far=(17,23,31):
  actual H(1)=4331573/239667820  (positive)
  A+B+C-D-E-F+G=-159950761/12223058820  (negative)

consec8, far=(15,16,17):
  actual H(1)=57167/17143140  (positive)
  A+B+C-D-E-F+G=-41563/2285752  (negative)

third-pocket active, far=(30,33,35):
  singles are all zero
  actual H(1)=198629/22185240
  pair-tax shadow=2841/7395080
  S_omega-P_omega norm=1483/85377600
```

The small bank for far `(15,16,17)` over `3003` primitive bounded cores gives:

```text
actual residual signs:   +2821 / -182
pair-tax shadow signs:   +1250 / -1753
triple packet G signs:   +1999 / -1004
pair layer signs:        +2651 / -352
actual and recursion have same nonzero sign in only 1368 rows
```

The four main lenses pick different leaders:

```text
actual residual leader:
  (0,5,10,11,12,13,14,15,16,17)

pair-tax shadow leader:
  (0,4,6,8,10,12,14,15,16,17)

Eisenstein imbalance leader:
  (0,4,6,8,10,12,14,15,16,17)

direct p0 leader:
  (0,9,10,11,12,13,14,15,16,17)
```

Tournament Analysis on these proof lenses is transitive:

```text
direct_p0 > eisenstein_imbalance > pair_tax_shadow > actual_residual
```

by the chosen observable "leader with larger direct `p0`, tie by actual
residual."  This is a warning: cube-root imbalance is structurally meaningful,
but it is not automatically direct cap risk.

## Revised Target

For a signed three-far proof, keep the packet triple

```text
(actual residual H(1), pair-tax shadow H(1)-2R_2, Eisenstein imbalance)
```

until the finite atlas or signed Abel estimate is applied.  Collapsing to only
one of these scalars loses information: the pair-tax shadow often reverses sign
relative to the actual residual, and direct risk can be led by a different row.

No proof of LRC(14) is claimed.

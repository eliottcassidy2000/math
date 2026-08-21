---
id: THM-3626
title: "AMM R=8192 adjoint horizon and finite-scale phase rebound"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + FINITE-EXACT + VERIFIED-EXACT /
  PENDING INDEPENDENT HOSTILE AUDIT.  At the archived last-failing Rule-A
  epoch (R,D0)=(8192,191), the exact trace dies at row 2045.  Its positive
  truncated-adjoint wall changes sign between cuts 773 and 774, forcing every
  admissible continuation avoiding that fatal inequality to depart by row
  773.  The exact ratio 773/8192 lies below (3-sqrt(5))/8.  Its error is
  9/8192 smaller than at R=4096 but 7/8192 larger than at R=2048, proving a
  finite-scale phase rebound rather than monotone sharpening.  No alternative
  feasibility, eventual limit, uniform extractor, or AMM bound is claimed.
source: kps-s189 / THM-3616 hostile dyadic continuation, 2026-08-21
audit: >
  PENDING -- author exact reconstruction and optimization replay passed 27
  gates, including the inherited R=4096 control, complete degree profile,
  fatal row, contiguous adjoint wall, boundary and active-ledger digests,
  exact radical comparisons, dyadic defects, AST assertion exclusion, and
  parent/source pins.  Independent hostile reconstruction remains required.
depends_on:
  - THM-3616-amm-R4096-adjoint-horizon-and-golden-finite-scale-obstruction
related:
  - THM-3601-amm-r2048-terminal-failure-adjoint-horizon
  - THM-3597-amm-dyadic-rule-a-transition-and-adjoint-horizons-through-R1024
script: 04-computation/amm_binary_sturmian_R8192_adjoint_horizon_thm3626.py
output: 05-knowledge/results/amm_binary_sturmian_R8192_adjoint_horizon_thm3626.out
script_sha256: 548e5ddaad0b7220687468552a448f3021ad16b3162a2c11409d139d31ee3be5
output_sha256: 4319f792b1c69b84ea30496ef0500d3666da73575d4de267a60428d5c757300b
semantic_sha256: 226de0151e7fe6f9ded73b11041bb83d259f8d4d5d44c738de8f9b706d052d59
hash_basis: LF-normalized bytes
---

# THM-3626 -- AMM `R=8192` adjoint horizon and phase rebound

**RESERVED / PROVISIONAL PROOF CANDIDATE + FINITE-EXACT + VERIFIED-EXACT /
PENDING INDEPENDENT HOSTILE AUDIT.**  This is the next hostile dyadic epoch
after THM-3616.  It certifies one archived failed Rule-A prefix and sharpens
the finite-scale phase diagnosis; it does not construct a repair.

## 1. Exact fatal row and adjoint wall

For `0<=i<8192`, use the inherited exact Fibonacci--Lucas comparison to set

```text
d_i(D0)=floor(log_5(phi^(2(8192+i))))+D0.               (1)
```

At the archived last-failing Rule-A offset `D0=191`, the exact profile has

```text
first death j=2045,
(d_0,d_2045,d_8191)=(5089,6312,9987).                   (2)
```

The fatal integer is negative and has `4937` bits.  Use the same halved junk
states, Lucas bounds, Pascal transitions, and positive truncated adjoint as
THM-3616.  If an admissible alternative agrees with this Rule-A trace through
row `s-1`, transposed propagation of the fatal-top evaluation supplies the
necessary inequality

```text
0<=B_s.                                                  (3)
```

The complete exact sweep gives one contiguous negative wall:

```text
B_773>0>B_774,
B_s<0 exactly for s=774,775,...,2044.                   (4)
```

Thus agreement through row `773` would force the false inequality
`0<=B_774`.  Every admissible continuation that avoids this fatal inequality
must depart from Rule A at or before

```text
q=773.                                                   (5)
```

The two boundary integers in `(4)` have bit lengths

```text
bits(B_773)=4939,                  bits(B_774)=4930,      (6)
```

and their ordered exact digest is

```text
d42558aa87d6c8e57f739abdb54758ec192d0c6081d3ab850ce04c746cedecc3. (7)
```

At the first negative cut, the positive multiplier ledger contains `808356`
active cells and has digest

```text
dcbbf6f02e3c4b2d01ad0d8e35a6e8b46a73fa8d08c07969306adc544c4a5ba5. (8)
```

These compact digests pin the enormous exact integers without printing them.

## 2. Exact phase rebound

The normalized departure coordinate is

```text
q/R=773/8192=0.0943603515625.                            (9)
```

Compare it with the candidate golden coordinate from THM-3601,

```text
theta_gold=(3-sqrt(5))/8=0.09549150281252628795... .    (10)
```

The new signed discrepancy is

```text
q/R-theta_gold=(-2299+1024sqrt(5))/8192<0.              (11)
```

The sign is exact because

```text
2299^2-5*1024^2=42521>0.                                (12)
```

Putting the absolute discrepancies at all three audited dyadic scales over
the common denominator `8192` gives

```text
E_2048=(2292-1024sqrt(5))/8192,
E_4096=(2308-1024sqrt(5))/8192,
E_8192=(2299-1024sqrt(5))/8192.                         (13)
```

Consequently

```text
E_4096-E_8192=9/8192>0,
E_8192-E_2048=7/8192>0.                                 (14)
```

The `R=8192` phase has rebounded toward the golden coordinate after the
`R=4096` deterioration, but has not recovered the `R=2048` accuracy.  This
three-scale order is exact and nonmonotone.  It neither proves nor disproves
eventual convergence, a convergent subsequence, or a phase-corrected law.

Two integer dyadic defects record the same finite-scale displacement:

```text
q_8192-2q_4096=773-2*382=9,
j_8192-2j_4096=2045-2*1014=17.                          (15)
```

Their relation to the local Sturmian phase, normalized offset `D0/R`, and
clamp-saturation profile remains open.

## 3. Exact recurrence and hostile controls

The proof uses the compressed top-distance recurrence from THM-3616.  Starting
from the fatal row, the backward multipliers satisfy

```text
flat degree step:   M_(i-1)=(1+z)M_i,
rising step:        M_(i-1)=Pi_(>=0)[z^(-1)(1+z)^2M_i]. (16)
```

All coefficients are nonnegative integers.  A shrinking forward pass pairs
them with the exact nonnegative clamp charges, and a suffix sum produces every
`B_s`.  The initial `2731` top-distance cells cover the entire fatal scan.

Before computing `(2)--(15)`, the companion independently replays the parent
control

```text
(R,D0,j,(d_0,d_j,d_(R-1)))
  =(4096,88,1014,(2537,3143,4986)).                     (17)
```

It then checks every degree step is binary, the forward and adjoint fatal
integers agree, the wall is contiguous, all boundary signs and metadata are
exact, and the radical comparisons in `(14)` are integer identities.

## 4. Scope and reproduction

Reproduce with

```bash
python3 04-computation/amm_binary_sturmian_R8192_adjoint_horizon_thm3626.py
python3 -O 04-computation/amm_binary_sturmian_R8192_adjoint_horizon_thm3626.py
```

The normal and optimized runs execute the same `27` optimization-safe gates,
contain no Python assertion statement, and agree with the stored transcript
after LF normalization.

This candidate does **not** replay the adjacent archived survivor or prove the
full entry polytope infeasible.  The adjoint wall constrains only continuations
agreeing with the archived failed Rule-A prefix through a cut.  It proves no
alternative feasibility, uniform extractor, asymptotic horizon law, value of
the AMM constant, or improvement to its bound.

Independent hostile audit remains the promotion gate.

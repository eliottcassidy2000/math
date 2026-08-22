---
id: THM-3689
title: "Fully transverse two-by-two sparse-support Keller gate"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT.  In the normalized planar chart
  P=x plus two nonlinear monomials and Q=y plus two nonlinear monomials,
  suppose both P monomials have nonzero x exponent, both Q monomials have
  nonzero y exponent, and every P/Q exponent pair has nonzero determinant.
  Then Jac(P,Q) is not constant.  The exact degree-unbounded proof enumerates
  all 715 possible no-singleton partitions of the eight active Jacobian
  contributions; every exponent-equality system forces repeated support or a
  linear monomial.  Thus any two-by-two sparse candidate must enter an axis
  or parallel-exponent boundary.  Those boundaries, larger supports, and
  JC(2) remain open; no counterexample is claimed.
source: jc-sparse-direct-search / 2026-08-22
depends_on: []
related:
  - THM-3418-one-monomial-nonlinear-fiber-keller-classification
  - THM-3544-planar-keller-target-pencil-total-degree-six-floor
  - THM-3551-one-ray-planar-jacobian-mate-no-go
  - THM-3587-consecutive-keller-fibre-factor-toric-and-coefficient-span-gates
script: 04-computation/jc2_fully_transverse_two_by_two_support_gate_thm3689.py
output: 05-knowledge/results/jc2_fully_transverse_two_by_two_support_gate_thm3689.out
script_sha256: 28676a46cd24f67e73d4fa9e748de20289fa53d31184a4d9971daf0d37206887
output_sha256: d9aa099dd13ca85839b7d9f4254f74da99c27fe002eee7283d1ee9e188e8cf2f
semantic_sha256: 1bde0faea0db3139bada3f6d73d35b9eb772be369431bbf532259fe81667bccf
hash_basis: raw LF bytes for files; ordered partition/solution/exit ledger for semantic hash
---

# THM-3689 -- fully transverse two-by-two sparse-support Keller gate

**PROVED + FINITE-EXACT + VERIFIED-EXACT.**  This theorem attacks the
smallest genuinely two-sided sparse chart directly in the source polynomial
ring.  It does not use a degree bound, Newton-area proxy, or failure to find a
coefficient solution: its obstruction occurs before the coefficient equations.

Work over a characteristic-zero field `K`.

## 1. Statement

Let

```text
A=(alpha_1,beta_1),       C=(alpha_2,beta_2),
B=(gamma_1,delta_1),      D=(gamma_2,delta_2)             (1)
```

be distinct pairs of nonnegative integer exponent vectors on each side, with

```text
|A|,|C|,|B|,|D| >= 2.                                    (2)
```

Take nonzero coefficients `a,c,b,d in K*` and put

```text
P=x+a x^alpha_1 y^beta_1+c x^alpha_2 y^beta_2,
Q=y+b x^gamma_1 y^delta_1+d x^gamma_2 y^delta_2.         (3)
```

Assume the **fully transverse** conditions

```text
alpha_1 alpha_2 delta_1 delta_2 !=0,
det(U,V)!=0 for every U in {A,C}, V in {B,D}.             (4)
```

Then

```text
Jac(P,Q) is not constant.                                 (5)
```

In particular, a planar Keller counterexample in the normalized two-by-two
support chart cannot be fully transverse.  It must lie on at least one of the
following boundary types:

```text
axis:       alpha_i=0 or delta_j=0;
parallel:   det(U,V)=0 for some cross pair;
support:    one side has fewer or repeated nonlinear monomials;
linear:     a displayed exponent has total degree at most one.              (6)
```

These are surviving **search cells**, not Keller solutions.

## 2. The eight labelled contributions

Subtract the constant contribution of `Jac(x,y)=1`.  Under `(2)--(4)`,
`Jac(P,Q)-1` has eight labelled nonzero contributions.  Their exponent
buckets are

```text
P_A: A-(1,0),                 P_C: C-(1,0),
Q_B: B-(0,1),                 Q_D: D-(0,1),

J_AB: A+B-(1,1),             J_AD: A+D-(1,1),
J_CB: C+B-(1,1),             J_CD: C+D-(1,1).           (7)
```

The four cross coefficients are respectively the nonzero products of
`a,c,b,d` with the determinants in `(4)`.  Characteristic zero and the first
part of `(4)` similarly make all four divergence coefficients nonzero.

If the Jacobian were constant, no bucket in `(7)` could contain exactly one
label: its nonzero coefficient would have nothing with which to cancel.
Consequently equality of exponent buckets would partition the eight labels
into blocks of size at least two.

## 3. Complete partition obstruction

There are exactly `715` set partitions of eight labelled objects with no
singleton block.  Their block-size profiles are

| profile | partitions |
|---|---:|
| `2+2+2+2` | 105 |
| `2+2+4` | 210 |
| `2+3+3` | 280 |
| `2+6` | 28 |
| `3+5` | 56 |
| `4+4` | 35 |
| `8` | 1 |

For each partition, equate the two coordinates in `(7)` inside every block
and solve the resulting affine-linear system over `Q` in the eight exponent
coordinates.  Every system exits the chart, with the exact exhaustive ledger

```text
A=C                                      676 partitions,
B=D                                       23 partitions,
|A|=1                                      8 partitions,
|C|=1                                      8 partitions.                 (8)
```

The four alternatives in `(8)` contradict respectively distinctness or
nonlinearity in `(1)--(2)`.  No sign, bounded exponent box, or coefficient
genericity is used.  This proves `(5)`.

The finite enumeration is a proof because the labels are fixed and an actual
constant-Jacobian identity produces exactly one of these set partitions.
Different blocks becoming accidentally equal merely gives a coarser
partition, which is already among the 715 cases.

## 4. Why the boundary is real

Dropping support cardinality immediately permits tame Keller maps.  For
arbitrary `alpha,delta in K`,

```text
P=x+alpha y^2,
Q=y+delta P^2
 =y+delta x^2+2 alpha delta x y^2+alpha^2 delta y^4       (9)
```

has `Jac(P,Q)=1`.  It is a triangular shear followed by a triangular shear,
not a counterexample.  Formula `(9)` is a one-by-three support control and
explains why a raw monomial count is not itself a nonproperness certificate.

For the direct counterexample search, `(6)` is the useful consequence.  A
next ansatz should retain an explicit collision or nonproper curve while
entering an axis or parallel boundary, or should add a third nonlinear
monomial on one side.  The theorem neither proves those boundary cells empty
nor supplies a collision.

## 5. Exact companion

The deterministic companion:

1. generates all set partitions by a restricted-growth recursion;
2. checks the profile count `715` independently of the exit ledger;
3. solves every equality system over exact rationals;
4. prints the actual consequence counts `(8)`;
5. verifies the positive tame control `(9)` and a fully transverse hostile;
6. rejects inactive Python `assert` statements.

Reproduce with

```bash
python3 -B 04-computation/jc2_fully_transverse_two_by_two_support_gate_thm3689.py
python3 -O -B 04-computation/jc2_fully_transverse_two_by_two_support_gate_thm3689.py
```

Both streams must agree byte-for-byte with the stored output.  This is an
all-exponent finite collision-pattern proof, not a bounded support census.
It proves no planar Jacobian counterexample and no case of `JC(2)` beyond the
displayed chart.  **QED.**

---
id: THM-4438
title: "JC2 row-fifteen relative response on the boundary Gm"
status: >
  PROVED FINITE-ROW RELATIVE TO THM-4426 + VERIFIED-EXACT + AUDITED.
  On THM-4426's rational boundary G_m over any characteristic-zero field,
  the prefix-preserving exact-valuation-fifteen source packet has rank-one
  image in the six-dimensional row-fifteen bracket cokernel. Its seven-free-
  parameter response kills all sixteen bracket coefficients, after which the
  only projected-depth debt is the inherited conic equation Q=0 and the
  terminal fibre is A^10. Within this packet the least-weight visible payer
  is p^3*y^6 of weight 24; weight 23 is invisible. This is not a statement
  about arbitrary lower-weight deformations, full source or termination,
  a Keller pair, JC(2), or DC(2).
source: row15_response_canon + root / next-sharp continuation, 2026-09-06
depends_on:
  - THM-4426-source-normal-row-fourteen-weight-eighteen-memory-repair
related:
  - THM-4410-source-normal-least-weight-twenty-row-thirteen-affine-continuation
  - THM-4415-source-normal-row-fourteen-same-row-response-rank-obstruction
script: 04-computation/jc2_source_normal_row15_relative_response_boundary_gm_referee.py
output: 05-knowledge/results/jc2_source_normal_row15_relative_response_boundary_gm_referee.out
script_sha256: be992bd9f3550c7a23334ed6d12f770d6422426412f17d3db2863ac55ce10374
output_sha256: a86ba2e50407989f7bc5784880c3267bbdb5d43a9710aadae25d1b119957b1cc
dependency_sha256: c647bdf2c894d60413b041197cff446884191bf1469372e6b4a6a55ae96fd52a
hash_basis: raw LF bytes
audit: >
  PASS. The hash-pinned extension retains all six raw left-nullspace
  compatibilities, computes exact rational-function ranks and nonzero constant
  pivot minors, and substitutes into every one of the sixteen bracket and
  ninety-one projected-depth coefficients. Normal, optimized, and fixed-hash-
  seed modes reproduce the frozen LF transcript in 259 live checks. The audit
  distinguishes the characteristic-zero positive lift from the rational-only
  irreducible-quartic hostile and does not infer arbitrary-weight leastness.
---

# THM-4438 -- JC2 row-fifteen relative response on the boundary Gm

**PROVED FINITE-ROW RELATIVE TO THM-4426 + VERIFIED-EXACT + AUDITED.**
The rational row-fourteen boundary continues through the row-fifteen bracket
and projected depth after one exact-valuation source response. `JC(2)` and
`DC(2)` remain **OPEN**.

## 1. Inherited boundary and statement

Work in THM-4426's source-normal chart. On

```text
Phi=eta=0,       alpha11=1,       c51=1087/135,        (1)
```

its simultaneous row-fourteen bracket/depth locus is a rational `G_m`,
parametrized by `s!=0`, with source coordinates `z=z(s), h=h(s)` satisfying
the exact conic equation

```text
Q(z,h)=0.                                               (2)
```

The row-fourteen terminal has ten free tangent coordinates. Adjoin the
complete prefix-preserving exact-valuation-fifteen packet

```text
R_15=sum_(b=0)^7 r_b p^(15-2b)y^b.                     (3)
```

For every characteristic-zero field `K`, every `s in K*`, and every
`r_1,...,r_7 in K`, there is a unique `r_0` for which all row-fifteen bracket
coefficients vanish. After that response all row-fifteen projected-depth
conditions vanish, and the terminal tangent fibre is `A^10`.

This statement is relative to the retained partial source and to the exact
boundary `(1)-(2)`. It does not assert that arbitrary source-normal data
enter this boundary.

## 2. Bracket quotient and response observable

Before `(3)` is added, the row-fifteen bracket matrix in the ten inherited
tangent variables has exact data

```text
shape 16 x 10,       rank 10,
pivots 0,...,9,      pivot minor 56056/390625.          (4)
```

Its raw cokernel therefore has dimension six. The inherited partial source
has two primitive debt representatives:

```text
P(z,h),       4Q(z,h).                                 (5)
```

The audit freezes the full quartic `P`; equation `(2)` removes only the
second debt. The response map of `(3)` into the six raw compatibilities has
rank one. Its single observable is

```text
Lambda=145r_0+30r_2+20r_4+24r_6,                      (6)
```

while `r_1,r_3,r_5,r_7` are bracket-invisible. Put

```text
C=10852621164972710686787843667734315747451565056000000000000000.
```

The exact response graph is

```text
r_0=-P/C-(6/29)r_2-(4/29)r_4-(24/145)r_6.             (7)
```

Substitution of `(7)` and the rational `G_m` parameterization kills all
sixteen raw bracket coefficients identically in
`Q(s,r_1,...,r_7)`. The constant pivot in `(4)` makes the conclusion stable
under base change to any characteristic-zero field with `s!=0`.

The use of all six raw compatibilities is essential. Conditions that were
proportional before the new variables were introduced can split after a
response is added; proportionality-deduplicated diagnostics are therefore
recorded but are not used to construct `(7)`.

## 3. Least weight inside the frozen-prefix packet

For a monomial in `(3)`,

```text
a+2b=15,       weight=2a+3b=30-b,       0<=b<=7.       (8)
```

Thus weights 15 through 22 do not occur in this exact-valuation diagonal.
The unique weight-23 term is `p*y^7`, namely the invisible variable `r_7`.
The first visible term is

```text
p^3*y^6,       valuation 15,       weight 24.          (9)
```

Setting every response except `r_6` to zero and taking

```text
r_6=-145P(z,h)/(24C)                                   (10)
```

makes `r_0=0` in `(7)` and kills the full bracket. This is a distinguished
one-coordinate section with constant rational denominator before the
`G_m` substitution.

The word "least" in `(9)` is strictly relative to the prefix-preserving
packet `(3)`. A lower-weight deformation entering an earlier row changes the
frozen prefix and must be re-solved; this theorem neither excludes nor
classifies such deformations.

## 4. Projected depth

Append the sixteen row-fifteen tangent coordinates after the bracket
response. The projected Hasse-depth selector has

```text
shape 91 x 16,       rank 6,
pivots 10,...,15,    pivot minor 1/64.                 (11)
```

Its sole primitive source condition is exactly

```text
-4Q(z,h)=0.                                             (12)
```

It is independent of all seven relative response parameters. Equation `(2)`
therefore kills all ninety-one raw depth residuals coefficientwise over
`Q(s,r_1,...,r_7)`, and the ten unpivoted tangent coordinates give the
terminal fibre `A^10`. The section `(10)` was separately substituted into
all sixteen bracket and all ninety-one depth coefficients.

## 5. A no-response rational hostile

Without `(3)`, pull back `P` to the rational parameter. It is a nonzero
rational multiple of

```text
N(s)/s^2,                                               (13)
```

where `N` is the primitive degree-four polynomial frozen in the audit output.
Exact factorization over `Q` leaves one irreducible quartic factor. Hence

```text
every s in Q* fails the row-fifteen bracket without a new response. (14)
```

This is only a rational hostile. An extension field can contain a root of
`N`; no extension-field no-response claim is made. Conversely, the positive
response `(7)` is a rational identity and works over every characteristic-
zero field.

## 6. Reproduction and scope

```powershell
python -B 04-computation/jc2_source_normal_row15_relative_response_boundary_gm_referee.py
python -B -O 04-computation/jc2_source_normal_row15_relative_response_boundary_gm_referee.py
$env:PYTHONHASHSEED='161803'; python -B 04-computation/jc2_source_normal_row15_relative_response_boundary_gm_referee.py
```

The script is hash-pinned to THM-4426's two-memory companion. It performs
259 optimization-live exact checks and emits an 84-line frozen LF transcript.
The full quartic debts, pivot ledgers, response graph, three rational controls,
and coefficientwise zero tests are documented in
`05-knowledge/results/jc2_source_normal_row15_relative_response_boundary_gm_audit.md`.

This theorem closes only the row-fifteen bracket and projected depth on the
named partial-source boundary. Complete lower-weight source deformations,
the global row-fourteen depth cover, later rows, finite termination, chart
entry, full `B_2`, a Keller pair, `JC(2)`, and `DC(2)` remain **OPEN**.

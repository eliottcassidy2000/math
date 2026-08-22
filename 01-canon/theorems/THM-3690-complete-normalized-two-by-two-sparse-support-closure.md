---
id: THM-3690
title: "Complete normalized two-by-two sparse-support closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In
  the displayed normalized planar chart P=x plus exactly two distinct
  nonlinear monomials and Q=y plus exactly two distinct nonlinear monomials,
  Jac(P,Q) is never constant.  There are no axis or parallel-exponent
  exceptions: an exact activity-mask, partition, and saturated Groebner
  computation closes all of them in arbitrary degree.  This is a support-chart
  obstruction, not a proof of JC(2), and no counterexample is claimed.
source: jc-sparse-direct-search / 2026-08-22
audit: >
  PASS -- root independently checked exact label activity, the no-singleton
  necessity, the four witness charts, every saturated nonzero condition, the
  affine extraction, and the final nonnegative-integral lattice contradiction.
depends_on: []
related:
  - THM-3418-one-monomial-nonlinear-fiber-keller-classification
  - THM-3544-planar-keller-target-pencil-total-degree-six-floor
  - THM-3551-one-ray-planar-jacobian-mate-no-go
  - THM-3587-consecutive-keller-fibre-factor-toric-and-coefficient-span-gates
  - THM-3689-fully-transverse-two-by-two-sparse-support-gate
script: 04-computation/jc2_complete_two_by_two_support_gate_thm3690.py
output: 05-knowledge/results/jc2_complete_two_by_two_support_gate_thm3690.out
script_sha256: bc36bcb6f24e7876c492f188647ac1e1a877e21875311c2c2d49f55225f39260
output_sha256: 8d7bb6239c8e98261c1192df0f724b3819a0cf29c549cf348bcc670a2d3b9af8
semantic_sha256: 24868754ee5dd4281469f6a195617624d58dabb15fb32acc4caba9b1f0be6e59
hash_basis: raw LF bytes for files; affine-exit ledger, reduced-basis certificates, and invalid exponent packets for semantic hash
---

# THM-3690 -- complete normalized two-by-two sparse-support closure

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This theorem closes the axis and parallel-exponent boundary left open by
THM-3689.  The obstruction is again in the source support itself, before any
coefficient-cancellation equations need to be solved.

Work over a characteristic-zero field `K`.  Write `X^(r,s)=x^r y^s`,
`e_1=(1,0)`, `e_2=(0,1)`, and `|(r,s)|=r+s`.

## 1. Statement

Let `A,C,B,D` be nonnegative integer exponent vectors satisfying

```text
A != C,                 B != D,
|A|,|C|,|B|,|D| >= 2.                                  (1)
```

For arbitrary nonzero `lambda,mu,nu,omega in K*`, put

```text
P=x+lambda X^A+mu X^C,
Q=y+nu X^B+omega X^D.                                  (2)
```

Then

```text
Jac(P,Q) is not constant.                               (3)
```

No transversality is assumed.  In particular, coordinates of the exponent
vectors may vanish and any subset of the four cross determinants may vanish.

The normalization in `(2)` is part of the statement.  General linear source
or target changes can enlarge or merge displayed supports, so `(3)` must not
be advertised as a coordinate-invariant monomial lower bound.  It proves no
case of the full planar Jacobian conjecture beyond this exact chart.

## 2. Eight potential contributions and exact activity

Subtract the contribution `Jac(x,y)=1`.  The eight labelled potential terms of
`Jac(P,Q)-1` have exponent buckets

```text
P_A: A-e_1,                       P_C: C-e_1,
Q_B: B-e_2,                       Q_D: D-e_2,

J_AB: A+B-e_1-e_2,               J_AD: A+D-e_1-e_2,
J_CB: C+B-e_1-e_2,               J_CD: C+D-e_1-e_2.    (4)
```

Their coefficient factors, after removing the nonzero displayed coefficients
in `(2)`, are respectively

```text
A_1, C_1, B_2, D_2,
det(A,B), det(A,D), det(C,B), det(C,D).                  (5)
```

Thus a label is inactive **if and only if** its corresponding factor in `(5)`
vanishes.  This uses both characteristic zero and the hypothesis that all four
source coefficients are nonzero; no generic-coefficient assumption is hidden.

Every active divergence term in `(4)` has total degree at least one and every
active cross term has total degree at least two.  Hence no nonlinear term can
alter the scalar contribution `1`.  If the Jacobian were constant, it would be
exactly `1`, and every occupied bucket in `(4)` would contain at least two
active labels.  The equality classes of active buckets therefore form a set
partition with no singleton block.

## 3. Complete activity/partition pass

There are `16` possible activity masks for the four divergence labels and `16`
for the four cross labels.  Enumerating every no-singleton partition of each
active label set gives exactly

```text
4,140 activity/partition systems.                       (6)
```

For each system, impose the inactive divergence equations and all two-coordinate
bucket equalities inside each partition block.  An exact affine-linear solve
eliminates `4,043` systems.  The ordered exit ledger is

```text
inconsistent                                      1,698
A=C                                                1,805
B=D                                                  284
|A| forced below 2                                   128
|C| forced below 2                                    96
|B| forced below 2                                    16
|D| forced below 2                                    16.               (7)
```

Every exit in `(7)` contradicts `(1)`.  The ordering only assigns a unique
diagnostic to a system that may have more than one defect; it is not used as a
logical exclusivity claim.  Exactly `97` affine residual systems remain.

## 4. Boundary saturation without unsafe division

For each residual system, impose every inactive cross determinant as an
equation.  The required open conditions are encoded by one product `S`:

```text
- every declared-active divergence coordinate;
- every declared-active cross determinant;
- |U|(|U|-1) for U in {A,C,B,D};
- one nonzero coordinate witness for A!=C;
- one nonzero coordinate witness for B!=D.             (8)
```

The two choices of witness on each side give four principal-open charts.  Their
union is exactly the distinct-support condition in `(1)`.  For nonnegative
integer exponents, `|U|>=2` is equivalent to `|U|(|U|-1)!=0`, so the third line
of `(8)` loses no admissible integer point.  It deliberately admits extra
complex points; those are removed only at the final lattice check.

No division by a possibly zero factor occurs.  A new variable `s` and the
polynomial equation

```text
s S - 1 = 0                                             (9)
```

encode each principal open.  An identically zero `S` makes that open empty and
is skipped.  The exact ledger is

```text
97 residual systems x 4 witness charts = 388 charts,
19 identically empty charts,
369 exact saturated Groebner computations,
64 nonunit reduced bases in 48 basis shapes.            (10)
```

Each of the `64` nonunit reduced bases is zero-dimensional and, more strongly,
is a square affine-linear system in its remaining parameters and `s`.  The
companion verifies this property and solves each basis with exact linear
algebra.  It then substitutes the point back into `(9)`, every inactive
determinant equation, the activity mask, every bucket equality, the selected
distinctness witnesses, and all nonlinear-total conditions.  Thus no
heuristic nonlinear root extraction or unverified saturation division enters
the proof.

The `64` opens collapse, with uniform multiplicity four from the witness
charts, to the following `16` exponent packets.  Coordinates are ordered as
`(A_1,A_2,C_1,C_2,B_1,B_2,D_1,D_2)`:

```text
(-1,2/3, 0,1/3, -2,5/3, -3,2)
(-1,2/3, 0,1/3, -3,2,   -2,5/3)
(-1,3,  -1/3,2, -2/3,2,  2/3,0)
(-1,3,  -1/3,2,  2/3,0, -2/3,2)
(-1/3,2, -1,3,  -2/3,2,  2/3,0)
(-1/3,2, -1,3,   2/3,0, -2/3,2)
(0,1/3, -1,2/3, -2,5/3, -3,2)
(0,1/3, -1,2/3, -3,2,   -2,5/3)
(0,2/3,  2,-2/3, 2,-1/3, 3,-1)
(0,2/3,  2,-2/3, 3,-1,   2,-1/3)
(2,-2/3, 0,2/3,  2,-1/3, 3,-1)
(2,-2/3, 0,2/3,  3,-1,   2,-1/3)
(2,-3,   5/3,-2, 1/3,0,  2/3,-1)
(2,-3,   5/3,-2, 2/3,-1, 1/3,0)
(5/3,-2, 2,-3,   1/3,0,  2/3,-1)
(5/3,-2, 2,-3,   2/3,-1, 1/3,0).                       (11)
```

Every row of `(11)` has a negative or nonintegral coordinate.  Consequently
none lies in `(Z_{>=0}^2)^4`, and all `97` residual systems are incompatible
with the exponent hypotheses `(1)`.  Together with `(7)`, this proves `(3)`.

## 5. Independent finite census and controls

As an independent representation check, the companion directly enumerates all
`42` nonlinear monomials of total degree at most eight, then all

```text
binom(42,2)^2 = 741,321
```

ordered left/right support-pair pairs.  It constructs the active buckets in
`(4)` directly, without masks, partitions, affine solving, Groebner bases, or
saturation.  No pair has all bucket multiplicities at least two.  This bounded
census is a check on the degree-unbounded proof, not its source.

Support drop is a real boundary.  For arbitrary `alpha,delta in K`,

```text
P=x+alpha y^2,
Q=y+delta P^2
```

has `Jac(P,Q)=1`; it is a tame one-by-three control.  Conversely, an explicit
collision alone does not discharge the Jacobian debt.  The two-by-two map

```text
P=x+a x^2-x^3,
Q=y+b x^2+d x^4
```

maps `(1,0)` and `(-1,0)` to the same point, but

```text
Jac(P,Q)=1+2 a x-3x^2,
```

so it is not Keller.  These controls keep source-support cancellation and
noninjectivity as separate obligations.

## 6. Reproduction and consequence

Run

```bash
python3 -B 04-computation/jc2_complete_two_by_two_support_gate_thm3690.py
python3 -O -B 04-computation/jc2_complete_two_by_two_support_gate_thm3690.py
```

Both streams agree byte-for-byte with the stored output.  The script uses exact
rational arithmetic, checks all count ledgers and all `16` final lattice
contradictions, and rejects inactive Python `assert` statements.

The useful search consequence is precise: in the normalized two-sided branch
where each component has at least two nonlinear monomials, a candidate must
leave the two-by-two cell, so at least one side needs a third nonlinear support.
One-sided support drops remain separate branches.  No nonproper Keller map is
constructed here.  **QED.**

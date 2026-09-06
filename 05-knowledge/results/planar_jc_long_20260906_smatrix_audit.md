# Independent audit of raw dual testers and integral Student responses

**Status: ACCEPTED analytically and by an independent exact implementation.**
This audits the scoped theorems in
[planar_jc_long_20260906_smatrix.md](planar_jc_long_20260906_smatrix.md).
No characteristic-zero Keller obstruction or full depth-compatible source
lift is inferred from the integral lattice calculation.

## 1. Analytic audit

**The generic-dual coefficient ideal is correct over every commutative
ring with identity.** Cauchy--Binet differentiated with respect to each
entry of `A` gives the displayed adjugate identity as an integer polynomial
identity. For fixed increasing row set `I`, a row outside `I` of
`det(A_I)B-A adj(A_I)B_I` is its ordered bordered determinant; rows inside
`I` vanish by the adjugate identity. Different `I` give disjoint monomial
supports in the independent entries of `U`. Each such monomial has
coefficient `+1` or `-1`, so no factorial, multiplicity, cancellation, or
division changes the coefficient ideal. Repeated rows give the expected
zero bordered determinant. The ideal is exactly the mixed maximal-minor
ideal for the complete displayed bank.

The membership qualifications are essential and correctly stated. Over a
field, full column rank of `A` makes vanishing mixed minors equivalent to
membership of every bank column in the image. On a localization with an
invertible maximal minor, the pivot inverse gives a unique solution and
every remaining raw row is one of these residual equations. Over a
general ring there need not be such a chart, and over `Z` saturation can
be strictly larger than the image. The examples `A=(2,0)^T, B=(1,0)^T`
and `A=0, B!=0` retain these two different failures. The identity does
not enforce a shared parameter omitted when assembling the raw operator.

**The all-index Student Smith theorem is correct for its declared
lattices.** I independently recovered the response from
[THM-4308 / source-normal-bracket-hasse-truncation-through-row-eight](../../01-canon/theorems/THM-4308-source-normal-bracket-hasse-truncation-through-row-eight.md),
equation (13): with `q=-(x^2+6)/2`, its action is
`-x theta+(x^2+6)theta'/(2m)`. Its monomial coefficients are exactly
`6j` and `j-2m` after multiplying the target coordinates by `2m`.
This is the unrestricted tangent-response operator at arbitrary index;
it supplies no existence of all the original source rows or satisfaction
of their depth restrictions.

For each parity, the nonzero matrix support is the stated bipartite path.
All edge weights are nonzero for `0<=j<m`. A square submatrix has at most
one perfect matching because two different perfect matchings would make
an alternating cycle in a forest. Conversely any matching picks distinct
rows and columns and is the unique perfect matching of that submatrix.
Thus its product is an actual nonzero minor, including its possible sign;
there are no hidden sums with cancellations. The two-path gcd recurrence
therefore computes every determinantal divisor. A full-column matching
exists, and ordinary integer Smith theory gives the displayed invariant
factors and one free cokernel coordinate. The operation bound is `O(m^2)`
arithmetic operations, as claimed, without a bit-complexity assertion.

**The saturation and common-denominator statements are also correct.**
The primitive left covector makes its kernel the saturation of the image.
The quotient is exactly the torsion part of the integer Smith cokernel.
Injectivity of the response makes the rational preimage unique. The order
of a compatible integral target class is its least rational-solution
denominator; a bank requires the least common multiple of these orders.
Consequently the sharp universal denominator is the largest invariant
factor, and integral membership is equivalent to vanishing of the target
classes. This interpretation uses the declared target lattice
`(1/(2m))Z[x]_(<=m)`; silently changing that lattice would change the
arithmetic question. Tensoring with a characteristic-zero field removes
all of this torsion and leaves the inherited single Student equation.

The row-14 seven-column bank really attains the full exponent `720720`.
It is compatible in the unrestricted response and does not repair a
restricted terminal-row or memory condition. The prime-support bound is
also valid: one selected maximal matching product has only primes at most
`max(m,3)`, and every determinantal divisor divides that product.

## 2. Cited-paper scope repair

The initial draft described the S-matrix budget immediately after the
condition that entries of `B` lie in `[0,1]`. I checked the actual primary
[v1 paper, Sections 2, 4.2, and 6.2--6.3](https://arxiv.org/html/2608.29750v1).
The budget also assumes even `n>2`, nonsingularity, and the contradiction
bound `||B^(-1)||_F<=2n/(n+1)`. The author added these hypotheses before
audit acceptance and corrected the section references. The constants
and the shared-column interpretation agree with those sections. This
audits the stated transfer and its scope, not the paper's complete proof
or its Lean development. No positive real energy is transferred to `C`.

## 3. Independent exact engine and finite universe

The standalone [audit source](../../04-computation/planar_jc_long_20260906_smatrix_audit.py)
imports no research producer. Its main integer engine uses elementary
Euclidean row and column operations. After clearing each pivot row and
column, any trailing entry not divisible by the pivot is moved into that
row, causing a subsequent Euclidean pivot decrease. This constructs the
entire Smith diagonal without calling the producer's path-product gcd
algorithm or a Smith-normal-form library. Connected support components
are extracted directly from the actual matrix to reduce intermediate
integer growth; their diagonal forms are combined by elementary `gcd/lcm`
operations. An undecomposed matrix calculation independently agrees for
every `2<=m<=14`.

For every `2<=m<=48`, it independently constructs the actual monomial
matrix and compares this integer diagonal with the theorem's formula
evaluated by prime valuations and minimum-cost nonadjacent selections.
For every `2<=m<=10`, it additionally enumerates every support matching
of the undecomposed matrix by choosing source columns and distinct target
rows. This third path uses neither parity paths nor their recurrence.
The forest proof identifies the resulting gcds with actual minors.

Four generic-dual size/seed controls recover coefficients directly from
ordered bordered determinants. The source checks every such coefficient
and excludes extra monomial support. It includes explicit rank-drop,
integral, zero-divisor-ring, and complex isotropic hostiles. Dense rational
nullspace reduction independently recovers the four inherited Student
covectors. Dense Gaussian solves verify all equations of the actual
two-column `m=2` and seven-column `m=14` banks, including the individual
denominators, their least common multiple, the full row-14 torsion order,
and the localized exponent `715`.

No finite control bank is being promoted to an all-index proof: the
all-index conclusion rests on the audited forest and Smith arguments.

## 4. Reproduction and manifest

```bash
python3 -B 04-computation/planar_jc_long_20260906_smatrix_audit.py > /tmp/smatrix_audit_normal.out
python3 -B -O 04-computation/planar_jc_long_20260906_smatrix_audit.py > /tmp/smatrix_audit_optimized.out
cmp /tmp/smatrix_audit_normal.out /tmp/smatrix_audit_optimized.out
cmp /tmp/smatrix_audit_normal.out 05-knowledge/results/planar_jc_long_20260906_smatrix_audit.out
```

Normal and optimized runs pass **13,868 always-active exact gates** and
agree byte for byte with the [frozen audit output](planar_jc_long_20260906_smatrix_audit.out).
The mathematical audit accepts all three new algebraic claims; the only
requested repair concerned the external citation hypotheses.

Raw LF-byte SHA-256 manifest:

```
audited producer source
17751f9596a6be43a7fd8757350fee679347eb513a0b318f7b52cd5867128f6d
audited producer output
40b41ace5ffc917124bbe2995ad8c5a8bc0a85ec2faf8461e50b97609870d7e2
independent audit source
250fcc4636ec73baaa7ab8c029fd10e267efe84f07f1867cf8e98cf068c14d51
independent audit output
49475a11995b66eb99d50d7dcd261fab6c2e65949e175b60f84e60b2b8817352
```

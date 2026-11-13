# Independent review of the observer and Clifford boundaries

**Verdict: PASS after the two scope repairs recorded below.** This review
certifies the elementary arguments at the exact hashes below; it does not
promote the root's RESERVED statements or any finite census. No root producer
file was edited by the reviewer. Base commit:
`a8b6d8dbf0c6b601813785efb0b5904b46c23d43`.

| Reviewed packet file | SHA256 |
|---|---|
| `ROOT_OBSERVERS.md` | `b0b4508aa28b76bd05d4f9e006d4139ad44267afbad98a38f8f16313b644d6a1` |
| `CLIFFORD_BOUNDARY.md` | `ee27ca8a5bd08b7db309286052c63e9a41ef215e5820e4bbe402c5a79b77c1f4` |

## Observer arguments O1–O3

**O1 passes with its corrected equality hypothesis.** A real unital subalgebra
of functions on a finite set separates its own equivalence classes. Finite
interpolation then puts every class indicator in the algebra. This proves the
partition/algebra identification, and invariance of the indicator basis is
equivalent to constancy of all target-block row sums. Intersections of the
specified invariant algebras remain unital, pointwise multiplicative, and
invariant. Hence the minimum exists and equals the coarsest stable refinement
of the initial block algebra. Closure and the poset-reflector language follow;
no extra site structure follows from that universal property.

The first reviewed draft merely required `V0` to **contain** the initial block
indicators. That did not justify its conclusion about closure of `Pi0` alone.
For example, take three points, the one-block partition, the identity operator,
and `V0=span{1,(0,1,2)}`. The partition is already stable, but the algebra
containing `V0` is the full three-point algebra by interpolation. The repaired
text sets `V0` **equal to** the initial block-constant algebra. Another valid
repair would have begun with the joint level-set partition of all of `V0`.
This is a hypothesis repair, not a failure of the finite-algebra mechanism.

The distinction from invariant linear spans is correct: the pointwise square
of `(0,1,2)` lies outside its span with the constant function. A minimum number
of allowed measurement columns is a third optimization problem. None of these
three minima should be substituted for another.

**O2 passes.** With nonempty partition blocks and positive weights, `M` is
invertible, `LE` is the weighted orthogonal projection, and `R` is its
complement. The exact defect is `EPRPL`. Multiplication by `M`, reversibility,
and `R^T D=DR` identify it with `B^T D B`; positivity makes vanishing equivalent
to `RPL=0`, which is strong lumpability. For `G=P-I`, the identities `RL=0`
and `ER=0` give the same defect at the second derivative. There is no Taylor
factor of one half in that derivative formula.

For the path example, the degree measure is reversible and positive. After
`k` strict rounds, only the first vertex of the remaining tail can reach the
previous singleton in one step. This proves the displayed partition induction
and exactly `m-1` strict rounds for every `m>=2`. The actual two-step return
is `1/2`, while the degree-weighted quotient gives `1/(2m-1)`, so the
second-order defect is positive already at the initial partition. Detection
depth and repair depth are different predicates.

**O3 passes with the narrowed consequence.** The displayed word-intertwining
identity follows by induction from each individual matrix intertwiner, and
the formal resolvent closes fixed-order walks. It says nothing by itself
about construction at the next tournament order. During review the empirical
discreteness consequence was narrowed to **strongly lumpable partition
quotients refining the specified initial observer**. That is precisely what
the coarsest stable-partition theorem rules out. A different linear realization
of selected word statistics is not excluded. The repaired prior-recurrence
link resolves to the final section of
[the preceding packet's README](../../joint_metagraph_dyadic_20261113/README.md).
The present review does not re-audit that inherited earlier-values recurrence.

## Clifford arguments C1–C3

**C1 passes.** An involution provides an ungraded representation of the
algebra with one positive-square generator. It also provides a grading, but
the involution itself is even relative to that grading. The text correctly
avoids claiming a graded Clifford module with that involution as an odd
generator. An invertible anticommuting map exchanges the two eigenspaces,
forcing equal dimensions. Conversely the displayed off-diagonal matrix has
square `-I` and supplies the `Cl_(1,1)` extension when dimensions agree.

For the ordinary-class permutation representation the difference is exactly
the number `s` of fixed classes, and a transitive tournament supplies at least
one. Adding `k_+` positive and `k_-` negative dimensions balances exactly when
`k_- - k_+ = s`; minimum total addition is `s`, attained by adding only negative
dimensions. Trivial rank-sixteen amplification multiplies the imbalance.
These are finite-dimensional linear statements, not new tournament classes,
a class-pairing operation, or a periodicity theorem. The general Clifford
facts are classical elementary algebra; the packet's contribution is the
explicit application and the retained incidence obstruction, not priority
for the algebraic criterion.

**C2 passes.** In balanced coordinates any invertible odd map has two
invertible off-diagonal blocks. Requiring its square to be `-I` gives the
stated form with `B` and `-B^-1`. Commutation with the block-diagonal incidence
operator is equivalent to `B C_+ = D B`; the other block equation follows
by inversion. Thus similarity is both necessary and sufficient. On the
symmetric real sector operators it could equivalently be tested by spectra
with multiplicity, but no such computational test is assumed here. Merely
adding the right number of dimensions does not supply the needed operator.
For example, `C_+=diag(1,2)` and `D=(3)` cannot be made similar by adding only
one negative-sector operator: the eigenvalue 3 remains. The displayed balanced
`diag(1,2)` two-sector hostile is already enough for the claimed boundary.

**C3 passes.** In the original ordinary-class basis, both products with the
nonnegative converse permutation are nonnegative. Their sum can vanish only
if both vanish, and invertibility of the converse permutation forces `J=0`.
This basis matters; diagonalizing the converse changes the nonnegative cone.
The incidence-preserving involutions counted in the signed-lift lane commute
with converse and are not these anticommuting Clifford operators. The shared
number `s` is an exact fixed-point/character dimension in the stated
representation, without a consequence about shifts in tournament order.

## Independent exact controls

[check_root.py](check_root.py) uses rational matrices and imports no producer
code. It independently tests:

- 2,560 pairs of a binary three-by-three matrix and a three-point partition,
  class-indicator interpolation, the row-sum/invariant-algebra equivalence,
  and word intertwiners whenever the partition is stable;
- path walks at `m=2,...,14`: exact Gram defects, derivative normalization,
  rejection of a mistaken Taylor factor, both return probabilities, and exact
  strict refinement counts;
- sector dimensions zero through four, the negative-square odd-map criterion,
  minimum added dimensions for `s=0,...,8`, an explicit nontrivial incidence
  intertwiner, the balanced spectral hostile, and all 512 nonnegative binary
  three-by-three matrices for a converse permutation with one fixed point.

All 3,605 gates pass in normal and optimized Python. These controls are
**FINITE-EXACT**; the all-order linear statements rely on the proofs above.
Both receipts bind the reviewed root-note hashes and are byte-identical.

```sh
python3 05-knowledge/results/incidence_closure_20261113/audit_bridges/check_root.py > 05-knowledge/results/incidence_closure_20261113/audit_bridges/root.normal.json
python3 -O 05-knowledge/results/incidence_closure_20261113/audit_bridges/check_root.py > 05-knowledge/results/incidence_closure_20261113/audit_bridges/root.optimized.json
```

| Independent artifact | SHA256 |
|---|---|
| `check_root.py` | `6840b0d54349278001436d7e6d216c1b0b5085f54271c9a3145f1bbdfe5a646c` |
| `root.normal.json` | `42113913c47ce7eee3d70fd5a9037379ab593179dcd0ad4fc7e0dfb7c12dfc5c` |
| `root.optimized.json` | `42113913c47ce7eee3d70fd5a9037379ab593179dcd0ad4fc7e0dfb7c12dfc5c` |

These are audit snapshots. Any later producer edit, including a status-only
promotion, needs an explicit hash rebind before the guard in the checker will
accept it as the reviewed input. This prevents a replay from silently
certifying unreviewed text.

## Status-only promotion rebind, 2026-11-13

The preceding report is the original prepromotion snapshot. Root has promoted
O1–O3 and C1–C3 to **PROVED IN NOTES**. Independently reversing each note's one
declared status edit in `../PROMOTIONS.json` recovers its exact reviewed hash:
`b0b4508aa28b76bd05d4f9e006d4139ad44267afbad98a38f8f16313b644d6a1` for
`ROOT_OBSERVERS.md`, and
`ee27ca8a5bd08b7db309286052c63e9a41ef215e5820e4bbe402c5a79b77c1f4` for
`CLIFFORD_BOUNDARY.md`. No mathematical text changed. Thus this explicit
rebind carries the existing PASS verdict onto the promoted bytes; references
to RESERVED status above describe the earlier review stage.

| Final promoted note | SHA256 |
|---|---|
| `ROOT_OBSERVERS.md` | `bbe31d87ba4a7d328add13194af9dde80f68a2d6d03bcadb5de3e47bce4063c7` |
| `CLIFFORD_BOUNDARY.md` | `6615b63e6e7c1bb692facbe99ade8a925b96cebb963ecdd8305d56f13ca7fa10` |

The manifest used has SHA256
`dbccf7ced4fe064dbe908d196ebdf507542018b2c8acf74ff57fdcc84db1a898`.
[promotion_verification.json](promotion_verification.json), SHA256
`d69e60772899562bb8c53f3ea2f5789077d375952e36d3cd493fea6097ac703c`,
records both exact inverse checks and the separate bridge status rebind.

Original receipts are preserved as `root.normal.prepromotion.json` and
`root.optimized.prepromotion.json`, both SHA256
`42113913c47ce7eee3d70fd5a9037379ab593179dcd0ad4fc7e0dfb7c12dfc5c`.
The earlier artifact table refers to that preserved snapshot. Only the two
expected root-note hashes changed in `check_root.py`. Normal and optimized
reruns again pass all 3,605 gates and agree byte for byte. Removing just the
`reviewed_sha256` and `source_sha256` provenance fields leaves exactly the
same JSON mathematical payload as the archived prepromotion receipts.

| Current rebound artifact | SHA256 |
|---|---|
| `check_root.py` | `b1d8ff1934dfa3d2b981ef53c95e9e10284307de6ba95ff2ea8a246346a9419e` |
| `root.normal.json` | `356721a1ff6aace1b71f7517b246b122cc3b664214be9b7ce2763bea58cd892d` |
| `root.optimized.json` | `356721a1ff6aace1b71f7517b246b122cc3b664214be9b7ce2763bea58cd892d` |

# Independent audit: a bounded Frobenius-cube interpretation

**Status:** FINITE-EXACT point counts and classifications; elementary deductions
PROVED below, with the two standard elliptic-curve facts clearly CITED. This is
a bounded Session143 follow-up, not a new all-extension investigation.

**Main outcome:** there is a genuine elliptic-curve realization of the surviving
numerical cubic `t^3−15t`. It is neither a magic square of squares nor evidence
that the user's unspecified “lambda cube” meant Frobenius. That referent remains
**UNRESOLVED**. No octonion, tournament, or Hecke-equivariant map is supplied.

## 1. Exact universe and direct count

Let `E_ab/F5` be `y²=x³+ax+b`, with `a,b=0,...,4` and
`4a³+27b² != 0 (mod 5)`. There are exactly 20 such models. The standard
short-Weierstrass theorem says every elliptic curve over F5 has such a model,
and all origin-preserving isomorphisms between these models are the scalings
`(a,b) -> (u^4 a,u^6 b)`, `u in F5*`.
[CITED: Milne, *Elliptic Curves*, II Theorem 2.1, pp. 50–51](https://jmilne.org/math/Books/EC2.pdf#page=55).

The audit counts directly in

- `F5`;
- `F25=F5[z]/(z²+2)`;
- `F125=F5[z]/(z³+2z+1)`.

The cubic has evaluations `1,4,3,4,3` on F5, hence has no root and is
irreducible. The quadratic has no root because 3 is not a square in F5.
All field elements are coefficient tuples. The checker independently verifies
`x^q=x` and a multiplicative inverse for every nonzero element. It forms the
histogram `R(c)=#{y:y²=c}` and computes

`N(E/Fq)=1+sum_x R(x³+ax+b)`.

Thus N3 is obtained by direct finite-field enumeration, **not inferred from
a recurrence**. The point at infinity is included explicitly. F5 counts have
a separate 25-pair `(x,y)` check. This is a different cubic extension polynomial
from the producer's implementation.

The following table is the independently obtained count data. Entries with
`±b` name two separate models; `b=0` names one.

| a | b | t=6−N1 | N1 | N2 | N3 |
|---:|:---:|---:|---:|---:|---:|
| 0 | ±1, ±2 | 0 | 6 | 36 | 126 |
| 1 | 0, ±2 | 2 | 4 | 32 | 148 |
| 1 | ±1 | −3 | 9 | 27 | 108 |
| 2 | 0 | 4 | 2 | 20 | 122 |
| 2 | ±1 | −1 | 7 | 35 | 112 |
| 3 | 0 | −4 | 10 | 20 | 130 |
| 3 | ±2 | 1 | 5 | 35 | 140 |
| 4 | 0, ±1 | −2 | 8 | 32 | 104 |
| 4 | ±2 | 3 | 3 | 27 | 144 |

In particular all nine labels `N1=2,...,10` occur. These labels, or the centered
traces `t=-4,...,4`, are scalar quotients: they do not identify the curve class.

## 2. What the cube means in this construction

For Frobenius on the two-dimensional rational Tate module, let its eigenvalues
be `α,β`. The standard point-count theorem gives
`α+β=t`, `αβ=5`, and `N_n=5^n+1−α^n−β^n`.
[CITED: Snowden, Lecture 3, “Point counting”](https://public.websites.umich.edu/asnowden/teaching/2013/679/L03.html).
The same formula appears in [Milne IV Remark 9.11(a), p. 156](https://jmilne.org/math/Books/EC2.pdf#page=161).

An elementary expansion now yields

`α³+β³=(α+β)³−3αβ(α+β)=t³−15t`,

hence `N3=126−t³+15t`. All 20 direct counts agree. This is a cube of a
rank-two operator, or cubes of its two eigenvalues; it is **not** simply the
cube `t³` of the trace.

For comparison,

`tr(Sym³ F)=α³+α²β+αβ²+β³=t³−10t`.

The checker computes both traces from the explicit companion matrix
`F=[[0,−5],[1,t]]`; the symmetric cube is built by substitution on the four
cubic monomials, not by inserting the desired formula. Their difference is
`5t`, so the traces coincide in this packet only at `t=0`. If a normalized
weight-two good-prime Euler factor is separately supplied, the coefficient of
`T³` in `(1−tT+5T²)^−1` is also `t³−10t`. This conditional formal-series
observation does not identify the point-count power trace with a Hecke
coefficient, construct a global modular form, or assert a Hecke map. The
normalization and good-prime Euler factor are independently checked against
[Sutherland, Lecture 24, Theorem 24.27](https://math.mit.edu/classes/18.783/2025/LectureNotes24.pdf).

The producer also states a valid formal-family corollary. For a common real
determinant q and the nine traces `t=-4,...,4`, the power-cube polynomial is
`t³−3qt` and the symmetric-cube polynomial is `t³−2qt`. Substitution into the
pinned all-real cubic classification gives magic-arrangeability exactly at
`q=5` and `q=15/2`, respectively. This is an immediate conditional-on-the-
inherited-theorem deduction, not a search over fields or a realization claim
for every q. The F5 model census supplies the separate arithmetic realization
of the first case.

## 3. Actual magic witness and its failure boundaries

Select one of the models for each trace using this coefficient grid:

```
(1,1) (1,0) (2,0)
(3,2) (0,1) (2,1)
(3,0) (4,0) (4,2)
```

Their traces are

```
−3   2   4
 1   0  −1
−4  −2   3
```

and their F125 point counts are

```
108  148  122
140  126  112
130  104  144
```

These are nine distinct integers. All three rows, all three columns, and both
diagonals sum to 378. Independently enumerating all 8! placements of the other
entries after forcing the center gives exactly eight arrangements. The same
complete check gives zero arrangements for the distinct values of `t³−10t`
and of `t³` on `t=-4,...,4`. This explicitly realizes the affine cubic family proved
in Session143's `magic_cycle/FULL_CUBIC_CLASSIFICATION.md` (pinned copy under
`sources/`). The present finite witness does not depend on rerunning that full
classification.

The center is 126, not an integer square. Therefore this is **not** a magic
square of squares. The original N1 labels in these same cells are

```
9  4  2
5  6  7
10 8  3
```

with row sums `15,18,21` and column sums `24,18,12`: even magic incidence is
not preserved by the entrywise replacement `N1 -> N3` on arbitrary layouts.
The construction provides a specified new arrangement of point counts; it is
not an automatic magic-preserving operation on a nine-cycle's label layout.

The nine F25 counts, one for each `t=-4,...,4`, are `36−t²`. Their sum is
264. In any real 3×3 magic square its center equals the mean of all nine
entries; this would be `264/9=88/3`, not an integer entry. Hence this entire
nine-entry multiset admits no magic arrangement. This short obstruction is
independent of the inherited general quadratic theorem.

## 4. Twelve classes is not a twelve-object correspondence

Since every nonzero u in F5 satisfies `u⁴=1` and `u⁶=u² in {1,4}`, the scaling
orbits are `(a,b)~(a,−b)`. Among the 20 nonsingular models there are four
singletons with `b=0`, and eight doubletons with `b!=0`, giving 12 classes.
Counting the traces of the representatives gives fiber sizes

`(2,2,2,1,1,1,1,1,1)`;

the double fibers are at traces `−2,0,2`.

The d=2 quadratic twist sends `(a,b)` to `(4a,3b)`. It induces an involution
on classes: applying it twice sends `(a,b)` to `(a,4b)`, a scaling-equivalent
model. It has no fixed class. If `a!=0`, the first coefficient changes from a
to 4a, whereas scalings keep a; if `a=0`, nonsingularity forces `b!=0`, and
`3b` is neither b nor −b. Thus there are exactly six twist pairs, explicitly
listed in `result.json`. Direct counts verify the expected trace sign flip.

For a separate control, the auditor enumerated all `2^10=1024` labeled
five-vertex tournaments, canonicalized under all 120 permutations, and
computed exact characteristic polynomials by integer Newton identities.
There are 12 isomorphism classes, but their adjacency-spectrum fibers are
`(3,2,1,1,1,1,1,1,1)`, with eight fixed classes under converse and ten converse
orbits. These agree with the pinned Session143 tournament certificate.

Consequently no bijection of the 12-object sets can preserve these quotient
fiber relations up to relabeling; the fiber partitions differ. Nor can one
intertwine elliptic quadratic twist with tournament converse; the involutions
have different fixed-point counts. A bare 12-element set bijection exists
but is mathematically uninformative. No spectral or duality bridge is claimed.

## 5. Verification and scope

Reproduce from any directory after copying this entire audit folder:

```
python3 check_independent.py > result.json
python3 -O check_independent.py > result.optimized.json
cmp result.json result.optimized.json
python3 test_negative_controls.py > negative_controls.json
python3 verify_manifest.py
```

There are no `assert` gates. Normal and optimized outputs are byte-identical.
Twenty-two negative controls run in both modes: corrupted certificates and
code mutations to reducible field, missing infinity point, wrong cubic
coefficient, square instead of nonsquare twist, and power-versus-symmetric
trace confusion. Two unchanged-certificate controls pass. The manifest pins
all audit inputs and output bytes; source pins name the original repository
commit and relative paths. `audit_producer.py` independently replays the
producer and compares all twenty point-count tuples, magic arrangement count,
trace fibers and tournament polynomial fibers. Its `producer_review.json`
pins the reviewed producer proof and code bytes. The producer's older lexical
“lambda cube” hit was inspected only as provenance: for `φ²=φ+1`,
`α=1+3φ` has `Norm(α)=−5` and `α²=5φ⁴`; it is not the determinant +5 pair.
No interpretation of the user's intended referent follows from that hit. No repository edits, commits, or pushes were made
by this independent auditor.

**Stopping condition reached:** the requested F5/F25/F125 finite bridge and
its representation/quotient boundaries are verified. No all-extension search,
new conjecture, resolution of the lambda-cube referent, or square-of-squares
existence claim is part of this audit.

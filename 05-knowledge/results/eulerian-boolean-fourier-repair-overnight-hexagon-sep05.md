# Conditional Fourier repair into the actual Boolean kernel

**PROVED CONDITIONAL FORMULA / FINITE-EXACT n=3,...,7 / INDEPENDENTLY
WRITTEN-AUDITED.** The formula below explicitly retains
the extra invertibility obligation. It is not a proved all-order Boolean
nullity formula, nor a Laplacian-gap theorem.

## 1. The missing matrix rather than another scalar weight

The source is the complete parity basis in the
[self-complementary trace dictionary](overnight_hexagon_sep05_boolean_selfdual.md).
Its fixed Fourier modes are joint weighted zero modes but fail Booleanity
already at n=4. The [four-cycle gauge obstruction](eulerian-boolean-kernel-overnight-hexagon-sep05.md)
rules out every two-sided diagonal repair for n>=5. The least-used
sidecar is the full crossblock on the complementary-pair basis.

Order the Eulerian classes by even/odd edge parity and write the native
Boolean adjacency as `B=[[0,C^T],[C,0]]`, where C has q_- rows and q_+
columns. Let Delta=q_+-q_-. In the even parity subspace, let S contain the
Delta self-complementary Fourier orbit sums and P the q_- sums of free
complementary orbit pairs. The proved trace dictionary says `[S P]` is
an invertible q_+-square matrix. Define

```text
A=C P,                  D=C S.                     (1)
```

Here D is the actual Boolean defect of the weighted zero modes. No
multiplicity-weighted eigenvalue replaces either matrix in(1).

## 2. Exact conditional compiler

**If det A is nonzero**, set

```text
K=S-P A^(-1) D.                                    (2)
```

Then C K=0. The columns of K are independent: a relation K v=0 gives
`[S P](v,-A^(-1)Dv)=0`, hence v=0. Since A is invertible, C has full row
rank, C^T is injective, and all of ker B lies in the even parity subspace.
Consequently (2), extended by zero on odd classes, is a complete native
kernel basis and

```text
nullity(B)=Delta.                                  (3)
```

All entries are rational when the canonical integer Fourier orbit sums
are used. The operation is a genuine full-matrix correction, not a
rescaling of vertices or identification of weighted and Boolean operators.

The hypothesis is equivalent to `ker C intersect image(P)={0}`. It is
stronger than mere index saturation. For example the abstract matrices
`C=[[1,0,0],[0,1,0]]`, `S=e_2`, `P=[e_1 e_3]` have `[S P]` invertible
and C of full row rank, but C P is singular. This is a linear-algebra
hostile to that implication, not an Eulerian Fourier counterexample.
Whether the specific canonical P satisfies the hypothesis at all orders
remains **OPEN**.

## 3. Complete finite controls and a genuinely nontrivial repair

The exact program materializes both complete spaces through n=7, quotients
by adjacent vertex transpositions, constructs every Fourier orbit sum on
every Eulerian class, and verifies the square full/parity basis matrices.
It checks weighted annihilation of S, computes the actual Boolean A,D,
and verifies both the corrected kernel and its dimension.

| n | Delta | size of A | det A | rank D | denominator lcm in A^(-1)D |
|---|---:|---:|---:|---:|---:|
|3|0|1|2|0|1|
|4|1|1|4|1|1|
|5|1|3|5120|1|5|
|6|4|6|-53687091200|3|1600|
|7|0|27|nonzero, retained exactly in output|0|1|

At n=4, S=(6,-2)^T on the even classes and P=(2,2)^T. Formula(2)
gives K=(4,-4)^T, exactly the Boolean mode up to scale. At n=6 every
individual fixed mode fails Booleanity, but rank D=3 for four columns,
so one linear combination already survives. The full six-by-four rational
correction is retained in the output; these are four independent repaired
native modes, not just a dimension count.

The n=7 case is balanced and has no modes to repair; it verifies the
free block but does not extend the nontrivial correction test beyond n=6.
The completed native census there has54 classes and32,768 labelled states.
All values are exact. Domain-based rational row reduction avoids an
unnecessary expression-growth cost; it changes no matrix or predicate.

```bash
python3 -B 04-computation/eulerian_boolean_fourier_repair_overnight_hexagon_sep05.py
python3 -B -O 04-computation/eulerian_boolean_fourier_repair_overnight_hexagon_sep05.py
```

[Source](../../04-computation/eulerian_boolean_fourier_repair_overnight_hexagon_sep05.py)
raw LF SHA256:
`7c056c2a1384ef26deaa275f1822003d64c56cd3b6a2cf263a266113e482dfe9`.
[Output](eulerian_boolean_fourier_repair_overnight_hexagon_sep05.out) SHA256:
`8932b4352fbf1feee46852b8b6040243681ee8f1379f66b34201eb01e46625a4`.
These are finite controls for the explicitly conditional compiler, not
an unbounded transversality proof.

The independent observer audited the complete formula, hypothesis
directions, dimensions, source construction and all named boundaries.
The n=7 computation is not represented as a separate independent census.

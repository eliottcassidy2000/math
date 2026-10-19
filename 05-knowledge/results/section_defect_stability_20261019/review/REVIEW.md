# Independent audit: the 29-adic section-defect bound of 445

Date: 2026-10-19. Verdict: **PASS. No mathematical correction is required.**

The argument excludes 88 failures on the three direction families `(1,0), (0,1), (3,9)` and therefore proves at least 89 failures on that triple. Rotational averaging then gives at least **445 failures among 12,615 undirected residue edges** for every section from precision 29 to precision 29² in the unramified quadratic norm graph.

This audit found no omitted incidence case, modular alias, normalization counterexample under the stated hypotheses, or gap between compact stress compatibility and membership in the full carry coset. The last distinction is essential: many explicit 88-error arrays satisfy every cycle sum and the compact stress but are still rejected by the global period. Neither 89 nor 445 is shown to be sharp. No chromatic bound or Euclidean realization follows.

## Inputs and dependency

The audited source is `../README.md` before its documented publication-presentation update, SHA-256:

`c3b05ada18143db426cceeafa8cbd4a51e5acbdb969a2eeaa4e35dc356231272`

The audited verifier is `../verify_one_extra.py`, SHA-256:

`ee4650ead4ec7d6f22cb56afa65f2aa8f18751f3fbbe4a76eba266327a831b35`

The submitted exact-results JSON has SHA-256:

`7914d80eb4139f33bdd925f2f3208594c6984dc15a27ba36b78d5d657001625b`

The prior published checkpoint is `../../balanced_digit_section_defect_20261019/README.md`. Its affine carry identity, cycle sums, binary seam classification, period, and rotation argument were read together with the existing independent audit. The new independent controls also recompute the relevant carries, actual precision-two errors, all seam offsets, and the rotations; they do not merely trust the previous numerical period list.

## 1. Dipole normal form is complete

Every one of the 87 cycles in the three families has error sum −2 in F29 and therefore has at least one failed edge. Exactly 88 failed edges leave precisely one cycle with two failed edges and every other cycle with one. Failures belong to unique cycles because the three families partition their own edges and are distinct modulo sign.

On the exceptional cycle, exchanging the names of its two failed origins changes separation `k` into `29−k`, so one can always arrange `1 ≤ k ≤ 14`. Selecting `P` for the binary field gives exceptional errors `−2+t` at P and `−t` at Q. Thus the displayed dipole form is exact. In addition to the stated `t ≠ 0`, an actual two-failure candidate has `t ≠ 2`; that extra restriction follows automatically from the nonzero error at P and is not needed in the exclusion proof. Its omission from the statement is harmless.

For the stated convention `D_e f(x)=f(x+e)−f(x)`, the segment `S=Σ(j=1..k) δ_(P+je)` satisfies `D_e S=δ_P−δ_Q`, with the stated sign. Translation and swapping do not introduce a missing orientation.

Substitution into the affine stress gives `−2 Σ s_d H_d+t s_e C=0`, where `C=D_uD_vD_wS`. Off the cube, bounded-digit uniqueness really forces all three *integer* second differences to vanish. The congruence has no nontrivial bounded digit solution: its integer absolute value is at most 26, and reduction modulo 3 finishes the argument.

The independent controls check all 1,134 combinations of exceptional family, short separation, and allowable nonzero exceptional-error pair. They also check all 125 digit triples.

## 2. All cube incidence cases are covered

The pair determinants are 1, 9, and −3. In particular the three relevant nonzero magnitudes are distinct modulo sign. A collision between two cube vertices would express a nonzero multiple of the exceptional direction as a signed sum of the other two directions. Taking its determinant with the exceptional direction excludes both the one-term and two-term cases. Hence the cube has eight distinct points for every nonzero separation, not only the short separations used later.

For a fixed derivative pair, the cube consists of two disjoint four-edge rectangles in its bipartite line-incidence graph.

- If the pair contains the exceptional direction, the third direction's coordinate along the other member of the pair is neither 0 nor ±1. The two faces therefore have disjoint line sets on one side. Their two-element sets on the other side can share at most one vertex: equality of both sets under a nonzero translation would require an element of order two in F29.
- If the pair excludes the exceptional direction, the face displacement has two nonzero coordinates. A shared line vertex on both sides would require both coordinates to be ±1. Their ratio would then be ±1, excluded by the same determinant condition.

There are consequently two four-cycles, either disjoint or joined at one vertex. Zero line sums make the values alternate on each face. If the faces share a vertex, all other degree-two vertices already enforce alternation, and the shared-vertex equation introduces no extra relation. The integer kernel is the free rank-two lattice generated by those face cycles. No nonintegral coefficient or finite-characteristic exception is needed.

The independent implementation examines **84 geometries**, all three exceptional directions and all 28 nonzero separations. This strictly extends the submitted 42 short geometries. For each of the 252 direction-pair cases, it verifies the exact incidence graph, rational and F29 kernel dimensions, and an integral saturated kernel basis. The coupled matrix has rank four over both Q and F29 in every geometry.

## 3. Coupling over F29 upgrades to equality over the integers

Dividing each field's value by the cube sign converts its two face coefficients into a function of the single omitted coordinate. The stress is a sum of three such one-coordinate functions equal to a constant. Varying one cube coordinate while keeping the other two fixed shows its function is constant in F29. This uses only that every stress coefficient is nonzero in the field.

Thus `H_d=λ_d C` modulo 29. Choose `λ_d` from one actual signed value of `H_d`; it is an integer in [−2,2]. At any other vertex, the discrepancy between this integer value and `λ_d C` has absolute value at most four. A discrepancy divisible by 29 must therefore be zero. Off the cube both sides vanish. This proves the required integer identity, not merely a field-valued substitute.

As a hostile exact control, for every geometry the independent script enumerates every six-tuple of bounded integer face amplitudes, checks the pointwise modular coupling, and requires integer scalar proportionality. There are exactly **125** solutions in each geometry, one for each triple `(λ_u,λ_v,λ_w)` in [−2,2]³ before imposing a particular right-hand amplitude. There are no additional modular-only solutions. Free coordinates of the integral kernel are actual cube entries, so bounding them by two does not omit any bounded field.

## 4. The short-column normalization lemma is valid, and its hypotheses matter

The lemma is an integer statement. Here it applies to an integer `J_d=I_d−λ_d D_d S`, whose other mixed derivative vanishes exactly over the integers by the preceding step. Summing that mixed difference on grid rectangles gives an actual integer representation `J(s,t)=F(s)+G(t)`.

For clarity, the key estimate in the lemma is `(p−2)m ≤ p`. If F is nonconstant, normalizing with an entirely unmodified column makes F binary with `m ≥ 1` ones. Every unmodified column then has shift zero. Since `p ≥ 5`, the estimate forces `m=1`.

In either exceptional column at most `p−2` entries can change. There are `p−1` zero positions of F, so an unmodified zero position remains and forces that column's shift to be 0 or 1. Equality of total weights gives `ΣG=0`, forcing both shifts to be zero. If F is constant, each column has an unmodified position and its constant value is binary immediately. The argument works for at most two exceptional columns, including zero or one.

The support hypothesis is correctly verified in the application:

- For `d ≠ e`, choose the e-axis as the row direction in the derivative grid. The support of `D_d S` lies in two distinct columns, with exactly k segment positions in each before any zero scalar factor. Since `k ≤ 14 < 27`, the bound is stronger than necessary.
- For `d=e`, the two dipole points have different column coordinates in either derivative grid, since e is independent of either derivative direction. Each column has one affected position.
- Multiplication by λ cannot enlarge this support. Every derivative has total sum zero and zero sum on each of its own d-cycles.

The independent controls use the actual change-of-basis coordinates for all 126 short-geometry/family cases. No visual axis assumption is substituted for this calculation.

An independent exhaustive search through nonconstant binary F on grid sizes 5 through 9, all possible exceptional shifts, and all admissible binary exceptional columns finds no counterexample. Constant F is settled directly by the unchanged-entry argument. The audit also verifies three genuine counterexamples when hypotheses are weakened:

1. Allowing `p−1` changes per column can hide a nonbinary J, with equal weight p. This independently reproduces the submitted boundary phenomenon at p=5.
2. Keeping the `p−2` support cap but dropping the weight-p restriction fails already at p=5 with equal weights 10.
3. Allowing p=4 fails with equal weight 4 and two changes in each exceptional column.

These hostiles are controls of the proof's boundaries; none satisfies all its hypotheses.

## 5. The correction is an actual local motion

The normalization preserves the sum on every d-cycle, not just total weight. Thus each now-binary `J_d` has exactly one point on every d-cycle. The inherited integer separability lemma makes it one straight seam in one of the two allowed orientations.

The coefficient equation yields `L_w=3L_u+9L_v` in F29. Since the projection matrix for u and v has determinant 3, there is a vector H with `ell_u(H)=L_u` and `ell_v(H)=L_v`. Linearity then gives the w-equation. Hence the error change is exactly `R(HS)`.

If the putative z were in the carry coset, subtracting this explicit motion would place `−2J` in that same coset. This step does not infer sufficiency from compact stresses or cycle sums. It transports an actual assumed section by a specified vertex motion, which is precisely what is required.

To challenge this distinction, the independent script enumerates all normalized candidates at origin P=0 for every short geometry: all seam offsets and both orientations in every family, every λ in [−2,2], binary corrected fields, the selected exceptional origin condition, and both exceptional errors nonzero. There are **42,576** candidate parameter tuples. Each has a valid local-motion coefficient vector and a nonzero global period. These are tuples within the proved normal form, not a claimed enumeration of arbitrary sections or necessarily distinct arrays across different geometries.

For one nontrivial example in each of the 42 short geometries, the script expands the entire error array and directly checks 88 nonzero entries, compact-stress vanishing, the pointwise local-motion identity, and period preservation. The submitted checker separately supplies 84 explicit full-array examples.

## 6. The period and rotational conclusion remain valid

The global period is a necessary condition for membership in the full carry coset. Integrating the u- and v-equations introduces arbitrary functions of y and x. In the horizontal period their contributions vanish: the former is multiplied by 29 and the latter telescopes. Thus its nonzero value cannot be repaired by a different integration constant or a motion missed by the chosen integration convention.

The independent script computes literal carries from integer representative differences. It directly integrates each of the 58 seams in each family and combines the resulting independent contributions for all **195,112** orientation/offset triples and all 29 horizontal rows. It does not use the predicted period formula to populate those tables. The only values are 2, 4, 13, 14, 15, 16, 25, and 27, each occurring 29³ times; none is zero.

As positive controls, three deterministic actual sections are evaluated directly modulo 29² at every edge of the triple, giving 7,569 direct comparisons with the affine errors. Their periods vanish, their cycle sums are −2, and their full affine stresses vanish. These controls guard against a period test that accidentally rejects every input.

Norm-one residue rotations have norm-one lifts modulo 29² by the elementary scalar correction in the prior proof. The independent controls verify all 30 lifts and the 15 opposite-direction pairs. Rotating the triple by each pair counts each family exactly three times. Therefore `3M ≥ 15·89 = 1335`, giving `M ≥ 445`. There is no need for a disjoint partition or compatible failure locations across triples.

## Reproduction, checker limitations, and files

Run:

`python3 -B independent_controls.py`

The independent checker uses only the standard library and explicit `require` checks, not Python assertions. It was also run under `python3 -B -O independent_controls.py`; normal and optimized receipts and stdout agree byte for byte. Its checks therefore remain active with assertion optimization enabled.

The submitted verifier was copied unchanged to `submitted-replay/` and run with assertions enabled. It passed and reproduced the submitted results JSON byte for byte. **The submitted verifier uses assertions and must be run without `-O`.** This is a reproducibility caution, not a defect in its mathematical proof or normal-mode verification. A publication reproduction note should retain that caution.

Files:

- `REVIEW.md`: this audit
- `independent_controls.py`: independent exact implementation, importing no submitted code
- `independent-results.json`: latest full receipt
- `independent-normal.json`, `independent-optimized.json`: matching normal/optimized receipts
- `independent-normal.out`, `independent-optimized.out`: matching stdout summaries
- `submitted-replay/`: unchanged submitted verifier and its exact normal-mode reproduction
- `input-sha256.txt`: initial input hashes
- `review-manifest.json`: final hashes, scope, and verdict

The finite controls support, but do not replace, the universal proof above. No literature-priority claim was evaluated or added.

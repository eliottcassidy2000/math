# Independent review: the 29-adic section-defect bound of 440

Date: 2026-10-19. Verdict: **PASS for the 440 theorem. PASS for the general corollary with the supplied scope clarification.**

The argument proves that every set-theoretic section of the norm-one graph reduction from precision 29² to precision 29 fails on at least 440 of the 12,615 residue edges. I found no gap in the universal argument. It improves the straight-cycle lower bound of 435 by five. Neither attainability nor a chromatic consequence follows, and the draft correctly makes neither claim.

The optional general corollary needs its intended anisotropic/unramified setting stated explicitly at the start. A separate, reviewed Section 6 overlay makes that clarification and explains the lifted coordinates. It does not change the 29-adic theorem, numerical bounds, or verifier.

## Main findings

1. The affine stress annihilates the carry as well as the linearized motion. There is no substitution of a homogeneous rigidity equation for the actual mixed-characteristic equation.
2. Equality in the three independent cycle bounds really forces binary, single-error cycle supports. The digit separation and binary separability lemma classify all possible equality supports as the eight stated seam orientations; no nonseam equality case is omitted.
3. The period test is necessary for every choice of integrated solutions. Its dependence on integration constants vanishes in characteristic 29. Every seam offset gives one of the eight nonzero values 16, 14, 4, 2, 27, 25, 15, 13.
4. Norm-one rotations transport actual sections at precision 29². The fifteen direction-pair rotations count every family three times, so the passage from 88 to 440 is valid without any compatibility assumption or disjoint-triple hypothesis.
5. The original verifier passes unchanged and reproduces its recorded JSON byte for byte. A separate implementation also passes substantially broader direct controls, detailed below.

## Audit of the universal implications

### 1. Graph model, edge errors, and cycle lower bound

The discriminant −3 is nonsquare modulo 29, so the displayed norm is anisotropic and represents the unramified quadratic model. There are 30 nonzero norm-one directions, 15 opposite pairs, 841 vertices, and 841 × 15 = 12,615 undirected edges.

Every section has a unique expression X(x) + p h(x), with h valued in the residue plane. Expanding the quadratic norm modulo p² gives precisely

z_d = c_d + ell_d(D_d h).

Thus z_d = 0 is equivalent to preservation of the edge. Reversing an edge replaces the actual difference by its negative, so its norm error is unchanged. For odd p, indexing by a single representative of each opposite pair counts every undirected edge once.

For a fixed integer representative d, the wrap vector satisfies

X(x+d) − X(x) = d − p k_d(x).

The quadratic expansion gives c_d = K_d − ell_d(k_d) modulo p. On each length-p cycle, the sum of the actual representative differences is zero, hence the sum of k_d is d. Therefore the cycle error sum is −ell_d(d) = −2Q(d) = −2 modulo p; the motion term telescopes. Every one of the 29 cycles in each family has a nonzero error.

If three distinct families together have exactly 87 failures, every individual cycle has exactly one. Its error must then equal −2. This justifies the binary form z_d = −2I_d, with one selected edge origin per d-cycle. There is no unproved implication from an average error count to pointwise equality.

### 2. The mixed affine stress

Since w = 3u + 9v, the linear forms satisfy ell_w = 3ell_u + 9ell_v. Translation differences commute. Consequently the motion terms in the twelve-edge expression cancel identically.

For the carry, the integer identity D_e k_d = D_d k_e follows directly from the four representative values at a parallelogram. Applying the remaining difference makes all three terms involve D_u D_v k_w, whose coefficient is −ell_w + 3ell_u + 9ell_v = 0. Constant K_d terms disappear. This proves the full affine identity for every section; it is not a statement checked only on the particular motion used by the verifier.

The independent operator check additionally verifies exact cancellation over the integers of the two vector coefficients of the common translation operator. The finite tests of carry values supplement the algebraic argument.

### 3. Equality supports and seam completeness

Binary second differences lie in {−2,−1,0,1,2}. Hence |A − 3B − 9C| ≤ 26 < 29. A vanishing residue is a vanishing integer. Reduction modulo 3 forces A = 0, and the remaining equality B = −3C with |B| ≤ 2 forces B = C = 0. All 125 triples were checked independently as well.

On a cyclic grid in independent directions d,e, a vanishing integer mixed difference implies

I(s,t) = I(s,0) + I(0,t) − I(0,0).

This follows by summing mixed differences over grid rectangles; periodicity creates no exception. If both one-variable summands vary, each integer range has width at least one, and the range width of their sum is at least two. A binary-valued function cannot have that range. Therefore I is invariant in at least one of the two directions.

An invariant binary support is a union of whole parallel lines. Every line parallel to either of the other two directions meets each d-cycle exactly once. The condition of one selected origin per d-cycle therefore permits exactly one such line, not several. All three pairs of directions here are independent, since 3 and 9 are nonzero modulo 29. The resulting eight orientations exhaust the equality cases, with arbitrary offsets retained.

### 4. The horizontal-period obstruction

The projections H_u = ell_u(h), H_v = ell_v(h) determine h because their matrix has determinant 3. The target residuals r_u and r_v have zero cycle sums, so their cyclic difference equations are solvable. All solutions differ by f(y) and g(x), respectively.

In the horizontal sum of the w-residual, a change f(y) contributes a multiple of p and a change g(x) contributes a translated cyclic sum minus itself. Both are zero in the residue field. Thus the period obtained using any convenient integrated solutions is a necessary condition for an actual h, independent of their arbitrary constants.

The summation-by-parts identity is valid modulo p:

sum_x H_u(x,y) = −sum_x (x+1)r_u(x,y) = 2(A_y+1).

The boundary carry lies at x = p−1, whose weight x+1 is zero modulo p. A w-parallel u-seam satisfies A_(y+b) − A_y = a, while a v-parallel seam has constant A_y. This gives the first contribution 2a eta_u.

For the v-equation, summing the b successive vertical differences after summing over x gives −2b eta_v. A horizontal seam contributes either zero or p selected origins to a row, hence zero in the field; a w-parallel seam contributes one. Its boundary carry is also a whole row and vanishes in the horizontal sum. The w-seam similarly contributes −2 eta_w.

Finally, direct quadratic expansion, including simultaneous coordinate wraps, gives horizontal carry sum −a(2a+b). The remaining quadratic wrap term contains an additional p and vanishes after reduction. The complete period is therefore

P = −2 eta_w + a(2a+b) − 2a² eta_u + 2b² eta_v.

For a = 3, b = 9, the eight values are exactly those stated in the draft, all nonzero. The offsets do not occur in this formula. Therefore 87 failures on the triple are impossible, and the integral failure count is at least 88. No stronger local bound or equality construction is claimed or established by this audit.

### 5. Rotation and the final count

Every residue norm-one multiplier g has a norm-one lift modulo p². If a lift G has norm 1 + pk, multiplying it by the scalar 1 − pk/2 fixes its norm modulo p². Conjugating a section by this lifted multiplication gives another genuine section, so the triple theorem holds for each rotated triple.

Multiplication by each of u, v, w is a bijection on the norm-one group modulo sign. Thus summing the fifteen rotated inequalities counts each family exactly three times, regardless of repeated edge locations across different inequalities. Hence 3M ≥ 15 × 88 = 1320 and M ≥ 440. This reasoning does not require the optional disjoint-triple partition found in the verifier.

## General corollary and the reviewed correction

Within an unramified quadratic extension of Q_p, the generalized proof is correct. Independent residue unit vectors can be lifted to a basis U,V with norms one modulo p². Their coordinate norm is X² + 2TXY + Y² modulo p², with T reducing to t. The same carry calculation then yields

eta_w = a²(1−eta_u) + tab + b² eta_v.

The four possible values of the right side before subtracting eta_w are a²+tab, 1−tab, tab, and b²+tab, using a²+2tab+b² = 1. Requiring all four to avoid 0 and 1 is exactly the stated period exclusion. The small-digit hypothesis supplies the same seam classification. Averaging over the (p+1)/2 opposite pairs gives

M ≥ p(p+1)/2 + ceil((p+1)/6).

The number p+1 is a property of the anisotropic norm-one group; it must not be inferred from the coefficient hypotheses alone. For example, p = 31, a = 3, b = 9 and Q = x²−xy+y² satisfy the digit and period hypotheses, but this form is split: its discriminant −3 = 28 is a square modulo 31. It has p−1 unit directions. This example diagnoses the missing explicit scope, not a demonstrated counterexample to the numerical bound by some other argument.

The original ending already mentions an anisotropic plane, so this is a clarification of its intended hypothesis. The supplied overlay states the unramified setting at the beginning, explains the coordinate lifts, and states which group count is used. It is the recommended publication version of Section 6.

## Exact-control audit

The original verifier is exact, standard-library-only code. Its rigidity, carry, deterministic-motion, cycle, seam-period, and rotation checks are meaningful and pass. Its modular-23 hostile digit triple is a valid control for the no-aliasing step.

One coverage distinction is worth retaining: the original loop over 195,112 offset/orientation tuples checks analytic support-sum identities at selected rows. It does not independently integrate and evaluate all those systems. The draft calls these “offset-independence identities,” which is defensible, but the number should not be described as a full-system brute-force integration count.

The independent implementation imports no submitted code. It uses y-major arrays and literal integer representative differences, and verifies:

- 37,845 direct comparisons of actual norm errors modulo 29² against affine errors, across all 15 direction families and three deterministic sections
- All 435 direction-cycle carry sums, and all corresponding cycle sums for those sections
- Exact translation-operator cancellation and all 841 carry stresses
- All 125 small-digit triples, with the modulus-23 hostile example
- Exhaustive binary-separability controls on all 2×2, 3×3, and 4×4 binary grids; the general lemma is proved above, not inferred from this sample
- Every seam offset in every family: 58 seams per family, their required invariances, and 5,046 one-per-cycle checks
- Direct integration of all 58 possible u-seams and all 58 possible v-seams; independent computation of all w-seam row sums
- Combination of those directly integrated contributions for all 195,112 full offset/orientation systems and all 5,658,248 horizontal row periods. This reuse is exact because the residual is additive in the three independently chosen seams; no predicted period formula is used to build the tables
- All 30 norm-one lifts, group closure, and the threefold rotation incidence
- Arithmetic equivalence of the general period criterion in every prime field from 5 through 101, including the split-form diagnostic; separate direct integrations for all eight orientations in one qualifying anisotropic example at each of 16 primes

Every check passed. The finite controls are not an exhaustive enumeration of arbitrary sections and are not presented as one. The universal conclusion rests on the proof audited above.

## Files, reproducibility, and scope

Run `python3 independent_checks.py` from this review directory. The script writes `independent-results.json` and prints the same JSON. The saved run is `independent-checks.out`.

The unchanged author verifier and its receipts are in the package's parent directory. `overlay-manifest.json` records the original and reviewed scope-clarification hashes; `SOURCE_PROVENANCE.json` also records the final presentation-only cleanup. The published README retains the full reviewed Section 6.

Original audited README SHA-256: `79603eb3bc14893c45b9080d3828bdeb6af544403b0b4c374cc8154e9cc83111`.

Original verifier SHA-256: `63a73900f8d10783c830bfbdf74d33c8b44cea9fbaae0a2d345b6aa59e729240`.

This review does not assess exhaustive literature priority or the external chromatic-status statements. They are not used in the proof. It supplies no speculative strengthening, no new chromatic bound, no Euclidean realization, and no claim that 440 is sharp.

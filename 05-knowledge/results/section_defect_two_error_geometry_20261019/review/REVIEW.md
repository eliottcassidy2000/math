# Independent audit: two extra errors at p = 29

Date: 2026-10-19. Verdict: **PASS for the frozen original packet. No mathematical correction required.**

The theorem excludes an 89-failure section with one three-failure cycle and one failure on each other cycle. The two-double-cycle result in this frozen packet is a conditional exclusion criterion requiring disjoint mixed-cube supports, coupled homogeneous dimension four, and the stated normalization hypotheses. The disjoint-support dimension-six hostile is correct. Numerical lower bounds remain **89 on the mixed triple and 445 on the full norm graph**.

The separate path-complement addendum has its own **PASS** review in `../path-complement/review/REVIEW.md`. That addendum removes the original criterion's extra normalization condition; it does not change the original files, remove the rank-four/disjointness conditions, or improve the numerical bounds. This review preserves the historical scope of the frozen original.

## Frozen inputs and dependency boundary

The original audited input hashes are retained below; published presentation hashes appear in `review-manifest.json`.

Frozen README SHA-256:

`2ea9c2a21fec044d5424443998ef6b6710c281a1dc214d456372e448772e1960`

The publication list contains exactly `README.md`, `verify_two_extra.py`, `verify_two_extra.out`, and `two-extra-exact-results.json`. Initial hashes of these and the publication list are retained in `input-sha256.txt`; all were checked unchanged after replay. The exploratory probes are not used.

The inherited carry identity, balanced-digit stress, seam classification, nonzero periods, and global 445 bound were read against the prior independent reviews in `../../balanced_digit_section_defect_20261019/review/REVIEW.md` and `../../section_defect_stability_20261019/review/REVIEW.md`. This audit does not claim a new independent reconstruction of every earlier carry calculation: those previously reviewed results are explicit dependencies. The new finite lemma, normalization proofs, integer lift, and local-motion reduction are audited here independently.

## 1. The 1,134-case universe is complete

An exceptional e-cycle with three distinct failed origins has a distinguished baseline P and two distinct other origins. Translation sends P to zero, leaving Q = q e and R = r e with q,r nonzero and unequal. Interchanging Q and R leaves the support and conclusion unchanged, so exactly three choices of e and 378 unordered pairs remain. No scaling of directions or unverified symmetry reduction is used.

The independent implementation evaluates every derivative pointwise on the entire 841-point grid. It imports no submitted code, uses y-major coordinates, includes all 29 line labels whether occupied or empty, and computes ranks by forward elimination with multiplicative-table inverses. This is separate from the submitted sparse derivative and Gauss-Jordan/nullspace implementation.

For every geometry it verifies:

- Both cubes have eight points and coefficients ±1.
- The shared face and each private face have exactly four points.
- The cube signs agree on the shared face.
- The full 36-variable line/coupling constraint matrix has rank 32 over F29.
- The four displayed model vectors are independent and lie in that kernel.
- Each individual W-space dimension agrees with the submitted certificate. The coupled rank also agrees, computed independently as the sum of domain dimensions minus the full homogeneous nullity.

The reported histogram is exact: 372 generic cases per exceptional direction, plus two groups of three exceptional cases per direction. All 18 nongeneric cases retain homogeneous dimension four. Every one of the 1,134 submitted per-geometry records agrees with the independent calculation.

The independent single-cube checks cover all **84** nonzero separations, including the 42 short separations used by the criterion. Their W-dimensions are [2,2,2] and homogeneous dimension is two. This supports the stated disjoint-line sufficient condition without a missing long-separation case.

## 2. The inhomogeneous reduction and integer lift are valid

The section's cycle sums force the displayed dipole normal form. For an actual three-failure cycle the coefficients t1 and t2 are nonzero, and the selected-origin error −2+t1+t2 is nonzero as well. The latter restriction is automatic from the original assumption and is not needed later.

Outside the union T, the compact stress has right side zero. Each mixed second difference is an integer in [−2,2]; the balanced-digit implication therefore forces all three components to vanish there over the integers. The independent control checks all 125 possible bounded digit triples.

Every component has zero line sums in its two derivative directions. The particular inhomogeneous solution is valid because each cube has zero sums in all three line families: in the exceptional family the original dipole has zero total on its line; in either other family one derivative telescopes. The homogeneous rank lemma then yields the claimed coefficient form modulo 29.

Each coefficient can be read from a signed private-face entry, so it has an integer representative in [−2,2]. On the common face the combined value has absolute value at most four, while the actual second difference has absolute value at most two. Their discrepancy has absolute value at most six, strictly less than 29, so the field equality is an integer equality everywhere.

As an additional exact control, for each geometry the independent script enumerates all **29² coefficient pairs** whose cube combination takes residues in {−2,−1,0,1,2}. Exactly 19 survive: the integer pairs a,b in [−2,2] satisfying |a+b| ≤ 2. All lift without any modular-only alias. This directly challenges the private-face step rather than relying only on numerical rank.

## 3. Both integer binary-normalization lemmas are sound

For the generalized column lemma, an untouched column normalizes a nonconstant row factor F to a binary vector. Every entirely untouched column then has zero shift. If F has m positive entries, the untouched columns contribute at least (p−r)m to the binary array I. Since p > 2r, m must equal one. Every exceptional column retains an unchanged zero of F, because there are p−1 such zeros and at most p−2 changes. Its shift is 0 or 1. Equality of total weights forces the sum of all shifts to zero, hence all shifts vanish. If F is constant, one unchanged position in each column settles the result directly.

The proof does not assume nonnegativity of the unknown additive J where it may differ from I. The only lower bounds on weight are taken on genuinely unchanged binary positions.

For the crossing lemma, a binary additive array on the clean rectangle must have at least one constant factor. The clean rectangle cannot be all ones because its area exceeds I's total weight. If both factors are constant there, normalize both to zero; unchanged strip positions make every remaining factor binary, and the weight equation makes the sum of all factor values one. If one factor varies, normalize it to 0 and 1, with the other factor zero on clean columns. The untouched-column weight bound forces exactly one clean one-row. The strip caps leave an unchanged clean zero-row position in each exceptional column. Again all remaining factors are 0 or 1, and the total-weight equation excludes any additional one. The symmetric case is identical.

The independent controls supplement these proofs with exhaustive finite checks:

- The column lemma for n = 5,7,9, every allowed number of designated columns, and all normalized nonconstant binary row factors and possible exceptional shifts. The exact change-budget/weight test examines 3,603 candidate additive arrays.
- The crossing lemma for (n,r,s) = (5,1,1), (7,1,2), (7,2,1), (7,2,2), and (9,2,2). A clean-rectangle gauge and the existence of an unchanged strip entry bound every potentially feasible exceptional integer factor, so these finite cases cover all additive integer candidates, not random samples. Exactly 70 feasible arrays across those cases remain, and all are binary.
- The explicit n = 29 boundary hostile has binary I of weight 29, nonbinary additive J of the same weight, and exactly 28 changed positions in each of two columns. Thus the p−2 cap cannot be weakened to p−1 on the lemma's information alone.

These small-grid controls are corroboration; the stated universal lemmas rest on the audited arguments above.

## 4. Largest-gap paths satisfy the hypotheses in every case

The three positive cyclic gaps sum to 29, so the largest is at least ten. Starting after that gap puts both other origins on a single arc of length at most 19. Both paths are nested initial segments, with the sign convention D_e S = δ_P−δ_Q.

For d unequal to e, the supports of the two D_d paths lie in two e-parallel columns. Nestedness limits each union to at most 19 points. For d = e, the correction lies at P,Q,R. Since e is independent of either derivative-grid axis, those points are in three distinct columns, with one point per column. All 1,134 geometries are checked independently in dense coordinates, including the rebasing sign and the sums on each individual d-cycle.

Subtracting the resulting correction from I gives an integral additive J. Its total weight remains 29, so the column lemma makes it binary. Its individual d-cycle sums remain one, so it is exactly the one-per-cycle object needed by the inherited seam obstruction.

## 5. The motion preserves the full carry coset

The coefficient relation is precisely L_w = 3L_u+9L_v. Since the u,v projection matrix has determinant three in F29, a constant vector h_j realizes all three L_dj. With the verified dipole sign, the algebra gives

z = −2J + R(h_1 S_1+h_2 S_2).

Thus an assumed actual section puts −2J in the same affine carry coset. This is a genuine explicit motion, not an inference that compact stress compatibility is sufficient for a section. The inherited nonzero global period then gives the contradiction.

The independent controls expand the full 841-point arrays for both potentials and all nine pairs of exceptional directions. They verify the projection compatibility and the complete pointwise identity in all three components, giving 22,707 edgewise comparisons. The two potentials have different lengths and translated origins, guarding against a hidden same-origin or equal-length simplification.

Since each of the 87 cycles already contains a failure, the only distributions at total 89 are one triple-failure cycle or two distinct double-failure cycles. Excluding the former is therefore the claimed complete distributional reduction.

## 6. The original conditional criterion and hostile have the stated scope

Disjoint cube supports make the four endpoint-cube modes independent and give exclusive entries for each coefficient. If the homogeneous dimension is four, these modes exhaust it. With the normalization hypotheses, the same integer lift and coset motion apply. For distinct exceptional directions the correction is supported in two derivative-grid rows and two columns: a same-component endpoint dipole can be assigned to the perpendicular strips. The crossing lemma applies at p = 29. For equal exceptional directions the original frozen statement correctly identifies a possible shared-column union of 28 positions; it makes only a conditional claim at that point.

The concrete separated example has support size 16, W-dimensions [4,4,4], direct full-constraint rank 44, and homogeneous dimension four. Hence the criterion is nonvacuous.

The hostile's endpoint cubes have disjoint eight-point supports, while the intermediate cube lies in their union. The two outer layers force any expression of the middle cube as a combination of the endpoints to have both coefficients zero, proving independence. Independent exact elimination yields W-dimensions [6,4,6], full-constraint rank 42 on 48 variables, coupled rank ten, and homogeneous dimension six. All six stated model modes are independently verified in the kernel and have rank six.

The endpoints are actual short u-dipoles on distinct u-cycles, so this is a relevant incidence geometry. It is correctly described as a linear-incidence hostile, not a binary error array or an attaining section. The displayed w-amplitude ±3 itself warns against that stronger interpretation. Nothing in this example proves that additional local-motion correction or another exclusion is impossible.

## Reproduction and artifacts

Run `python3 -B independent_controls.py` in this directory. All independent checks use explicit exceptions rather than Python assertions. Normal and optimized (`python3 -B -O`) runs give identical JSON and stdout, retained separately.

The submitted verifier was copied unchanged to `submitted-replay/` and run there. Its output JSON and stdout reproduce the submitted files byte for byte. **The submitted verifier uses assertions and should be run without `-O`.** That is a reproduction caution, not a gap in its normal-mode verification.

Five deliberate fault controls detect a changed cube sign, duplicated failed origin, incomplete separation universe, the false claim that disjointness forces dimension four, and a wrong stress sign.

Main files:

- `REVIEW.md`: this audit of the frozen original
- `independent_controls.py`: independent dense-grid exact implementation
- `independent-results.json`: complete independent certificate receipt
- `independent-normal.json`, `independent-optimized.json`, and their `.out` files: matching normal/optimized runs
- `submitted-replay/`: unmodified submitted checker and exact reproduced outputs
- `input-sha256.txt`, `review-manifest.json`: frozen-input and review provenance
- `../path-complement/review/`: separate audit and controls for the separately frozen addendum

Neither sharpness, a 90/450 bound, a chromatic consequence, nor literature priority is established or implied by this review.

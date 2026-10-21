# Independent final audit of the Session 23 extension synthesis

Verdict: **PASS for the exact synthesis bytes below. No mathematical correction or source edit is needed.** Reviewed 2026-10-21. This is a bounded synthesis audit, not a new research lane, a replacement for the nine independent packet reviews, or a fresh literature-priority assessment. No repository modification or Git operation was performed.

Audited synthesis: `session23b-root/SESSION23_EXTENSION_SYNTHESIS.md`.

SHA-256: `6701313d62bd1b2930ed5ecf064580110dfb0994af4aac3f54c6cc075417b824`.

## Evidence and exact source binding

I read the synthesis against the final primary sources and independent reviews for capacity, reciprocity, explicit pentagon flows, flow-decorated SCC words, tetrahedral duality, geometry, Fano moves, signed stabilizers, and acyclic palettes. The pentagon lane includes both the raw triangulation addendum and the separate five-case anti-fixed derivation. I inspected their proof arguments, numerical formulas, finite-control scopes, and inheritance disclaimers; I did not rerun the large censuses.

`source_bindings.json` pins 127 exact files covering the nine publication packets, their independent reviews and receipts, the ten linked wrappers, and this synthesis. Its keys retain original workspace paths, with a separate explicit published path per file. At audit time every one of the nine extension publication manifests was verified for file size and SHA-256, and each bundled producer/reviewer artifact was byte-identical to its original workspace counterpart. The first-phase wrapper was verified against `SESSION23_SYNTHESIS_PUBLICATION_RECEIPT.json` instead. Its prior results are context, not newly proved here.

The only transient issue was a wording-only duality status promotion after the first review hash. The duality reviewer has explicitly checked and rebound that promotion. Its final proof hash is `af297f77570dc82643e1453eaccc3f7b6976001e4dc016996a07a01d96d43729`, and the final review, manifest, receipts and published bundle now agree. No unresolved binding issue remains.

## 1. Capacity transport and the two different Ehrhart laws

The incidence sign is consistent: with inflow-minus-outflow incidence, complementation on reversed arcs gives `B' f' = Bf - M B_D 1`. For nonpath reversal from a regular tournament, `B_D 1 = -(e_last-e_first)`, so the target demand is positive `M(e_last-e_first)`. The first endpoint has output, the last has input. Target outdegrees are respectively `r+1,r,...,r,r-1`; the partner is not regular. The 24/24/14 capacity-one hostile is correctly typed. The invariant `M s+b`, centered divergence, and fixed-edge centered quadratic energy follow algebraically. No electrical constitutive law is inferred.

The labelled-score complement law is explicitly inherited from THM-646, with the path-direction convention noted in the primary source. The affine integral flow lift is the added construction. Equal labelled scores imply equivalence of unweighted bounded-circulation polytopes, not equality of automorphism quotients or total-use weights. The cyclic7/Paley7 counts 175/189, shared capacity-one count 2640, and different fourth-use coefficients 28/21 are correctly used to demonstrate information loss. The synthesis never turns this into an A000568 enumeration result.

For the full regular polytope, dimension `(n-1)(n-2)/2`, interior shift two, exact numerator degree `d-1`, palindromicity, and the separate order-one exception are correct. The full order-five polynomial agrees with the six-degree Ehrhart count and the numerator `1,17,72,72,17,1`. Since its value at one is positive, the stated denominator exponent is reduced and exact.

For the zero-on-the-same-path face, strong connectivity after path deletion follows from the strict cut bound for odd `n>=5`. Dimension is `(n-1)(n-4)/2`, not the full dimension. The exact codegree is three, giving numerator degree `delta-2` and no general palindromicity. At order seven, dimension nine, degree seven, constant term one and leading coefficient 34 agree with the independently checked endpoint-demand count. The score interpretation has `r-1` first and `r+1` last, correctly opposite to the nonpath-complement partner's score shift.

The added linear term is elementary and valid: the coefficient of `w` in `h*(w)/(1-w)^10` is `h*_1+10`, so `h*_1=91-10=81`. This term is also present in the final bound capacity proof. It does not compute any other missing coefficient. The order-five simplex count and capacity/use series are correct; all three displayed capacity-one use polynomials have three unweighted objects but different weights.

## 2. Reciprocity, rarity, and explicit pentagon series

The formal-variable reciprocity theorem retains regularity, total arc-use grading, and the lattice of the fixed space. The anti action reverses arc endpoints as well as permuting vertices. The all-ones translation supplies shift `m`, and the odd-order/anti-involution argument supplies the common dimension parity. Anti-involution dimension `r^2` is not incorrectly assigned to every anti element. The denominator-factor parity qualification is retained, without pretending an arbitrary presentation is reduced.

The general symmetry-rarity theorem is cumulative and assumes a fixed self-converse regular tournament of odd order at least five. Its exponent is `-r(r-1)`, hence `-6` at order seven. Capacity boxes have the same exponent but need not have the same constant. The cumulative R5 constant `3150/17` is correctly separated from the pointwise `2100/17`. The synthesis uses M as an explicitly named cumulative threshold in this paragraph; it does not identify that cutoff with capacity or substitute it into a capacity generating function.

All three Paley fixed-action formulas match their independently derived orbit-coordinate sums and satisfy sign minus and shift 21. They do not supply the missing raw Paley or anti-involution numerator. Strong nonregular T4 has common shift ten but opposite raw/anti-fixed signs, exactly as stated.

The R5 raw numerator, ordinary Burnside average and converse-fixed numerator are correct. Their proof scopes remain the certified 25-cone unimodular fan and the independent five-case lattice decomposition, not fitting. Exact rational reduction independently confirms minimal recurrence orders **22,26,22,32**, minimal quasiperiods **60,60,120,120**, and the stated three pointwise leading coefficients. The required primitive-root poles are genuine. Dividing the leading coefficients yields `2100/17`; cumulative summation multiplies this ratio by `6/4`, yielding `3150/17`.

## 3. SCC decorations and actual invariants

The T4 cone is genuinely free on the two triangles and four-cycle; the anti map swaps just the triangles. Cross-SCC circulation vanishing and ordered-word factorization are correctly retained. All-order rationality is explicitly conditional on the inherited finite strong-atom and center bounds and concerns **underlying** h and underlying fundamental black degree, including degree zero.

Independent opposite halves require the square of the ordinary decoration series; mirrored halves require the grading substitution `q -> q^2`. The distinct neutral denominators in the synthesis reflect this exactly. The unequal C3 multiplicity hostile is valid. No decorated stabilizer is substituted into the original h, no arbitrary single anti-map replaces an anti-coset average, and no circulation is confused with an Euler traversal.

## 4. Duality and metric geometry

The face rule and the omitted-opposite-pair rule give the canonical cycle-point-to-edge-line duality. The orbit-size obstruction `4+3` versus `1+6` forbids an equivariant point bijection, not the duality. Neither construction automatically identifies the metric features with Fano points.

The rational tetrahedron counterexample uses three finite distinct collinear labels and is independent of A's interpretation. The complete exceptional-locus statement is restricted to the second-intersection interpretation and preserves the finite-definition and rank-loss caveats. Its valid locus has precisely the two repeated pairs, leaving five distinct noncoplanar points, not seven concyclic points.

The positive generic theorem states the independence condition on v,s,w. It has a different sphere, three inverse pairs `(D,E),(B,F),(C,4C)`, and strictly exterior A and H'. The synthesis does not silently replace the false original construction by this one or assert six distinct points at every special independent configuration. The source's uniqueness caveat where both original pairs coincide remains consistent with the summary.

## 5. Fano moves and colors

The eight coherent-boundary completions, H=171 obstruction, and two Paley completions after one facial reversal match the exact source. The 16/112 census and Hamming-code locus are expressly inherited. The Q7 inclusion is isometric for triangle reversal; quotient return retains a nonzero weight-three residual codeword, so it is not full-presentation return. Adding parity alone does not retain that residual state.

The full quotient matrix, class sizes, eigenvalues, H-weighted move-word function, and hitting means are correct. Small exact matrix algebra independently reproduces them. The synthesis specifies the two different transition laws and explicitly separates triangle reversal from the original fixed-path nonpath-complement operation. All seven-vertex anti permutations are odd, whereas Fano collineations are even; exactly two contained Paley designs are established in the independent Fano review and are exchanged by the anti actions. No higher-order equitability is assumed.

Gluck's imported odd-order regular-subset theorem is correctly credited and the uploaded paper is not misrepresented as a geometry source. Binary ordinary rigidity, the possible residual anti-involution, ternary signed rigidity for n>=2, and empty/singleton kernel exceptions are accurately summarized. The three regular-color polynomials match the independent inclusion-exclusion calculation and give least palette sizes 2,3,4. Equivariance cannot choose a single symmetry-breaking subset of an unchanged symmetric object.

The acyclic palette law is a separate predicate. Its rational budget series has reduced denominator degree four; lexicographic substitution yields 7 rather than 9 for the double Paley product. The depth recurrence has the displayed initial values. The non-C-finiteness proof uses bounded carry plus eventual periodicity forced by a hypothetical recurrence, then a divisibility contradiction. It retains a reduced nonintegral multiplier and positive initial state, and does not infer periodicity from a finite alphabet alone. Distinguishing, acyclic and unit-distance coloring are not conflated.

## 6. Boundaries and inherited work

The session chronology is consistent with the existing first-phase synthesis and session status record. Incoming commit `1e25cb70aec303cb6e8515bad610c9e99ddb5f7f` is acknowledged as preserved external work, not re-proved or used as a premise for the new flow or rounded-carry theorems. This audit performs no Git ancestry or publication verification and does not independently certify that incoming theorem.

The synthesis leaves complete original metagraph wiring, uniform153N and all-horizon native-clock obligations, a general A000568 speedup, and the named flagship conjectures open. It also correctly marks the remaining order-seven path-face coefficients and two large Paley numerators as uncomputed. The final proposed tests are identified as future targets rather than conclusions. No claim of literature priority is made.

## Reproduction

The checker performs **404 always-active checks**, including 127 exact-file bindings, ten wrapper destinations, 66 prose/formula/scope anchors, and the small algebra controls above. It does not rerun the expensive independent enumeration controls. It requires SymPy.

Workspace layout:

    python3 session23b-synthesis-review/check_synthesis.py --root /path/to/workspace
    python3 -O session23b-synthesis-review/check_synthesis.py --root /path/to/workspace

Published layout:

    python3 session23_extension_synthesis_20261021/session23b-synthesis-review/check_synthesis.py --results-root /path/to/repo/05-knowledge/results

`normal.json`, `optimized.json`, and `published-layout.json` are byte-identical. The third run tested the full explicit published-path mapping in a temporary filesystem layout; it does not assert that any push has occurred. Any change to a bound source, review, receipt, wrapper, or synthesis requires an explicit new binding and audit rather than silently updating the manifest. Work stops at this scoped final audit.

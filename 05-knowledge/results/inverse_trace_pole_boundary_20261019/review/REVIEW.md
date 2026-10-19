# Independent audit: inverse trace poles and alternating monodromy

Date: 19 October 2026.

## Verdict

**PASS, scoped proof audit. No blocking mathematical defect found.**

The main detector, exact local pole formula, maximal-cancellation packet, and even-degree alternating-monodromy gate follow under the stated hypotheses. The local discrepancy hostile is valid. The arguments do not prove regularity of any inverse trace for an arbitrary Keller map, exclude the residual quartic A4 possibility, or resolve the Jacobian conjecture.

Reviewed proof: `TRACE_POLE_BOUNDARY.md` before the documented status-only publication update.
SHA-256: `e8f95b7c28e7ba84c1e0e9cf64aa5aa994ca8f70491799e8d4cf1f9f820dc9bd`.

The submitted control script writes its neighboring report file. Its reviewed code is included here unchanged.

## Proof audit

### 1. Finite normalization and actual finite sheets

The factorization through the finite normalization is justified. The relevant finite-generation property is that a finitely generated C-algebra is Nagata/excellent; normality alone would not imply finite normalization in every setting. Here the base is C[u,v], so there is no missing finiteness hypothesis. Elements of its integral closure in L are integral over C[x,y] as well, hence belong to C[x,y]. The algebraic version of Zariski's Main Theorem then identifies the entire quasi-finite source with an open subset of that normalization.

Lemma 2.1 is correct. If an affine prime D=V(g) were omitted, g(P,Q) would be a nowhere-vanishing complex polynomial, hence a nonzero constant. Dominance excludes that. Étale openness then promotes one intersection to a dense open intersection with D. At a source generic point above D, the map is unramified; after a geometric transverse strict henselization, an actual inverse branch has both coordinates integral and inertia fixes it.

This uses the polynomial affine source and étaleness, rather than a bounded value of one arbitrarily chosen rational observable. No global section along D is claimed.

Lemma 2.2 correctly removes a possible codimension-two loophole. If both x and y are regular at z in the normalization, they give a local map to the affine source, inverse to the source open immersion on a dense open. Separatedness extends that identity and puts z in the source open. Normality makes every failure of coordinate regularity visible at a height-one valuation. A finite map preserves the dimension of boundary divisors, so their images are target divisors. Properness on a neighborhood meeting that image would force the source open immersion to be proper, hence closed and surjective onto the irreducible finite normalization over that neighborhood, a contradiction.

**Counting convention:** d counts geometric generic branches, including ramification multiplicity over a boundary valuation. The finite branches are those having affine source centers. These are not, in general, the number of reduced points in the normalization fiber.

### 2. Exact pole detection and maximal cancellation

For maximum inverse pole order lambda>0, the r most singular values have nonzero scaled residues alpha_i. The inequality ord(s_k) >= -k lambda is automatic. If their first r power sums all vanished, characteristic-zero Newton identities would force their nonzero product to vanish. Therefore at least one k <= r attains equality, proving the exact maximum formula.

Consequently a surviving bounded source sheet leaves at most d-1 escaping values of every source polynomial h. Repeated conjugate values when h is nonprimitive cause no problem: the trace counts them with multiplicity, and the lemma is explicitly a multiset statement.

If exactly m values escape and the first m-1 traces are regular, r<m is impossible. All m have common pole order, their first m-1 leading power sums vanish, and their leading polynomial is Z^m-c with c nonzero. It is separable in characteristic zero. The m-th trace has exact pole m lambda, which is an integer because the trace lies in the base discretely valued field. This proves the root-of-unity packet without imposing a common ramification denominator in advance.

The vector addendum is also sound when the phrase “leading pole vectors” means the m nonzero vectors under consideration. For generic linear forms, every projected vector is nonzero. Ratios of their projected values belong to the finite set of m-th roots of unity, so those rational functions on the irreducible space of linear forms are constant. This gives collinearity and the complete root-of-unity packet. Ordinary collinearity alone is insufficient, as the draft states.

### 3. Generic-coordinate and d-1 hypotheses

The simultaneous generic-coordinate choice exists over C. Choosing a transverse direction outside the zeros of the top homogeneous part of P makes the leading coefficient in y a nonzero constant, which can be normalized to one. A generic affine linear x separates the finitely many generic source points, hence is primitive. The finitely many forbidden primitive directions and forbidden transverse directions can be avoided simultaneously.

The monic P(x,y)-u relation is load-bearing: integral x implies integral y at every valuation over a finite target divisor. It therefore upgrades “x is bounded” to “the entire inverse source point is bounded.” Primitivity ensures that the minimal polynomial used later has degree d and distinct generic roots. The general trace-pole formula for an arbitrary h does not require primitivity.

The separate Newton argument in Section 5 is a valid independent proof of the one-moment saving: regular T_1 through T_(d-1) make a_1 through a_(d-1) integral; one bounded root forces a_d integral; the entire monic polynomial then has integral roots. The monic transverse equation bounds y as well.

Thus the pole supports of the first d-1 moments equal the complete reduced nonproper locus, not merely its generic subset. Polynomiality of those moments implies finite F. A finite complex Keller map is a covering of C^2, and connectedness and simple connectedness force degree one. Finite birationality onto the normal target yields the polynomial inverse.

The d=1 empty-list discussion is consistent: finite-sheet survival leaves no possible escaping branch, and the boundary lemma supplies properness. No bound on unknown generic degree is inferred from finitely many moments of an unrelated germ.

### 4. Even-degree alternating gate

Assume the first d-2 traces are regular at every affine prime. At a nonproper divisor, the detector forces exactly d-1 escaping x-values. Monicity then makes the remaining bounded x-value precisely the unique actual finite source sheet.

Writing m=d-1 and a=-ord(T_m)>0 gives lambda=a/m. Distinct root-of-unity leading residues imply that every difference between two escaping values has valuation -lambda. The difference between an escaping value and the unique bounded value has the same valuation. There are binom(m,2)+m=d(d-1)/2 such pairs, so the discriminant order is exactly -da. There is no possible cancellation in this multiplicative discriminant calculation.

At an affine prime outside the nonproper set, the full finite map is étale. After an unramified splitting extension, all roots and their differences have integral valuation. Each difference is squared, so the discriminant order is even. The discriminant here is that of the primitive-coordinate polynomial; projection coincidences can produce even zeros without ramification. The draft correctly distinguishes these zeros from branching.

When d is even, every affine prime order of the nonzero rational discriminant is even. Since C[u,v] is a UFD and every nonzero complex constant is a square, this implies a square in C(u,v). No separate check at projective infinity is required: factorization of a rational function already proves the square globally. The Vandermonde square root then shows that every element of the geometric monodromy group is even.

The sharper inertia calculation is correct over a geometric transverse trait. Continuation multiplies leading residues by exp(-2 pi i a/m). Distinct leading residues label the escaping roots uniquely, so no further permutation can hide in higher terms. The cycle type is 1 times (m/g)^g, g=gcd(a,m). This statement concerns geometric inertia, not the decomposition group or unramified residue-field monodromy.

For d=4, polynomial T_1,T_2 excludes S4. The only transitive subgroups of A4 on four letters are A4 and the regular V4. The inherited finite-sheet/normal-generation theorem excludes V4; its proof uses purity and trivial affine étale fundamental group. That extra inherited input is correctly named rather than smuggled into the square-discriminant proof. The A4 lane is left open.

### 5. Local discrepancy hostile

The displayed local model has

- x=z^(-1), y=z^4 w, with dx wedge dy=-z^2 dz wedge dw;
- t=z^3, s=-w/3, with dt wedge ds=-z^2 dz wedge dw.

Thus ramification index three and different/relative canonical order two match exactly. Starting on P1 x P1, the successive charts y=z^n w_n for n=1,2,3,4 come from four blowups with centers over x=infinity, outside the original affine plane. The pulled-back source form has order n-2 in the nth chart, giving order two in the final chart. These boundary blowups preserve a tree dual graph.

The resulting map t=x^(-3), s=-x^4 y/3 is rational and has Jacobian one, but has no retained finite sheet over t=0. It is not a polynomial Keller map.

The f4 cotangent hostile genuinely supplies the missing retained finite sheet while preserving the same escaping local type. In the boundary coordinates q=1/z and p=z^4 w, its exact target coordinates are

Q=z^3/(1+z^4),
P=w(1+z^4)^2/(z^4-3).

Their Jacobian in (z,w) is -z^2, Q/z^3 tends to one, and P tends to -w/3. The finite branch q approximately Q, p approximately P is separate. This confirms the draft's claim of the same cubic escaping behavior, including its limitation to a rational global map.

The hostile establishes that individual valuation realizability, local volume preservation, and the numerical discrepancy equality do not exclude the packet. It does not establish simultaneous global polynomial realizability or A4 global monodromy.

### 6. Differential escalation

The dual derivations have polynomial coefficients under the Keller hypothesis, extend the target partial derivatives, and commute with separable field trace. If T=h^(-a)b with b a unit at the divisor and partial_w h nonzero modulo h, the leading term of partial_w T has order -a-1 and cannot cancel with the derivative of b. Induction proves the exact pole tower.

The local-finiteness consequence is properly conditional: a finite-dimensional span of polynomial iterates has traces with a common bounded denominator. The increasing pole tower contradicts that. The text does not assert local finiteness of arbitrary Keller dual derivations.

## Computational audit

### Submitted controls

Both copied runs passed **808 exact checks**:

- `python3 submitted_controls_copy.py`
- `python3 -O submitted_controls_copy.py`

Normal and optimized stdout and JSON are byte-identical. All **7** packet-file checksums in the original manifest and all **12** repository-input checksums in SOURCE-DEPENDENCIES.json match. The file `dependency-hash-check.json` records the individual checks. Matching a dependency hash identifies the reviewed input; it is not an independent proof of all its contents.

### Independent controls

The independent script uses direct binomial sums over explicitly specified root packets, then compares with coefficients of the logarithm of the root polynomial. It does not import or call the submitted Newton or companion-matrix implementations.

It passed **575 exact checks**, covering:

- 209 multisets of nonzero leading residues from {1,-1,i,-i}, through cardinality six, including multiplicities;
- 34 valued-packet cases, including mixed slopes, fractional pole orders, translated roots, leading units, a merged quartic packet, and repeated nonprimitive observables;
- exact detection depths r, m, and d-1;
- retained-root discriminants via resultants, and explicit inertia permutations;
- even discriminant zeros caused by bounded projection collisions;
- all-escaping no-survivor examples requiring the d-th moment;
- the four boundary blowup charts and both rational local models.

Run `python3 independent_checks.py` or `python3 -O independent_checks.py` in this directory. Both outputs and JSON are byte-identical. This finite universe supplements the proof audit; it does not certify a global Keller realization or an all-degree theorem by enumeration.

## Citation check and scope

The four external sources were opened and checked directly:

- [Stacks 05K0](https://stacks.math.columbia.edu/tag/05K0): quasi-finite separated factorization as an open immersion followed by a finite morphism.
- [Stacks 03GT](https://stacks.math.columbia.edu/tag/03GT): the affine relative-normalization formulation of Zariski's Main Theorem.
- [Stacks 0AVB](https://stacks.math.columbia.edu/tag/0AVB): the height-one intersection statement over a Noetherian normal domain.
- [Stacks 02GL](https://stacks.math.columbia.edu/tag/02GL): geometric structure of étale fibers.

I also read the relevant finite-sheet/normal-generation proof in THM-2633, `01-canon/theorems/THM-2633-derangement-character-obstruction-and-d4-keller-exclusion.md`, and the boundary-tree construction in THM-3951, `01-canon/theorems/THM-3951-affine-plane-boundary-incidence-forest-and-equianharmonic-survivor-nonentry.md`.

No historical-priority claim was tested or is made. No repair is required for the audited statements. Optional exposition improvements would be to state the geometric-branch counting convention near the main findings and to display the exact f4 boundary-chart formula above. Neither changes a theorem or closes an open lane.

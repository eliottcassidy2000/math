# Overnight prime arithmetic: current checkpoint

**ONGOING SESSION, 20 October 2026.** This is a research-navigation checkpoint, not a final session close-out. The [publication ledger](PUBLICATION_LEDGER.md) pins seventeen reviewed packets at baseline `0ad127453c8ad4b2e31b96bb29a9b7e3d5ae56e7`. Later approved continuations may supersede the frontier below; work still under review is not a proved dependency.

**Main distinction:** several precise families are now closed, not general Collatz, complex JC, arbitrary odd-group symmetry refinement, or all prime-field definability. “Sign-even” means trivial edge-reversal character on every graph automorphism, not even vertex degree. Labelled maps, normalizer orbits and full isomorphism types are kept separate.

## 1. The Boolean-affine obstruction and its sharp degree successor

Read [the affine ideal theorem](../affine_circulant_ideal_obstruction_20261020/AFFINE_CIRCULANT_IDEAL_OBSTRUCTION.md), [validation](../affine_circulant_ideal_obstruction_20261020/VALIDATION.md), then [the Boolean-degree theorem](../circulant_boolean_degree_20261020/BOOLEAN_DEGREE_BOUNDARY.md) and [independent review](../circulant_boolean_degree_20261020/review/REVIEW.md).

**PROVED.** For an odd prime p=2m+1, a Boolean-affine multiplier-equivariant rule from labelled circulant tournaments to labelled sign-even circulant graphs inducing a type bijection exists exactly when m is odd or m=2: p=3 mod4 or p=5. This forbids all affine Boolean coefficients in the obstructed primes, even if coefficients use field addition or a chosen primitive root. It does not forbid nonlinear affine-group-equivariant algorithms: the [all-prime derivative-path construction](../circulant_necklace_bijection_20261020/CIRCULANT_NECKLACE_BIJECTION.md) already exists.

Every affine equivariant rule is one cyclic convolution h(R) after the antiperiodic derivative. Write m=e q, q odd, e=2^v2(m). In F2[x]/(x^m-1), the ideal U=(x+1)^e is the largest cyclic linear code inside the even-primitive-root target. A rule reaching the empty graph is trapped in U; the adjacent-pair word 1+x is outside U when m is even and m>2. This is a necessary-nonlinearity theorem, distinct from the earlier abstract multiplicative-group naturality obstruction.

The exact optimal affine type coverage is also proved. If A(m) is the number of target necklace types, the optimum is A(m) for odd m, [A(m)+A(m/2)]/2 for v2(m)=1, and A(m)/2 for 4|m. At p=13 this is 4 of 6 types. At p=17 it is 8 of 16, or only 1 of 16 if the empty graph must be retained.

**PROVED, stronger successor.** If e=2^a>=4, every equivariant type bijection has Boolean degree at least ceil((e-1)/a). The bound is unbounded across actual primes: prime divisors of Fermat numbers suffice, with no assumption of infinitely many Fermat primes. Exact minima are two at p=13 and three at p=17, attained by explicit quadratic and seven-monomial cubic witnesses. At p=17 the successful labelled fibers have sizes 2, 4 and 16; it is not a labelled bijection.

**FINITE-EXACT controls:** the affine packet exhausts all h and input words for m<=12; the degree packet includes all 2^29 quadratic equivariant maps at m=8 and independent witness/degree checks. The all-prime lower bound is proved independently of that finite enumeration. General sharp degree minima, cubic sparsity optimality and first-order ring definability of some other rule remain **OPEN**; no general circuit-size or first-order lower bound follows merely from unbounded Boolean degree.

## 2. Complete mixed-prime symmetry census through twelve

Read [the complete census through ten](../odd_mixed_prime_census_20261020/README.md), [its audit](../odd_mixed_prime_census_20261020/AUDIT_AND_NEXT_TARGET.md), and [the eleven/twelve extension](../mixed_prime_census_extension_20261020/README.md) with [verification scope](../mixed_prime_census_extension_20261020/VERIFICATION.md).

**PROVED classification; FINITE-EXACT numerical census.** The universe is odd-order subgroups H<=S_n with at least two distinct prime divisors, up to ambient S_n-conjugacy. Embeddings and fixed vertices are retained. There are no classes for n<=6, then respectively 1,2,2,5,8,10 classes for n=7,8,9,10,11,12: 28 degree/class rows in total. This is a complete bounded census, not a selected-group list or a census of every graph.

Combining it with the inherited all-p-group theorem, ten is the least failure degree for any odd subgroup. The only failing classes through twelve are C15 on 5+3+1^f:

| n | f | tournament types | sign-even graph types | difference |
|---:|---:|---:|---:|---:|
| 10 | 2 | 32 | 31 | 1 |
| 11 | 3 | 172 | 168 | 4 |
| 12 | 4 | 1376 | 1360 | 16 |

The exact defect identity is tau_T-tau_E=F_E+Z-F_T, where F counts normalizer-presentation fusion and Z counts relatively-even normalizer orbits that become globally sign-odd. At n=12 it reads 16=34+14-32: the deficit is not just the 14 escaping graph types. In general Z counts restricted relative orbits, not automatically distinct global graph types.

**PROVED family mechanism.** For p=1 mod3, the diagonal D=Cp semidirect C3 on p+3+f points and the independent product P=(Cp semidirect C3) times C3 have the same fixed cube but different normalizers. For f=0,1,2, exact source fusion compensates sign escape: tau_T=tau_E=c_f a_p with c_f=2,8,32 and a_p the number of internal affine-multiplier tournament types. The all-f relative identities do not force equality, but the pronormal-base theorem separately proves it for 0<=f<=p-4.

**Selected FINITE-EXACT probes, not a larger census:** at p=7, n=13 both embeddings have 172/172 types; at n=14 each has 1376/1375. The [scope note](../mixed_prime_census_extension_20261020/FAMILY_SCOPE.md) records the weaker replay scope of these four probes. Complete classification at n>=13 and uniform structural treatment of further embeddings remain **OPEN** at this checkpoint.

**PROVED selected degree-thirteen successor.** The [double-triple boundary](../double_triple_symmetry_boundary_20261020/README.md), with [completed independent review](../double_triple_symmetry_boundary_20261020/REVIEW_STATUS.md), structurally proves equality at n=11,12 for H_d=C5 times synchronized C3 on 5+3+3 and H_i=C5 times C3^2 with independent triple cycles. At n=13, adding two fixed vertices gives exact deficit six for each embedding. The selected FINITE-EXACT totals are 768/762 and 256/250 respectively; these are not all degree-thirteen subgroups.

The exact proof retains the extension predicate on eligible C5 copies: four globally even fusion types plus two escaping relative types give 4+2=6, with no source C5-copy fusion. These C5-based defects must not be confused with the synchronized H-normalizer counts F_T=256,F_E=4,Z=258. A useful proved side result is the sharp regular-tournament module bound 3|A|<=d for a proper nontrivial homogeneous module A; equality is a directed triangle of three regular modules. It rules out a single 5-cycle in a transitive odd degree-thirteen action without a full group classification. A complete degree-thirteen census remains open.

## 3. Square fields and weighted twin defects

**PROVED square-field family; FINITE-EXACT full q=25 census.** For p=3 or 5 mod8, odd s>=1, q=p^(2s), and H=F_q^+ semidirect C_((q-1)/8), the [square-field theorem](../affine_square_field_census_20261020/AFFINE_SQUARE_FIELD_CENSUS.md) gives tau_T(H)=tau_E(H)=1. Only the empty graph is sign-even. These are mixed-prime symmetries when q>9; q=9 is the p-group boundary. Semilinear fusion joins the two tournament multiplier types, while explicit odd-sign automorphisms eliminate every nonempty graph type. At q=25 all 16 fixed tournaments form one global type with automorphism order 75; the 16 fixed graphs form six global types, only one sign-even. This closes the prior named F25 test, not all odd non-p subgroups at prime-power degree, and says nothing analogous for cyclic rings Z/p^r.

**PROVED local formulas and sharp first boundary.** The [twin-quotient packet](../twin_quotient_defect_weights_20261020/README.md) expresses the triple-packing defect by a signed weighted false-twin cycle index. Every graph of order <=11 has defect -1,0,or 1; sign-odd graphs there have no negative defect. At order twelve exactly K_(6,6) and the independent triple blow-up of P4 have magnitude greater than one, with defects -2 and +2. The unique first sign-odd negative example is 2K_(3,3), with defect -1. Odd graphs have both positive and negative defects of unbounded magnitude.

Thus the naive all-order unit-weight split/join cancellation and “negative means even, positive means odd” partition are **REFUTED**. The inherited p-group theorem already makes the global defect sum zero. An explicit structural cancellation on a suitably decorated carrier remains **OPEN**; the results do not rule out every decorated pairing.

## 4. Collatz: the complete integer one-swap mechanical classification

Read [local carry and prime-power gates](../collatz_local_prime_power_gates_20261020/README.md), [sharp sparse norms](../collatz_one_swap_sparse_norms_20261020/SHARP_ONE_SWAP_CONTENT_BOUND.md), [the critical square-root-seven gap](../collatz_one_swap_critical_gap_20261020/CRITICAL_CLOCK_SQRT7_GAP.md), and [the positive-family closure with its audits](../collatz_one_swap_exclusion_20261020/README.md).

For a chronological parity word with length q, p ones and carry B, set D=2^q-3^p and g=gcd(|D|,B). The reduced denominator is |D|/g. Integer-cycle recognition is g=|D|; positivity additionally requires D>0.

- **PROVED local recognition:** any single chronological neighboring-rotation carry difference, together with D, recovers g exactly. A nonchronological difference can retain an extra subclock factor. Full prime-power valuations and the same word across prime tests must be kept: prime support alone and different local witnesses do not suffice.
- **PROVED sharp norm gate:** for a coprime mechanical shape and one cyclic adjacent unequal swap, g<=3^p-2 for p>=2, attained infinitely often at each p. This leaves only q=ceil(p log2(3)) for a positive integer cycle. At that critical shape, legal edit positions force the stronger g<7^(p/2).
- **PROVED positive-family exclusion:** for every q>=2 and 0<p<q, one such swap of the zero-intercept lower mechanical word gives no positive integer Collatz cycle except (q,p)=(2,1). The primitive proof combines the norm gap, a **CITED** explicit Matveev bound, exact rational continued fractions and a complete residual certificate. A separately audited repetition-factor obstruction removes coprimality. Crossing-cut swaps, rotations, duplicate occurrences and outputs returning to balanced necklaces are included.
- **FINITE-EXACT certificate, with its actual universe:** the closure has 18 large convergents with 100<=p<10^12 and 59 small coprime critical shapes / 3,548 swap occurrences. The earlier 6,492,980 occurrences through q=400 remain finite evidence, not an all-necklace or trajectory theorem.

Historical “one-swap OPEN” statements in the predecessor proofs and receipts are superseded only for this specified positive mechanical family. They are retained as evidence lineage, not rewritten as if the finite scans had proved the closure. The positive theorem alone does not close negative clocks, several-edit bounds or arbitrary parity words. **General Collatz convergence and general nontrivial positive-cycle exclusion remain OPEN.**

**PROVED signed continuation.** The [negative-clock theorem](../collatz_one_swap_negative_exclusion_20261020/NEGATIVE_ONE_SWAP_EXCLUSION.md), with [independent audit](../collatz_one_swap_negative_exclusion_20261020/audit/NEGATIVE_ONE_SWAP_AUDIT.md), classifies the same nonconstant mechanical one-swap family for D<0. The only integral shape is (q,p)=(3,2), D=-1: the two edited words 101 and 110 start at -7 and -5 and belong to the same cycle -5 -> -7 -> -10 -> -5. Every other negative instance is excluded, including noncoprime shapes via the separate repetition-factor proof. The constant parity cycle {-1} is outside the unequal-swap family and is not denied.

The negative proof retains its own signs and bounds: legal edits give g<6^(p/2), hence q=floor(p log2(3)) and a signed logarithmic gap; Matveev and exact continued fractions close the large range. Its small residual has 59 shapes / 3,348 swaps. This is not obtained by silently reusing the positive logarithm inequality. The [minus-seventeen boundary](../collatz_one_swap_negative_exclusion_20261020/MINUS17_THREE_SWAP_BOUNDARY.md) proves that the genuine -17 cycle is exactly three cyclic adjacent swaps from the mechanical necklace, by an explicit path and an integer-flow lower bound, independently checked on all 330 words at (q,p)=(11,7). It is a three-edit hostile, not a two-edit counterexample. Primitive two-edit negative classification and nontrivial positive two-edit exclusion remain open.

### Two-edit continuation: structural progress and an honest method barrier

The [two-swap packet](../collatz_two_swap_norm_barrier_20261020/TWO_EDIT_FRONTIER.md) proves a complete ordered gap test, moving-anchor carry formula and six canonical integral sparse norm forms for the stated primitive shapes. It retains undo paths, cut crossings, overlapping moves and distinct moved-one labels. Path length two is not necessarily minimum necklace distance two. The norm gate preserves full prime-power divisibility but loses the clock/common-root relation.

**FINITE-EXACT:** the exploratory lower/upper-critical slice through p=100 covers 199 shapes and 919,272 ordered paths, with no integral output of exact necklace distance two. This is not a complete residual universe: a two-edit critical-clock reduction has not been proved. The separate all-clock window 1<=p<=10, p<q<=80 covers 745 shapes / 102,684 paths. Together with the elementary q<=2p necessary bound for positive cycles and q<p log2(3) for negative clocks, it proves all-length exclusion of integral exact-distance-two outputs **only at weights p<=10**.

**PROVED method barriers, not cycle counterexamples.** A legal critical path at (p,q)=(4,7) has stripped norm 343>3^4 but actual content 1. A legal endpoint-chain family has raw norm pth-root limit 3 at unbounded coprime critical shapes, without resolving its 3-adic stripping there. The separately [reviewed stripped-norm addendum](../collatz_two_swap_stripped_norm_barrier_20261020/NEAR_CRITICAL_STRIPPED_NORM_BARRIER.md) makes a different explicit infinite choice: p=2*3^a, a>=2, and q is the smallest integer at least ceil(p log2(3)) with q=1 mod6. These are legal coprime positive-clock two-swap shapes with 0<q-p log2(3)<6. Their resultant is prime to 6, so the stripped modulus M equals its raw absolute norm and M^(1/p) tends to 3. Thus every proposed eventual bound M<=C c^p with fixed C>0 and c<3 fails on this nearby-clock family.

The last theorem is **not exact-critical-only**, not a lower bound on actual carry content, and not an obstruction to every content-sensitive exclusion method. The exact-critical stripped-norm asymptotic and unbounded two-edit integer-cycle exclusion remain open. The first nearby member has actual content 17, while the six other tested members through a=8 have content 1; those finite observations are not an infinite content theorem.

## 5. All-prime branches, codes and ring lifts: precise bridges only

Read [Prime-local branches, codes and lifts](../prime_local_branch_code_lifts_20261020/PRIME_LOCAL_BRANCHES_CODES_AND_LIFTS.md) and [validation](../prime_local_branch_code_lifts_20261020/VALIDATION.md).

**PROVED:** for oriented linear or affine gain-switching over F_(p^s) or Z/p^s, the p-regular restricted Burnside branch and its complement are nonnegative integers, except p=2,n=2 where both are 1/2. The fixed-space formula and p-integral partition terms combine with the **CITED** Frobenius root-count theorem: Z[1/p] intersect Z_(p)=Z. This supplies numbers, not a canonical partition or matching of types. The F2 affine carrier is tournament switching; its linear carrier is graph switching, whose complementary branch is different.

The same packet gives a proved classical Pascal/repeated-root-code identification and exact minimum-distance staircase, and the exact order p^(s+ceil(log_p F)-1) of a size-F unipotent Jordan block over Z/p^s for F>=2. It also gives all state-cycle counts at F=p+1 modulo p^2. Equal reduction modulo p can lose period; finite fields and finite rings are not interchangeable. CRT is lossless here and selected local cycles have pairwise coprime periods, so total cycle counts multiply. None of these facts gives Collatz descent or glues unrelated one-prime symmetry matches.

## 6. Source conductor and the nodal-cylinder exclusion

Read [Source conductor and collision multiplicity](../source_conductor_collision_multiplicity_20261020/README.md) and its [independent review](../source_conductor_collision_multiplicity_20261020/review/INDEPENDENT_REVIEW.md). The original proof's dated review-pending text is superseded by that completed audit.

**PROVED in the stated actual primitive monic Keller graph gauge:** the classical divided-difference D-resultant equals the normalization-conductor generator, the collision module's zeroth Fitting determinant, and R_Z(P,Q,x)/J after the specified monic normalization. Its norm is the square of the order-index determinant. Each source-line specialization has conductor degree twice its total affine delta invariant; no blanket conductor base-change principle is assumed. The minimum nonzero generic conductor degree two reduces the graph image ring to a polynomial cylinder over an ordinary nodal curve.

The resultant identities are explicitly credited to classical work, not claimed as new priority. Exact conductor controls cover 30 producer families / 293 checks and 14 independently constructed families / 125 checks; collision-bearing controls are explicitly non-Keller.

**PROVED successor, with CITED embedding input.** The [nodal-cylinder period obstruction](../nodal_cylinder_period_obstruction_20261020/README.md) and [independent audit](../nodal_cylinder_period_obstruction_20261020/review/AUDIT.md) exclude a generating Keller pair for E=C[x]+(y^2-1)C[x,y]. If C[x,P,Q]=E, the conormal determinant lambda is a nonzero constant and the exact polynomial period integral from -1 to 1 of P Q_y is (8/15)lambda. Its derivative gives zero integral of the Jacobian, contradicting a constant nonzero Jacobian. Generation is essential: membership, generic generation, or an invertible conormal matrix alone is insufficient.

The proof uses the [source-verified Abhyankar–Sathaye nodal-cubic embedding theorem](../nodal_cylinder_period_obstruction_20261020/SOURCE_VERIFICATION.md) fiber by fiber, not an unproved relative ambient automorphism over C[x]. Combined with the conductor classification, it proves that a noninvertible Keller map under the exact primitive monic source-gauge hypotheses must have generic source-conductor degree at least four, equivalently generic affine delta at least two. This is **not** a generic sheet-number bound, not a claim that every exceptional slice has delta>=2, and not complex JC. The next open sector is conductor degree four, not a reopened generic sheet-degree-four problem. Nodal finite controls have 285 checks across 17 producer generating families and 491 checks across 36 independent families; they do not prove the imported all-embedding theorem.

## 7. Active frontier and inherited routes

The next work is deliberately separated from the proved record:

- **OUTSIDE THIS CHECKPOINT / NOT IMPORTED:** later cylinder, nonprimitive mechanical two-edit and low-delta continuations belong to the follow-on queue. Subsequent approvals do not automatically enlarge this dated seventeen-packet checkpoint.
- **INCOMING, PRE-AUDIT / NOT A CHECKPOINT DEPENDENCY:** concurrent [THM-4724](../../../01-canon/theorems/THM-4724-affine-fixed-point-criterion-restricted-twisted-burnside-for-every-submodule-even-n-switching-half-and-the-dc-obstruction-to-uniform-circulant-rules.md) and its [odd-prime note](../odd_primes_20261020_affine_criterion_switching_even_n_uniform_rules.md) now have author status **PROVED (pre-audit)**, not a reserved stub. Their fixed-point, switching, derivative-family and power-sum claims are preserved as incoming work, not independently audited or imported into this seventeen-packet checkpoint. [HYP-9341](../../hypotheses/HYP-9341-what-remains-open-after-thm-4722.md) retains the incoming update alongside the reviewed reconciliation.
- **OPEN:** full n>=13 mixed-prime subgroup census; structural decorated defect cancellation; sharp general Boolean degree; first-order definability beyond the particular parity projection; unbounded two-edit Collatz exclusion, actual-content bounds and exact-critical stripped-norm asymptotics; p-singular gain formulas and higher ring state-period strata.
- **OPEN flagship boundaries:** general Collatz, complex JC and LRC(14). The results above do not change those statuses.

The common useful lesson is to retain the missing coordinate: primary height for affine rules, embedding multiplicity and sign for symmetry descent, legal edit position and a common full-valuation witness for Collatz, ring lift for exact periods, and multiplication/conductor multiplicity for source geometry. These are proved local connections with different carriers, not a universal equivalence between the problems.

Prior route: [Prime/symmetry session map](../prime_power_bijections_session_20261020/README.md), [geometry/arithmetic map](../geometric_arithmetic_session_20261020/README.md), and [earlier map chain](../prime_power_bijections_session_20261020/PUBLICATION_LEDGER.md#prior-session-maps). No theorem identifiers are reserved by this checkpoint.

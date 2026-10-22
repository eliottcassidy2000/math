# Session29 independent geometric audit

Audit date: 2026-10-22. Verdict: **ACCEPT mathematical claims within the stated scopes**, with the minor definition and testing qualifications below. No repository files changed by this audit. Independent checker normal and optimized results are included; producer checkers were copied into this directory before execution so their output did not mutate producer artifacts.

## 1. Direction-enriched Dehn tensor

The claimed common-frame dissection identity is correct. A precise elementary domain is compact three-dimensional polyhedral bodies admitting finite face-to-face convex-cell dissections, with disjoint interiors. More general finite polyhedral dissections can be refined into this setting. At an open subedge of the common refinement, intersect with a small perpendicular disk. The sectors of the piece interiors partition (a) the original dihedral sector, (b) a half disk at a new boundary seam, or (c) the whole disk at an interior seam. Their angles therefore sum respectively to the original angle, pi, or 2pi. All contributions use the identical direction u of that subedge and its identical length. Hence the tensor relation is an exact Q-bilinear identity modulo pi Q. If the subedge lies in the interior of a face of a piece, its formal angle is pi and adds zero; a subedge through a piece interior has formal angle 2pi and also adds zero. Introducing those zero terms handles arbitrary common refinements. Length subdivision is independent and additive. Vertices are measure-zero endpoints and create no edge term.

For nonconvex bodies the phrase “interior dihedral angle” must allow reflex angles and, for nonmanifold situations, a sum of local sectors. The present proof is unambiguous on convex cells and their finite unions; this explicit convention is recommended rather than silently including arbitrary degenerate polyhedra.

The regular tetrahedron nonvanishing witness checks: the six normalized directions evaluate under the harmonic p=x^4-6x^2 y^2+y^4 to a sum of -1; their common interior angle is arccos(1/3). Its ratio to pi is irrational: if it were rational, twice its cosine would be both a rational number 2/3 and an algebraic integer, which is impossible. Evaluation gives the nonzero simple Q-tensor -length tensor [angle]. No continuity of angle coordinates is assumed. The tensor is not an invariant under arbitrary independently rotated reassemblies.

## 2. Strengthened full-rotation no-go and quantifiers

The real-linear no-go in the report is correct and can safely be strengthened to Q-linear. Write H for the underlying real harmonic-quartic representation, and W_Q for the Q-span of all gh-h. W_Q is automatically stable under every real scalar a, since a(gh-h)=g(ah)-ah and ah is again a permitted vector. Thus W_Q equals W_R. The invariant inner-product argument gives W_R=H because an O(n)-invariant harmonic quartic is zero. Consequently every Q-linear O(n)-invariant map out of H is zero, not just every real-linear map.

For the actual direction-enriched target H tensor_Q B with B=R/pi Q and trivial rotation action on B, the same conclusion holds: express any h as a finite sum of differences gh_i-h_i, then tensor that identity with b. All pure tensors, hence the whole target, lie in the Q-span of rotation differences. This avoids an implicit field mismatch between the first-factor real representation and the Q tensor product. The target also has a natural real-vector-space structure through H if desired.

The no-go concerns linear/additive observable targets and their rotation coinvariants. It does not show that all framed geometry is zero, that spectra vanish, or that arbitrary nonlinear invariants are impossible. The continuous additive angle-factor obstruction is separately correct by density of pi Q. The scaling obstruction is valid for a map required to intertwine degree-one and degree-four dilation actions; a fixed normalization is outside that hypothesis.

## 3. Gram-cycle identities and sharp bounds

Independently expanded the rank-one cyclic trace identity and the mixed cross-Gram formula. With v_x=xx^T-(m/n)I, K=(x dot y)^2-m^2/n and c=2Nm^2/[n(n+2)], the report's constants in tr(A^2) and tr(A_X A_Y) are correct. The pair formula uses spherical second moments to sum (x dot y)^2. The identity tr A=0 itself follows from equal norms and the chosen c; the stronger identification as a harmonic quartic defect uses the second-design hypothesis. Antipodality is correctly retained when inferring a full spherical 4-design from even-moment vanishing.

The inequalities

-sqrt((d-1)V2/d) <= lambda_min <= -sqrt(V2/[d(d-1)])

are correct for d>1, including V2=0. The extreme patterns are (-(d-1),1,...,1) and (d-1,-1,...,-1), up to scaling. Thus they are sharp for arbitrary trace-zero symmetric operators, not necessarily for geometric shell operators. The Hessian sign corollaries additionally require b(t)>0 in the pencil b(t)(r(t)I+A), which holds in the stated Gaussian modular setting. They are not valid for an arbitrary pencil with negative b.

The independent checker tests these squared inequalities on rational spectra in dimensions2 through14 and both extremal patterns. This is corroboration, not a replacement for the Cauchy/nonnegative-coordinate proofs. The claimed Gram PSD lower bound lambda_min>=-c is immediate and correct.

## 4. Rank16 spectra and uniform saddle proof

The two classical root-system spectra are correct, including all dimensions and the trace-free tangent convention. Independently constructed doubled roots, then formed the 480 by480 integer matrix C=16K. This is a different calculation from the producer's 135 by135 basis-form comparison. The nonzero Gram/operator eigenvalues coincide, so tr M^k=tr C^k/16^k for k>=1. For the binomial shift use tr I=135, not480. All intermediate integer bounds are explicitly checked against int64 capacity; both ordinary and python -O executions give identical JSON.

The resulting angular traces are exactly:
- E8 square: 0, 30720, 1146880, 132382720.
- D16 plus: 0, 30720, 1146880, 49807360.

Their fourth gap is82575360. This calculation does not require enumerating480^4 quadruples or trusting a numerical eigensolver.

The all-temperature proof is sound. Differentiating theta(t)=t^-4 theta(1/t) gives g(t)+g(1/t)=-4; negativity of both derivatives implies -4<g<0. The cross-block coefficient (t theta'/32)(t theta'+4theta) is therefore strictly negative. The single-E8 block coefficient h8=t^2 theta''/40+t theta'/8 is positive for t>=1 by minimum norm2 and 2pi/5>1. For t<=1 its transformation h8(t)=t^-4 h8(1/t) follows by applying Poisson summation to exp(epsilon H/2)E8 and differentiating twice; duality sends epsilon to -epsilon and second derivatives preserve the sign. This t^-4 is the single-factor identity; the full rank16 Hessian transforms with t^-8. There is no inconsistency in using the single-factor identity before multiplying by the positive other theta factor.

The relative scalar direction is larger by96b>0, using the valid common scalar-plus-angular factorization. This proves exact inertia(71,64,0) for all positive t. For D16 plus the high diagonal branch is always positive; its low branch is negative at1 and positive in each far tail. Hence at least one reciprocal pair of degeneracies exists. Uniqueness or monotonicity is not established or implied.

The exact self-dual ratio was independently reduced from the stated Ramanujan derivative identities to40/3-216/y. The elementary bound9<y<27 suffices for all signs claimed. The producer checker separately verifies the rational upper estimate and Taylor estimate.

## 5. Tuple-theta interpretation and primary-source pins

The cyclic polynomial in the ordered k-tuple Gram entries, restricted to diagonal2, is a finite coefficient functional of genus-k theta. Thus a fourth-trace discrepancy is a genuine genus-four separating functional. Equality of lower angular traces alone would not prove equality of the entire lower-genus series; the report correctly cites that stronger classical input separately.

Verified primary sources on2026-10-22:
- Heimendahl et al., arXiv:2105.07868v2, dated24May2022, https://arxiv.org/html/2105.07868v2 : section5.2 explicitly states the two M spectra and equal ordinary theta series. Sections1.1.2 and5.2 describe the all-temperature E8-square saddle assertion as numerical evidence in that paper. Its existence there does not establish that the elementary proof is new, and the report appropriately makes no novelty claim.
- Oura et al., arXiv:0811.2259v2, https://arxiv.org/pdf/0811.2259 : PDF page16, discussion beginning with Witt's problem, explicitly gives equality through genus3 and separation in genus4. The genus-four difference is the Schottky form and vanishes on the Jacobian locus. Full Siegel-space separation, not separation on that locus, is the asserted result.

## 6. Reproduction and limits

The copied producer checkers both passed. The Dehn producer uses Python assert, so optimized execution would remove many tests; no -O robustness claim should be attached to that checker without replacing asserts. Its stdout also precedes the final two checks, while its final JSON contains them; this is a reporting-order quirk, not a failed test. The rank16 producer uses require checks and avoids that issue.

Our independent_check.py uses explicit require throughout and gives byte-identical ordinary and -O JSON. Hashes of audited inputs and generated audit artifacts are in MANIFEST.json. These checks do not prove general theorem statements by finite enumeration; the analytic arguments above supply those proofs. No actual rank88 angular eigenvalues, complete tournament invariants, Collatz convergence theorem, or Hadwiger–Nelson result follows from this audit.

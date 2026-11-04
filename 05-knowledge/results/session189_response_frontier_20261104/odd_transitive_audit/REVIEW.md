# Independent review of mixed-quartic tomography and conversed-child integration

**PASS. The complete frozen proofs, named banks, producer replays, independent normal/optimized receipts, and physical-label controls are approved within the scope below.** This review covers both frozen producer packets:

- Odd-transitive proof: SHA-256 `4f0123e3244b4515d80122d23e8bb62cd671c569a411d3ad9602043ea2a25a6a`.
- Conversed-child anchor proof: SHA-256 `492fe4d4f81e332437f70ce9db9b69b785745e779c01dbbc19afd78360892111`.

No mathematical correction to either producer was required. All statements remain source-relative and retain their stated hypotheses. No arbitrary-prime closure or unrestricted reconstruction follows.

## 1. General two-channel mixed quartic

The identity is valid for arbitrary internal tournaments, not only transitive ones. The eight disjoint three-path-cover terms have positive quartic signs: each port supplies one negative outgoing-incidence factor, so their product has positive sign. Both port orders and both orders of each marked pair are retained.

Eliminating exterior runs gives the stated internal path matrix M(tau), endpoint factors alpha C and beta R, cross entries R_q rho_i and ell_i C_q, and exterior start/end corrections beta H ell_i and alpha H rho_i. Every support multiplication is squarefree and joint; independent endpoint marginals cannot replace these expressions.

The derivative constant and orientation are correct. In the squarefree algebra, M'=M J M, hence M'_uv=R_u C_v, C'_j=H C_j, and R'_j=H R_j. Differentiating the four terms defining X gives exactly

    X' = [A]{H X_before_extraction
             +(C_p R_q+C_q R_p)(C_r-R_r)}.

This is precisely the second exterior channel in the converse difference. It is not multiplied by an additional factorial or an excursion-order factor. The first channel has sign alpha S_i-beta E_i; the second has sign alpha beta(ell_i-rho_i).

## 2. Transitive kernel, odd cancellation, and even boundary

For source order 0,1,...,m-1, the source transitive simple-path matrices and the rank-one expansion are correct. Extracting the three marked variables leaves the stated expression with P_left, P_right, and D. The resulting geometric sum is an identity of formal power series; the finite-difference transform makes it a polynomial.

The conversion constant is exactly two. Newton's formula and the identity for the geometric sum of binomial(s,j) produce t^j, not t^j/j! or j! t^j. Equivalently, expanding the ordinary polynomial in powers s^a gives transform coefficient j! S(a,j), which our independent implementation uses. The degree drops to at most m-3 because the two leading ordinary powers cancel.

For odd L=m-2, the binomial pairing gives the constant ordinary polynomial one. Therefore the combined X is exactly two and its derivative is exactly zero. This cancels every power of the exterior excursion tau over the integers. The exterior capacity interpretation is correct: C_zi=[R]alpha S_i and C_iz=[R]beta E_i, with ordered nonempty cover slots and the gap arc untested.

For even L, the combined ordinary polynomial is 2s+1, giving X=2+4t and X'=4. The producer retains the resulting second-channel term and does not assert its general nonvanishing. The m=3 odd boundary and m=4 even boundary have the stated constants.

A nonzero signed binomial combination implies at least one nonzero named coefficient and then at least one actual finite-difference corner. Negative combination weights do not require a new physical operation. The count k=(m-1)/2 of odd-order quartic terms and the stated conservative actual-query bounds are correct.

## 3. Inherited ordered-witness positivity and the root-family consequence

We read the immutable S188 [TT3 exterior-skew proof](https://github.com/eliottcassidy2000/math/blob/8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7/05-knowledge/results/session188_structural_bridges_20261104/tt3_degree_tomography/PROOF.md) and its [final independent review](https://github.com/eliottcassidy2000/math/blob/8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7/05-knowledge/results/session188_structural_bridges_20261104/tt3_degree_tomography_audit/REVIEW.md). Its skew theorem concerns a virtual singleton root and the exterior, so it has no distinguished-child-size restriction. Equal COMPLETE witness profiles and log-concavity are used; they are not omitted.

Grouping the outside c-contact uniformly over X or Y gives exactly the sum of its individual saturated quartics. There is no extra Boolean alias contribution: the two a contacts and two c contacts already exhaust bidegree (2,2). Prescribed fixed exterior bits cannot contaminate that saturated coefficient either.

The prime-root original-module argument is valid, including crossing descendants. If the core Hamilton count increases, the core-zero query, with the prescribed exterior bits, separates by hereditary profile monotonicity and a positive first ambient substitution coefficient. Otherwise every relevant fixed child has the same complete profile. The target distinguished child has Hamilton count one and is transitive.

The previously reviewed distance theorem, including its exact fixed-ambient all-cover transfer, applies to the internal quartics. The frozen distance proof has SHA-256 `dcaad9182ddea6a9520a1fd15deaa9f32bb60b1f18e9e54dae52b2dc3619d824`. It leaves only the original physical transitive order or its whole converse. The whole-converse branch is excluded by the compressed anchor lemma reviewed below. In the unchanged branch, whole-core saturated-quartic invariance and exact equal-profile restoration of the other uniform children reduce to the local distinguished-child converse. The binomial identity and ordered-witness inequality then separate.

The witness union stays a module after the distinguished child is removed, even after adjoining the fixed ambient, because ambient names see the whole core uniformly. No exterior operation is absorbed into the restricted root-family theorem. The result applies to original source modules contained in the core, with its prime root selected oddly. It is not pairwise reconstruction among arbitrary outputs.

## 4. General conversed-child anchor integration

We checked the precise inherited [U2–U4 uniform-compression theorem](https://github.com/eliottcassidy2000/math/blob/8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7/05-knowledge/results/session188_structural_bridges_20261104/returns_next/UNIFORM.md). Leaving every physical name of the distinguished child as an individual singleton is a legal partition, regardless of that child's shape. The other whole children remain compressed module blocks. The chosen anchor is uncovered in that partition; it is not required to lie outside every original module.

The target distinguished child already equals the source child's literal named converse. At a query uniform on every other child, restoring those other target children to their source converses preserves every augmented-cover coordinate exactly over the integers. Only after that restoration is the whole-converse parity theorem used. Thus no quartic of an arbitrary target comparison is discarded.

The compressed count q=|D_z|+k-1 and battery bound q+2=|D_z|+k+1 are correct. The parity coefficient d counts blocks, not physical vertices. The same physical exterior masks are used for the target and restored whole-core converse. The singleton-child case is included. The total-order-two boundary is expressly separated from the order-at-least-three U2 identity, rather than incorrectly applying E1 there.

Under full named-response equality for an original-module XOR target, other highest changed nodes are restored using exact profiles and uniform exposure to the chosen ports. Consequently no child at a highest changed node can be its source whole converse. Combining this with distance tomography for transitive children of order at least five removes the converse ambiguity and leaves the literal source order. The quotient-only candidate remains unresolved by this lemma; its existence or invisibility is not asserted.

## 5. Complete independent named-bank recount

`bank_audit.py` imports no producer implementation. `count.cpp` independently counts actual augmented Hamilton paths by a predecessor subset recurrence, whereas the producer uses a forward recurrence. Every source and target physical name and port orientation is retained.

All 101 graphs across the two producer packets were fully recounted. All 173,831 named coefficients and supporting small-support query values match exactly, including all zeros:

- Degree zero: 101.
- Degree one: 1,570.
- Degree two: 12,277.
- Cubics: 51,368.
- Quartics: 108,515.

The independently reconstructed data files are byte-identical to all five odd-transitive bank parts and the complete anchor bank. This is not merely selected-coordinate verification or reuse of the producer's answer as the count. Only the physical rows and specified query basis are supplied to the independent counter.

## 6. Independent kernel, physical-query, and family controls

`kernel_audit.py` builds the complete squarefree simple-path matrix by subset dynamic programming with polynomial edge weights [arc]+t. It obtains X from two disjoint path factors and obtains the derivative channel independently from the three-factor expression. It checks 159 internal kernels, including all labeled three- and four-name internal tournaments and transitive controls through eight names. Ordinary monomial expansion followed by Stirling transformation independently verifies the closed transitive kernel and binomial combinations through order 61. These controls support the structural proof; no all-order conclusion is inferred from the finite range.

The same program checks 55 additional mixed quartic contrasts with genuinely nontransitive internal children and arbitrary declared exterior incidences using 1,760 actual physical queries. Both exterior channels are retained. These controls test the general-module lemma, not just its transitive specialization.

`queries_audit.py` independently recounts all 672 retained grouped physical queries and all 84 anchor source/target/restored-converse query counts. It checks physical masks, fixed exterior bits, full graph labels, all source/local-converse relations, the actual fully covered original-module quotient targets, the source ordered-witness modules and incidences, complete matching log-concave witness profiles, and complete matching child profiles for the anchor examples. The local and quotient-only grouped contrasts agree exactly. The positive differences are 176,176,240 in the three producer witness controls, including a nontrivial strongly connected ambient.

`family_audit.py` exhausts all 4,096 root-containing original-module families of the named source C3[TT5,TT2,TT2], in both the local and a cyclic fixed-ambient context. The ten original TT5 intervals and the two pair modules are verified to be source modules. Interval XOR covers all 1,024 possible distinguished-child orientations; both pair selections are retained independently. All 8,192 family/context comparisons are detected by one fixed deduplicated 180-query battery.

Those controls exercise every branch: 7,232 scalar separations, 944 internal-order separations, eight unchanged-child mixed separations, and eight conversed-child anchor separations. Four additional complete physical named banks retain 13,210 coefficients, and 15,258 literal queries support them and the scalar cases. Every physical transport of a transitive target is explicitly checked on the entire target graph, keeps both ports and exterior names fixed, and transports the query with the names. The two pair selections leave the selected query uniform, as required for exact restoration. This is not an independent unlabeled relabeling of the comparison.

## 7. Reproduction and boundaries

All checks use exact integers and explicit conditions that remain active under optimized Python. At most fifteen vertices occur in the producer augmented banks, and at most thirteen in the additional family banks; Hamilton counts and sixteen-term finite differences fit signed 64-bit arithmetic. No software is installed.

Normal and optimized runs of each independent program must produce identical receipts. `replay_audit.py` copies only each producer's frozen allowlist to temporary directories, regenerates the controls, runs both producer verification modes, and checks every resulting allowlisted byte against its pin. It never changes the producer packet. The complete source and artifact pins are retained with the receipts.

Approval includes the general two-channel identity, exact odd/even kernels, the conditional odd-transitive original-module root-family theorem, and the general conversed-child integration. It does not remove the ordered-witness/equal-profile/log-concavity hypotheses, assert the even residual is nonzero, prove arbitrary-prime closure, exhibit a collision, or extend source-relative detection to arbitrary pairs of orbit outputs.

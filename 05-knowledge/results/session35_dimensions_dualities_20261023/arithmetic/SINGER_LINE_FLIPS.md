# Singer57 and the parity boundary of projective line flips

**PROVED:** all-order projective-plane parity obstruction and exact quadratic triangle-observer continuation theorem. **FINITE-EXACT:** explicit Singer planes q=2,3,7 and independent full Walsh transform at q=3. Classical Singer construction is not claimed as new. No Hamiltonian-count extrapolation or metagraph connectivity conclusion.

## 1. Inheritance and explicit connection contract

The recovered source is `05-knowledge/results/session33_continuation_states_20261023/FANO_CONTINUATION_SPECTRUM.md`: seven disjoint triangle-edge blocks give128 regular tournaments, with a special Hamiltonian observer supported on seven fourth-order Walsh characters. Its all-r coding extension expressly does NOT construct further tournaments. We now extend the actual projective edge-block operation, not its special Hamiltonian identity.

Least-used prior route: **THM-3234, `01-canon/theorems/THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate.md`**. Its Singer action is affine on F_p² and its added Heisenberg action generates AGL2, losing a pointed head. The present Singer action is projective on PG(2,q), a DIFFERENT carrier. What transfers is the need to check enlarged symmetry before preserving a distinguished orientation.

Canonical hostile: **THM-1131, `01-canon/theorems/THM-1131-pg213-one-twelfth-gap-refuted.md`**, whose final section separates a Singer difference set from a tournament skew half-set. Here the same mismatch is8 versus28, not14 versus91. The other hostile is parity: a K8 cannot be a regular tournament, so the Fano local directed-triangle mechanism cannot extend unchanged.

Board: finite-field trace / additive CRT / independent line reversals / score parity / triangle Walsh spectrum / allowed symmetry. The new result follows by preserving the operations and discovering the parity boundary rather than identifying all objects of size57.

## 2. A concrete (57,8,1) difference set

Use F343=F7[t]/(t³+t²+1) and alpha=t+5t². The checker verifies irreducibility (no roots, since cubic), and primitivity by the proper prime-divisor order tests for342. Projective points are the multiplicative cosets F343*/F7*, identified with Z57 by i->[alpha^i]. The trace-zero projective line is

D={8,15,19,38,48,50,51,56}.

Every nonzero residue of Z57 appears exactly once as an ordered difference of two distinct elements of D. Thus its57 translates are the lines, each of size8; every vertex pair lies on exactly one line. Their57 copies of K8 partition all1596 edges of K57 because57*28=1596.

Classically this is the Singer trace construction. For a prime power q, the trace pairing on F_(q³)/F_q is nondegenerate. Trace-zero is a2-dimensional F_q subspace, giving q+1 projective points. For a nontrivial multiplicative coset beta, the functionals Tr(x) and Tr(beta*x) are independent, so their common kernel has dimension1. This gives exactly one point common to the two lines, proving the difference-set property. The cited finite-field trace facts are standard; the explicit q2,3,7 instances are reconstructed independently in code.

Historical primary reference: James Singer, *A theorem in finite projective geometry and some applications to number theory*, Trans. AMS43 (1938),377–385, DOI https://doi.org/10.1090/S0002-9947-1938-1501951-4. The publisher PDF fetch returned403; no claim of inspecting its full text is made. A verified modern primary paper, Mészáros–Rónyai–Szabó, https://arxiv.org/abs/1908.05591, studies planar Singer difference sets through norm-one groups. None of its projective-norm-graph conclusions is imported here.

### The exact CRT bridge and its loss ledger

The chain of group isomorphisms is

F343*/F7* --chosen alpha exponent--> Z57 --CRT--> Z3 x Z19.

It preserves cyclic multiplication as addition of exponents, translation of D, and ordered-difference uniqueness. It does NOT preserve addition in F343, projective incidence without retaining D, or multiplication in the ring Z57. The generator and trace-zero subset are essential sidecars. This is not an isomorphism F343≈Z57, nor a field with57 elements.

Frobenius x->x^7 induces exponent multiplication by7, which is (a,b)->(a,7b) in the CRT coordinates. It preserves D; its fixed members in D are19 and38, and the remaining members form the cycles(8,56,50) and(15,48,51). Thus the factor3 has an actual role in the order-three Frobenius action, while the factor19 is the other cyclic coordinate. This does not identify CRT square-root sign flips with projective line flips: the former form a two-generator group; the latter have57 independent generators.

## 3. The odd-order parity obstruction

Let a projective plane of order q have v=q²+q+1 points and lines. Fix any tournament on its points. For each line L reverse all edges in its K_(q+1) block. The blocks partition edges, so these are v commuting independent involutions, giving an injective F2^v family of labelled tournaments. Reversing every line is tournament converse.

Let B be the point-by-line incidence matrix, reduced modulo2, and m the line-flip mask. At a point x on line L, reversal changes the score by q−2*outdeg_L(x). Hence for **odd q** the vector of score-parity changes is exactly Bm.

**Theorem P (PROVED).** For odd q, ker_F2 B is generated by the all-ones vector. Therefore exactly the empty mask and the full mask preserve every point's score parity.

Proof. Each point is on q+1 lines (even), and two distinct points share one line. Thus BBᵀ=I+J overF2. Since v is odd, I+J has kernel precisely the constants and rank v−1. Therefore rank B>=v−1; B1=0 gives the reverse bound. Its kernel has dimension1 and contains1. QED.

**Regularity corollary.** If the baseline is regular, exactly TWO tournaments in the full2^v line-flip family are regular: the baseline and its converse. Any other mask changes at least one score parity and therefore cannot preserve the common score(q²+q)/2. At q7 this is an all2^57 theorem, not an enumeration.

For even q, each line has odd cardinality, so choose any regular tournament on each line independently. Combining them gives a global regular tournament, and EVERY line reversal preserves regularity. There are2^v members in each such reversal family. This includes the q2 Fano construction but asserts no common Hamiltonian-count formula. For odd q regular local blocks are impossible because their sizeq+1 is even.

## 4. A universal triangle-observer spectrum and its minimal continuation

Write h(x,L)=2*outdeg_L(x)−q for x in L and0 otherwise. For line signs z_L in{+1,−1}, score is

d_x(z)=(v−1)/2+(1/2) sum_L h(x,L)z_L.

Let C(z) be the number of directed3-cycles. Each transitive triple has exactly one source, so

C(z)=binom(v,3)−sum_x binom(d_x(z),2).

The linear terms vanish because sum_x h(x,L)=0 in every line block. Each pair of distinct lines intersects at exactly one point w, giving

C(z)=Cbar−(1/4) sum_{L<M} h(w,L)h(w,M) z_L z_M,

where

Cbar=binom(v,3)−v*binom((v−1)/2,2)−(1/8) sum_{x,L} h(x,L)^2.

These identities hold for every baseline, regular or otherwise. The triangle observer always has Walsh degree at most2. For **odd q**, every h(x,L) on an incidence is odd and nonzero. Consequently **every one of the binom(v,2) quadratic coefficients is nonzero**, for every baseline.

**Theorem C (PROVED).** For odd q, two masks have identical directed-triangle counts after every word of specified line flips iff they are equal or global complements. Thus the minimal deterministic continuation state has exactly2^(v−1) classes. The real linear span of all translated observation functions has dimension1+binom(v,2).

Proof. A displacement a preserves all observations iff translation by a fixes the complete function C. Walsh characters are linearly independent, and every weight-two character occurs, so this requires a_L+a_M=0 modulo2 for every pair L,M. Therefore a is constant. Both constant displacements do preserve C because converse preserves3-cycles. The constant Fourier coefficient is positive: C is nonnegative and nonconstant (the displayed nonzero coefficients), hence its mean cannot vanish. Every nonzero Fourier coefficient contributes one independent character to the translate span. QED.

An independent mean formula is Cbar=[binom(v,3)−v*binom(q+1,3)]/4+sum_L C(T restricted to L). Noncollinear triples have edges in three different line blocks and become uniformly oriented, hence cyclic with probability1/4; collinear triples are only globally reversed and retain cyclicity. This also proves Cbar>0 directly.

At q7 this means2^56 continuation states, requiring56 bits for a lossless deterministic state representation, while the translation span dimension is1597. A scalar current count can have far fewer distinct values; that does not reduce the controlled continuation state. It is the triangle observable, not an asserted Hamiltonian observable, which proves this boundary.

At the Fano cyclic-line baseline every local h is0, so the directed-triangle count is constant14. Its special Hamiltonian fourth-order spectrum remains valid. Thus there is no conflict: changing q parity destroys local balance, and changing observer changes the continuation problem.

For any q, Hamiltonian-path count is invariant under full reversal, by reversing vertex words. Thus its line-flip Walsh expansion has no odd-weight terms. Nothing here shows that all even terms occur, that its continuation quotient has2^(v−1) states, or that its spectrum is quadratic. Those are separate questions.

## 5. Tournament symmetry guardrails

The full projective group acts2-transitively on distinct points: two independent representatives can be extended to bases and mapped to each other. Therefore some projectivity swaps any specified ordered pair, and no tournament orientation is invariant under the full group. In fact PSL3(q) already contains such a swap: exchange the first two basis vectors and negate the third, obtaining determinant1.

The cyclic Singer subgroup alone has odd orderv, so invariant circulant tournaments do exist: choose one of each pair{d,−d}, giving2^((v−1)/2) choices. But D has8 elements at q7, whereas the connection set needs28; D itself is not such an orientation. The Jacobi character modulo57 is symmetric under d->−d, as established in ARITHMETIC.md, and does not fix this.

## 6. Exact controls and reproduction

Run `python3 session35-arithmetic/check_singer.py` and also with `python3 -O`. Checks are explicit exceptions, not Python asserts. They construct all finite-field arithmetic, difference sets, line partitions, binary incidence ranks, and quadratic coefficients from scratch.

- q2:128/128 masks regular; triangle count14; binary incidence rank4; zero nonconstant triangle coefficients.
- q3: all8192 masks checked, exactly2 regular; triangle counts range16..91. An independent fast Walsh transform of the enumerated values gives exactly the constant plus all78 degree-two terms, and its support independently gives only the two claimed translation stabilizers.
- q7: all57 line blocks,1596 edge assignments, and1596 nonzero quadratic coefficients checked; binary incidence rank56. Sixty specified masks tested (empty, full,57 singleton flips, one mixed mask); exactly2 of these are regular. The all-mask assertion is Theorem P, not this sample. The tested triangle range7046..7714 is not claimed to be the full range.

The script produces `singer_receipt.json` and stdout stored as `singer_check.out`. Hypothesis/test loops: common57 cardinality -> explicit field/CRT map; regular Fano extrapolation -> even-block obstruction -> rank theorem; fourth-order-blindness extrapolation -> full quadratic triangle support -> exact continuation-state theorem; full projective symmetry ->2-transitive obstruction -> Singer-only orientation survivor.

Independent audit: `session35-singer-audit/AUDIT.md` reports PASS, using an independently constructed coordinate plane and separate Fourier/parity checks. Its q3 baseline is a different orientation gauge, so its triangle range55..91 need not equal this script's16..91. The theorems are baseline-independent, not the extrema.

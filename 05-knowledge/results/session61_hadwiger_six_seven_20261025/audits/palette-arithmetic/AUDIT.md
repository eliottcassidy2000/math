# Session61 independent audit: palette/deck and residue/CM bridges

**Verdict: PASS for the two producer packets, with their stated imported results and scope boundaries.** No mathematical correction to either producer `RESULT.md` was required. One consequential scope clarification was requested and verified corrected in the parent synthesis (§5): the ramified-prime-three three-color bound is for **K=Q(sqrt(-3))**, not for every larger field containing K. In particular Q(zeta165) has no conjugation-fixed prime above three.

This is an independent proof audit, not a new HN bound. No repository files were edited or pushed by this auditor. The parent controls publication and concurrent changes. Input byte pins are in `INPUT_PINS.json`.

## 1. What was actually audited

Inputs:

- `/tmp/session61-palette-deck/RESULT.md` and `check.py`;
- `/tmp/session61-residue-k3/RESULT.md` and `check_arithmetic.py`;
- parent `/tmp/session61-parent/SYNTHESIS.md`, especially §§3–5;
- the pinned predecessor statements, including actual Session60 two-torsion, canonical-character, genus-two and seven-line definitions, and corrected THM-4628's proof rather than its title.

Inheritance board: fixed-label bit lift; full versus projected gluing obstruction; perfect matchings and labeled Fano systems; norm-one integrality; CM splitting; quotient direction. Canonical hostiles are the sectionless bipartite triangle double cover, even nonzero vector holonomy, three-color source with K7 quotient, and even-order but non-fixed cyclotomic primes. The audit checked each bridge's map, target predicate, lost coordinate and stopping boundary.

The complex-surface classifications remain the explicitly inherited Session60 geometric theorems. Here their input integer d and its entire labeled-kernel census are independently rederived. This audit does not newly reprove surface minimality or export those classifications to characteristic three. The external six-color lower-bound theorem, its source and any Lean build are outside this audit's assigned claims.

## 2. Palette/deck proof audit

### 2.1 Exact fixed-ternary obstruction

For a finite simple graph and fixed q:V→{0,1,2}, an edge with unequal q-values cannot cause a collision in c=2q+b. An edge with equal q-values is proper exactly when b(u)+b(v)=1 in F2. These are precisely the equations δb=1 on the spanning subgraph H_q. A graph has no 2-cells, so every edge cochain is a cocycle; solvability is precisely [1]=0 in H¹(H_q;F2). The kernel is the functions constant on each connected component, giving 2^c solutions whenever nonempty. Isolated vertices and the empty graph obey this formula. The existential q must be retained when deciding six-colorability.

For arbitrary binary gain a, a graph-homomorphic section of the double cover is the same vertex potential. On a connected base, a nonzero cycle voltage joins the two sheets; if all cycle voltages vanish, the cover separates into two copies. In particular the all-one cover is always bipartite by its sheet coordinate, whether or not a section exists. C3→C6 is the correct hostile. No false implication from cover bipartiteness enters the packet.

### 2.2 The 15 matchings and the completion torsor

The even-subset quotient of F2^6 by the all-one vector is four-dimensional, and the dot product descends to a nondegenerate alternating form. Each of its 15 nonzero elements has a unique two-element representative. Distinct such representatives are orthogonal exactly when disjoint. Their sum is the class of the remaining complementary pair. Consequently a two-dimensional isotropic subgroup gives precisely one perfect matching, and every perfect matching gives one such subgroup.

With a marked seventh point 0, the three lines through 0 partition the other six points into three pairs. Any line not through 0 must meet each pair once. Once endpoints are framed, requiring each cross-pair pair to lie in exactly one transversal forces the four triples of a single parity. Both parities satisfy all Fano axioms and are distinct. Swapping one pair's endpoints fixes the matching, even its pair ordering, but exchanges these two completions. The stabilizer therefore fixes no point of the fiber, proving the asserted nonexistence of an S6-equivariant section. Ordering the pairs does not remove the obstruction; framing endpoints does.

The interpretation as Richelot kernels is legitimate for the actual marked genus-two Jacobian in characteristic zero. It does not assert that all targets are smooth genus-two Jacobians, that they are pairwise nonisomorphic, or that every permutation is an automorphism of a fixed curve. Bruin–Doerksen §4 explicitly supplies the pair-partition description of maximal isotropic subgroups. [Primary source](https://arxiv.org/html/0902.3480v2#S4).

### 2.3 What parity forgets

For fixed channel ordering, A=F2³ and p:A→F2 is surjective with two-dimensional kernel. A graph gain omega admits compatible full frames exactly when it is a coboundary; completion frames test the projected coboundary equation. The constant coefficient sequence splits as vector spaces, hence also after applying the graph cochain complex and cohomology. Equivalently, choose a spanning forest: one A-valued obstruction remains on each of b1 fundamental cycles. Projection retains one bit and discards two on each cycle. This proves the short exact sequence and dimension 2b1, including disconnected graphs.

The triangle with constant gain (1,1,0) has that same nonzero total gain, while every projected gain is zero. It is a genuine counterexample to replacing the vector obstruction by parity. It says nothing against recoloring the graph with a different ternary assignment. The producer correctly excludes channel permutations from this constant-coefficient statement.

Selecting one endpoint/root in each pair gives an injective A→F2^6/<1>: no nonzero word supported on the three selected positions equals the all-six word. Total parity of the image equals p. Thus the marked canonical-character restriction is exact; it is not a claim that palette permutations and root sign flips are the same geometric action.

### 2.4 First, second and third moments

For each chosen code K, a uniform element has coordinates s_i=a·v_i for uniform a∈F2³ and seven distinct nonzero v_i. Character orthogonality yields E((-1)^(a·w))=1 exactly when w=0. No single v_i or sum of two distinct v_i vanishes; a sum of three vanishes exactly on a Fano line. This proves all stated moment identities. Binary first and second moments determine the one-/two-coordinate distributions, so every statistic determined solely by these marginals is blind across all 30 kernels. This statement is about uniform subgroup distributions, not arbitrary tilted measures or affine cosets.

The three actual crossings 012,034,056 are part of the geometric input. Their cubic sum counts d. The 15 central matchings overlap the actual matching in 0,1,3 pairs: respectively 8,6,1 matchings. Multiplication by the two completion choices gives 16,12,2. No two-pair overlap is possible because the third pair is then forced. The numbers count labeled kernels/covers, not distinct abstract surface isomorphism classes.

The product of three canonical characters has complementary four-element support: each canonical support is 1+e_i, so three sum to 1+e_i+e_j+e_k. It descends through the projective sign quotient. The producer's projective-character guard is correct.

### 2.5 Boundary-relation and geometric guards

For fixed q, elimination of interior bit variables gives an affine parity relation or the empty relation. Varying q unions these relations and can destroy affine or pairwise completeness. The exact six-point Eisenstein graph has the stated inner triangle and three pendant ports; exhaustive enumeration yields 48 three-colorings, 24 NAE boundary assignments and two lifts each. The abstract K4 join has six-palette boundary relation 'at most two colors', hence 96 assignments with unrestricted pair shadows. Adding the port triangle is K7. Its K4 forbids an injective planar unit-distance realization, so it supplies no planar lower bound.

The six-equality-gadget equivalence is valid with nonempty six-coloring set and unit-separated omitted-edge ports. Adding their edge forbids six colors and needs at most one additional color; conversely deletion from an edge-minimal finite non-six-colorable unit graph yields exactly the required equality relation. These qualifiers are retained.

## 3. Arithmetic proof audit

### 3.1 Fixed primes exactly characterize norm-one integrality

For cP=P, v_P(cu)=v_P(u), so ucu=1 forces 2v_P(u)=0. Conversely P and cP are distinct maximal ideals, so select an integer a∈P\cP. Then c(a) is a P-unit, and a/c(a) has positive P-valuation. Its inverse is norm-one and has a pole. This proves both the integrality and unit iff statements, including the inverse-direction quantifier. An edge difference outside the local ring prevents both endpoints from lying in any one translated local ring, so changing the component basepoint does not repair this mechanism.

At a fixed unramified prime whose residue involution is q-Frobenius, reduction of w/c(w), with w a unit lifting any nonzero residue, is w^(1-q). Its image is exactly the norm-one subgroup of F_(q²)^*. Every finite residue edge therefore lifts individually. On a source component, the residue image is the additive subgroup generated by those norm directions; the finite norm graph is a disjoint union of translates of the resulting Cayley component. All its components have equal chromatic number. Thus requiring a proper coloring to factor through this reduction costs exactly χ(G_q), without any unproved connectivity assumption. The source graph itself need not require that many colors.

This is compatible with Madore's valuation/coset pullback mechanism and general quadratic-form extension. [Primary source, Proposition 3.2 and §6.6](https://arxiv.org/pdf/1509.07023). It does not import the obsolete 2015 plane lower-bound status.

### 3.2 Eisenstein splits, triangles and misleading quotient cliques

For p≠3, X²−X+1 splits exactly at p≡1 mod3; such primes are exchanged by conjugation. A fixed prime in a larger conjugation-stable field would contract to a fixed prime in K, contradiction. The exclusion persists under extension. For a fixed unramified prime in characteristic not three, the surviving order-three root is inverted by conjugation, forcing q≡−1 mod3.

Normalize a norm-graph triangle to (0,1,z). N(z)=N(z−1)=1 forces Tr(z)=1 and z²−z+1=0. In characteristic three there is one solution z=−1. Otherwise the two roots contribute precisely when q≡2 mod3, including characteristic two. This confirms that G7 is triangle-free and is not the lattice's split-prime K7 quotient.

In the omega-sixth-root convention, alpha=2+omega and its conjugate 3−omega both have norm seven. At omega=5 mod7 the first vanishes and the second does not; at omega=3 the reverse holds. The norm-one ratio therefore gives the explicit pole witness. The lattice's six unit directions map onto all nonzero residues mod2 and mod(alpha), producing K4 and K7 from its exact three-coloring a+2b mod3. Vertex/edge surjectivity gives no coherent clique lift and does not reverse χ(source)≤χ(target).

For **K alone**, the ramified prime above three is fixed and its residue conjugation is trivial. Norm-one edges reduce to ±1, giving the exact three-color field-plane bound in conjunction with the triangle. No such ramified-prime claim was made in the producer for every containing field.

### 3.3 Complete fixed-prime classification in Q(zeta165)

For an unramified rational prime ell, the decomposition group is the cyclic subgroup generated by ell mod165 and c corresponds to −1. Thus c fixes a prime iff −1 belongs to that cyclic subgroup, not merely iff its order is even. CRT decomposes the units as C2×C4×C10. Projecting to mod55, the only fixed classes are −1 of order two and 19,24,29,39 of order ten. Their half-orders are odd. The additional mod3 projection must therefore be −1. CRT gives exactly 164 (degree two), and 29,74,134,149 (degree ten).

At ell=3,5,11, inertia permits arbitrary action on the ell-factor, while the prime-to-ell conductor part has decomposition subgroup generated by ell. Hence the correct ramified test is −1∈<ell> modulo55,33,15. It fails in all three: half-order powers are respectively34,23,11 rather than−1. Thus there is no unexamined ramified fixed prime.

The degree-two class first reaches a prime at659: 164 and494 are even; 329=7·47; trial division through25 verifies659. Every degree-ten case has ell≥29, so q=ell^5>659. This is an all-prime argument, not an inference from the finite prime sieve. The old candidate109 fails the extra mod3 projection.

### 3.4 Spectral consequence and its exact scope

For an odd-q norm graph, an additive character indexed by a≠0 has eigenvalue obtained by summing psi(Tr w) over N(w)=N(a)=c. The exact number of w with trace t and norm c is 1−eta(t²−4c). The exact number of x≠0 satisfying x+c/x=t is 1+eta(t²−4c). Since a nontrivial additive character sums to zero over Fq, the two exponential sums are negatives. This derives the needed Kloosterman identity for every odd prime power, not merely the split form x²+y² at selected primes.

The remaining deep input is the explicitly CITED Weil bound |Kl(q,c)|≤2sqrt(q), c≠0, valid for finite fields, not residue rings. Conrad's primary article states precisely this scope in its introduction/formula(1). [Primary article](https://kconrad.math.uconn.edu/articles/kloosterman.pdf). No claim of reproving Weil is made here.

For completeness, Hoffman's implication is elementary: an independent-set indicator f with mean alpha has zero f^T A f. Writing f=alpha·1+g and using the least eigenvalue lambda<0 gives 0≥N[d alpha²+lambda alpha(1−alpha)], so alpha≤−lambda/(d−lambda). A k-coloring has some class of proportion≥1/k, hence k≥1−d/lambda. Here d=q+1 and |lambda|≤2sqrt(q), giving the stated lower bound. At q=659 the exact inequalities are

    4·659·12² < 660² < 4·659·13².

Therefore the real lower bound lies strictly between13 and14, and integrality forces at least14 colors. It is increasing for q>1. This is a lower bound for the specified **single-prime residue-factor method only**, while the same embedded field graph admits a seven-coloring by restriction from the plane. Neither a source lower bound nor an upper bound14 for the method follows.

### 3.5 CM and model scope

For p≥5 the displayed elliptic curves y²=x³−1 and y²=1−x³ have good reduction. In the inert class p≡2 mod3, cubing permutes Fp, so the character trace vanishes exactly. In the split class, expansion of (x³−1)^((p−1)/2) has a unique contributing positive exponent p−1 to the power-sum congruence; its coefficient is the stated signed binomial coefficient and is nonzero modp. The constant term sums to p=0 modp and the degree is below2(p−1), so no omitted term contributes. The second curve is the quadratic twist and has the same ordinary/supersingular status.

Katsura–Oort explicitly states the j=0 supersingularity criterion for p≠2,3 in §6, printed p.135. Only this elliptic criterion is imported. [Primary source](https://numdam.org/item/CM_1987__62_2_107_0.pdf). The geometric genus-two maps pull back the elliptic differentials to 2u du/v and −2 du/v, independent for p≥5. These computations do not identify the original rational genus-two twist with C6 over Q, do not transfer rational point counts between twists, and do not assert any positive-characteristic K3 Picard-rank classification. The producer correctly retains those boundaries.

## 4. Independent computation and hostile controls

`check_independent.py` imports no producer code. Its alternative methods are:

- exact-cover enumeration of the 21 pairs by seven triples, yielding all30 Steiner triple systems;
- complement-of-line construction of the codes and direct character sums;
- incidence-matrix elimination versus truth-table fibers on every simple graph of orders0..5 (1,100 graphs,59,810 binary-gain instances), plus5,422 fixed ternary cases;
- direct cyclotomic orbit tests, prime sieve and ramified decomposition tests;
- explicit finite-field arithmetic for q=2,3,4,5,7,8,9,11,13,16,17,25,27,49, checking irreducibility via the multiplicative law, circle sizes, triangle counts, and exact odd-characteristic trace-fiber formulas;
- elliptic point counts by square-multiplicity tables for all93 primes5..499, plus the genus-two quotient trace identity.

There are **79,859 always-active checks**. Normal and optimized outputs are byte-identical. Both original producer scripts were separately rerun in normal and optimized modes: palette11,508 checks; arithmetic1,853 checks. Original stored receipts agree with these reruns.

Three targeted mutants are preserved with their exact errors in `NEGATIVE_CONTROLS.json`: replacing fixed-prime membership by even order; replacing30 Fano systems by15; and deleting one of the two lost bits per cycle. All three fail normally and with `python -O`. The tests are finite controls, not substitutes for the arbitrary-graph/all-prime proofs above.

Reproduction:

    python3 -B check_independent.py > independent.normal.json
    python3 -B -O check_independent.py > independent.optimized.json
    cmp independent.normal.json independent.optimized.json

## 5. Parent synthesis disposition and operational scope

The parent synthesis was reviewed, not modified by this auditor. The parent applied the following changes, and the final §§3–5 were read back successfully:

1. Scope the prime-three field-plane bound to K=Q(sqrt(-3)). This avoids contradicting the later no-fixed-ramified-prime result for Q(zeta165).
2. Say the central K3 quotient's associated genus-two curve is **geometrically** C6; preserve the rational-twist distinction.
3. Scope the pairwise observer statement to the uniform one-/two-character marginals, and identify16/12/2 as labeled-kernel counts.

No mathematical blocker remains for the two producer packets or the corrected parent synthesis §§3–5. The final synthesis hash is recorded in `SYNTHESIS_READBACK.json`. Coordination command `python3 agents/processor.py --status --peek --limit 8` did not run its queue check because `.machine-id` was absent; no registration or repository file was invented. Team messages used the functioning collaboration tools. An optional PDF screenshot fetch failed with a cache-miss error; the CM criterion was read from the primary PDF's extracted text and no geometric visual claim depends on that screenshot. These are not theorem-validation failures.

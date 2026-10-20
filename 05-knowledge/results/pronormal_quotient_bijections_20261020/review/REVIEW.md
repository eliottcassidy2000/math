# Independent review: pronormal bases and normalizer-quotient localization

**Verdict: PASS.** Theorems 1 and 3, the fusion/sign descent Lemma 2, and the added pronormal-lifting Lemma 4 are correct as stated. No mathematical correction is required. The examples respect the distinction between normalizer-relative and global equivalence. The claimed finite counts pass independent bounded controls.

Reviewed source: `session12-prime-local/PRONORMAL_QUOTIENT_THEOREM.md`, SHA256
`c33157c00a065b2e470d59e67fcbe21e91084fb799343943186cd398fe136684`.

This review reads the repository and producer packet without modifying either. All review artifacts are staged separately.

## 1. Odd-orbit Fourier pairing

Let the odd group K act on the unordered pairs, with orbits O_i. Each orbit size divides |K| and is odd. A K-element cannot reverse a pair: a permutation of odd order cannot have a 2-cycle on its two endpoints. Thus each pair orbit admits two consistent orientations, and Fix_T(K) is an affine d-cube.

Choose one K-invariant tournament as the origin. If a K-invariant graph X is the union of the orbits indexed by a subset I, then on this affine cube

    phi_X(t) = c_X (-1)^(sum_{i in I} t_i),  c_X in {+1,-1}.

Indeed reversing one orbit changes the product by (-1)^|O_i| = -1. Therefore the restrictions indexed specifically by K-invariant graphs give every cube character exactly once, up to signs. Their inverse is 2^(-d) times the transpose. This is the required perfect pairing; the warning about arbitrary restricted characters is justified.

N(K) acts by signed permutation on this character basis. K fixes every tournament in the cube. It also acts trivially on each K-invariant graph's character: its graph-stabilizer sign is a homomorphism from an odd group to {+1,-1}. Hence both actions factor through Gamma=N(K)/K.

Over Z[1/2], taking invariants does not require dividing by |Gamma|. For an ordinary permutation basis, orbit indicators form a basis of invariant vectors. For a signed permutation basis, a sign-odd stabilizer forces a coefficient to equal its negative and hence to vanish; a sign-even orbit contributes one signed orbit-sum basis vector. This directly establishes the two invariant bases and the invertibility of their integer transition matrix M over Z[1/2].

## 2. Correct quotient weights and localization

Write D_T and D_X for the diagonal matrices of Gamma-orbit sizes. Orthogonality of the cube characters gives

    M^t D_T M = 2^d D_X.

Both matrices are square and invertible. Inverting/rearranging and using |orbit(x)|=|Gamma|/|Gamma_x| gives exactly

    M diag(|Gamma_X|) M^t = 2^d diag(|Gamma_T|).

Thus the weights in the source are the normalizer-quotient stabilizer orders, not full automorphism orders and not full normalizer-stabilizer orders with an unaccounted |K| factor.

For completeness, let a_i=|Gamma_Ti| and b_j=|Gamma_Xj|. The inverse relation is

    (M^-1)[j,i] = b_j M[i,j] / (2^d a_i).

At an odd prime q, both M and M^-1 are q-integral, while det(M) is a q-unit. A nonzero determinant term modulo q supplies a perfect matching of q-unit entries. On such an entry the inverse relation gives v_q(a_i)<=v_q(b_j). Taking determinants in the weighted identity shows equality of the sums of these valuations. Consequently every matched pair has equal valuation. This rederives the inherited matching step without assuming that every q-unit entry has equal endpoint weights.

For a Sylow q-subgroup Q of Gamma_T, use arbitrary lifts in N(K) of its elements. The lifts need not form a q-group, and no splitting of a quotient extension is needed. K acts trivially, so they induce a well-defined Q-action on the graph copies. The transporter cocycle and the orientation-character cocycle cancel, making each signed summand constant on a Q-orbit. Every non-singleton orbit has q-divisible size. A q-unit coefficient therefore forces a Q-fixed copy, and equal matched valuations make Q a Sylow subgroup of that copy's Gamma-stabilizer. Transport back gives Gamma-conjugacy.

A q-subgroup-presence predicate is determined by the Gamma-conjugacy class of a stabilizer's Sylow q-subgroup. Thus one matching preserves all the specified q-subgroup predicates and every Boolean combination. This argument concerns one quotient prime at a time.

**Optional exposition improvement only:** add the preceding sentence about arbitrary lifts explicitly to Theorem 1. The proof is already valid; this prevents a reader from mistakenly inserting a splitting assumption.

## 3. Pronormal descent and ambient embedding predicates

Both parts of Lemma 2 check out.

For fusion, suppose x'=g.x and both objects are K-fixed. Then K and gKg^(-1) lie in Aut(x'). A pronormal conjugator a lies in their generated subgroup and hence in Aut(x'). It can be chosen with ag in N(K), and (ag).x=x'. Therefore the intersection of a full type with Fix(K) is exactly one N(K)-orbit.

For signs, if X has an odd automorphism g, its sign-kernel A^+ contains K and gKg^(-1), because both are odd groups. The pronormal conjugator a lies in their generated subgroup inside A^+. Consequently ag normalizes K, remains an automorphism of X, and is still odd. Thus every globally odd K-fixed graph is already odd relative to N(K). The reverse implication is immediate.

Theorem 3 does not confuse ambient conjugacy with quotient conjugacy. If gHg^(-1)<=Aut(x) for a K-fixed representative x, pronormality provides a in Aut(x) for which n=ag normalizes K. Then nHn^(-1)<=N(K) intersect Aut(x). Conversely any such normalizer embedding is an ambient embedding. Since both subgroups contain K, passage to N(K)/K is exact. This proves the equivalence of the ambient H-presence predicate and the quotient H/K-presence predicate, including when H itself is nonpronormal.

Any type admitting H also admits K, so restriction of the sector bijection gives the stated equality tau_T(H)=tau_E(H). No classes outside the K-admitting sector are lost from that count.

## 4. Added lifting lemma

Lemma 4 is valid for arbitrary finite G; oddness is only needed for the subsequent tournament application.

With L=<H,gHg^(-1)>, pronormality of K gives a in L such that ag normalizes K. Then J=(ag)H(ag)^(-1) lies in L intersect N(K), and H,J both contain K. Pronormality of H/K in N(K)/K gives a conjugator in <H/K,J/K>, which lifts to b in <H,J><=L. Equal quotient subgroups give equality of H and bJb^(-1), since both contain K. Hence ba in L conjugates gHg^(-1) to H.

The proposed iterative construction therefore works. It must use full preimages of pronormal subgroups in the ambient normalizer quotient at each stage; the stated warning against an arbitrary abstract Sylow tower is appropriate.

## 5. Examples and hostile cases

### Six vertices

An independent implementation uses affine pair-orbit coordinates and full S6 relabellings. It reproduces:

| K action | Fixed objects each | Relative types each | Global types each |
|---|---:|---:|---:|
| One 3-cycle | 128 | 12 | 11 |
| Two synchronized 3-cycles | 32 | 4 | 3 |
| Independent C3 x C3 | 8 | 1 | 1 |

The first two global predicates intersect in the unique third type on both sides. Their tournament fusion multiplicities are respectively ten ones plus one two, and two ones plus one two. Each has exactly one relative-even/global-odd graph, with 9 and 6 edges respectively; the producer's explicit masks identify these as K_(3,3) and two disjoint triangles. All 720 vertex permutations are used for global equivalence and sign tests.

This is correctly presented as compensation beyond pronormal descent, not evidence that the nonpronormal C3 base satisfies Lemma 2.

### C15 at orders eight, nine, ten

At n=8,9 the displayed C5 is an ambient Sylow subgroup; at n=10 it is not. The explicit two-disjoint-5-cycle witness generates C5 x C5 of order 25 and admits no internal conjugator between the factors. The K_(5,5) part swap reverses 25 edges and does not normalize the chosen factor.

The inherited prime-module argument also gives the tournament counts directly. The C3 and regular C5 modules are prime and strongly connected, hence intrinsically unique when at most two other vertices remain. With r=0,1,2 roots the quotient presentations number 2,8,64, and only the last case has a free two-root swap, yielding 2,8,32 types.

Independently checking all 74 independent-block graph quotient presentations across n=8,9,10, using intrinsic full twin classes and their complete weighted automorphism groups, gives 2,8,31 sign-even graph types. Nonempty internal C3/C5 graphs have explicit odd automorphisms extending by the identity, so this reduction loses no even graphs. Thus the boundary counts are independently supported without an unrestricted S10 census.

### C7 on ten vertices

The independent affine-cube control obtains d=9, 512 fixed objects per side, 24 relative types, and |det(M)|=2^88. It checks the weighted identity entry by entry. More strongly, it verifies the signed localization congruence and existence of a fixed graph copy for all 300 coefficients nonzero modulo 3, and finds a full unit matching inside equal quotient-Sylow-action blocks.

The five extension counts are exactly 24,12,4,2,2 on both sides. The joint line-predicate histogram is:

- (false,false,false): 10
- (false,true,false): 2
- (true,false,false): 10
- (true,true,true): 2

The nonpronormal diagonal F21 witness independently has subgroup order 21 and generated-group order 63, with no conjugator in the latter. Conceptually its two diagonal quotient lines are distinct in the abelian C3 x C3 quotient, so this is also an elementary proof that Theorem 3 genuinely extends beyond the pronormal-H corollary.

The count's global interpretation uses the proved Sylow/pronormal descent. Neither source nor review claims an independent full S10 census for this C7 example.

### C5 wr C3 on fifteen vertices

The claimed value 1=1 is correct. P=C5^3 is an ambient Sylow 5-subgroup; each invariant tournament has circulant five-blocks and uniform cross-block arcs. All four internal choices are multiplier-isomorphic, and the block 3-cycle forces the quotient directed triangle.

Any nonempty internal circulant graph has an odd automorphism extending by identity, so even graphs have three independent blocks. The two block-cycle-invariant quotients are empty and complete. The complete quotient yields K_(5,5,5), with a two-part swap of sign (-1)^25. Only the empty graph survives. For the P-sector alone, the quotient has the ordinary three-vertex tournament/even-graph counts 2=2; the two graph types are the empty graph and K_(5,10).

## 6. Reproduction and scope

Run from the workspace root:

    python3 -B session12-prime-local-review/independent_controls.py > session12-prime-local-review/independent_controls.json
    python3 -O -B session12-prime-local-review/independent_controls.py > session12-prime-local-review/independent_controls.optimized.json
    cmp session12-prime-local-review/independent_controls.json session12-prime-local-review/independent_controls.optimized.json

The independent script imports no producer code or repository atlas and uses only Python's standard library. Its normal and optimized outputs agree byte for byte.

The producer's two scripts were also replayed in normal and optimized modes. All four replay files agree with the frozen producer outputs. This replay is separate from the independent implementation.

The all-order proof is accepted; finite controls are bounded checks, not substitutes for that proof. The result makes no claim of simultaneous control across quotient primes, equality of full automorphism groups, a natural labelled equivariant bijection, or computational efficiency in n. No broader census or literature-priority claim was attempted.

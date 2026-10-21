# Signed primitive tournament exceptions: independent verification and scope repairs

Status: PROVED relative to the cited Gluck theorem and Feit–Thompson; FINITE-EXACT controls. This is an independent audit, not a claim of external literature priority. The research and review were prepared separately before publication.

Integration status: the independently audited companion [structural proof](../session24-structure/STRUCTURE.md) now closes HYP-9353 by a new argument. Statements below about the inherited OPEN frontier describe the earlier source state, not the current result; they do not retroactively validate its premature completeness prose.

## Results first

1. The user's three-exception statement is correct when primitive means the faithful vertex permutation action: for a finite tournament of order at least three, precisely `C3`, `C5{1,2}`, and `P7` need three named colors after converse symmetries are admitted.
2. For tournaments of order at least three, ordinary vertex-primitivity and signed vertex-primitivity are equivalent. Every vertex-primitive tournament is self-converse. This clarification is absent from the incoming theorem.
3. The Paley conclusion holds for every finite-field prime power `q ≡ 3 mod 4`, not merely prime q: only `q=3,7` are exceptions.
4. Independent confirmation of the companion structural-reduction result: **every exceptional tournament is vertex-transitive under its ordinary automorphism group**. Consequently exceptions have odd order and are strongly connected. Intransitive, even-order, non-strong, and unequal canonical-module cases cannot supply missing exceptions.
5. Converse-fixed regular-subset orbits do not imply converse-fixed Hamilton-path orbits. For `P7`, all four ordinary regular-subset orbits are fixed, but only three of nine Hamilton-path orbits are fixed. For `C3[C3]`, the corresponding counts are two of two versus eleven of thirty-nine.

## 1. Inheritance and fresh-input audit

The prior reviewed packet is `signed_tournament_palette_rigidity_20261021` (initially located in `session23-publication/final-remote-layout`, then available in the refreshed math checkout). It already proves the universal three-color upper bound, the full-vertex faithfulness boundary, residual-path kernel warnings, the Paley7 obstruction, and the distinction between an equivariant family and an equivariant choice.

The refreshed incoming theorem is:

`01-canon/theorems/THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md`

Its frontmatter claims statements 1–5 as proved but also says **pre-audit**. These already include the orbit-pairing lemma, the exact primitive classification, the Paley prime-power deduction, strong equal-block substitution, and the triangle-tower family. They are inherited claims being independently audited, not new discoveries in this report. I independently confirm the primitive and Paley proofs. The present report adds the ordinary/signed primitivity equivalence and independently verifies the vertex-transitivity reduction developed in the companion structural-reduction report. The latter is one shared result, not a second novelty claim.

At inheritance, `HYP-9353` was explicitly OPEN in the only-if direction. The source `seven_point_pencil_fano_and_converse_gluck_20261021.md`, section 4, nevertheless says that THM-4743 decides the exceptions “only in the Gluck circulants and their triangle towers.” That prose incorrectly promotes the conjectural completeness assertion. The justified replacement at inheritance was: the primitive exceptions are classified, and the tower family is an infinite proved family; completeness was still open. The companion structural proof supplies the missing all-order necessity argument in this session.

The same incoming documents cite THM-2450 as an all-order cyclic-ternary grammar classification. That file itself is a historical finite-scope source:

`01-canon/theorems/THM-2450-rigid-self-converse-classes-are-cyclic-ternary-towers.md`

explicitly limits its own identification to an exhaustive census through order ten. Its old OPEN language is superseded by [THM-2477](../../../../01-canon/theorems/THM-2477-rigid-strong-classification-through-sixty-two.md) and the [all-order Hamiltonian-bypass proof](../../merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md#3-rigid-strong-tournaments-are-completely-classified): the strong path-rigid atoms are exactly `{1,C3}`, all path-rigid tournaments are ordered stacks of those atoms, and the self-converse ones are palindromes. Here “rigid” means `H(T)=|Aut(T)|`, one ordinary Hamilton-path orbit, not `Aut(T)=1`. Separately, current [THM-2454](../../../../01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md) gives the strong pure-blue atoms `{1,C3,T5}`. Neither classification is the signed-exception triangle-substitution family: `C3[C3]` is exceptional but has `h=39`, not one. Current repaired canon outranks THM-2450's historical frontier.

The incoming `G1–G5` script was absent from this sparse checkout during this audit. Therefore the claimed exhaustive order-nine and circulant order-fifteen censuses were not reproduced here and are not certified by this report. The independent finite controls below are smaller and explicitly scoped.

## 2. Definitions and faithful carrier

Let T be a finite tournament on V, with n at least two. Write `G=Aut(T)` and let `E=Aut±(T)` be the union of its automorphisms and anti-automorphisms, viewed as actual permutations of V. An anti-automorphism reverses every arc. A permutation cannot be both positive and negative when n is at least two, so sign is well defined. If anti-automorphisms exist, G is normal of index two in E; otherwise E=G.

A subset S is regular when its **setwise** stabilizer is trivial. This is a distinguishing binary coloring with two named colors, not a pointwise base, proper coloring, acyclic coloring, or an unlabeled partition allowing color interchange. The signed distinguishing number is the least number of named colors whose stabilizer in E is trivial.

For a faithful finite permutation action, primitive means transitive with no invariant partition other than singletons and the whole carrier. It does not mean that T has no nontrivial module. A module is a set on which every outside vertex has uniform orientation. Tournament indecomposability is a different condition.

An explicit indecomposable but not vertex-primitive five-vertex tournament has outneighborhoods

- 0: {2,4}
- 1: {0}
- 2: {1}
- 3: {0,1,2}
- 4: {1,2,3}.

It has no nontrivial module, ordinary automorphism group 1, and signed group of order two. Its degree classes are `{0}`, `{1,2}`, and `{3,4}`; an automorphism cannot swap either two-element class because it would reverse that class's arc, so the ordinary group is visibly trivial. The independent exhaustive module test verifies indecomposability. Thus one may not substitute modular primeness for the group-action hypothesis in Gluck's theorem. This is a distinction between hypotheses, not a counterexample to the modular-prime classification, which the companion structural proof now closes.

For n=0 or 1, use the faithful permutation image if discussing vertex colors: it is trivial. A formal sign extension that retains both `(id,+)` and `(id,−)` has an invisible kernel and cannot be distinguished by any vertex coloring. For n=2 the signed group is C2, primitive and binary distinguishable, while the ordinary group is not transitive. This is the small exception to the equivalence below.

## 3. Ordinary and signed primitivity are equivalent for n≥3

**Proposition.** For a finite tournament T of order n≥3, `G` is primitive on V if and only if `E` is primitive. In either case T is self-converse, n is odd, and in affine coordinates its signed point stabilizer is the ordinary point stabilizer times the central scalar involution `−I`.

**Proof.** A tournament has no nonidentity automorphism of order two: a transposed vertex pair would have its arc reversed. Thus G has odd order and is solvable by Feit–Thompson; E is solvable as an extension by a subgroup of C2.

If G is primitive, the usual structure theorem for primitive solvable groups gives a regular elementary-abelian normal subgroup V0, identified with the additive group of `F_p^d`. The tournament is a Cayley tournament for V0. Its connection set S satisfies `S∩(−S)=∅` and `S∪(−S)=V0\{0}`. Therefore the map `x↦−x` is anti. Since E contains the primitive group G, E is primitive as well.

Conversely suppose E primitive. Its normal subgroup G cannot be trivial: otherwise a transitive group of order at most two would have degree n≥3. The G-orbits are E-blocks, so G is transitive, and n divides its odd order. Write the primitive solvable affine group as `E=V0⋊M`, where `V0` is regular elementary abelian of odd prime-power order. If E=G the claim is immediate (and the previous paragraph would produce an anti-map, so this case does not actually occur). Otherwise `V0≤G` and M has even order. Choose an involution t in M. It is negative and fixes the affine origin. No anti-map can fix two distinct vertices, so t has no nonzero fixed vector. Since t²=I in odd characteristic, its eigenvalues are ±1 and its +1 eigenspace is zero. Hence `t=−I`. Thus t is central in M and

`M=(M∩G) × ⟨−I⟩`.

Any `M∩G`-invariant subspace is invariant under scalar −I, hence under M. Irreducibility of M therefore implies irreducibility of `M∩G`, which is the affine criterion for G to be primitive. ∎

The elementary-abelian affine structure and inversion mechanism are classical, not a novelty claim. For a primary-source occurrence of the primitive odd-group Cayley-tournament/inversion argument, see Cameron–Gadouleau, *Remoteness of permutation codes*, European Journal of Combinatorics 33 (2012), 1273–1285, proof of their odd-order transitive-group result: https://arxiv.org/abs/1110.2028 and https://doi.org/10.1016/j.ejc.2012.03.027 .

## 4. Exact primitive classification

**Theorem.** If T is finite of order at least three and its ordinary, equivalently signed, vertex permutation action is primitive, then

`D±(T)=3` exactly for `T ≅ C3, C5{1,2}, P7`,

and `D±(T)=2` otherwise.

**Proof.** Gluck's Theorem 1 classifies primitive solvable permutation groups with no regular subset. Its exceptional degree and point-stabilizer data are:

- degree 3: stabilizer order 2;
- degree 4: orders 3,6;
- degree 5: orders 2,4;
- degree 7: order 6;
- degree 8: order 21;
- degree 9: orders 16,24,48, or the order-eight dihedral group.

The actual signed tournament group has order either odd or twice odd. The only entries in this list having 2-part at most two are `S3` on 3 points, `D10` on 5 points, and `AGL(1,7)` on 7 points. Equivalently the surviving group orders are 6,10,42. The ordinary index-two subgroups have orders 3,5,21.

In degree three, invariance under the regular C3 action forces the directed triangle. In degree five, invariance under C5 gives a skew connection set in `Z5\{0}`; its four possibilities are all multiplier-isomorphic to `{1,2}`. In degree seven, affine group structure identifies the ordinary subgroup with translations and the order-three scalar subgroup. Its two nonzero scalar orbits are `{1,2,4}` and its negative, so the connection set is one of them. These give P7 and its isomorphic converse.

Conversely, C3 and P7 have signed groups of orders 6 and 42, larger than their largest constant-cardinality subset layers (3 and 35), so no regular subset can exist. For C5, complement to a set of size at most two. Every singleton is fixed by a reflection, and every pair `{u,v}` by `x↦u+v−x`; hence no binary coloring is regular. Three colors suffice in each example by separately coloring two vertices and giving the rest a third color (for C3 all vertices are separate). The universal three-color theorem was already proved in the prior packet. Finally, the primitive action is nontrivial, so every nonexception has distinguishing number exactly two, not one. ∎

The deep classification input is David Gluck, *Trivial Set-Stabilizers in Finite Permutation Groups*, Canadian Journal of Mathematics 35(1) (1983), 59–67, Theorem 1 on printed page 65; Corollary 1 on page 66 gives the odd-order regular-subset theorem. DOI: https://doi.org/10.4153/CJM-1983-005-2 . The user's uploaded paper text and the publisher's primary PDF were both checked. The primitive statement is a deduction from that classification, not a theorem stated about tournaments in the original paper.

## 5. Paley prime powers, without an automorphism-classification dependency

Define Pq on the finite field Fq, for a prime power `q≡3 mod4`, by `x→y` exactly when `y−x` is a nonzero square. Here `q=p^f`, with `p≡3 mod4` and f odd.

Every map `x↦ax+b`, `a≠0`, is positive when a is a square and negative when a is a nonsquare. Thus E contains the sharply two-transitive group `AGL(1,q)`, so E is primitive. The theorem above leaves exactly q=3 or 7. This proof needs neither the exact full automorphism group nor Carlitz's theorem; in particular it does not accidentally omit composite prime powers such as 27.

An alternative direct ordinary-primitivity proof uses translations and square multipliers. A nonzero additive invariant subspace containing w contains all `(q−1)/2` elements of `w·(Fq*)²`, plus zero. A proper Fp-subspace has size at most `q/p≤q/3`, too small. Thus the ordinary affine subgroup is primitive.

The notation P5 must not be used for the five-vertex exception under this tournament convention: q=5 gives a Paley graph, not a Paley tournament. The five-vertex tournament is `C5{1,2}`.

Literature attribution boundary: ordinary tournament two-color distinguishability is classical. Primary later accounts include Meslem–Sopena, https://arxiv.org/abs/1608.04866 (published DOI https://doi.org/10.1016/j.dam.2018.10.021), and Lozano, https://arxiv.org/abs/1707.05549 . Targeted searches for converse/anti-automorphism distinguishing classifications did not locate an earlier exact signed three-exception statement. This is not a comprehensive priority search and does not establish novelty; moreover the exact statement is already incoming THM-4743 in this repository.

## 6. Independent review evidence for the companion transitivity reduction

The following proof was independently derived here and passed to the companion structural-reduction report. It is included as review evidence, not as a separate theorem-publication claim.

**Theorem.** Let T have at least two vertices. If `G=Aut(T)` is intransitive on its vertex set, then `E=Aut±(T)` has a regular subset.

**Proof.** Gluck gives a G-regular subset S. If T is not self-converse, this is already E-regular. Otherwise E has twice odd order and contains an involution τ outside G. Because G is normal, τ permutes the G-orbits. Every G-orbit has odd cardinality. On a τ-invariant G-orbit, the involution τ has a fixed vertex; since an anti-map fixes at most one vertex in the entire tournament, at most one G-orbit can be τ-invariant. There are at least two G-orbits, so τ swaps some pair O,O′. Their equal cardinality m is odd.

Every anti-map has the same induced permutation of the G-orbits, since all anti-maps lie in Gτ. If `|S∩O|≠|S∩O′|`, no anti-map can stabilize S, so S is E-regular. If the occupancies agree at k, replace S with `S△O`. Because O is G-invariant, symmetric difference with O commutes with G and preserves the G-stabilizer exactly. The two occupancies become m−k and k, which differ because m is odd. Again every anti-map is excluded. ∎

Consequences:

- An exceptional tournament is ordinary vertex-transitive, so its order is odd.
- It is strong: the SCCs of a tournament have a canonical linear order, which every ordinary automorphism fixes; vertex-transitivity forces just one SCC.
- A canonical modular decomposition of an exceptional tournament has equal isomorphic maximal strong modules, since automorphisms act transitively on vertices and permute those modules. Unequal canonical-block constructions cannot produce missing exceptions.
- At this intermediate stage the remaining frontier was vertex-transitive modular-prime tournaments whose vertex action is imprimitive, and the palette behavior of homogeneous modular decompositions. This transitivity theorem alone does not close HYP-9353; the companion structural proof subsequently closes it using minimal blocks and the full ternary kernel.

The closest inherited mechanism is Gluck's regular subset. The added operation is symmetric difference with a G-orbit. It preserves ordinary rigidity exactly while changing a negative-channel occupancy. The cheapest hostile to the necessity of transitivity is order two, and the construction correctly gives a distinguishing singleton there.

## 7. Hamilton paths and metagraphs: the right carrier is essential

Let `mathcal H(T)` be the set of ordered directed Hamilton paths and `H(T)=|mathcal H(T)|`. Ordinary G acts freely, so

`h(T)=|mathcal H(T)/G|=H(T)/|G|`.

If τ is anti, its action on ordered paths must reverse their positions:

`(v1,...,vn) ↦ (τ(vn),...,τ(v1))`.

This induces an involution on the ordinary path-orbit set. If it fixes fH orbits, then the signed path-orbit count is `(h+fH)/2`. This is a different carrier from the regular-subset orbit palette. Exceptionality says the latter involution fixes every point; it does not say the path involution does.

Independent exact controls, with order listed as `(H, |G|, h, fH, signed path orbits)`:

- C3: `(3,3,1,1,1)`.
- C5{1,2}: `(15,5,3,3,3)`.
- P7: `(189,21,9,3,6)`.
- C3[C3]: `(3159,81,39,11,25)`.
- Nonexceptional C7{1,2,3}: `(175,7,25,9,17)`.

In particular P7 is exceptional on binary vertex subsets but has six of nine Hamilton-path orbits moved by converse. The triangle and its two-level tower have the same two ordinary regular-subset orbits but one and thirty-nine Hamilton-path orbits. No equality or preservation of metagraph adjacency follows from their equal palette sidecar.

For an individually marked ordered full Hamilton path, the ordinary residual stabilizer is trivial. The signed residual stabilizer has order at most two, and any surviving element sends the ith vertex to the `(n+1−i)`th. A binary color on the whole vertex carrier can always kill that residual involution by marking one moved vertex. Thus global signed exceptionality is not an obstruction after a full path is already individually marked.

If colors are restricted instead to omitted vertices of a partial path, retain the prior packet's kernel warning. The signed residual action is faithful for at least two omitted vertices. With zero or one omitted vertex an anti-map can be invisible; then no coloring of that omitted carrier kills the kernel. Binary failure can still occur for a faithful omitted action: put a directed triangle U between a source a and sink b, and mark the directed path `(a,b)`. The residual action on U is S3, so two omitted-vertex colors fail. These statements concern the actual residual subgroup, not the full induced subtournament group by default.

## 8. Independent finite receipts

`verify_small_actions.py` imports no producer modules. It enumerates all labeled tournaments of orders 2 through 5 (2+8+64+1024=1098), enumerates their automorphisms and anti-automorphisms directly, tests all vertex subsets and all potential permutation blocks, and checks:

- ordinary binary regular subsets always exist;
- ordinary/signed primitivity equivalence for n≥3;
- primitive tournaments are self-converse;
- the stated primitive exceptions in this range;
- intransitive ordinary automorphism groups, and in particular non-strong tournaments, are not exceptional.

It separately examines C3, C5{1,2}, P7, C7{1,2,3}, C3[C3], and the explicit modular-prime but action-intransitive example. Hamilton-path counts are checked both by permutation enumeration and independent subset dynamic programming; path-orbit counts use explicit group orbits with reversed negative transport.

The degree-nine tower is an important missing-hypothesis control: its signed group has order 162 and zero regular binary subsets, but is imprimitive and decomposable. It does not contradict the three-exception theorem.

Reproduce from the workspace root:

    python3 session24-signed/verify_small_actions.py --output session24-signed/receipt.json
    python3 -O session24-signed/verify_small_actions.py --output session24-signed/receipt-O.json
    cmp session24-signed/receipt.json session24-signed/receipt-O.json

Normal and optimized receipts match byte-for-byte. Explicit exceptions, not removable assertions, enforce the checks. The all-order statements are proved above; the finite checks are controls, not substitutes for the classification or a proof of HYP-9353.

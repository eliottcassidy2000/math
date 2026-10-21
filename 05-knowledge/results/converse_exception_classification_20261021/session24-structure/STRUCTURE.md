# Binary signed-distinguishing tournaments: structural classification

Session 24, 2026-10-21. Status: **PROVED relative to the cited Gluck and Feit–Thompson theorems; independently audited, with FINITE-EXACT controls**.
The elementary structural steps below are new derivations in this session,
not assertions of literature priority. The primitive seed classification is
inherited from THM-4743 and independently verified in the companion signed
report. The proof and both-mode finite replay passed the companion review.

## 1. Scope, inheritance, and notation

All tournaments are finite, with at least two vertices unless stated otherwise.
Write G(T)=Aut(T), E(T)=Aut±(T), and call T exceptional if no vertex subset has
trivial setwise stabilizer in E(T). Let R(T) be the set of G(T)-orbits of
G(T)-regular subsets; q(T)=|R(T)|. For self-converse T, any anti-automorphism
induces the same involution alpha on R(T); f(T) is its number of fixed points.
Then T is exceptional exactly when f(T)=q(T).

Inherited dependencies and boundaries:

1. **CITED:** Gluck, *Trivial Set-Stabilizers in Finite Permutation Groups*,
   Canadian Journal of Mathematics 35 (1983), 59–67, Corollary 1: every odd-order
   finite permutation group has a regular subset. DOI
   https://doi.org/10.4153/CJM-1983-005-2 . Primary PDF:
   https://www.cambridge.org/core/services/aop-cambridge-core/content/view/9667B85D9DDD0FE911C8382287B86167/S0008414X00022550a.pdf/div-class-title-trivial-set-stabilizers-in-finite-permutation-groups-div.pdf
2. **INHERITED, INDEPENDENTLY VERIFIED:** THM-4743,
   `01-canon/theorems/THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md`,
   primitive signed-action exceptions exactly C3, R5=Z5{1,2}, P7, with
   respective q=2,6,4. This uses Gluck's primitive solvable classification and
   Feit–Thompson. Its regular-palette and tower results are inherited, not new.
3. **INHERITED:** the universal signed distinguishing bound 3, published in
   `05-knowledge/results/signed_tournament_palette_rigidity_20261021/`.
   A self-contained proof is recalled below.
4. **OPEN AT INHERITANCE:** HYP-9353 claimed that the entire exceptional family
   is X[C3^[k]], X in {C3,R5,P7}, k>=0. The work here targets its outstanding
   unequal-module and nonmodule-imprimitivity cases. No bounded census is used
   as an all-order premise.

Prior-repository search covered the exact THM/HYP files, their source synthesis,
the signed-palette publication packet, and searches for intransitivity,
vertex-transitivity, regular subset orbits, and palette recursions. The inherited
HYP explicitly leaves unequal blocks, non-strong composites, and nonmodule
imprimitivity open. No claim of external novelty is made.

Live concepts: vertex-orbit parity; regular-subset palettes; majority quotients;
primitive block actions; ternary marker kernels; substitution recurrences.
Hostile controls: strong nontransitive C3[C3,1,C3]; proper ternary base kernels;
C3[R5] versus R5[C3]; twisted palettes with nontrivial converse action.

## 2. Elementary preliminaries

The group G(T) has odd order: an involution would exchange some pair of vertices,
reversing their arc. If T is self-converse, E(T) has twice odd order, G(T) is its
normal subgroup of index 2, and all anti-automorphisms form one coset.

Every anti-automorphism has at most one fixed vertex, because two fixed vertices
would leave their arc unchanged. Every anti-involution therefore has exactly
one fixed vertex when the vertex set has odd size.

Complementation is fixed-point-free on R(T). If gS=S^c for an odd-order g, then
iterating an odd number of times gives S=S^c, impossible. In particular q(T) is
positive and even. An anti-automorphism acts on R(T) because it normalizes G(T);
its square is in G(T), so the induced action is an involution. The union of its
moved ordinary orbits is exactly the set of E(T)-regular subsets, and each pair of moved ordinary
orbits gives one regular E(T)-orbit.

**Universal bound 3.** Choose a G(T)-regular subset S by Gluck. Its E(T)
stabilizer is trivial or has order 2: two distinct anti-elements would have a
nonidentity ordinary product stabilizing S. In the latter case its unique
nonidentity element is an anti-involution tau. Choose a vertex moved by tau and
split that vertex off from its one of the two color classes. The resulting
three-coloring has trivial E(T)-stabilizer.

## 3. Every exception is ordinarily vertex-transitive

**Proposition 1.** If G(T) is intransitive, T has an E(T)-regular subset.

Proof. If T is not self-converse, use Gluck directly. Otherwise let tau be an
anti-involution, which exists since |E(T)| is twice odd. Each G(T)-orbit is odd,
and tau permutes these orbits. It can preserve at most one individually:
otherwise its restrictions to two odd orbits would each have a fixed vertex,
contradicting the preceding fixed-vertex observation. Thus when G(T) is
intransitive, some pair O,O' of its orbits is exchanged by tau. Every anti-map
induces this same exchange, since all differ from tau by an element of G(T).

Choose a G(T)-regular subset S. If |S intersect O| differs from |S intersect O'|,
no anti-map stabilizes S. If both occupancies equal k, replace S by S symmetric
difference O. This preserves its G(T)-stabilizer exactly because O is G(T)-
invariant. The two occupancies become |O|-k and k, unequal since |O| is odd.
Again no anti-map stabilizes the resulting subset. QED.

Consequences: every exceptional T is vertex-transitive, has odd order, and is
strongly connected. The latter follows because automorphisms preserve each
component in the linearly ordered strong-component decomposition.

**Palette lower bound.** If G(T) has k vertex orbits, then q(T)>=2^k. On each
orbit O, its faithful induced odd-order group H_O has a regular subset and its
complement in two different H_O-orbits. Any of the 2^k choices of local orbit
types gives a globally G(T)-regular subset, and distinct choices cannot be
G(T)-equivalent. More precisely, if r_O is the number of regular H_O-orbits,
then q(T)>= product_O(r_O |H_O|)/|G(T)|. Hence q(T)=2 itself implies ordinary
vertex-transitivity. If T is also self-converse, the two complementary regular
orbits have different cardinalities (|T| odd), so alpha fixes both and T is
exceptional.

## 4. Majority quotient for arbitrary group blocks

This step does not assume that permutation blocks are tournament modules.

Let T be exceptional, G=G(T), E=E(T), and let B_1,...,B_m be any nontrivial
E-invariant system of blocks. Because G is vertex-transitive, it is transitive
on the blocks, each block has the same odd size b, and m is odd and >=3. The
induced tournaments on all blocks are ordinarily isomorphic; denote one by B.
They are self-converse, since an anti-map takes any block to another block by
an anti-isomorphism, followed by an ordinary identification back to B.

Define A, the **majority quotient tournament**, on {1,...,m}: orient i→j exactly
when more than half of the b^2 arcs between B_i and B_j point from B_i to B_j.
There is no tie because b is odd. Ordinary elements of E induce automorphisms
of A and anti-elements induce anti-automorphisms of A. In particular the kernel
K of E's action on the blocks has no anti-element: an anti-element fixing all
blocks would induce the identity anti-map on the nontrivial tournament A.
Equivalently, reversing the odd number b^2 of arcs between two fixed blocks is
impossible. Therefore Q=G/K has odd order and index 2 in Ebar=E/K.

Choose ordinary identifications B_i≅B. Each element of E acts between blocks by
ordinary maps or anti-maps of B, all of the same sign as its permutation in
Ebar. Thus E embeds in the full sign-matched wreath permutation group

  W = { (pi; beta_1,...,beta_m): pi in Aut±(A),
        beta_i in Aut±(B), all beta_i of the sign of pi }.

W acts by automorphisms and anti-automorphisms of the lexicographic product A[B].
We do NOT assume W is the full signed automorphism group of A[B], and do NOT
infer that T=A[B]. Only the containment E<=W as permutation groups is used.

**Proposition 2 (palette forcing).** In this situation B is exceptional,
q(B)=2, and A is exceptional.

Proof. First suppose B is nonexceptional. It is self-converse, so its converse
involution moves some pair x,alpha(x) of ordinary regular-subset orbits. Color
the m blocks with these two types according to an ordinary distinguishing
binary coloring of A (Gluck). Choose a representative local regular subset of
each specified type. An ordinary stabilizer must induce the identity on A and
then be identity in every block. An anti-stabilizer would exchange the two
types. That would force equal counts of the two types, impossible for odd m.
Hence this subset is W-regular and therefore E-regular, a contradiction.

Now B is exceptional, so alpha fixes every one of its q(B) regular types.
If q(B)>=3, choose a signed distinguishing three-coloring of A and realize its
colors by three distinct regular types of B. The same blockwise argument gives
a W-regular subset, a contradiction. Because q(B) is positive and even,
q(B)=2. Finally if A were nonexceptional, its signed distinguishing binary
coloring, realized by these two fixed types of B, would again be W-regular.
Thus A is exceptional. QED.

This also proves the earlier canonical-module reduction: in any decomposition
T=A[B] with the copies a nontrivial signed-invariant block system, exception forces both
factors exceptional and q(B)=2. For canonical maximal strong modules,
vertex-transitivity supplies equal isomorphic factors, so unequal canonical
blocks and non-strong composites are ruled out.

## 5. Minimal blocks are triangles

Assume E is imprimitive, and choose a non-singleton E-block of minimum possible
size, together with its translates. Its setwise stabilizer induces a primitive
permutation group on the block: a nontrivial smaller block for that induced
action would also be an E-block, contrary to minimality. This induced group is
a subgroup of Aut±(B), so Aut±(B) is primitive too (a block system invariant
under the larger group would be invariant under the primitive subgroup).

By Proposition 2 B is exceptional with q(B)=2. The inherited primitive
classification leaves B=C3,R5,P7; their palette counts 2,6,4 force B=C3.
Thus every minimal non-singleton block has three vertices inducing a directed
triangle. The block count m is odd and >=3.

## 6. The ternary kernel is full, forcing genuine modules

Fix cyclic coordinates Z3 on each triangle block, so its ordinary automorphisms
are translations and its anti-automorphisms are x↦−x+a. The block kernel K is
therefore a vector subspace of V=F3^m, acting by coordinate translations.
Conjugation by an element of G permutes the coordinates of K; the translations
in that element cancel in conjugation. Since G is transitive on blocks, this
coordinate-permutation action is transitive.

Choose a binary block coloring c regular for Q=G/K, using Gluck on Q's faithful
action on the blocks. For EVERY marker vector x=(x_1,...,x_m) in F3^m form the
vertex subset S_x as follows:

- in a block colored 0, take the singleton with vertex x_i;
- in a block colored 1, take the doubleton whose OMITTED vertex is x_i.

Every S_x is G-regular. Indeed an ordinary stabilizer must preserve singleton
versus doubleton cardinalities, hence must fix c, hence lies in K; no nonzero
coordinate translation fixes the singleton or omitted vertex in that block.

Since T is exceptional, each S_x is fixed by an anti-element. Its block
permutation must preserve c. There is at most one anti-element tau in Ebar
preserving c: the quotient of two would lie in Q and preserve c. There is at
least one, by taking x=0. Therefore the SAME tau is the only possible block
permutation of an anti-stabilizer for EVERY x. Moreover tau^2 belongs to Q and
preserves c, so tau^2=1. It is a nonidentity anti-element of the majority
tournament, and as an involution of odd m blocks it has a fixed block i.

Choose a lift h in E of tau. On marker vectors h acts by

  h(x) = -P_tau x + a

for some fixed a in V. All lifts of tau are exactly k h, k in K. Thus the
existence of an anti-stabilizer for every S_x states exactly

  (I+P_tau)x-a belongs to K for every x in V.

Taking x=0 gives a in K. Subtracting gives (I+P_tau)V contained in K. For the
fixed block i, (I+P_tau)e_i=2e_i, so e_i belongs to K. Conjugating by G gives
every coordinate vector e_j in K. Therefore K=F3^m.

Each individual triangle rotation consequently is an ACTUAL automorphism of
T fixing all vertices outside that triangle. Its transitivity on the three
vertices forces every outside vertex to have uniform orientation to the whole
triangle. Thus every triangle block is a tournament module, and

  T ≅ A[C3].

This is the precise missing implication from arbitrary permutation blocks to
modules. Majority arcs alone lose all minority-arc data. The full-kernel
argument, not the majority quotient itself, proves that no minority arcs can
remain in an exceptional tournament.

## 7. Classification and exact scope

**Theorem (PROVED relative to the cited classical group theorems and the
independently verified primitive seed classification).** A finite tournament T on at least two vertices has no
regular subset for Aut±(T) if and only if

  T ≅ X[C3^[k]],  X in {C3,R5,P7},  k>=0.

Proof of necessity. Proposition 1 makes T ordinary-vertex-transitive. If E(T)
is primitive, apply the primitive seed classification. Otherwise Sections 4–6
give T=A[C3] with A exceptional and fewer vertices. Induct on |T|. The
lexicographic product is associative, so the resulting factors assemble into
the stated family.

Sufficiency is inherited from THM-4743's tower substitution theorem. It can
also be seen directly: C3 has exactly two ordinary regular subset types,
singletons and doubletons, both fixed by converse. In X[C3^[k]], nonregular
local subsets admit a local nonidentity automorphism; when all local subsets
are regular, their two types reduce the stabilizer question to the previous
exceptional factor. Canonical triangle/tower modules ensure the required full
wreath actions. QED.

This closes the old nonmodule-imprimitivity route after independent audit;
it is not merely a reduction to prime vertex-transitive exceptions. In
particular no prime tournament with imprimitive signed action is exceptional.
The word "prime" here means indecomposable as a tournament, NOT primitive as
a permutation action.

## 8. General twisted-palette recursion under C3 substitution

Let B be a strong self-converse tournament, q=q(B), f=f(B), and assume the three
copies of B in C3[B] are its canonical maximal strong modules. This holds for
strong B and the prime quotient C3. Set q'=q(C3[B]), f'=f(C3[B]). Then

  q'=(q^3-q)/3,          f'=f(q-1).

Proof. Ordinary regular subsets correspond to nonconstant triples of regular
B-palette types, modulo cyclic rotation. Nonconstant triples have free C3
action, giving q'. Converse acts by a reflection of the three positions
combined with alpha in every coordinate. A specified twisted reflection fixes
f q triples: f choices at its fixed position, and q choices for one of the
exchanged coordinates, the other being its alpha-image. Exactly f of these
are constant, leaving f(q-1). A nonconstant triple is fixed by at most one of
the three twisted reflections (two would yield a nontrivial rotation).
Counting all three reflections and dividing by the three representatives in
each ordinary orbit gives f'=f(q-1).

Consequently the number of regular signed subset orbits of C3[B] is

  ((q^3-q)/3 - f(q-1))/2,

and the actual number of signed-regular vertex subsets is

  3 |Aut(B)|^3 ((q^3-q)/3 - f(q-1)).

If B is exceptional, f=q, so the number of signed-regular orbits simplifies to
binom(q,3). Thus C3[B] is exceptional exactly when q=2. In particular:

- B=C3: (q,f)=(2,2), (q',f')=(2,2), no signed-regular subset;
- B=R5: (6,6)→(70,30), 20 signed-regular orbits and 15,000 actual subsets;
- B=P7: (4,4)→(20,12), 4 signed-regular orbits and 222,264 actual subsets.

The tower fixed point (2,2) and the 15,000 count were inherited numerical
facts; the two-variable recursion supplies their common mechanism and retains
the converse-fixed palette coordinate f that q alone forgets.

## 9. Sparse family generating function

The proved-family order support contains exactly one class at each of
3^(k+1), 5*3^k, 7*3^k, k>=0, since these three sets of orders are disjoint.
Its ordinary generating function is

  F(z)=sum_{k>=0}(z^(3*3^k)+z^(5*3^k)+z^(7*3^k)),
  F(z)=z^3+z^5+z^7+F(z^3).

It is not D-finite over characteristic zero. Its coefficient sequence has
arbitrarily long consecutive runs of zeros followed by a nonzero coefficient.
A P-recursive sequence cannot do this: beyond all integer zeros of the leading
polynomial coefficient of a fixed finite-order recurrence, a run of as many
zeros as the recurrence order propagates forever. For order zero the sequence
would be eventually zero immediately. D-finite formal power series have
P-recursive coefficients, so F is not D-finite.

This statement about the known family does not require the new classification.
The independently audited classification now identifies F with the generating
function of ALL exceptions. No natural-boundary theorem is asserted.

## 10. Reproducibility and loss ledger

Run `python3 verify_structure.py --output receipt.json` and
`python3 -O verify_structure.py --output receipt-O.json`. Assertions are not
used for validity gates, so optimization does not disable checks. The finite
universe and controls are printed in the receipts. These checks supplement,
not replace, the all-order proofs.

Typed bridges:

- T→G-vertex-orbits preserves ordinary stabilizer constraints; forgetting
  within-orbit data is repaired by a G-regular subset and orbit toggling.
- Arbitrary E-blocks→majority tournament preserves global signed action but
  loses minority arcs. The full ternary kernel restores module status.
- Regular block subsets→palette colors preserves stabilizer existence under
  sign-matched wreath actions; actual orbit representatives and group-order
  factors restore subset counts.
- C3 substitution→(q,f) recursion preserves ordinary and converse-fixed regular
  orbit counts; q alone loses the twisted-reflection obstruction.

No Hamiltonian-path-count classification, metagraph equivalence, or LRC result
follows merely from this signed-distinguishing classification.

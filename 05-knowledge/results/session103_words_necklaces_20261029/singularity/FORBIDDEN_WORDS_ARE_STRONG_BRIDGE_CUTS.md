# The actual forbidden words have an intrinsic strong-bridge meaning

**PROVED, independently audited with the explicit block-size correction.** This note proves a structural equivalence for the recovered insertion grammar. It does **not** prove that the native two-run observer is singular exactly on this structure. The latter remains the spectral obligation in this note. The independent audit is `../audit/STRONG_BRIDGE_PROOF_AUDIT.md`; its verdict explicitly requires nontrivial blocks, because singleton strong cuts already occur at the seed.

Pinned definition source: commit `d3132eafbfa9c43c6d9f1110566b1e9c3a3286e1`, `05-knowledge/results/operation_closed_synthesis_20261029/THEORY.md` and `tournament_probe.py`. Earlier work at `75bd8b27b6d` did not yet contain these definitions; that snapshot is not overwritten.

## 1. Grammar and structural predicate

Start with the strong four-tournament, with unique cycle 0→1→2→3→0 and diagonals 0→2,1→3. A new vertex x loses only to its selected parent p and beats every other old vertex. Insert x immediately after p in the marked unique cycle.

A type-A state has one maximum-outdegree vertex a (indegree one); legal parents are pred(a), labelled a, and a, labelled b. A type-B state has two such vertices a0→a1, adjacent in the marked cycle; legal parents are pred(a0), a0, a1, labelled c,d,e respectively. The seed is type B. These are exactly the recovered operations. Statements here are relative to the generated family and do not require reproving the cited completeness classification of all unique-cycle tournaments.

Call (R,Q) a **strong bridge cut** if R,Q partition the vertex set, both blocks have size at least three and their induced tournaments are strongly connected, and there is exactly one arc directed from R into Q, say q→r. All other cross arcs point Q→R. The size condition is explicit: singleton tournaments are strongly connected by a common convention and must be excluded.

**Claim.** A generated tournament has a strong bridge cut if and only if its generating word contains `be` or `ee`.

This is an intrinsic graph predicate. The proof does not identify cyclic rotations of construction histories.

## 2. Elementary facts used

**F1: the newest maximum.** After every insertion, the new vertex is the unique indegree-one vertex in type A, or is the second one in the type-B pair. Indeed every old vertex except the selected parent gains one incoming edge, whereas the new vertex has indegree one. If the parent had indegree one, precisely it and the new vertex remain indegree one, with parent→new. Otherwise only the new vertex has indegree one. This also recovers the stated child types.

**F2: cycle order across a bridge cut.** In any Hamiltonian cycle of a tournament with only the cross arc q→r from R to Q, all vertices of R occur consecutively, ending at q, and all vertices of Q occur consecutively, starting at r. The cycle crosses from R to Q exactly once; any second such crossing would repeat q→r. Hence any strong bridge cut in this family is visible as two complementary intervals of the marked unique cycle.

**F3: SCC order.** The strong components of a tournament are totally ordered C1→C2→…→Ct. Any spanning directed path visits them in this order. A singleton vertex in Cj has at least the total size of its preceding components as indegree inside the tournament; a vertex in a nontrivial strong component has at least one additional incoming edge internally. Nontrivial strong components have size at least three.

## 3. Every legal insertion preserves an existing strong bridge cut

Suppose (R,Q) is such a cut in T. Every vertex in R has at least |Q|−1≥2 incoming edges from Q. Therefore every maximum-outdegree (indegree-one) vertex lies in Q.

Its cycle predecessor also lies in Q. The only Q vertex whose predecessor in the marked cycle lies in R is the entry vertex r (F2). But r has the incoming bridge q→r and at least one incoming edge from inside strongly connected Q, so r has indegree at least two and is not a maximum.

Every legal parent is a maximum or the predecessor of a maximum. Thus p∈Q. Adjoining x to Q preserves strong connectivity: old Q is strong; every old vertex reaches p and then x; x reaches every old vertex other than p directly, and reaches p through any other old vertex and a path in Q. All new edges to R point x→R. Consequently (R,Q∪{x}) remains a strong bridge cut. This establishes structural persistence, not pointwise persistence of individual observer-kernel vectors.

## 4. A new cut can first appear only at `be` or `ee`

Let T+x have a strong bridge cut (R,Q), but T have none. Write p for the parent of x.

First x cannot lie in R: any vertex of R receives at least |Q|−1≥2 edges from Q, while x has indegree one. Hence x∈Q. The parent p also lies in Q, because x needs an incoming edge within strongly connected Q, and p is its only incoming neighbour. Put Q0=Q−{x}.

If Q0 were strongly connected, (R,Q0) would already be a strong bridge cut of T, contrary to hypothesis. Therefore Q0 is nonstrong. Since Q is strong, |Q0|≥2. Let its SCCs be C1→…→Ct, t≥2.

**The parent p lies in Ct.** Otherwise Ct has no route to p in Q0, hence no route to x, because the only edge into x is p→x. It also has no route to an earlier component without using x. This contradicts strong connectivity of Q.

Assume for contradiction |Q0|≥3. Every old maximum lies in Q0: a vertex in R receives at least |Q0|−1≥2 incoming edges from Q0. No vertex of Ct can have indegree one in T. If Ct is a singleton, at least two vertices precede it unless all of Q0 has size two. If Ct is nontrivial, it has an internal incoming edge and at least one vertex precedes it. Thus p is not a maximum.

Nor can p be the marked-cycle predecessor of a maximum m. An indegree-one vertex of Q0 can lie only in C1, or in a singleton C2 preceded by a singleton C1. In the first case, its cycle predecessor lies in C1 or R, not Ct. In the second case its predecessor is the vertex of C1, again not Ct. Here the marked cycle visits R, then C1,…,Ct consecutively by F2–F3, which apply to the same cross partition of T. This exhausts the possible legal parents and gives a contradiction.

Therefore |Q0|=2. Write its arc u→p. The bridge from R must enter u, not p, because otherwise T could not have a Hamiltonian cycle visiting both u and p (F2). Accordingly u and p each have indegree one in T: u receives only the bridge, and p receives only u→p. They are the type-B maximum pair. Choosing parent p is exactly operation e.

Since R is strong, T has at least five vertices, so it has a preceding insertion letter. Because T is type B, that letter is b,d,or e. It cannot be d. In a d insertion the parent u was the first old maximum, and the other old maximum v had u as its only incoming neighbour. Removing u and the newly inserted p therefore leaves v as a source in R, contradicting strong connectivity of R. Hence the preceding letter was b or e: the new cut first appears at `be` or `ee`.

This argument proves the direction “cut implies forbidden factor” by induction: either the cut already existed, or the last two letters are one of the forbidden pairs. The seed has no cut because it has only four vertices.

## 5. Each forbidden pair creates a cut

Suppose the next two operations are b,e or e,e. Let u be the selected maximum in the first of these operations, y its newly inserted child, and x the child inserted by the final e. Their arcs are

    u→y→x→u,

so Q={u,y,x} is a directed triangle. The old maximum u has just one incoming neighbour q in the older tournament. Both y and x beat every vertex outside Q. Hence all cross arcs point Q→R except q→u, where R is the older tournament with u deleted.

R is strongly connected. If the older tournament has at least five vertices, F1 says that the selected maximum u (the unique maximum for b, or the second maximum for e) is its newest vertex; deleting it recovers the preceding generated strong tournament. At the four-vertex boundary only e is possible; deleting the second seed maximum leaves a directed triangle, directly checked from the seed arcs. Thus (R,Q) is a strong bridge cut. Section 3 preserves it under every later legal letter.

Together Sections 3–5 establish the claimed all-order structural equivalence, with both blocks explicitly required to have size at least three.

## 6. What this closes—and what it does not

This replaces the formal forbidden-factor language by a graph-native predicate. The all-order observer conjecture can now be stated without construction histories:

> In this generated unique-Hamilton-cycle family, M2 is singular if and only if the tournament admits a strong bridge cut.

The combinatorial equivalence alone proves neither spectral implication. A bridge cut factors the Hamilton-cycle count as

    (#Hamilton paths of R ending at q) × (#Hamilton paths of Q starting at r).

Thus in the present family both rooted counts equal one. This sidecar matters: a single-bridge composition of two arbitrary strong tournaments need not have singular sign sector. An S4/S4 bridge choice with four Hamiltonian cycles supplies a hostile.

Individual kernels are not frozen by the legal grammar: the exact minus-sector nullities of ee,eec,eed,eee are 4,2,4,3. Consequently a proof of spectral persistence must transport some certificate or adapt its support, not assume every old kernel vector survives.

The earlier repair by a relay vertex is outside the legal grammar (the relay has indegree three). It also destroys strongness of the high block, replacing it by a strong block followed by a sink. Thus it is consistent with the intrinsic strong-bridge criterion and remains a useful guard against unrestricted insertion arguments.

## 7. Exact audit

    python3 -O session103_words_necklaces/singularity/check_actual_bridge_grammar.py

The checker reconstructs every pinned source word through order 12 with chronological names and an explicit cycle, verifies its rotation-canonical adjacency code against the pinned source, enumerates all strong bridge cuts, and compares their existence with forbidden-factor occurrence. For each discovered cut it checks every legal next insertion preserves the cut, including the order-12 to order-13 boundary. It also computes the four kernel nonpersistence controls from scratch. Exact counts are in `actual_bridge_grammar_receipt.json`. The source observer data are inherited finite-exact evidence; the combinatorial proof is not inferred merely from agreement with them.

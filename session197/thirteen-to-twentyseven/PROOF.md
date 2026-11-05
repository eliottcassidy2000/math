# No tournament with 13 Hamilton paths embeds induced in one with 27

Session 197. **All-order induced-extension obstruction.** This is not a nonexistence claim for either scalar count.

**Theorem.** There are no tournaments T,U with H(T)=13 and H(U)=27 such that T is an induced subtournament of U.

Consequently no append-only graph construction can implement the branch 13->27, even once, regardless of how many ports it uses, how they change, how many vertices it adds, or which graph representing the scalar 13 it chooses first. Rebuilding that branch must cease to preserve the old tournament as an induced graph. This does not obstruct unrelated target-dependent constructions of 27.

## 1. Reduction to strong components

Strong components of a tournament form a transitive chain, and Hamilton-path counts multiply along that chain. Every nontrivial strong tournament has at least three Hamilton paths: each of its at least three vertices can occur as a Hamilton endpoint. Every tournament has at least one Hamilton path.

Since 13 is prime, T has a strong component C with H(C)=13; every other component has count one and is a singleton. If T embeds induced in U, the vertices of the strong C lie in one strong component D of U. Hamilton-path counts are monotone under induced extension: insert each new vertex at a legal gap of every old Hamilton path, using a fixed rule; deletion of the new vertex recovers the old path. Thus H(D)>=13. But H(D) divides H(U)=27. Of the positive divisors 1,3,9,27, only27 is at least13, so H(D)=27.

It therefore suffices to classify all nontrivial strong tournaments of count at most27 and test whether an H27 core contains an H13 core. Section2 supplies a complete all-order certificate: the cap27 mixed-extension frontier is empty at order eight, and a strong-deletion theorem rules out every later order. No external strong-order inequality or graph-order cutoff is assumed. Arbitrary singleton padding on either side is harmless because the strong source core still lies in one target strong component.

## 2. A complete regenerable strong-core census

Every strong tournament of order at least four has a strong vertex deletion. Here is the elementary proof used by the inherited strong-core generator. Choose a proper strong vertex subset S of maximum size; a directed triangle exists. If its size is n-1, stop. Otherwise every outside vertex is uniform toward S, since a mixed vertex would make a larger proper strong subset. Strongness supplies a vertex l with S->l and a vertex u with u->S, with l->u. Then S+l+u is strong. If more than two outside vertices remain, this contradicts maximality. If precisely two remain, deleting any vertex of S gives a strong proper subset of order n-1, because every remaining S vertex lies on the directed triangle through l,u. This again contradicts maximality.

Therefore all nontrivial strong tournaments arise from the triangle by successive mixed one-vertex extensions. The cap H<=27 is safe throughout each deletion chain by induced monotonicity. The verifier enumerates every such mixed orientation and deduplicates exactly: order vertices by outdegree, and minimize adjacency rows over every permutation within equal-degree classes. Every isomorphism preserves these classes, so no graph can be lost by this canonicalization.

The complete cap27 inventories are:

- order3: one initial class;
- order4: 6 extension candidates,6 accepted,1 class;
- order5: 14 candidates,14 accepted,6 classes;
- order6: 180 candidates,48 accepted,14 classes;
- order7: 868 candidates,15 accepted,4 classes.

- order8: 504 candidates,0 accepted,0 classes.

This final empty frontier closes the induction. A hypothetical strong graph of larger order and H<=27 has a chain of strong vertex deletions down to order eight, with every deletion still under the cap by monotonicity. But the complete order-eight frontier is empty. Thus the preceding rows exhaust the nontrivial strong cap27 universe at every order. This is not a claim of stabilization from several empty-looking finite windows: the strong-deletion theorem supplies the explicit predecessor at each larger order.

The source regenerates these banks from scratch, without reading an inherited graph bank. A separate complete labelled enumeration through five reproduces the small banks. Every retained graph is emitted in the receipt with its exact Hamilton count.

## 3. The exact core table and containment check

There is one strong H13 core, of order five:

    [2,12,9,17,7].

The complete strong H27 list is:

    [2,32,11,19,7,29],
    [16,5,9,3,46,15],
    [2,4,17,7,43,79,31],
    [2,4,33,7,15,91,31],
    [2,16,3,7,45,79,31].

The first two have order six; the last three have order seven. The verifier counts every induced five-vertex subtournament of all five graphs. None has Hamilton count13. The two order-six subset-count multisets are both

    [3,3,9,9,9,15].

Each order-seven multiset contains only values from {1,3,5,9}; the complete multisets are in the receipt. Thus none contains the H13 core, completing the all-order proof.

## 4. Scope and positive controls

Both counts are physically realized by the displayed rows. The obstruction is specifically to induced containment. It does not claim that27 is absent, that no binary operation works at other counts, or that 13 cannot be extended to different larger counts.

For example the H13 core has a literal transitive-clone family 13,25,49,... at vertex0. Its first cloned graph has rows [34,12,41,49,39,2] and count25, and contains the original core induced. The obstruction distinguishes the opposite binary branch 13->27.

A further hostile control is important: the analogous statement for15->29 is false. The target rows [2,4,73,83,103,79,3] have H29, while the induced subset [0,1,2,3,4,6] has rows [2,4,41,51,39,3] and H15. Thus a forbidden all-order literal-clone return at15 can still have a one-time induced extension to29. The theorem here does not transfer clone obstructions indiscriminately to every graph operation.

In particular changing to another H13 representative before an append-only step cannot help: every such representative has the same strong H13 core, with only singleton components outside it. A successful 13->27 rebuilding operation must change or remove enough old structure that the original H13 tournament is not retained induced.

## 5. Reproduction

Run `python verify.py` and `python -O verify.py`; both must reproduce the frozen receipt. The verifier is standalone standard-library Python. Its endpoint-subset dynamic program counts exact Hamilton paths, while exact degree-class canonicalization and complete mixed extension generation prove census coverage. All failures raise explicit exceptions.

The finite search is sufficient because Section1 reduces induced containment to strong components and Section2 proves complete generation under the cap together with an empty frontier that excludes every larger strong order. There is no inference from an arbitrary bounded append search. The exploratory two-port search that suggested this question is not a premise of the proof.

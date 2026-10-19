# An affine cleanup from square sumsets to induced complete bipartite graphs

2026-10-19. **PROVED elementary bridge**, with a **CITED** unbounded sumset consequence and a separate **FINITE-EXACT** K_3,3 height threshold. No external-priority claim.

## 1. Exact cleanup lemma

Let A and B be finite sets of positive integers such that every element of A+B is an integer square. Define

    A' = {25a+1 : a in A},
    B' = {25b-1 : b in B}.

Then A' union B' induces K_(|A|,|B|) in the positive square-sum graph.

**Proof.** Both maps are injective and their values are positive. The two images are disjoint because their residues modulo five are respectively 1 and 4. For each cross pair,

    (25a+1)+(25b-1)=25(a+b),

which is a square. A sum within A' is 2 modulo five, and a sum within B' is 3 modulo five. Both are quadratic nonresidues modulo five, so there are no within-part edges. These are all possible pairs. QED.

More generally, the transformation preserves the cross-edge square predicate in both directions, even without assuming that every cross sum is a square: 25(a+b) is an integer square exactly when a+b is. It also removes every possible within-part edge.

If A and B both lie in [1,M], the resulting labels lie in [1,25M+1]. A bound only on B does not supply this height bound unless the size of the elements of A is also controlled.

## 2. Cited consequence: induced K_3,m for every m

Dujella–Elsholtz, *Sumsets being squares*, Acta Mathematica Hungarica 141 (2013), 353–357, Theorem 1.1, provides positive three-element sets A and corresponding B contained in [0,N], with A+B consisting of squares and |B| growing without bound. See [the original paper, Theorem 1.1](https://www.math.tugraz.at/~elsholtz/WWW/papers/papers44.pdf) and [DOI 10.1007/s10474-013-0334-8](https://doi.org/10.1007/s10474-013-0334-8).

Discard 0 from B if present, losing at most one element. For any prescribed positive integer m, take m elements of a sufficiently large remaining B and apply the cleanup lemma. Therefore:

**Corollary.** Every K_3,m has an induced realization by distinct positive integer square-sum labels.

This is an existence consequence with anchors allowed to vary. It does not assert a common fixed three-anchor set for every m, a bound for arbitrary prescribed anchors, or a quantitative height estimate derived solely from the displayed bound on B. It therefore does not conflict with the separate prescribed-seed obstruction: a universal extension rule for three specified neighbors can fail even after common square scaling.

No conclusion conditional on Bombieri–Lang is used here.

## 3. A direct K_3,3 witness and its exact finite threshold

The two parts

    A = {1,193,673},
    B = {3,168,483}

already induce K_3,3 without applying the affine cleanup. Their cross-square roots, in the displayed row and column orders, are

    2   13   22
    14  19   26
    26  29   34.

The three within-A sums are 194,674,866; the three within-B sums are 171,486,651. None is a square. All six labels are distinct and positive.

**FINITE-EXACT threshold.** The least positive integer ceiling supporting K_3,3 is **673**, both for a subgraph realization and for an induced realization.

The lower bound is certified by the standalone `verify_bipartite_bridge.py`, using a complete pair-intersection scan in the square-sum graph on [1,672]. It checks all 225,456 unordered anchor pairs. For each pair u,v, put S=N(u) intersect N(v). A third left vertex w can complete a K_3,3 exactly when at least three members of S are adjacent to w. The implementation counts these incidences through the neighbor sets of vertices in S. No third vertex outside {u,v} has count at least three for any pair. This excludes every K_3,3 subgraph, regardless of possible within-part edges. In particular it excludes every induced copy. Monotonicity excludes all smaller ceilings as well.

At ceiling 673 the same independent scan finds the displayed witness, whose complete adjacency check supplies the induced upper bound. This finite exhaustive result does not establish any general formula for the minimum height of K_3,m.

## 4. Exact controls and scope

The verifier uses only exact integer arithmetic and no imported producer code. Normal and optimized Python give byte-identical receipts with **269,721 explicit checks**. It checks the affine cross-square equivalence for all source pairs 1<=a,b<=300, every within-part pair in that range, the six-label witness, and the complete pair-intersection scans at ceilings 672 and 673.

The affine lemma is proved for all positive integer inputs. Arbitrary-m existence uses the cited sumset theorem. The height 673 is a bounded exhaustive certificate with an explicit matching witness. These three scopes are kept separate.

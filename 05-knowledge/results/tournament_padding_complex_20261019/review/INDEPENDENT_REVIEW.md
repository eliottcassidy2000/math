# Independent review: exact tournament padding complex

2026-10-19. **PASS.** All displayed identities, exactness assertions, rank formulas, and parity-position conclusions in [the proof](../PADDING_COMPLEX.md) are valid. This review certifies the mathematics, not external priority or a classification of good trees.

Final reviewed proof SHA256: `5d864fcad7b47b836065adec1d4b97e8daae11951a310e4727c5323daddd8424`. The author incorporated the explicit d_(-1)=0 notation during review.

## Quantifiers and edge cases

1. The spaces are finite-dimensional F_2-valued isomorphism-invariant functions on all tournaments of each order. Both deletion summation and universal-source evaluation preserve invariance. Multiplicities of isomorphic vertex deletions must be counted, as they are in the proof and checker.
2. The double deletion identity counts each unordered pair twice, including when the resulting subtournaments are isomorphic for additional reasons. No freeness assumption is needed.
3. In the homotopy identity, deletion of the newly added source contributes the original value once; all old-vertex deletions match the terms of d h. At n=0 the unique empty tournament and its unique one-vertex augmentation give h_1 d_0=id. The maps h_0 and d_(-1) are zero. The author now explicitly names d_(-1)=0, removing the initial minor notation concern.
4. Exactness holds in every degree including zero: d f=0 implies f=d h f, while d squared zero supplies the reverse containment. At degree zero d_0 is injective, so its kernel and the preceding zero image both vanish. Source adjunction is an explicit preimage certificate; it does not assert that the preimage has a single tree representative.
5. For an n-vertex oriented graph, embeddings are injective into an n-vertex tournament and hence bijective. The isolated-vertex formula is obtained by choosing its image. Under source adjunction, the vertex mapped to the new source must have indegree zero, and that condition is sufficient because all its outgoing edges are available. Isolated vertices are included; all possible preimage vertices count with their actual multiplicities. For the empty graph, d_0 gives the one-vertex graph and the h_0 statement is the stipulated zero map.
6. A spanning tournament embedding is an isomorphism. Its number at a tournament of the same type is the automorphism-group order. A nonidentity involution has a transposed vertex pair and would reverse the arc between them, impossible. Cauchy's theorem therefore gives odd automorphism-group order, also for the empty tournament. The tournament embedding parities are exactly the class-indicator basis, so they span the whole function space.
7. Rank-nullity and exactness give t_n=r_(n-1)+r_n. Solving this recurrence gives the stated alternating integer sum. With t_0,...,t_5=1,1,1,2,4,12, the listed ranks are 1,0,1,1,3 and r_5=9. The distinction between direct matrix verification through d_4 and theorem-derived r_5 is correctly retained.
8. Constants satisfy d_n 1=(n+1 mod 2)1. For odd n, an image equal to one would contradict d_(n+1)d_n=0, since n+2 is odd. Thus every even-order oriented graph with at least one isolated vertex fails constant odd embedding parity. No restriction on its other components is needed.
9. If n is even and d_n f=c times one, applying the contracting identity yields f=c times one+d_(n-1)h_n f. This also works at n=0. It does not imply f itself is constant. The kernel of padding can contain nonconstant functions.
10. The order obstruction is sharp: a directed path on positive even n has constant odd spanning parity, and its one-isolate padding evaluates to n+1 odd contributions. The one-vertex graph comes from the empty graph and handles the remaining odd order. Multiple isolated vertices still force zero parity by d squared zero.

## Rédei source check

The only non-elementary named input in the sharpness examples is the classical odd-Hamiltonian-path theorem. A directly inspected primary modern treatment is Thomas Schweser, Michael Stiebitz, and Bjarne Toft, *The Tournament Theorem of Rédei revisited*, arXiv:2510.10659v2 (18 February 2026), https://arxiv.org/html/2510.10659v2 . Section 1 explicitly derives odd Hamiltonian-path parity from Rédei's stronger theorem. Paths there are vertex sequences, which matches injective embeddings of a fixed directed path. Its Section 3 supplies a proof of the stronger theorem. The use and attribution in this packet are correct.

## Independent exact controls

`independent_controls.py` uses directed edge sets, literal bijection enumeration, and explicit images of linear maps. It imports no producer code. Its normal and optimized runs agree byte-for-byte on **8,832 checks**:

- Every oriented graph through four vertices: counts 1,1,3,27,729.
- The source-deletion embedding identity through graph order four, against every labelled host of the required lower order.
- The padding embedding identity through source order three, against every labelled host of the required higher order.
- The even-order isolated-vertex obstruction for all two- and four-vertex oriented graphs with an isolate.
- Independent tournament classes through order four; class-indicator embedding basis, matrix ranks, contracting identity on every function, square-zero identity, explicit kernel primitives, and constant formulas, including degree zero.
- Direct padded-path parity for every labelled tournament at orders one, three, and five.

The producer checker was rerun normally and with optimization. Both reproduce its saved 35-check receipt exactly, including classes through order five and the direct ranks through d_4. These bounded computations support the proof; none is used to extrapolate its universal quantifiers.

No files inside `math/` were edited, and no Git write was performed.

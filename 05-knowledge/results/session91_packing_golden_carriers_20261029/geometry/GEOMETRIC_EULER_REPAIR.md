# Sparse carriers that retain actual unit-distance geometry

Repository source pin: `90e704a8f35614d2085778d1b518f03ff0d750ca`. Named packets below are under `05-knowledge/results/` unless a different path is given.

**Status:** PROVED construction; exact finite controls below. This extends Session90's disjoint Euler repair with an explicit geometric padding family. The graph-theoretic mechanisms are elementary; no external-priority claim is made.

## The representation contract

A faithful unit-distance realization assigns distinct planar points to vertices and has an edge exactly when two points have Euclidean distance one. Crossings of drawn segments are allowed. A unit-disk contact realization additionally requires every distinct pair of points to be at distance at least one. The second condition is stronger: an ordinary unit-distance graph can have arbitrarily high degree, whereas the additional separation imposes a kissing constraint.

These are different again from the extremal number of unit distances among n points and the Hadwiger–Nelson chromatic number. Nothing below solves either extremal problem. The realization convention is that of Alon–Kupavskii, [Two notions of unit distance graphs](https://arxiv.org/abs/1306.3916), also used in the inherited `odd_sign_spectra_trijections_20261028/LRC_UNIT_DISTANCE.md`.

**Theorem.** Given a connected finite simple graph Q of cycle rank beta=e(Q)−v(Q)+1, there is an explicit graph P(Q) with:

1. e(P)=v(P)−1;
2. the **abstract graph** group Aut(P) naturally isomorphic to Aut(Q);
3. the same Royle endpoint-reversal sign character as Q;
4. P faithfully unit-distance realizable in the plane iff Q is;
5. P realizable with all point separations at least one iff Q is.

If a realization of Q is supplied, the extension is explicit and keeps all its points unchanged. The added vertex count is O(beta log(beta+1)). This does not preserve a containing square, packing density, total diameter, or a specified connected contact pattern.

## An integer-coordinate rigid tree family

For beta>=1 put b=ceil(log2 beta), L=3b+12. For each integer code k with 0<=k<beta, define T_k by the points

- spine: (i,0), i=0,...,L;
- mandatory tooth: (3,1);
- optional tooth: (6+3j,1) when bit j of k is one, j=0,...,b−1.

Use distance-one pairs as edges. Each tooth is joined only to its spine point, and spine edges join consecutive integers. Distinct teeth are at least three horizontal units apart. Therefore these are trees, every pair of points has distance at least one, and the description is a faithful unit-distance realization with integer coordinates.

The full spine is the unique diameter path: each tooth lies at least three spine steps from the left endpoint and at least nine from the right endpoint, so no tooth can replace a spine endpoint in a maximum-distance pair. An automorphism must therefore preserve the spine and either fix it or reverse it. Reversal fails because there is a tooth at position3 and no tooth at L−3. Thus each T_k is rigid. The orientation of its spine is intrinsic, and the optional teeth recover k, so all T_k are pairwise nonisomorphic.

The number of vertices is L+2+popcount(k)<=4b+14. This deliberately simple family prioritizes exact recognition and geometry, not minimum constants.

## Proof of the theorem

For beta=0, set P(Q)=Q. For beta>=1, let P(Q) be the disjoint union of Q and T_0,...,T_(beta−1). Each added tree has one fewer edge than vertices, so

    e(P)−v(P) = (beta−1)−beta = −1.

Q is the unique cyclic component, so every automorphism preserves it. The added trees are rigid, distinct, and acyclic, giving Aut(P)=Aut(Q), including their actual actions on the original vertices. Every added edge is fixed pointwise, so its contribution to Royle endpoint-reversal parity is zero. This proves the sign-character statement, not merely preservation of the numerical automorphism order.

For a supplied planar realization of Q, compute its coordinate bounding box. Translate the first padding tree so that its leftmost x-coordinate lies at least two units to the right of Q's maximum x-coordinate. Place each later tree at least two units to the right of the preceding tree's maximum x-coordinate. The y-coordinates of each padding tree remain0 or1. Any points in distinct components then differ in x by at least two, so their Euclidean distance exceeds one. Thus faithful unit distances are preserved, and so is the optional minimum-separation condition. Rational input coordinates allow rational, indeed integer-increment translations. Conversely any realization of P restricts to a realization of its Q component.

This proof explains why adding disconnected components is a mathematically substantive freedom. The scalar relation e=v−1 can be achieved without changing any difficult local geometry inside Q. The known independence of trijectivity and faithful unit-distance realization is not contradicted: this theorem does not supply a trijection for P(Q).

## Abstract automorphisms are not spatial isometries

The theorem preserves the graph automorphism group, not its realization by Euclidean motions. In fact, when beta>=1, the displayed point realization of P(Q) has **trivial global Euclidean isometry group**. An isometry preserving the point set induces a graph automorphism. That automorphism fixes each distinct rigid padding component pointwise. One such component contains three noncollinear points, forcing the Euclidean isometry to be identity.

For example, a square C4 core has abstract automorphism group D8 and a symmetric point realization. The repaired graph still has abstract group D8, but its enlarged point set has no nonidentity Euclidean symmetry. Its surviving abstract automorphisms preserve unit distances, not all pairwise distances. This is a concrete geometry/combinatorics distinction, and agrees with the separate packing lane's metric versus arithmetic symmetry boundary.

## A precise square-packing convention hostile

Place axis-parallel unit squares at the same centers as T_k. Their interiors are disjoint. Positive-length side contacts correspond exactly to the tree edges. But every tooth also makes corner contacts with the two neighboring spine squares.

Consequently, if t is the number of teeth, the graph of **all boundary contacts** has 2t additional edges. For T_0 at beta=1, L=12 and t=1: there are14 squares,13 positive-length side contacts, but15 boundary-contact edges. The latter graph has e=v+1, not e=v−1.

Thus the same point configuration can represent a faithful unit-distance tree or a square contact graph with cycles, depending on which contact predicate is retained. This is an exact version of the repository's warning that an attachment carrier must preserve the actual observable. It does not relate the enclosing-square optimum of eleven rotated squares to a tree count.

## Links to prior work and limits

- Session90's `trijection/SPARSE_EVEN_CARRIERS.md` supplied faithful uncoloured graph encodings and exact automorphism preservation, but not a unit-distance realization of arbitrary inputs. This theorem only pads a core whose metric realization is already available; it cannot turn a non-realizable core into a realizable one.
- Session90's `capacity/CAPACITY_AND_EULER_REPAIR.md` proves that, for an odd-automorphism cyclic core retained as a component, Omega(beta log beta) padding vertices are required. The present geometric construction attains that order while respecting the extra minimum-distance condition. The lower bound is scoped to disjoint core-preserving repair, not all encodings.
- The latest `odd_sign_spectra_trijections_20261028/LRC_UNIT_DISTANCE.md` gives faithful theta carriers without weak trijections and a rigid trijective carrier containing K_(2,3), which has no planar unit-distance realization. Those are inherited separation examples.
- Arbitrary high-degree unit-distance stars are not counterexamples to packing kissing bounds because their neighboring points need not be separated from one another.

The useful positive bridge is exact: graph Euler repair can be made geometric without losing the source's symmetry or distance relation. The missing quantities are the containing region and its density/optimality. No stronger packing theorem follows.

Run `python3 check_geometric_repair.py` and `python3 -O check_geometric_repair.py`. The checks use integer squared distances, independent tree-center automorphism signatures, component counts and a corner-contact hostile. No floating tolerance or external package is used.

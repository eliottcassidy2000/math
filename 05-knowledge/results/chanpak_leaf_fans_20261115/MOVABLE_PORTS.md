# Moving unit ports: exact bounded tests and recursive obstruction

Use `L_q` for appending a vertex that beats exactly the old vertex q,
and write H for the number of directed Hamiltonian paths and
delta(q) = H(L_q T) − 2H(T). The finite tests concern the raw signs
(delta = ±1), not just repeat-stable signs.

## Checked bounded result

The standalone census runs a subset DP for Hamilton paths starting at
each vertex of each induced subset. If E_q is the sum, over even-order
paths starting at q, of H of the complement, then delta(q) = 2E_q − H.
This is the one-new-vertex term in the odd-cycle-packing expansion.
For each source carrying both signs the census really constructs every
unit child and recomputes the child's entire port menu.

* Among all labelled tournaments through order 6, the mixed-sign totals
  at orders 2–6 are 2, 6, 48, 360, 2880. The order-6 examples have
  H = 1, 5, 9, 37 (one rigid isomorphism class apiece; the census
  records the complete ordered SCC forms).
* Among chronological histories from a singleton, with arbitrary old
  parent at each birth, the mixed-sign totals at orders 1–8 are
  0, 1, 1, 5, 9, 9, 9, 9. From order 5 through order 8 the nine
  histories comprise the transitive singleton stack (one history),
  a padded strong order-4 H5 core (four), and a padded strong
  order-5 H9 core (four). These are not claimed to exhaust larger
  orders.
* In either census, not even one positive child of a mixed source still
  has both signs. The negative child sometimes does; for example
  (6,12,8,1) encodes an H5 source of defect vector (1,3,−1,−1);
  appending at parent 2 gives H9 with both signs. Appending at
  parent 0 instead gives H11 with defect vector (1,15,7,3,−3).

There is a scalar reason the first two H values cannot start a full
depth-two binary tree: H1 would produce 7, and H5 would produce
21. The failures at H9 and H37 are additional finite structural
obstructions, not consequences of those scalar exclusions.

## Absolute-target binary skeleton

The two maps h ↦ 2h±1 have a unique odd inverse when the target N
is odd. For N ≡ 1 mod 4 that inverse is (N+1)/2 and requires a
negative port; for N ≡ 3 mod 4 it is (N−1)/2 and requires a positive
port. The exclusion of 7 and 21 therefore forces independent entry
at their children 13, 15, 41, 43: no sequence of unit moves can
produce those targets without using an unattainable parent. A
seed at each of these children and at 1, plus proof of directed
closure along the admissible odd-parent edges, would suffice for
coverage. Mere coverage of parent values is not directed closure:
the chosen physical witness has to carry the child's required
continuation.

A full binary tree from a single seed h has depth-k interval
[(h−1)2^k+1,(h+1)2^k−1] in odd steps. If h>3 it misses open
intervals between depths. A nontrivial full binary seed, even
if found, would not by itself prove spectrum completeness.

## Vertex extension check through order eight

The C++ replay picks the least unseen order-7 binary orientation
code and marks its entire relabelling orbit. There are 456
representatives. At order 8, extend each of those by every allowed
outneighborhood of the eighth vertex (58,368 extensions); this
covers every isomorphism class by deleting one vertex, even though
some classes occur several times. This also has no positive child
that keeps both signs. The order-8 mixed sources reach H as high
as 473, and throughout orders 2–8 only mixed sources of H1 or H5
have a negative child that keeps both signs. None of these are
claimed to forbid a larger unrestricted example.

## Repeat-stable re-selection has zero binary entropy

The cross-port stability criterion says that after appending at a
repeat-stable q, a distinct old p is repeat-stable exactly when
q's old outneighborhood is {p}; the newborn would be repeat-stable
only if q lay outside the old terminal SCC. But a repeat-stable
q must lie inside the terminal SCC: otherwise the endpoints of
all Hamilton paths of T−q lie in the unchanged terminal SCC, and
q beats every endpoint. So the newborn cannot enter the stable
menu after a stable append.

If a move ever switches from stable q to the unique possible
distinct stable p, then p has outdegree at least two by that time.
Indeed p could not have been a sink before the first move at q:
then all paths of T−q would end there, in q's outneighborhood
and contradicting q's stability. In addition p beats each vertex
newly attached at q. The first append at p therefore leaves
no other stable port, and successive appends at p keep its
outdegree unchanged. In particular, any chain of repeat-stable
moves switches identities at most once.

From a seed with r initially stable ports, there are at most
r k histories of exactly k stable attachments, by the first
root and the possible time of the switch. A run of a repeated
stable moves at a nonsink increases H by 2 E (2^a−1), with E>0
the initial even-prefix sum. Hence each nonsink run in an
output of count H≤X has length O(log X). A sink run cannot switch
and produces no new count. Any finite seed pool generates only
O(log² X) distinct totals up to X by repeat-stable attachments,
even allowing all defects, so this entire grammar alone cannot
be cofinite in the odd numbers.

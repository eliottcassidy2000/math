---
id: THM-4988
title: "The 7/21 gaps disappear at endpoint boundaries; faithful graph carriers preserve the specified action"
status: >
  PROVED, independently audited in Session232. FINITE-EXACT controls have
  their stated universes; inherited all-order theorem certificates remain
  explicit dependencies, not newly rerun claims. Chan-Pak completeness OPEN.
session: Session232-automorphisms-20261109
source: 05-knowledge/results/session232_automorphisms_20261109/boundary_graphs/PROOF.md
audit: 05-knowledge/results/session232_automorphisms_20261109/graph_audit/REPORT.md
---

# THM-4988

[Full proof](../../05-knowledge/results/session232_automorphisms_20261109/boundary_graphs/PROOF.md) · [Independent audit](../../05-knowledge/results/session232_automorphisms_20261109/graph_audit/REPORT.md)

For a finite nonempty tournament T, let G=Aut(T), e(v) count Hamilton paths ending at v,
and G_v fix v. The endpoint map gives an exact fibre decomposition

    P(T)/G -> V(T)/G,
    fibre over Gv = P_v(T)/G_v,
    q_T(v)=e(v)/|G_v|,   H(T)/|G|=sum_[v] q_T(v).

Transport an endpoint to v; two paths now ending at v are G-equivalent
exactly when they are G_v-equivalent. The latter action is free. This proves
the formulas. The vertex stabilizer cannot be replaced by the full group:
P7 has e(v)=27, |G|=21, |G_v|=3, giving9 rather than the nonintegral27/21.

## Sharp boundary witnesses

A bitmask lists unordered pairs i<j lexicographically; bit1 means i->j.
The six-vertex mask26 is asymmetric, has H=19, and terminal vector
(1,5,7,3,2,1). The seven-vertex mask290 is asymmetric, has H=55, and terminal
vector(3,9,21,9,7,3,3). In both cases v=2 gives the named endpoint quotient.
Thus 7 and21 occur at this intrinsic boundary despite their global exclusion.

These minimum orders,6 and7, are exact: deletion of the final vertex injects
paths ending there into Hamilton paths of T-v. Exhaustive maxima on orders
1,...,5 are1,1,3,5,15, excluding smaller witnesses for either raw endpoint
count or its stabilizer quotient. The displayed finite counts are reproduced
by literal permutation enumeration and separate subset DP.

## Infinite path-action graphs

Prepend k ordered C3 components to either witness B. The full group is C3^k,
fixes the designated endpoint, and both raw total and endpoint counts multiply
by3^k. With generators rotating the triangle components, the path-action
graph is respectively19 or55 copies of the Hamming graph H(k,3). Its endpoint
subgraph consists of exactly7 or21 such components. A free group action has
Cayley-graph components, and C3^k with these generators gives H(k,3).

This graph records a specified action. Its full abstract graph-automorphism
group need not be C3^k; generator labels and the action are the retained data.

## Faithful ordinary simple graph encoding

For n>=3 replace each arc u->v by a fresh path u--L--R--v, attach one leaf to
L and a pendant path of length2 to R, and leave original vertices untagged.
Call the graph Gamma(T). It has n+5*binom(n,2) vertices,6*binom(n,2) edges,
and2*binom(n,2) leaves.

Leaf peeling recovers the 2-core, a subdivided K_n. Core vertices with pendant
length0,1,2 are intrinsically the original,L,R roles. Thus Gamma(T) decodes T;
every graph automorphism restricts to a tournament automorphism and each
such automorphism extends uniquely. Hence Aut(Gamma(T)) is isomorphic to G.
Its odd order makes every vertex-permutation sign and the graph's
edge-reversal sign positive. This is not an even-degree assertion.

Expanding the adjacent arcs of a Hamilton word yields a bijection with
intrinsically specified admissible original-vertex traversal traces in Gamma(T).
It is equivariant, so complete and endpoint path-orbit assertions transport.
It does NOT preserve ordinary graph Hamilton-path counts: Gamma(T) has at
least6 leaves and no spanning graph path. The trace definition is indispensable.

Independent uncoloured graph enumeration covers all72 labelled order3/4
inputs and six larger controls, including the112-vertex encodings with
Aut orders7 and21. The structural theorem is all-order; the enumeration is
a control. Inherited global gaps remain THM-1370 and THM-4812 with their
full paths listed in the source proof. Chan–Pak remains OPEN.

# A global converse-forest obstruction at every order at least fifteen

## Result and scope

For every integer n >= 15, two distinct ordinary tournament isomorphism classes have the same ordered pair of exact incidence signatures

    (E(T), E(T^op)).

Consequently the signature-incidence multigraph for E33, or for any fixed quotient of its nine symbols, has parallel distinct class edges at every such order. It is therefore not a forest. The same obstruction applies to the unquotiented exact E signature itself.

This rules out an all-order **global forest certificate** using these signatures. It does **not** establish nonuniqueness of the actual endpoint-constrained port wiring: capacities, supported cells, parity, or other authenticated constraints can still resolve a cyclic signature-incidence graph. Nor does it contradict the order-eight three-color result. Padding below is not strong, so no all-order strong-tournament statement is made here.

## Exact seed

Let T_i be the circulant tournament on Z/15Z with u -> v exactly when v-u lies in S_i, where

    S_1 = {14,13,12,11,5,6,8},
    S_2 = {14,13,12,4,10,9,8}.

Write c_v for the number of cyclic triangles through v and c_uv for the number through the unordered pair {u,v}. Then

    E(T) = multiset_u [ multiset_{u -> v} (c_v,c_uv) ].

Direct cyclic-triple enumeration gives, for both seeds and both converses, fifteen copies of the same row:

    [(28,3),(28,3),(28,3),(28,4),(28,4),(28,5),(28,6)].

The classes are distinct because their exact numbers of directed Hamiltonian paths are respectively 196264905 and 196277235. The separate checker recomputes both counts by subset dynamic programming; it does not infer nonisomorphism from a numeric class identifier.

## Two exact padding operators

Every cyclic triangle through u uses exactly one outgoing edge of u. Thus

    c_u = sum_{u -> v} c_uv.

In particular E determines the multiset of root counts c_u by summing the second coordinates of each row.

Let P_+(T) add a dominating source s, so s -> v for every old vertex v. No new cyclic triangle appears. All old E rows are unchanged, and the new row is the multiset of (c_v,0) over all old vertices v. Since E determines the multiset of c_v, equality of E is preserved by P_+.

Let P_-(T) add a dominated sink z. Again no new cyclic triangle appears. Each old row gains one (0,0), and the new row is empty. Hence P_- also preserves equality of E.

Converse exchanges the operators, with the new vertex kept in the same position:

    P_+(T)^op = P_-(T^op).

Therefore equality of both E orientations is preserved whenever a new source is added. Iterating P_+ on the two seeds gives the required equal ordered pair at every order n = 15+r.

The classes remain distinct: a dominating source must be first in every directed Hamiltonian path, and deleting it gives a bijection with the Hamiltonian paths of the old tournament. Both original unequal path counts are therefore preserved through every iteration.

## Why the forest certificate fails

For a signature S, the converse signature-incidence multigraph has left vertices S(T), right vertices S(T^op), and one edge for each ordinary tournament class T. The two distinct padded classes constructed above join exactly the same pair of vertices whenever S is E, E33, or a quotient of E33. The edges must be retained separately. They are parallel edges, a two-edge cycle in a multigraph, so the graph cannot be a forest.

Equivalently, the two corresponding columns in the margin-incidence matrix are equal. Their difference is a nonzero kernel vector. This is a failure of the **global** margin-only forest criterion, not a certificate that the vector can be realized inside a particular nonnegative, endpoint-constrained port table.

## Replay and finite controls

`forest_obstruction_check.py` uses only the Python standard library and no imported producer code. It recomputes the exact seed signatures in both orientations, both Hamiltonian counts, the outgoing triangle-balance identity, and source/converse-sink signature transformations through order 25. The all-order conclusion follows from the proved induction above, not from those finite controls.

    python incidence/forest_obstruction_check.py
    python -O incidence/forest_obstruction_check.py

Both modes reproduce `forest_obstruction_receipt.json`. The checker raises explicit errors rather than using removable Python assertions. The receipt records the full common seed signatures and its explicit check count. This is a separate addendum; the earlier frozen incidence, converse, ternary, and peeling-depth packets are unchanged.

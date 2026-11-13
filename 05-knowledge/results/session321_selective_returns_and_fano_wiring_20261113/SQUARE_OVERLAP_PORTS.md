# The first square-sum return is a change of boundary channel

**Proved:** the ear/cut bijections and sharp surgery costs below, conditional only on the displayed finite exact census. **Finite-exact:** the paths and exposed covers on the stated ranges. This is neither a new eventual square-sum existence theorem nor a proof of Hadwiger–Nelson, BSD, or Chan–Pak.

## 1. Source transfers

- [Square-sum Hamiltonicity and the Gerbicz construction](../collatz_mod6_20260922_w6_square_sum_hamiltonicity.md) distinguishes endpoint-aware degree-two forcing from the invalid assumption that every degree-two vertex is internal. [The graceful and square exchange note](../graceful_square_exchange_20261019.md) gives the four-cycle requirement for a two-edge fixed-endpoint switch on a fixed vertex set. Here adding two new vertices gives a different native operation: substitute a three-edge ear for one exposed matching edge, and repair the other exposed pair.
- [Hadwiger–Nelson interface profile](../six_colour_interface_profile_20261015f.md) and [component palette matching](../six_colour_component_pairs_20261016a.md) retain the full labeled relation across a boundary, not just whether each port can occur separately. The finite transfer here uses exactly that distinction: count oriented path pieces jointly across prescribed ports, and keep their connectivity pairing. This is a shared proof discipline, not a reduction of chromatic number to square-sum paths.
- [BSD component recursion](../bsd_component_recursion_20261013c.md) proves the exact block identity \(\Phi(A(T_1\to T_2))=\operatorname{diag}(\Phi(A(T_1)),\Phi(A(T_2)))\) for \(\Phi(X)=X^2+X+1\) over \(\mathbb F_2\), before interpreting nullity in the stated Monsky arithmetic setting. Its strength is the explicit cancellation of the interface block, not extrapolation from numerical ranks. The square-sum channel calculation below takes the complementary route: **do not cancel** the boundary pairing, because a legal-degree pairing can close a cycle instead of joining paths. The same repository notes separate Selmer dimension from actual elliptic-curve rank. No elliptic interpretation of the square-sum counts is claimed.

## 2. Exposed-cover factorization

Let \(Q_n\) have labels \(1,\ldots,n\) and edges between distinct \(u,v\) iff \(u+v\) is a square. Paths are identified under reversal. Write \(E_G(a,b)\) for the number of pairs of nonempty disjoint path words in \(G\), each starting at its prescribed mark \(a\), respectively \(b\), whose supports partition \(V(G)\). Their other endpoints are unconstrained. This is an *ordered by marks* pair, so there is no factor of two.

If new vertices occur in one forced ear \(a-v_1-\cdots-v_\ell-b\), deleting the ear and retaining the two old pieces gives an exact bijection between Hamiltonian paths using the ear and \(E_G(a,b)\). The inverse reverses the piece starting at \(a\), traverses the ear, and appends the piece starting at \(b\). Equivalently, \(E_G(a,b)\) counts Hamiltonian paths using the distinguished virtual edge \(ab\) in \(G+ab\), but the two-piece count retains the interface even if this edge is not a legal square sum.

Finite enumeration for sizes \(15,\ldots,25\), in order, gives
\[
1,1,1,0,0,0,0,0,3,0,10.
\]
The first two moves are just endpoint appends: attach 16 to the endpoint 9, then attach 17 to the endpoint 8. Every old edge survives both steps. On the other hand the three Q23 paths lose \(6,6,5\) old edges from **each** of the unique Q15, Q16, and Q17 paths, when ordered by Q23's other endpoint \(2,9,22\) opposite 18. Thus the two endpoint appends are genuine nested overlaps but do not carry across the gap.

At size 24 the newcomer has neighbors 1 and 12. It cannot be an endpoint: deleting such an endpoint would give a Q23 Hamiltonian path ending at 1 or 12, neither of which occurs. Thus Q24 uses exactly the virtual channel \(E_{Q_{23}}(1,12)\). It has count zero.

At size 25 the newcomer has exactly neighbors 11 and 24, so it is internal, as its removal as endpoint would leave a Q24 path. The preceding vertex 24 cannot be an endpoint of the resulting full path: the endpoint segment \(24-25-11\) would leave a Q23 path with endpoint 11. Thus the only possible ears are
\[
11-25-24-1,\qquad 11-25-24-12.
\]
They are mutually exclusive. The count factors into the disjoint union of the **two** marked channels, not a scalar continuation of the size-24 channel.

| channel | exposed cover count | second endpoint opposite the compulsory 18 |
|---|---:|---|
| \(E_{Q_{23}}(1,12)\) | 0 | none |
| \(E_{Q_{23}}(11,1)\) | 4 | 3, 8, 10, 23 |
| \(E_{Q_{23}}(11,12)\) | 6 | 2, 4, 9, 11, 13, 22 |

The replay enumerates exposed covers directly on complementary old supports, independently of its full-path recursion, and tests that expanding the ears gives exactly the Q25 paths. The loss of Q24 and recovery at Q25 are explained by the change of required boundary relation, not by adding a constant count.

## 3. The exact overlap and the endpoint-conditioned obstruction

The number of deleted old edges from each Q23 path to each Q25 path is the following. Rows use the second endpoint of Q23 opposite 18; columns are the second endpoint of Q25 opposite 18.

| old / new endpoint | 2 | 3 | 4 | 8 | 9 | 10 | 11 | 13 | 22 | 23 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 2 | 2 | 2 | 2 | 3 | 3 | 3 | 2 | 3 | 3 |
| 9 | 3 | 3 | 3 | 3 | 2 | 4 | 4 | 3 | 4 | 4 |
| 22 | 3 | 2 | 2 | 2 | 4 | 3 | 3 | 2 | 3 | 3 |

**At least two cuts are necessary** even if the endpoints may move. The two old endpoints of the ear, 11 and \(z\in\{1,12\}\), are not endpoints of any Q23 path. Both must lose an incident edge, and there is no Q23 edge \(11-z\) to pay for both. The overlap table shows that two cuts are sufficient for each old path when the second endpoint may change.

**Same-endpoint sharp costs are \(2,2,3\)**. With exactly two cuts and the same old endpoints, one cut is incident to 11 and one to \(z\), and after adding the ear one residual old edge must join the other exposed endpoints. On all three old paths the possible opposite ends at 11 are 5 and 14. At \(z=1\) they are 3 and 8; at \(z=12\) they are 4 and 13. The only square-sum pairing across those sets is \(5+4=9\). Thus the two-cut surgery is forced: delete \(11-5\) and \(12-4\), add \(4-5\) and the \(z=12\) ear.

The final gate is connectivity. In general two separated cuts split an old path into left, middle, and right pieces. Label the exposed ends in old order \(a,b,c,d\). Pairing \((a,b),(c,d)\) restores old connectivity. Pairing \((a,c),(b,d)\) also reconnects a spanning path, reversing the middle piece. Pairing \((a,d),(b,c)\) instead joins the outer pieces and closes the middle piece, so it fails; degree data alone do not see this. Subdividing one pairing edge by a new ear does not change the verdict.

For endpoint 2 and endpoint 9 the forced operation uses the connected cross pairing. For endpoint 22, orient the old path from 18 to 22. The two cuts leave a middle piece with ends 11 and 12, so the ear closes a cycle there and the last pairing fails. The three-cut witness deletes
\[
3-22,\quad11-14,\quad12-13
\]
and inserts \(3-13,\ 14-22\) and the \(11-25-24-12\) ear. Its actual old and new paths, with common orientation, are

    18 7 9 16 20 5 11 14 2 23 13 12 4 21 15 10 6 19 17 8 1 3 22
    18 7 9 16 20 5 11 25 24 12 4 21 15 10 6 19 17 8 1 3 13 23 2 14 22

This is a finite model of a native repair rule: keep exposed port labels, which pieces they belong to, and which actual endpoints must survive. A move that balances an interface but fails the connectivity pairing cannot be promoted to a global return. That is the precise ingredient to test when carrying exact-target tournament repairs between families.

Replay: run python3 replay_square_ports.py or python3 -O replay_square_ports.py.

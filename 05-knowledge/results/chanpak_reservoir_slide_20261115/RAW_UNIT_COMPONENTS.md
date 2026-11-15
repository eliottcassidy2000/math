# Raw unit moves are confined to a boundary strong component

This conclusion does **not** require a repeat-stable root. It
describes all reducible tournaments and shows that ordinary
count-preserving prefix/suffix padding cannot rearm unit fans.

As in [the leaf law](../chanpak_leaf_fans_20261115/CROSS_PORTS.md),
write \(\delta_F(q)=2E_q-H(T)\), where \(E_q\) sums over
nonempty even-order paths starting at \(q\), weighted by the
Hamiltonian path counts of the complement. Reversal gives the
analogous reverse defect.

## Outside the terminal component: a large positive defect

Take components in ordinal order. Let \(A\) be all components
preceding \(q\)'s component \(C\), and \(B\) all following
components. Suppose \(B\ne\varnothing\). Choose Hamilton paths
\(P_A,P_C,P_B\) of the three blocks, using the empty path for
an empty \(A\). Write \(P_C=U\,q\,V\),
where \(U,V\) may be empty. Form the path
\(qV\) followed by either none or the first vertex of \(P_B\),
choosing the unique option with even order. Such an option exists
because \(B\) is nonempty, and the join from \(C\) to \(B\) is
forward. The complementary path is \(P_A\), then \(U\), then the
unused suffix of \(P_B\), again with valid forward joins.

The produced pair of paths remembers the original tuple: use
the intersection of each path with \(A,C,B\), and reunite the two
pieces of \(P_C\) and the two pieces of \(P_B\) in their known
order. Therefore
\[
E_q\ge H(A)H(C)H(B)=H(T),
\quad\text{and}\quad \delta_F(q)\ge H(T).
\]
In particular, if \(H(T)>1\), no such root can have defect
\(+1\) or \(-1\). No conditional or asymptotic estimate is used.

## Inside the terminal component: a factor law

If \(q\) belongs to the terminal component \(C\), no starting
path can leave \(C\). Every complement Hamiltonian path is a
concatenation of its path in the old prefix \(A\) and its path
in the remaining part of \(C\). Therefore
\[
E_q(T)=H(A)E_q(C),\qquad
\delta_F^T(q)=H(A)\delta_F^C(q).
\]
A forward unit defect forces \(H(A)=1\) and a true unit defect
in \(C\). For \(H(T)>1\) this \(C\) must be nontrivial.
Every nontrivial strong component has at least three paths by
rotating its Hamiltonian cycle, so \(H(A)=1\) is equivalent to
all preceding components being singleton. The converse holds
by the same displayed product. Reversing the tournament gives
the exact dual for reverse unit defects.

Thus with \(H(T)>1\):
* Two or more nontrivial components forbid **all** forward and
  reverse raw unit moves.
* A sole nontrivial component surrounded by singletons on both
  sides also forbids every raw unit move.
* If the sole nontrivial component is terminal, all forward
  unit moves are exactly that component's unit moves; with a
  nonempty prefix, no reverse unit move exists. The initial
  case is dual.
* A tournament having both a raw forward unit and a raw reverse
  unit **must already be strong**.

So inserting neutral sources/sinks at the extremes does not
create new unit moves even at unstable roots. Padding can
preserve the direction whose target remains at its proper
boundary, and disables the other orientation. It cannot fix
the loss of signs in the [two-ended examples](TWO_DIRECTION_RETURN.md).
This only addresses ordinary ordinal padding, not count-preserving
rewiring of the strong component itself.

For the transitive exception \(H=1\), a port with \(d\) old
outneighbors has forward defect \(-1\) for \(d=0\) and
\(2^d-1\) otherwise; reverse uses indegree. In particular
the forward \(+1\) is the next-to-last vertex whenever it exists.

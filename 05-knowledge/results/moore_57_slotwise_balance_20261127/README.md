# Slot routing strengthens the degree-57 tournament obstruction

2026-11-27. The [rooted balance comparison](../chan_pak_gap_and_moore_balance_20261127/README.md)
gave a sharp block-mean benchmark and an entrywise lower bound
of 158 for the dense encoding. Here the actual Moore route
assignment gives 160 for the Paley-59 deletion and every regular
57-vertex quotient. There are also score-sensitive bounds and
exact transport-group checks in the two known smaller examples.
All degree-57 results are conditional on a Moore-perfect rooted
matching table and use the indicated dense encoding only.

## 1. The exact block expansion

Let \(m=56\), \(S\) be the 57-vertex skew tournament quotient,
\(g_{ij}=(S^2)_{ij}\), and \(C_i\) the internal skew matrix of
branch \(i\). The dense cross block is \(S_{ij}(J_m-2P_{ij})\).
Form the defect digraph by retaining pairs \(|g_{ij}|=3\),
oriented \(i\to j\) when \(S_{ij}g_{ij}>0\).
The orientation is invariant under switching the quotient, since
the two factors gain the same sign. A pair with \(|g_{ij}|\ge5\)
already forces a maximum absolute entry at least 262 by
the earlier block-mean formula and even parity.

For a defect edge \(i\to j\), align the direct matching to
the identity by relabelling \(j\). Multiply the squared
cross block by \(\mathrm{sgn}(g_{ij})\), which also makes the
factor on the internal terms \(+1\). For \(x_a=(C_i\mathbf1)_a\),
\(y_b=(C_j\mathbf1)_b\), the resulting signed entries are
\[
 W_{ab}=156+4q_{ab}+x_a-y_b-2(C_i[a,b]+C_j[a,b]).
\]
Indeed \((J-2P)(J-2P')=(m-4)J+4PP'\); Moore balance says that
the direct and the \(55\) two-step routes partition the slots.
Thus \(q_{aa}=0\), \(q_{ab}=\pm1\) for \(a\ne b\), and
every row and column of \(q\) sums to three. In particular
there are 29 positive entries per row and per column.
The replay checks the underlying expansion against actual
matrix multiplication at all 756 unordered-orientation
cross-block positions for Hoffman–Singleton with
nonconstant internal row sums.

## 2. Five defect steps rule out 158

Suppose every cross entry is at most 158 in magnitude.
In the aligned slot frame, the upper bounds on the row and
column sums give, at each matched slot \(a\),
\[
 54x_a-2y_a\le100,\qquad 2x_a-54y_a\le100.
\]
The 100 is \(56\cdot158-3\cdot54^2\).
Row sums are odd and in \([-55,55]\). Thus a row sum at
least three at the tail of an edge forces at least 31
at the matched slot in its head. A second outgoing step
would force at least 787, impossible. Dually a row sum at
most minus three cannot be at the end of a two-step walk.
Matchings are permutations even if a branch is revisited.

On the middle defect edge of a directed walk of length five,
the tail has a two-step incoming walk and the head a two-step
outgoing walk. Hence \(x_a\ge-1\) and \(y_a\le1\) for every \(a\).
Put \(t_a=x_a-y_a\). Then \(t_a\ge-2\), and some \(t_a\ge0\)
because the total is zero. At such an \(a\), the sets of
positive row and column slots of \(q\) each have size 29
among 55, so share at least three \(b\)'s. On any one of them,
\[
 W_{ab}+W_{ba}=320+t_a+t_b\ge318 .
\]
The internal matrix terms cancel by skew symmetry. This
contradicts the upper bound \(316\). The leaf tournament
has even order, so its off-diagonal squared entries are
even, giving the asserted 160. The proof also works for
a single defect edge with \(x\ge-1,\ y\le1\) from the start.
In particular **any quotient** with all branch interiors
near-regular (row sums \(\pm1\)) forces 160: if there is no
defect edge, some pair has \(|g_{ij}|\ge5\), since a 57-vertex
quotient cannot be doubly regular even after switching.

For a two-vertex deletion of any doubly regular tournament
of order 59, \(g_{ij}=1+u_i u_j+v_i v_j\) off the diagonal,
where the sign columns \(u,v\) go to the deleted contacts.
The four pattern classes form defect cliques oriented
by \(S\). One has at least 15 vertices and therefore gives
a directed walk of length five. In the Paley example deleting
57 and 58, the replay checks \(16\to8\to6\to11\to0\to4\).

For a regular 57-vertex quotient, \(g_{ij}=55-4r_{ij}\) with
\(r_{ij}\) the common predecessor count. If no larger pair
already wins by the mean, the only off-diagonal entries
are \(-1,3\). Since \(S^2\mathbf1=0\), each branch has 28
defect neighbors. The identity \((S^3)_{ii}=0\), with
\(S\mathbf1=0\), gives zero sum of \(S_{ij}\) over these
neighbors, so each vertex has 14 outgoing and 14 incoming
defect edges. That supplies the walk (or a directed cycle).
The conclusion extends to any switching of a regular quotient.

## 3. Score-sensitive screening, including transitive interiors

Let \(R=\max_{i,a}|(C_i\mathbf1)_a|\).
If every branch has both an incoming and an outgoing defect
edge, the maximum entry magnitude \(B\) satisfies
\[
 56B\ge 3\cdot54^2+52R .
\]
If the extreme is positive, use its outgoing edge: the row
sum of the scaled block at the extreme is
\(3\cdot54^2+54R-2y_a\) with \(y_a\le R\). If negative, use
an incoming edge and the column bound. Thus \(R\ge5\) forces
at least 162. The Paley-59 deletion above has minimum
incoming and outgoing defect degree five, and all Paley
deletions are equivalent by affine coordinates up to
converse. If a Paley deletion ever attains 160, every branch
internal outdegree must therefore be between 26 and 29.
A regular quotient has the same restriction unless the
mean already gives 262.

A different bound uses a single entry. Across a defect
edge \(i\to j\), take \(a,b\) maximizing \(x_a\) and minimizing
\(y_b\). Since \(4q_{ab}\ge-4\) and the internal term is
at least \(-4\),
\[
  B\ge148+\max_a x_a-\min_b y_b.
\]
In particular with **transitive interiors on all branches**,
the maximum is at least 258 for *every* quotient: a
defect edge has extrema \(55\) and \(-55\), while if every
pair avoids magnitude three the mean already yields 262.
Thus taking the Hamiltonian-minimal transitive branch
interiors imposes a much larger signed-square defect than
the near-regular target. None of the bounds claims attainability
or rules out a Moore graph.

## 4. Transport groups in the known controls

Choose a reference branch and use its matchings to identify
all branch slot sets. The transport group on its slots is
generated by transports along the reference-to-branch-to-branch-to-reference
closed walks. The closure does not depend on which branch is used:
closed walks at a second base are carried to the first by
transport out and back along a direct matching. The triangles
through the reference generate all closed walks by insertion
of those path pairs.

The replay enumerates the permutation closure for all 30
root/reference choices of Petersen and all 350 choices of
Hoffman–Singleton. They are \(C_2\) on two slots and the full
\(S_6\) on six slots, respectively. These are generated transport
groups, **not** automorphism groups or edge stabilizers.
In particular the smaller examples do not justify requiring
a hypothetical 56-slot completion to use a regular group.
The [regular-group obstruction of Smith and Montemanni](https://www.mdpi.com/2227-7390/14/23/4318)
concerns that narrower construction.

Run `node replay.js` for the exact controls in
[replay.out](replay.out).

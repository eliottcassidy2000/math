# An elementary order-seven cutoff, without the Busch floor

The [rooted emptying proof](EMPTYING_EXCLUSIONS.md) checks orders at
most six. For all larger orders, its invocation of Busch can instead
be replaced by the following simpler inequality sufficient for 21:
every strong tournament on at least seven vertices has at least 23
Hamilton paths.

## Order seven

If every indegree and outdegree is at least two, the seven outdegrees
are in {2,3,4} and sum to 21. At most three are 4, and equally many are
2. The number of directed triangles is C(7,3) minus the sum of
C(outdegree,2), hence at least 35−24=11. The odd-cycle packing formula
alone gives H≥1+2·11=23.

Otherwise, taking the converse if needed, empty an outdegree-one
vertex. Its unique successor u is in the first strongly connected
component of the remaining tournament W on six vertices. Add back
the outlet with the known rule H_new=H(W)+2E_u(W).

If W is strong, the order-six deletion table and the balanced case of
the emptying proof give H(W)≥15. Also E_u(W)≥4, by the following short
argument. A Hamilton path starting at u gives contributions of lengths
2,4,6. (For existence, take a longest path starting at u: a vertex
outside it that lost to any path vertex could be inserted after some
path vertex or appended. Thus every omitted vertex would dominate
the whole path, contradicting strongness.) If u has two successors,
there is already an extra length-two contribution. If it has only
one successor w, inspect the complement of {u,w}. A triangle there
gives weight at least three for the length-two term. Otherwise the
complement is transitive; w must beat its first vertex a, since
neither u nor any other complementary vertex can reach a directly.
The paths starting u,w,a and continuing with any odd subset of the
three remaining vertices already give four distinct even terms.
Thus in the strong-host case H_new≥15+2·4=23.

If W is decomposable, its first component has size 1,3,4 or 5.
Write it as A→B and use
`H_new=H(A)H(B)+S_u(A)Q_2(B)`.
In size 1, the guaranteed 2^5 ordered subsets of B give H_new≥33.
In size 3, S_u(C3)=3 and |B|=3 give H_new≥3+3·8=27.
In size 4, H(A)=5 and S_u(A)≥5 by its four-root data, hence
H_new≥5+5·4=25. In size 5, H(A)≥9 and S_u(A)≥7 by the order-five
rooted table (indeed its minimum is 9), giving H_new≥9+7·2=23.
All lower bounds suffice.

## Move the cutoff to every higher order

An elementary extension gives strong induced sets of every order
from three to the ambient order. Start with a directed triangle.
Given a proper strong induced S, an outside vertex with both in-
and outneighbors in S extends S by one. Otherwise each outside
vertex is in D (dominates S) or L (loses to S). Strongness of the
ambient tournament forces both groups and an arc x→y from L to D:
if there were no such arc, D could not be entered from S∪L.
Now discard any one vertex of S and add x and y. The resulting
set is strong even if the remaining vertices of S no longer
induce a strong tournament, because every ordered pair of them
can be routed through x→y. Its size exceeds |S| by one.
Thus any ambient strong tournament of order at least seven
contains a strong induced seven.

Hamilton path count is nondecreasing under adding a vertex: every
old Hamilton path has at least one valid insertion of the extra
vertex somewhere, including endpoints, and removing it recovers
the old path. Thus the strong-seven bound 23 extends to all higher
orders. Together with the rooted table and multiplication across
strong components this proves the all-order exclusions without
the precise Moon–Busch minimum formula.

The same extension supplies the full Hamilton-cycle fact without a
separate citation, if needed: a longest cycle would admit a mixed
outside vertex inserted in one edge; if all outside vertices were
uniform, the loser→winner pair could be inserted in any edge.
Opening a full cycle just before the desired root supplies its
rooted Hamilton path. The longest-rooted-path argument above
avoids even this separate intermediary.

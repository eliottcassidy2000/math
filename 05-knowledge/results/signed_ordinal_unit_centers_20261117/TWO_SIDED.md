# A two-sided entropy bound for every unit triangle switch

Use the two-port triangle switch of the [construction](../three_visit_triangle_switch_20261117/README.md):
triangle b→c→y→b and port incidences x→b, a→b, c→x, a→c,
x→y, y→a. The two variants differ only on x–a. The quotient Q
is transitive, so put the triangle between ordered tapes L and R.
Every L vertex beats the triangle; every R vertex loses to the
triangle; all L→R. Each tape vertex sees both ports the same way,
but the contacts can vary arbitrarily on both tapes.

The number of distinct positive counts ≤ X across both variants,
all splits and all binary words is
$$
 O(X^{\log_3 2}\log(2X)).
$$
Thus every such pair has exact gap two, yet the entire family
misses infinitely many primes and cannot place the center anywhere
on demand. This covers every unit-gap case of the stated
triangle construction, not other two-port gadgets. Indeed its
gap is (2H(Q)), and (H(Q)=1) iff Q is transitive: a directed
triangle already contributes two beyond the empty pack.

## Trim to active tapes

Call a left vertex a gate when the ports beat it, and a right
vertex a gate when it beats the ports. These are returns in the
opposite direction from the tape order. Before the first left
gate, all left vertices are successive global sources; after the
last right gate, all right vertices are successive global sinks.
Delete those redundant vertices without changing H. Hence each
retained left word starts at a gate; each retained right word
ends at a gate. Reversing the indexing on L but not on R makes
both normalized words end at a gate for bookkeeping.

If there are at least two gates on a tape, write i and i+r
for the positions of its last two gates in this normalized order,
with i≥1 and r≥1. The intervening r−1 positions are nongates.
There are at most (2^{i-1}) words for each (i,r).
Assign the tape to three increasing runs. On R the last vertices
of the first two runs are the gates at i and i+r, respectively;
the third run may be empty. Then positions before i have three
choices and positions between i and i+r have two (second or
third run), giving (3^{i-1}2^{r-1}) fillings. On L reverse the
construction: runs 1 and 2 (of runs 0,1,2) start at the actual
first two gates, and run 0 may be empty. The same count applies.

## Physical traces

If both tapes have at least two gates, use the same trace
in both branches:
$$
 L_0\,b\,R_0\,x\,L_1\,c\,R_1\,a\,L_2\,y\,R_2 .
$$
The internal runs L_1,L_2,R_0,R_1 are nonempty by the gate
assignments; the outer runs L_0,R_2 can be empty. Every tape-to-tape
join is explained by an intervening port or triangle vertex.
The joins R_0→x, R_1→a, x→L_1, a→L_2 are guaranteed by the gates,
and all remaining joins across triangle/tape follow by the quotient
order. No port edge is used. The independent fillings of the two
tapes yield at least (3^{i+j-2}2^{r+s-2}) distinct physical paths.

The one-tape growth bounds are valid with the other tape still
present and arbitrary. Use the following traces, where L or R
without a subscript is the entire other tape, placed at the outside.
For two right gates use
$$
 L\,c\,R_0\,x\,y\,R_1\,a\,b\,R_2,
$$
and for two left gates use
$$
 L_0\,c\,x\,L_1\,y\,a\,L_2\,b\,R .
$$
The successive fixed-vertex arcs x→y, a→b, c→x and y→a
used there are present. These traces give (3^{i-1}2^{r-1})
paths on the high-gate tape. If a tape has exactly one gate,
its normalized shape is S^tN, where S means nongate.
On R use
$$
 L\,c\,x\,y\,R_0\,a\,b\,R_1,
$$
with last vertex of R_0 the gate; the t preceding vertices are
free for either R run. On L use
$$
 L_0\,c\,x\,L_1\,y\,a\,b\,R
$$
with first vertex of L_1 the gate; the t succeeding vertices
are free for either L run. Thus H ≥ 2^t, independently
of what the other tape does. A tape with zero gates trims
to empty. Again no displayed trace uses the switched arc.

## Summing the addresses

Set α=log_3 2. In the both-high-gate case, at fixed (r,s) the bound
 3^(i+j−2)·2^(r+s−2) ≤ X gives
$$
 \sum_{i,j}2^{i+j-2}
 \ll (1+\log X)(X/2^{r+s-2})^\alpha
$$
over those (i,j): at total t=i+j−2 there are t+1 pairs
and their combined count is (t+1)2^t. Sum geometrically
over r,s. If exactly one tape has at least two gates, its
(3^{i-1}2^{r-1}) bound gives (O(X^\alpha)) addresses;
the other tape has only (O(\log(2X))) possible lengths
by its independent one-gate bound and has exactly one word
per length. If neither does, there are (O(\log^2(2X)))
pairs of lengths. Multiplying by two for the variants concludes.

The mechanism is the joint use of the same two ports: each return
run can cross from L through the triangle to R, so three-way freedom
on the two sides multiplies. Moving the triangle among a transitive
quotient cannot avoid this entropy cost. Literal DP and trace
checks on short words are in verify_two_sided.py.

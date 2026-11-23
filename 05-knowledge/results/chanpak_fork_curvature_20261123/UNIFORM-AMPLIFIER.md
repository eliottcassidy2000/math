# The exact gap for uniform tournament attachments

Fix a nonempty tournament \(M\) and let \(D=H(M)\). In a uniform
extension \(G\), each exterior vertex either beats all of \(M\) or
loses to all of it. Let \(C=G/M\) contract \(M\) to the marked point
\(m\). Let \(P_k(M)\) count ordered covers by \(k\) nonempty directed
paths, and set \(Q_2(M)=2D+P_2(M)\).

For any fixed exterior, a Hamilton path's maximal module segments
give an exact expansion
\[
 H(G)=\sum_{k=1}^{|M|} b_k P_k(M),
\]
where \(b_k\) counts exterior skeletons with \(k\) ordered module
positions; necessarily \(b_1=H(C)\). Every gap between these positions
is a nonempty exterior path from a vertex beaten by \(m\) to one
beating \(m\).

## The reentry dichotomy

If \(m\) is a singleton SCC of \(C\), no exterior gap is possible:
with \(m\) such a gap would be a directed cycle. Then the formula
reduces to \(H(G)=H(C)D\).

Conversely, suppose \(m\) lies in a nonsingleton SCC and \(|M|\ge2\).
There is a directed simple cycle through \(m\). Replace its start
and end by two distinct copies \(m_1,m_2\), with identical exterior
incidences, to make a directed path. Insert each unused exterior
vertex into this path in turn. The elementary insertion rule
(prepend, append, or place at the first switch from a path vertex
beating the new vertex to one it beats) preserves the relative
order of the original vertices. Thus a nonempty exterior gap persists
between the two copies. This gives a valid two-position skeleton,
so \(b_2\ge1\). Further \(b_1=H(C)\ge3\), since \(C\) is not
transitive and its Hamilton count is odd. Therefore
\[
 H(G)\ge3D+P_2(M)=D+Q_2(M).
\]
Equality occurs for the cyclic sandwich: add two vertices \(b,c\)
with \(M\to b\to c\to M\). Its contracted triangle has three
single-position orders, and its only two-position skeleton is
\(m_1,b,c,m_2\). No third position is possible. The sandwich is
strong even if \(M\) is not.

In particular this is the exact minimum for any proper strong
uniform extension of \(M\) of order at least two. More generally if
a \(k\)-position skeleton occurs, preserve any \(j-1\) of its gaps,
insert all remaining exterior vertices between or around the
retained skeleton, and permute the resulting \(j-1\) bridges.
This gives \(b_j\ge(j-1)!\) for every \(2\le j\le k\). Cutting
a Hamilton path at \(j-1\) positions yields
\(P_j(M)\ge\binom{|M|-1}{j-1}D\), giving corresponding multi-visit
lower bounds.

## Seven is not a uniform amplification ratio

If there is no reentry, ratio 7 would require the impossible
\(H(C)=7\). If there is reentry and \(n=|M|\ge6\), cutting every
full Hamilton path once gives \(P_2(M)\ge(n-1)D\), and the ratio
is at least \(n+2\ge8\).

For a strong five-vertex \(M\), there is a directed four-cycle.
For example, take a directed Hamilton cycle numbered modulo five.
A forward distance-two chord skips a vertex; if all these chords
are backward, use \(0\to3\to1\to4\to0\). Let \(v\) be the
omitted vertex and \(u\) a cycle vertex. If \(u\to v\), put \(v\)
in the first slot and a rotation starting at \(u\) in the second;
otherwise use a rotation ending at \(u\) and put \(v\) last.
This is a two-slot cover that fails to concatenate, so is counted
in \(P_2\) in addition to the \(4D\) valid cuts. In the nonstrong
case \(Q_2/H\) multiplies across SCCs, because the vertices along
each slot occur in SCC order. Each SCC ratio is at least its order
plus one by cutting full paths (including empty ends); the product
is strictly greater than six when total order is five and there are
at least two SCCs. Thus all order-five modules satisfy \(P_2>4D\).

The three order-four profiles \((D,Q_2)\) are \((1,16),(3,24),(5,32)\),
all with \(P_2>4D\). In orders one through three, \(D\) is 1 or 3,
and ratio 7 would produce one of the known holes 7 or 21.

This rules out uniform attachment of any size as a physical
multiplier by 7. It does not apply to the two-read fork's cofactor
\(40Q_2+D\), which need not equal \(D\), nor to nonuniform attachments
or changing arcs of the input module.

## Several uniform modules

In a strong quotient \(C\) on at least three vertices, simultaneously
inflate its vertices to disjoint nonempty modules \(M_i\). Write
\(D_i=H(M_i)\), \(n_i=|M_i|\). Count all skeletons meeting every
module once; then for each \(n_i\ge2\) count one fixed skeleton
meeting only \(M_i\) twice. These families are disjoint, yielding
\[
 H(G)\ge
 \left(H(C)+\sum_{n_i\ge2}\frac{P_2(M_i)}{D_i}\right)\prod_i D_i
 \ge \left(3+\sum_i(n_i-1)\right)\prod_i D_i.
\]
For instance, two or more actual prime-count modules of height
\(p\ge5\) cannot be packaged into one strong host of height
\(7p^r\) when they are the \(r\) marked modules; each has order
at least four, so the lower ratio is at least nine.

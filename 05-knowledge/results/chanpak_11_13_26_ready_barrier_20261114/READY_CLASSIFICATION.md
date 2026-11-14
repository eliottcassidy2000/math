# Every ready root has one outgoing arc and a nonterminal predecessor

A ready pair \((a,b)\), with \(b>0\), means that the full
root signature of a physical tournament \((T,v)\) is
\[
                    (H,R_+,R_-)=(a+b,a+b,-a).
\]
We now characterize all such graphs, not just those generated
by the previous bridge. In particular **every** physical
ready pair satisfies \(b\ge2a\), so the terminal-seven
counterexample cannot be avoided by a ready graph outside
the arithmetic cone.

## 1. The signed subset-path profile is a literal endpoint count

For every tournament \(Q\) and root \(u\),
\[
                  R_-(Q,u)=-\,\#\{\text{Hamiltonian paths of }Q
                                             \text{ ending at }u\}.
\]
Here is an algebraic proof with no enumerative assumptions.
Work with commuting square-zero variables \(x_i\), put
\(X=\operatorname{diag}(x_i)\), and let \(A\) be the directed
adjacency matrix. The multiaffine path series (including the
empty path) is
\[
 F_A=1+\mathbf1^tX(I-AX)^{-1}\mathbf1 .
\]
Repeated vertices vanish because \(x_i^2=0\).
The corresponding root-start series evaluated at \(-x\) is
\(G_u(-x)=-x_u e_u^t(I+AX)^{-1}\mathbf1\).
The full-set coefficient of \(G_u(-x)F_A\) is the signed
weighted complementary-subset profile by definition.

We may replace \(A\) by \(A+I\) inside this root-start
expression: a diagonal edge repeats the already weighted
vertex and thus gives zero. Since \(A+A^t=J-I\), let
\(K=I-A^t X\); then \(I+(A+I)X=K+\mathbf1\mathbf1^t X\).
The rank-one inverse formula over this ring gives
\[
 (K+\mathbf1\mathbf1^t X)^{-1}\mathbf1
   =\frac{K^{-1}\mathbf1}{1+\mathbf1^t X K^{-1}\mathbf1}.
\]
The denominator is \(F_{A^t}=F_A\) by path reversal on each
subset, and is invertible because its constant term is 1.
So \(G_u(-x)F_A=-x_u e_u^tK^{-1}\mathbf1\).
Its full-set coefficient counts negatively the paths in the
converse starting at \(u\), as asserted.

Consequently \(R_-(Q,u)=0\) exactly when \(u\) is **not** in
the terminal strong component: all paths visit the strong
components in their unique ordinal order, and every specified
vertex of a strong component can be a path endpoint by rotating
a Hamiltonian cycle (or trivially for a singleton).

## 2. The unsigned equality forces one outgoing arc

The \(R_+\) profile counts complementary ordered path pairs
\((P,S)\) such that \(S\) starts at the root \(v\), with
\(P\) allowed to be empty. Concatenation is a full path unless
\(P\ne\varnothing\) and \(v\) beats its last vertex. Every
full path is obtained exactly once by splitting at \(v\).
Thus \(R_+\ge H\), and equality forbids every nonjoinable pair.

Suppose \(v\) has distinct outneighbors \(x,y\). Work in the
suffix \(W\) of vertices reachable from \(v\), including \(v\).
Let \(X_x\) be the vertices other than \(x\) not reachable
after deleting \(x\) from this suffix; define \(X_y\) similarly.
The two sets are disjoint: a path to a vertex in both would
visit both \(x,y\), but \(v\) can shortcut to the later one.
If both sets were nonempty, each arc between a point of \(X_x\)
and a point of \(X_y\) would have to point both ways. Indeed
each point is still reachable after deleting the other set's
named neighbor, and even one arc in the wrong direction restores
the forbidden reachability. Therefore one of these deletions,
say \(x\), leaves all of \(W-x\) reachable from \(v\). The first
strong component of \(W-x\) contains \(v\); rotate its Hamiltonian
cycle, then follow the other components, to obtain a spanning
path starting at \(v\). Take the complementary path through
the earlier components of \(T\), followed by \(x\).
These are complementary paths that cannot join, since \(v\to x\).
Thus \(R_+=H\) forces outdegree at most one. In a ready
signature with \(b>0\), it must be exactly one.

## 3. Complete structural grammar and the cone

Let the unique outneighbor of \(v\) be \(u\), and put \(Q=T-v\).
The zero-bridge identity, or the split definition directly, gives
\[
 H(Q)=a,\qquad R_+(Q,u)=b,\qquad R_-(Q,u)=0.
\]
For the last two, every nontrivial path starting at \(v\)
immediately uses \(u\); the singleton contributes \(H(Q)\) to
\(R_+\) and \(-H(Q)\) to \(R_-\).
By the signed identity \(Q=A_0\to B\to C_0\) with \(u\)
in a marked strong component \(B\) and **nonempty** later
part \(C_0\) (the earlier part may be empty). Exactly
\[
 a=H(A_0)H(B)H(C_0),\qquad
 b=H(A_0)R_+(B,u)F_2(C_0).
\]
The vertex-rooted paths cannot visit earlier components and
split the later components freely across their two color slots.
Since \(R_+(B,u)\ge H(B)\) and \(F_2(C_0)\ge2H(C_0)\),
the claimed cone \(b\ge2a\) follows. In fact \(b\) is even,
since swapping the two ordered slots of a nonempty vertex
set pairs the \(F_2\) summands.

Conversely, take any such \(A_0,B,C_0\), and add \(v\) beating
only \(u\). Its profiles satisfy the ready signature
\((a+b,a+b,-a)\), by the same zero-bridge law and \(R_-=0\)
for \(Q,u\). This is a characterization of physical ready
graphs, not a claim that every cone point has those factors.
The \(A_0\) factor disappears for primitive pairs (earlier
singleton padding may simply be deleted). Every later strong
component incurs at least a factor 2 in \(b/a\) and in \(b\).
Unlike the earlier positive-bridge compiler, the later part
need not be transitive or have any singleton followers.

For example take \(Q\) to be a singleton \(u\) followed by the
strong four-vertex core \((8,1,3,6)\), then add \(v\) beating
only \(u\). This gives the new ready pair \((5,32)\).
No last **positive single-contact** bridge realizes this pair:
it would require a rooted \(H=5\) source with \(R_+\) a power
of two below 32. Its only nontrivial strong component is the
unique strong four-vertex \(H=5\) class (Busch excludes larger
ones). Its root menu is \(5,6,9,7\); an exterior root before
this component has \(R_+\ge32\), one after it is a multiple
of 5, and trailing singleton factors cannot erase an odd
factor. So no smaller binary-rooted source supplies it.
The terminal length-seven gate gives the exact total
\(7\cdot2027=14189\), with prime cofactor 2027.

## 4. The terminal-seven obstruction is unconditional within ready geometry

For target \(7p=180131\), every physical ready prefix lies
in the above cone. For any terminal length-seven width power
\(y\), its ordinary count is
\[
                   (57+68y)a+(68+81y)b.
\]
Thus \(193+230y\le180131\), leaving powers \(2,\ldots,512\).
The exact arithmetic replay finds only \(y=4,(a,b)=(59,410)\)
in the full cone. The companion obstruction note verifies
that no rooted strong \(H=59\) core has \(R_+=205\) or 410,
and excludes the other placements of the root in \(Q\).
So no arbitrary ready source at **any** terminal-seven
positive width gives this total. Extra earlier corridor
steps do not help: they would still leave a ready source.

This is not a new missing count. The 15-vertex two-hub word
`1002303022110` gives exactly 180131 by independent literal
Hamiltonian path DP. The family fails, not the target.

Run `python ready_class.py` for the signed identity and the
new ready witness, and `python alternate_target.py` for an
independent count-180131 witness.

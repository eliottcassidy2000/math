# The cost of making a product of induced blocks strong

Let $V(T)=\bigsqcup_i B_i$ be any partition into nonempty blocks. A
cycle is mixed if it meets more than one block. The odd-cycle packing
polynomial gives the exact identity
$$
 Q_T(z)-\prod_i Q_{T[B_i]}(z)
 =\sum_{\substack{\mathcal M\ne\emptyset\\
            \mathcal M\text{ a packing of mixed odd cycles}}}
 (2z)^{|\mathcal M|}
 \prod_i Q_{T[B_i\setminus V(\mathcal M)]}(z).
$$
Indeed every packing chooses its mixed cycles uniquely; the remaining
packings are entirely internal and independent after those vertices
are removed. In particular the difference is coefficientwise
nonnegative. At $z=1$ equality holds precisely when every strongly
connected component of $T$ is contained within one block.

For the last characterization, if no mixed directed triangle exists,
each strong component $C$ of a single induced block has uniform
incidence to every vertex outside that block. Otherwise an outside
vertex $x$ beats some but loses to other vertices of $C$.
Strongness supplies an arc from the part beaten by $x$ to the part
that beats $x$, closing a mixed triangle. Components in the same
block already have uniform incidence. Thus contract all the blockwise
components. A directed triangle in the quotient could not be contained
in a single block, and a mixed one would lift to a mixed triangle.
The quotient is therefore transitive. It follows that all strong
components of the whole tournament are internal. The converse is
immediate since any directed cycle lies in one strong component.

One quantitative corollary (without counting any packings with more
than one mixed cycle) is
$$
 H(T)\ \ge\ \prod_i H(T[B_i])
     + 2\sum_{C\text{ mixed directed triangle}}
            \prod_i H(T[B_i\setminus V(C)]).
$$
Each summand is the entire packing channel with $C$ the only mixed
cycle, so there is no double counting. This is why cross-arc edits
alone cannot turn the ordinal $C_3\to T\to C_3^*$ sandwich into a
strong tournament while retaining its exact count $9H(T)$, even if
no symmetry is required. A repair must also remove some internal
packing contribution or change the factor construction.

For the specific [cyclic wrap](CYCLIC_WRAP.md) a simpler independent
obstruction suffices: $Q_T(2)\ge(n+1)H(T)$ for $n=|T|\ge1$,
because each full path gives $n-1$ distinct ordered nonempty
two-covers by cutting at an edge and the two empty-slot choices
contribute $2H(T)$. Thus $H(W(T))\ge(n+2)H(T)$; for $n\ge8$
this rules out an exact factor nine even without examining the ports.

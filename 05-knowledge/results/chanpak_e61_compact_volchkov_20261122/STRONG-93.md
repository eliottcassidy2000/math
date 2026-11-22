# A new strong lower class family at exponent 61, and its whole one-flip obstruction

**PROVED.** There are at least **237,918,897,216,000 pairwise nonisomorphic strong tournaments on 93 vertices** with \(E_{61}-1\) Hamilton paths, in the following explicit construction. None of them has a single-arc neighbor with \(E_{61}+1\). These are new lower witnesses outside the known order-124 three-port family. Their large multiplicity does not solve the stronger seam problem.

## Integer exact cover by ten small components

For a tournament \(C\), set \(h(C)=H(C)\), and \(q(C)=\sum_{S\subseteq C}H(C[S])H(C-S)\), the count of ordered two-path covers with empty slots allowed. Both reads multiply in ordinal sums: each path crosses component cuts in order, so split both ordered paths at the cut, and glue the two slot indices component by component. Use these strong factors, with outgoing masks indexed from zero:

| Factor rows | \(h\) | \(q\) |
|---|---:|---:|
| 256,289,67,23,295,141,59,95,236 | 1321 | 24800 |
| 6,72,42,49,71,83,141,63 | 331 | 5368 |
| 768,257,3,7,815,399,831,863,524,46 | 451 | 21136 |
| 256,257,3,23,39,399,319,95,156 | 305 | 9268 |
| 2,4,1 (two copies) | 3 | 12 |
| 2,8,3,5 | 5 | 32 |
| 256,257,515,263,783,31,831,127,164,427 | 1057 | 38656 |
| 8,5,17,6,11 | 13 | 104 |
| 4,33,18,7,11,29 | 31 | 328 |

Add 24 singleton factors \((h,q)=(1,2)\). Across the 34 ordered components
\[
         h(U)=2^{60}-1,\qquad q(U)=2^{61}(2^{60}-1).
\]
This is a paired multiplicative selector, not a numerical approximation; the 24 singleton factors supply the remaining power of two once the two odd-part products match. The factors total 91 vertices. Form \(T\) with cyclic structure \(U\to a\to b\to U\), where \(a,b\) are two singleton vertices. The completion is strong regardless of the order of the components of \(U\). Paths with \(a\) before \(b\) traverse \(a,b\) consecutively and read an arbitrary two-slot cover of \(U\). Those with \(b\) before \(a\) start at \(b\), end at \(a\), and read an ordinary spanning path in \(U\). Hence \(H(T)=h(U)+q(U)=(2^{60}-1)(2^{61}+1)=E_{61}-1\).

The set \(U\) is the unique maximal nontrivial proper module (a module is seen uniformly by each exterior vertex). A module including a singleton and any vertex of \(U\) must contain the other singleton because it distinguishes the two. A module containing both singletons must contain every vertex of \(U\), since every \(u\) distinguishes them. Thus an isomorphism between wrappers fixes \(U\) setwise and preserves its intrinsic **ordered** strongly connected component types. The ten non-singleton factors are all strong and pairwise nonisomorphic except for the two triangles, as is already witnessed by their distinct \(h\) otherwise. Therefore the number of classes in this particular reorderable family is exactly \(34!/(24!\,2!)\).

## Class-wide single-arc exclusion, without scanning all orders

We use \(q(C)\ge1\), \(h(C)\ge1\) including empty \(C\), and \(q(C)\equiv0\bmod4\) for every order \(n\ge2\). For the latter pair complementary subsets: there are \(2^{n-1}\) pairs, each made of two equal odd products by Rédei. Also \(h(C)\le m!\) and \(q(C)\le(m+1)m!\) on \(m\) vertices, by ordering each chosen subset and its complement freely. Set \(L=11\cdot10!\), a common bound for the components in the table, and \(R=q(U)/h(U)=2^{61}>6L^2\).

**A marked edge inside one component \(C\subset U\).** All ordinal cuts around \(C\) still exist after the flip. Put \(U_{\rm out}\) for the other factors, so the wrapper changes by \(h_{\rm out}\Delta h(C)+q_{\rm out}\Delta q(C)\). Because \(q_{\rm out}/h_{\rm out}=R\,h(C)/q(C)\ge R/L\), any nonzero \(\Delta q(C)\) dwarfs \(|\Delta h(C)|\le L\): the magnitude is at least \(h_{\rm out}(R/L-L)>2\). If \(\Delta q(C)=0\), the jump is a multiple of \(2h_{\rm out}\), and \(h_{\rm out}>1\) since several other nontrivial components remain. A singleton contains no internal edge.

**A marked edge spanning two different components of \(U\).** The old arc belongs to no directed cycle. On every induced subset containing the two endpoints the same is true. Reversing adds odd cycles but removes none, so the [odd-cycle packing formula](https://arxiv.org/abs/2307.05569) makes all induced counts nondecreasing. Thus \(q(U)\) cannot decrease, and its increase is at least \(2\Delta h(U)\) from the whole-set and empty-set terms. If \(\Delta h(U)>0\), the wrapper jumps by at least six because \(\Delta h(U)\) is even. Otherwise its jump is a multiple of four. Neither is two.

**The arc \(a\to b\).** Its reversal makes the three-block structure ordinal, with path count \(h(U)\), so it cannot move up by two.

**A marked contact with \(a\) or \(b\).** Write \(S_u,E_u\) for the full paths in \(U\) starting or ending at \(u\), and \(s_u,d_u\) for the pointed weak two-slot covers with one designated slot starting or ending there. On reversing \(u\to a\), precisely \(E_u+d_u\) old paths disappear: the full \(U\) path ending at \(u\) if \(b\) came first, or the first slot ending at \(u\) before consecutive \(a,b\). The gained paths use \(a\to u\); \(b\) must be first because its only in-arc is \(a\to b\), which cannot be used when \(a\) is immediately followed by \(u\). The prefix of the path between \(b\) and \(a\) is a nonempty \(U\)-path, so \(s_u-S_u\) such new paths occur. Therefore
\[
  \Delta_{u,a}=s_u-d_u-S_u-E_u,\qquad
  \Delta_{b,u}=d_u-s_u-S_u-E_u
\]
by converse symmetry for the second case. If \(U=A\to C\to B\), \(u\in C\), slot splitting gives
\[
 \begin{aligned}
 s_u&=h(A)q(B)s^C_u,&d_u&=q(A)h(B)d^C_u,\\
 S_u&=[A=\varnothing]h(B)S^C_u,&
 E_u&=[B=\varnothing]h(A)E^C_u.
 \end{aligned}
\]
If both sides have at least two vertices the full-path terms vanish and both other terms are divisible by four, excluding a jump of two. Otherwise assume by converse that \(|A|\le1\); \(B\) is nonempty because components have size at most ten. Now \(s^C_u\ge1\) (put \(u\) alone in the designated slot) and \(d^C_u,S^C_u\le L\). The non-dominant terms have combined magnitude at most \(2L\,h(B)\): \(q(A)=2\) for a singleton, in which case \(S_u=0\), or \(q(A)=1\) otherwise. But
\[
      q(B)/h(B)\ge R/(2L)>3L.
\]
Both signed differences then have magnitude greater than two. This argument is independent of the component permutation, and exhausts all arcs.

See the [subset-DP replay](check_strong_93.py) for the full factor table, integer products, class count, and literal small controls of the contact and wrapper identities. The construction tests a genuine alternative class grammar, not merely new labels in the known neutral component.

## All unmarked slot profiles collide, not just the two scalar factors

This is a concrete extension of the [pointed-monoid slot geometry](../session320_second_moments_and_native_contexts_20261113/F1_SLOT_GEOMETRY.md), not just its name. For any number \(j\) of ordered weak slots, \(Z_U(j)\) multiplies over the same ordered components, because restriction and concatenation at each ordinal cut are inverse in every slot. Therefore every reordering has the entire same weak profile. For any fixed quotient with one occurrence of the module \(U\), collapse maximal runs of \(U\) in the \(j\) exterior slots. The number of fillings of \(r\) prescribed nonempty runs is the finite difference
\[
        B_r(U)=\sum_{s=0}^r(-1)^{r-s}\binom{r}{s}Z_U(s).
\]
Thus \(Z_T(j)\) also agrees for all wraps in the family, whatever \(j\) is. This realizes an all-slot collision inside one exact strong count fiber of astronomically many classes.

In the inherited pointed monoid \(C_m\cup\{0\}\), each \(r\)-occupied-slot cover on \(n\) vertices has \(m^{n-r}\) relative unit phases. Setting the unit group to the trivial group keeps every cover; evaluating its finite-field polynomial in \(q-1\) at \(q=1\) generally does not. The collision above therefore also holds for the unmarked ordered slot-phase counting polynomials at every unit order. It does not hold as a claim about all directed marked incidences: one two-vertex slot and \(n-2\) singleton slots already give a read that recovers the arc. In an absolute seam inverse those marked reads, rather than more evaluations of the same unmarked slot profile, are the information to add.

## A second cyclic exact cover on 91 vertices

The same paired multiplicative search gives six distinct strong factors with triples \((h,\operatorname{oddpart}(q),v_2(q))\)
\[
 (3963,5889,4),\ (331,671,3),\ (1661,6355,3),\
 (1525,6605,2),\ (1599,6951,2),\ (217,1,12).
\]
Their outgoing masks are embedded in the checker. Add 35 singleton factors to make 89 vertices in \(U\), then add the two exterior vertices of the directed triangle. It has the same two target products and gives \(41!/35!=3,237,399,360\) pairwise nonisomorphic strong 91-vertex lower classes. The class-wide no-flip argument applies without change: all components have order at most ten, the same ratio \(R\), and plenty of nontrivial exterior factors remain when one component is marked.

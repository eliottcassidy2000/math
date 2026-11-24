# Fixed substitution templates cannot prove Chan–Pak by count coverage

**Proved, independent of Moore graph existence.** Fix any finite collection of tournament templates. Start with a singleton, and recursively substitute previously constructed nonempty tournaments into their vertices: all arcs between two blocks take the direction of the corresponding template arc. Allow arbitrary ordered joins, arbitrary fixed seeds (absorbed into the templates), and arbitrary recursion depth. Let \(\mathcal H\) be the set of resulting Hamilton-path counts. For every \(\epsilon>0\),
\[
  |\mathcal H\cap[1,X]|=O_{\epsilon,\mathrm{templates}}(X^\epsilon).
\]
Equivalently, the Dirichlet series of **distinct** realized heights converges at every real \(\epsilon>0\). In particular no fixed upper bound on the sizes of recursively substituted quotients is enough for a full-spectrum Chan–Pak construction. Even one degree-57 Moore tournament as a fixed template, all its fixed-order representations, and all their fixed-order edits together do not change this conclusion.

## Profiles and what substitution preserves

Let \(H(T)\) count directed Hamilton paths, with \(H(\varnothing)=1\). For \(k\ge1\), define
\[
 Z_k(T)=\sum_{f:V(T)\to[k]}\prod_{j=1}^k H(T[f^{-1}(j)]) .
\]
This counts ordered \(k\)-tuples of paths (empty allowed) that partition the vertices. Ordered joins multiply these numbers at every \(k\): each colored path must run through its share of the left block before its share of the right block. This profile is the principal specialization of the [Rédei–Berge function of Grinberg and Stanley](https://arxiv.org/abs/2307.05569). Their odd-cycle packing formula for \(H\), applied independently in each color, gives a useful elementary polynomial form:
\[
 Z_k(T)=\sum_{\mathcal C} 2^{|\mathcal C|}
 k^{\,|T|-|V(\mathcal C)|+|\mathcal C|},
\]
where \(\mathcal C\) ranges over packings of disjoint directed odd cycles of length at least three. A chosen cycle must have one color; unoccupied vertices choose colors independently. For nonempty \(T\) the polynomial has no constant term, positive leading coefficient 1 in degree \(|T|\), and \(Z_k\ge k^{|T|}\).

Substitution acts on **full profiles**, not on the scalar \(H\). For each output color, compress its path into a sequence of child labels, merging each maximal interval from the same child. Consecutive labels must follow an arc of the template. If label \(i\) has \(r_i\) occurrences over these ordered sequences, the local choices are exactly the ordered nonempty \(r_i\)-path covers of the \(i\)-th child. These are determined by inclusion-exclusion from its \(Z_j\) for \(0\le j\le r_i\). Sum the products over possible sequences; \(r_i\le |T_i|\), so this is finite. In particular if two children have identical full profiles, replacing one by the other does not change the output profile. This is needed below to sum over distinct profiles rather than over derivations.

## The key growth inequality

For every nontrivial *strongly connected* template \(B\) of size \(r\), and integer \(k\ge1\), we have
\[
 Z_k(B[T_1,\ldots,T_r])\ \ge\ c(k,r)\prod_i Z_{k+1}(T_i),\qquad
 c(k,r)=\left(\frac{2}{3(k+1)^2}\right)^r .
\]
Pick one of the \(k\) output colors as special. In each child of size at least two, choose an ordered \((k+1)\)-path cover in which the two special pieces are nonempty. Each singleton child goes wholly to the special color. Every ordinary output color concatenates its pieces in a fixed Hamilton order of the occupied child labels, using Rédei's insertion theorem.

For the special output color, use the directed Hamilton cycle of \(B\) (Camion's theorem). If there are nonsingleton children, pick a fixed Hamilton path of the template induced on them. Start the full cycle at the first label of that path, travel it once as a linear path using the first special pieces and singletons, and take its closing arc back to the start to follow the nonsingleton-only path using all second special pieces. At the second visit it is a different physical vertex. If all children are singletons, only the first lap is needed. This construction is injective: child labels identify the two ordered visits assigned to each nonsingleton, and the local word at each visit and color is recovered. In particular, omissions of ordinary colors do not create any overcount.

Put \(z=k+1\). The number of allowed covers in a nonsingleton \(T\) is
\[
 Z_z(T)-2Z_{z-1}(T)+Z_{z-2}(T).
\]
For a monomial \(z^a\) of degree \(a\ge2\), this second difference is at least \(2z^{a-2}\): it counts maps to \(z\) colors using both designated colors, and we may assign the first two inputs to the designated colors in either order. The degree-one coefficient of the profile polynomial is at most \(2H(T)/|T|\), because it is twice the number of directed Hamilton cycles when \(|T|\) is odd and is zero otherwise. Cutting the cycles into paths gives that bound, with at most one cycle closing a given path. Since \(Z_z(T)\ge zH(T)\) by monochromatic covers, at least a third of \(Z_z\) comes from degrees at least two when \(|T|\ge3\). For \(|T|=2\) the degree-one coefficient is zero. Thus the second difference is at least \(2Z_z/(3z^2)\). A singleton contributes one choice instead of its factor \(z\); \(1/z\ge 2/(3z^2)\) suffices too. This proves the displayed inequality.

## Dirichlet descent

A nonstrong template decomposes into its linearly ordered strongly connected components. Substituting into it is equivalent to substituting into those smaller components and taking an ordered join. Enlarge the original finite family by their constituent strong templates and absorb one-vertex components as identity operations. This leaves only finitely many nontrivial strong templates and ordered join (which can be binary).

Fix \(\sigma>0\). Let \(F_k\) sum \(Z_k^{-\sigma}\) over distinct full profiles in the resulting family, including the empty profile only as a bookkeeping identity. At sufficiently large integer \(k\), this is finite: the number of expression trees of order \(n\) is at most \(C^n\) for some fixed \(C\) (all operations have fixed finite arity at least two and finite labels), and \(Z_k\ge k^n\). For any \(k\), call a root generator a singleton or an output of a nontrivial strong template; the inputs may themselves be arbitrarily complicated. Summed over distinct root-generator profiles, their weights are at most
\[
  1+\sum_{B} c(k,|B|)^{-\sigma}F_{k+1}^{|B|},
\]
whenever \(F_{k+1}\) is finite. Indeed each root generator comes from some tuple of child profiles, and such a tuple determines its output. The constant \(1\) overestimates the singleton term.

Every full profile is a product, pointwise in \(k\), of root-generator profiles, by flattening the top ordered joins. Multiplication commutes at the profile level. For \(k\ge2\) every nonempty generator has \(Z_k\ge2\). Thus an Euler product over generator profiles bounds \(F_k\) above by \(\prod_g(1-Z_k(g)^{-\sigma})^{-1}<\infty\): sum over multisets of generators, overcounting coincident products. This descends to \(k=2\). At \(k=1\) use the same generator estimate with \(F_2\); sum over distinct **generator heights**, not full profiles, and omit height 1 from the Euler product since it changes no output height. All nontrivial height factors are at least 3 by Rédei's parity theorem. The output height series is finite. Finally every term indexed by a height \(\le X\) has weight at least \(X^{-\sigma}\), proving the claim.

## Implication for a Moore-based compiler

The obstruction covers recursive substitution, not arbitrary target-guided arc edits or inter-block incidence rules with a growing quotient. It also does not contradict the repo's [order-seven *rooted bank* finite coverage](../../../01-canon/theorems/THM-4885-one-apex-join-representation-mass-pi2-over-12.md): a bound on some reusable cores is not a uniform bound on all quotient templates or on nonuniform contacts of the full outputs.

To use Moore geometry for Chan–Pak, the candidate ingredient should be a **selector**: consistent but nonuniform choices of inter-block contacts that can expand with the target, allowing unbounded quotient complexity. A single degree-57 graph offers a large new seed or a blueprint for gluing; it cannot be a sufficient seed for any grammar covered by this theorem, nor is it known to be necessary. The pentagon already gives exact +2 families.

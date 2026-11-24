# A count-pressure obstruction even for certain unbounded catalogs

This strengthens the [finite-template theorem](../moore_chanpak_counts_20261124/bounded-templates.md) without making any assumption about a degree-57 Moore graph.

Let \(\mathcal B_r\) be a catalog of strongly connected tournament templates of order \(r\), with finitely many at each \(r\). Fixed seeds and ordered joins are allowed. Define
\[
 W_s(r)=\sum_{B\in\mathcal B_r} H(B)^{-s}.
\]
**Theorem.** If for each \(s>0\) and \(L<\infty\),
\[
 \sum_{r\ge3} W_s(r)L^r<\infty,
\]
then arbitrary recursive substitution from this catalog yields \(O_\epsilon(X^\epsilon)\) distinct Hamilton counts up to \(X\), for every \(\epsilon>0\). The condition is equivalent to \(\limsup_{r\to\infty}W_s(r)^{1/r}=0\) for every \(s>0\).

This includes any catalog of at most exponentially many templates per order whose \(\log H(B)/r\) tends uniformly to infinity. It also allows **superexponentially many**: for example \(e^{o(r\log r)}\) templates per order and \(H(B)\ge r^{c r}\) for some uniform \(c>0\) suffice. Thus unbounded quotient order by itself is not enough. A fixed finite collection of nonstrong templates can be split into its strong components and joins as before.

**Weighted high-color initialization.** In addition to the transfer bound \(Z_k(B[T_i])\ge c_k^r\prod_i Z_{k+1}(T_i)\) in the finite proof, the following weaker but multiplicative inequality is useful:
\[
 H(B[T_i])\ge H(B)\prod_i H(T_i).
\]
Indeed independently pick a Hamilton word in each child and a skeleton Hamilton word, and concatenate the children as whole blocks. The output word uniquely recovers all choices. Applying this along a recursion tree gives \(H(T)\ge\prod_v H(B_v)\) over its internal nodes. Ordered joins have template height one. Since \(Z_k(T)\ge k^n\) and \(Z_k(T)\ge H(T)\),
\[
 Z_k(T)^{-\sigma}\le k^{-n\sigma/2}\prod_v H(B_v)^{-\sigma/2}.
\]
Here finite seeds may be regarded as template nodes with singleton children. For any fixed \(s=\sigma/2>0\), the hypothesis implies \(W_s(r)\le A_s^r\) for some finite \(A_s\ge1\); the binary join adds only one fixed term. For a plane tree with \(n\) singleton leaves and all internal arities at least two, their sum is at most \(2n-2\). Summing over template labels therefore costs at most \(A_s^{2n}\); summing over the at most exponentially many plane tree shapes still costs only \(C_s^n\). Choose \(k\) large for summability. The sum over distinct full profiles is bounded by the sum over expressions.

**Descending in the color parameter.** For \(0<\delta<\sigma\) and a strong root template \(B\), combine \(Z_k\ge H(B)\) with the color transfer bound:
\[
 Z_k(B[T_i])^{-\sigma}\le
 H(B)^{-\delta}\,c_k^{-r(\sigma-\delta)}
 \prod_i Z_{k+1}(T_i)^{-(\sigma-\delta)}.
\]
Thus finiteness of the profile sum at \((k+1,\sigma-\delta)\) implies a finite weighted root-generator sum at \((k,\sigma)\) by the catalog hypothesis, since the child sum and \(c_k^{-(\sigma-\delta)}\) enter only as a finite exponential base. Full profile substitution is single valued, as proved in the finite case. An Euler product over generator profiles controls ordered joins for \(k\ge2\); at \(k=1\) take it over distinct generator heights greater than one. For a desired exponent \(\epsilon\), initialize at a sufficiently large \(k_0\) with \(\epsilon/2\), then give each descent the same increment \(\delta=\epsilon/(2(k_0-1))\). This yields convergence at \((1,\epsilon)\).

**Consequence for a constructive search.** Under an exponential cap on catalog size, a grammar that could cover all odd counts must, at a minimum, have arbitrarily large strong templates of height bounded exponentially in their order, with a fixed exponential base along some subsequence. Alternatively it must use more templates or an operation outside substitution. Large templates chosen only because their Hamilton counts are huge do not remove the obstruction.

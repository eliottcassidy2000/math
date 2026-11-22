# A split-quotient selector for two exact apex cuts

**Proved selector, not a uniform lifting theorem.** Let \(T_-,T_+\) be a labelled one-arc seam and \(W^\epsilon_{uv}\) the even positive two-slot weight in [CUT-GAUGE](../chanpak_mirimanoff_unit_lifts_20261122/CUT-GAUGE.md). Use normalized cuts \(q_\epsilon=Q_{T_\epsilon}/2\). To preserve the signed gap and reach \(E_{p+1}\pm1\) from \(E_p\pm1\), the exact target is
\[
(q_-,q_+)=(D,D),\qquad D=(E_{p+1}-E_p)/2=2^{p-2}(3\cdot2^p-1).
\]
Everything below also tests any specified integer pair of targets and any initial pins on contacts.

## 1. Split the arcs, not the absolute target

Put \(x_i=[i\in S]\in\{0,1\}\), \(w^\epsilon_{uv}=W^\epsilon_{uv}/2\). If \(h^\epsilon_u=\sum_{u\to v} w^\epsilon_{uv}\), the cut polynomial is
\[
 q_\epsilon(x)=\sum_u h^\epsilon_u x_u-\sum_{\{u,v\}}w^\epsilon_{uv}x_ux_v .
\]
For an arbitrary split \(L\sqcup R\) of the unpinned vertices, absorb contributions from pinned vertices into an affine term. The two cuts become
\[
 \mathbf q(x,y)=\mathbf k+\mathbf A(x)+\mathbf B(y)-
 (x^\top K_-y,\;x^\top K_+y),\qquad K_\epsilon=(w^\epsilon_{uv})_{u\in L,v\in R}.
\]
Here the local vectors \(\mathbf A,\mathbf B\) contain all within-side terms as well as the exact weighted outdegrees, rather than a guess from symmetric weights. Choose shared left and right feature maps (over \(\mathbb Q\)) so that the displayed two bilinear terms depend only on the feature values \(\alpha(x),\beta(y)\). For example, one can use bases of the column spaces of \([K_-\ K_+]\) and \([K_-^\mathsf T\ K_+^\mathsf T]\). Arbitrary rational bases are fine; all values must be exact.

For every \(\alpha\), store the distinct attained vectors \(\mathbf A(x)\) and one witnessing contact per vector; for every \(\beta\) do likewise with \(\mathbf B(y)\). For each group pair the cross value is fixed; hash the larger list and test complements of the smaller in its exact two-dimensional target. A solution exists if and only if a group pair answers yes, since every joined witness is a valid physical contact. For enumeration rather than existence, store lists of witnesses, not just one. The number of lookups, excluding the arithmetic to form the half profiles and cross values, is at most
\[
 2^{|L|}+2^{|R|}+\sum_{\alpha,\beta}\min\bigl(|A_\alpha|,|B_\beta|\bigr).
\]
A rank bound by itself does not bound the number of distinct integer feature values: even one dyadic projection can encode an entire contact tape. Thus this is a compression guarantee only when actual feature fibers are large; it is an unconditional exact test regardless.

## 2. Common modules give proved small interfaces

Let a nonempty proper common module \(M\) have the same induced orientation in both hosts. A module here means each outside vertex either beats every member or loses to every member, in each host. The cuts can differ arbitrarily inside the outside set; the hypotheses are automatic if a module is unchanged and the marked reversal is wholly outside. Take an ordered two-slot cover, collapse each maximal run inside \(M\) to one port symbol, and keep all outside vertices distinct. If one slot ends at \(u\in M\) and the other starts at \(v\notin M\), the skeleton prescribes a distinguished last run at \(u\), and the other direction prescribes a distinguished first run at \(u\). Feasibility at each boundary between a run and outside vertex is independent of its filling, by the module hypothesis.

For each total number of runs \(j\), let \(E_j^-(u)\) count ordered nonempty \(j\)-slot covers of \(M\) ending at \(u\) in one designated slot, and \(E_j^+(u)\) count the analogous covers starting at \(u\). The choice of designated index does not change either count, by permuting slots. Group the compressed skeletons by \(j\); the precise weight factorization for \(u\in M,v\notin M\) is
\[
 W^\epsilon_{uv}=\sum_{j=1}^{R}
 \left(E_j^-(u) O^\epsilon_{j,-}(v)+E_j^+(u) O^\epsilon_{j,+}(v)\right),
\]
where \(O^\epsilon_{j,-}(v)\) and \(O^\epsilon_{j,+}(v)\) count the appropriate two-slot skeletons ending/starting as just specified, and \(R\) bounds their port visits across both hosts. Every pair of compressed skeleton and filling reconstructs a unique two-slot cover, so there is no path-order or split-position multiplicity. The \(2R\) module features are shared, and the joint cross matrix \([K_-\ K_+]\) has rank at most \(2R\), without assuming that the full hosts are ordinal sums or nonstrong. Imposing pins merely makes their contributions offsets and leaves the matrix bound intact.

For a transitive \(M\) of size \(m\) and position \(i\) in increasing arc order, the shared features are explicit:
\[
 E_j^-(i)=\sum_{k=1}^{j}(-1)^{j-k}\binom{j-1}{k-1}
 k^{i-1}(k-1)^{m-i},\qquad E_j^+(i)=E_j^-(m+1-i).
\]
Fix \(i\) in the distinguished slot. When \(k-1\) other slots are available, a predecessor may use any of \(k\) slots and a successor any of the \(k-1\) others. Inclusion–exclusion makes those other slots nonempty; the distinguished slot is already occupied. The convention \(0^0=1\) covers endpoints and \(j=1\). Both cut matrices use these same transitive feature functions even when the two outside profiles are unrelated.

**Ordinal control.** In a common one-way split \(L\to R\) the bound sharpens to two shared interface functions on the unchanged side. Write \(H_L^-(u)\), \(H_R^+(v)\) for full paths ending/starting at the indicated vertex, and \(Z_L^+(u)\), \(Z_R^-(v)\) for ordered two-slot covers with one designated slot starting/ending at the indicated vertex, the other allowed empty. Then
\[
 W_{uv}=H_L^-(u)H_R^+(v)+ Z_L^+(u)Z_R^-(v).
\]
The first summand comes from ending the first slot in \(L\) and starting the other in \(R\), forcing one full path in each side. In the opposite order, split each slot at the ordinal seam; this gives the independent pointed covers counted by the second summand. Inside either side, every weight is multiplied by the other side's ordinary Hamilton count. A literal unit seam with the reversal in \(L\) then forces \(H(R)=1\): Hamilton counts multiply and their difference is two, while all counts are odd. Hence a large unchanged nontransitive ordinal suffix is not a new nontrivial interface for an existing unit seam; the module theorem genuinely goes beyond this ordinal case.

## 3. A middle-mark simplification and an exact parity census

If the flipped edge joins the two middle vertices \(a,b\), then \(W^-_{ab}=W^+_{ab}\), since the edge lies between slots in either defining cover. After pinning all other contacts the gap polynomial \(q_+-q_-\) is therefore affine in just \(x_a,x_b\):
\[
 q_+-q_-=g_0+g_a x_a+g_b x_b.
\]
Thus enumerate exactly the up to four middle patterns solving this equation, and only then test the lower absolute cut. In particular a necessary condition for the four patterns to survive is equality of all three coefficients with zero, not merely that the common middle weight is even. In the other direction the absolute cut may retain its mixed middle term, so this is not scalar linearization of both cuts.

For any fixed pins, reduce the normalized cuts modulo two and regard their two coordinates as Boolean quadratic forms on the \(m\) free variables. This gives an exact polynomial-work census of the simultaneous congruences \(Q_\epsilon\equiv 2D\bmod4\), not an estimate from independent edges. For each nonzero \(\lambda\in\mathbb F_2^2\), let \(c_\lambda\) be the constant of \(\lambda\cdot\mathbf q\) and \(f_\lambda\) its part with \(f_\lambda(0)=0\). Form its alternating polar matrix \(B_\lambda\) with off-diagonal entries the modulo-two combinations of normalized weights. Write \(R_\lambda\) for its radical and \(r_\lambda\) its (necessarily even) rank. Let
\[
 S_\lambda=\sum_{x\in\mathbb F_2^m}(-1)^{f_\lambda(x)}.
\]
If \(f_\lambda\) is nonzero on the radical, \(S_\lambda=0\): translation by a radical vector of value one negates the sum. Otherwise choose a symplectic basis \(e_i,f_i\) of a complement to the radical, where each chosen pair has polar product one and other products vanish. The two-variable sum on a pair, with linear values \(A=f_\lambda(e_i),B=f_\lambda(f_i)\), is \(2(-1)^{AB}\). Consequently
\[
 S_\lambda=2^{m-r_\lambda/2}(-1)^{\sum_i f_\lambda(e_i)f_\lambda(f_i)}.
\]
The number of assignments with \(\mathbf q\equiv \mathbf d\bmod2\) is exactly
\[
 \frac14\left(2^m+\sum_{\lambda\ne0}(-1)^{\lambda\cdot\mathbf d+c_\lambda}S_\lambda\right).
\]
Gaussian elimination supplies the radical, rank and symplectic pairs. Pinning changes the constant and linear terms, so the rank alone is not an exact test, but the full formula is. If all three ranks are at least four, every parity label has at least \(2^{m-4}\) witnesses. Conversely a zero exact census discards the entire subcube before either integer cut is evaluated. Running the census on either one-bit child and taking a nonempty child reconstructs a parity-compatible contact, with no truth table. For the perfect step \(D\) is even at \(p\ge3\), so the required parity label is precisely \((0,0)\). None of the rank conclusions ensures this parity witness is an integer lift.

## 4. Common skew symmetry genuinely folds the exact selector

Suppose the same involution \(\iota\) reverses the arcs of both hosts. Since only one edge changes, its unordered endpoint pair is \(\iota\)-invariant; an orientation-reversing involution cannot fix both endpoints, so it swaps them. For either host, reversing each path and exchanging the two slots gives \(C(u,v)=C(\iota v,\iota u)\), hence \(w_{uv}=w_{\iota u,\iota v}\). Restrict the contact to be constant on each orbit \(U\) of \(\iota\). In this chamber the two mirror members of the middle marked pair receive the same contact. The two *directed* cuts reduce simultaneously to the same undirected quotient form
\[
 q_\epsilon(x)=\sum_{U<V}J^\epsilon_{UV}\,[x_U\ne x_V],
 \qquad
 J^\epsilon_{UV}=\tfrac12\sum_{u\in U,\ v\in V}w^\epsilon_{uv}.
\]
Indeed the involution exchanges every arc from \(U\) to \(V\) with a same-weight arc from \(V\) to \(U\). Thus \(J^\epsilon_{UV}\) is integral, irrespective of individual weight parities after normalization. There is no contribution within an orbit. Any anti-automorphism of a tournament fixes at most one vertex; hence the quotient has at most \(\lceil |T|/2\rceil\) contacts. Complementing all of them changes neither cut, so one orbit can be pinned without loss in this chamber.

This also gives an exact *visibility criterion in the symmetric chamber*: every symmetric apex extension preserves the gap iff all the quotiented weights \(J^-_{UV}=J^+_{UV}\). Necessity follows by applying the same two-point polarization to the quotient cuts and dividing by two; sufficiency is immediate. This does not contradict the full-cube noninvisibility result: the off-chamber patterns are exactly where the flipped arc itself can have one contact on each side. It also gives an exact exponential-size fallback in genuinely symmetric hosts: test \(2^{\lceil |T|/2\rceil-1}\) orbit patterns rather than \(2^{|T|}\), or feed the smaller quotient into the split selector. Exhaustion here makes no claim about asymmetric contacts.

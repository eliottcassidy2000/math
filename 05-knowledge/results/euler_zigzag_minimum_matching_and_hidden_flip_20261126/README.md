# The Euler zigzag response at minimum order

The [previous matching bound](../fourth_simplex_parity_cube_and_matching_response_20261126/README.md) forces at least \(2r+1\) vertices for a nonzero full contrast on \(r\) disjoint marked arcs. Here is a uniform labeled family at equality, strong in every state for \(r\ge2\). A second calculation classifies **every** single-arc tie at the base state and separates true isomorphism moves from changes that the uncolored Hamiltonian count misses.

## Explicit minimal strong family

On \(0,\ldots,2r\), set every arc forward except \(2r\to 0\). Mark the matching \(i\to r+i\) for \(1\le i\le r\). Give each reversal a sign \(-1\). Removing the marked arcs still leaves a strongly connected graph for \(r\ge2\): from 0 all internal vertices are directly reachable, and any internal vertex except \(r\) reaches \(2r\) directly; \(r\) uses \(r+1\). The fixed back edge returns to 0.

By [Grinberg–Stanley's cycle expansion, Theorem 1.38](https://arxiv.org/abs/2307.05569) and the matching bound, only spanning cycles using the entire matching survive the contrast. Each gets weight 2. Contract the marked pairs within such a cycle; every cycle uses \(2r\to0\), so reading it from 0 puts block \(r\) last and forward. A reverse-to-reverse connector always exists, a forward-to-forward connector never exists, and either mixed connector exists exactly for increasing block labels. In particular, the free word of length \(k=r-1\) ends reverse when \(k>0\).

Weight each reversal by \(x\). The forward slots of the free word are isolated. Maximal runs of forced comparisons are increasing: there is a possibly empty initial block of length \(2j\) with \(j\) reversed slots, followed by a sequence of blocks of length \(2j+1\) with \(j+1\) reversed slots, including singletons. The boundaries between successive blocks are unconstrained reverse-to-reverse steps. With block lengths \(\ell_i\), the labelings number \(k!/\prod_i\ell_i!\). Thus the EGF of the weighted spanning cycles with the whole matching is
\[
\sum_{k\ge0}C_k(x)\frac{u^k}{k!}
=\frac{\cosh(u\sqrt{x})}{1-\sqrt{x}\sinh(u\sqrt{x})}.
\]
The series only uses integer powers of \(x\). At \(x=-1\) it becomes \(\sec u-\tan u\), so the Hamiltonian contrast is \(2(-1)^{r-1}E_{r-1}\), where \(E_k\) is the [Euler zigzag number](https://arxiv.org/abs/0912.4240). The replay enumerates every Hamiltonian count through \(r=6\): the contrasts start \(2,-2,2,-4,10,-32\).

This gives a sequence beyond figurate numbers with a precise role on both sides: in normalized height \(N=(H-1)/2\), the full contrast agrees up to sign with the degree-\(r\) simplex finite difference in directions \((1,\ldots,1,E_{r-1})\). This is an aggregate identity, not a pairing of individual state heights or a proof of the two spectra.

## Every one-arc response of the base tournament

Let \(T_n\) be the all-forward order on \(0,\ldots,n-1\) except for \(n-1\to0\), also defined at even \(n\). Write \(P_T(z)=\sum_\mathcal C2^{|\mathcal C|}z^{n-\sum_{C\in\mathcal C}(|C|-1)}\) over vertex-disjoint directed odd cycle packs, including the empty pack. For positive integer \(z\) this counts path covers by \(z\) labeled colors, each color carrying exactly one directed Hamiltonian path on its possibly empty vertex set. Write \(O_m(z)=((z+1)^m-(z-1)^m)/2\) and \(B_m(z)=((z+1)^m+(z-1)^m)/2\).

Reversing the distinguished back edge gives \(\Delta P=-2z O_{n-2}\). For a forward arc \(i\to n-1\) with \(0<i<n-1\), set \(a=i-1,\ b=n-i-2\). Then
\[
\Delta P(z)=2\big[z^{a+2}O_b(z)-z^{b+1}B_a(z)\big].
\]
The deleted cycles through the edited arc have an even subset below \(i\) besides the anchor 0; the added cycles through the new back edge have an odd subset between \(i\) and \(n-1\). The two backward arcs share a tail, so neither one cycle nor a disjoint pair can use both. Reflect labels and take the converse for \(0\to j\).

For an internal forward arc \(0<i<j<n-1\), use \(b=j-i-1,\ s=n-b-4\). Its response is
\[
\begin{split}
\Delta P(z)={}&2z^{s+3}O_b(z)-2z^{b+1}O_s(z)\\
 &+z\big[(z+1)^s(z+2)^b-(z-1)^s(z-2)^b\big]\\
 &+4z^2\sum_{\substack{1\le t\le b\\t\ \mathrm{odd}}}\binom bt O_{s+b-t}(z).
\end{split}
\]
In order, these are new cycles using only \(j\to i\), deleted cycles using \(i\to j\), new cycles using both backward arcs, and newly possible disjoint two-cycle packs. In a cycle using both backward arcs, a vertex strictly between \(i,j\) can be omitted, placed before \(j\), or used after \(i\), giving the \((z\pm2)^b\) factor; all \(s\) outside vertices have just one potential location. A cycle using only \(j\to i\) takes \(t\) middle vertices, leaving \(s+b-t\) candidates for the separate old-back-edge cycle.

**Classification for \(n\ge5\).** Interior adjacent flips preserve the entire polynomial, by exchanging the two consecutive labels. Every interior nonadjacent flip strictly increases \(H\): at \(z=1\) the two-back-edge cycles include an injection of the deleted cycles by placing the same outside subset in the only available locations, and the term \(2O_b(1)\) is positive. At an endpoint, for \(a,b\ge1\), the response at 1 is \(2^b-2^a\). When \(a=0\) or \(b=0\) there is no tie for \(n\ge5\). Thus the only other ties are the two midpoint-to-endpoint flips when \(n\) is odd. The replay checks all pairs for \(3\le n\le10\), at two or three evaluations (411 independent subset-DP comparisons).

For odd \(n=2r+1\), the flip \(r\to2r\) therefore has
\[
H_{\rm base}=H_{\rm edit}=2^{2r-1}+1,\qquad
\Delta P(z)=z^r(z^2-1)\big[(z+1)^{r-2}-(z-1)^{r-2}\big].
\]
For \(r\ge3\) this tie is between different isomorphism classes: the difference is strictly positive at every integer \(z\ge2\). At \(z=2\), the difference is \(3\cdot2^r(3^{r-2}-1)\). This is an explicit mechanism for nontrivial equal-height branches: what is conserved is a single evaluation, not the path-cover polynomial.

## Transfer and the absolute-height test

The [three-switch](../pollock_chanpak_literal_parity_cube_20261126/README.md) and [four-switch](../fourth_simplex_parity_cube_and_matching_response_20261126/README.md) insertion gaps both equal \(P_M(2)\) for inserted module \(M\). At \(r=3\), the midpoint flip holds the module's \(H=33\) fixed while raising \(P_M(2)\) from 612 to 660. The three-switch endpoints move from \((209097,209709)\) to \((226569,227229)\); the four-switch endpoints move from \((26200533,26201145)\) to \((27861429,27862089)\). In particular, fixing the module's \(H\) does not independently fix the absolute starting height.

In Pollock units the corresponding gap parameters are 306 and 330. For each of the four exact starts the replay exhausts the positive 3-cube plus common fifth tetrahedron stencil: factor the gap into three positive directions and enumerate the positive lower corner until its low sum exceeds the start. There are 190, 304, 2142 and 3168 candidate placements, respectively, and none leaves a tetrahedral remainder. This excludes only that marked stencil at those endpoints. An exact signed response can transfer across the research problems without transferring either absolute endpoint.

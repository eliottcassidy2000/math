# Count the output fibers: a uniform Moore bound and exact one-apex inversion

Second overnight transfer: auxiliary matching completions on the Moore side and removable padding on the Chan–Pak side give misleading multiplicities. Passing to the actual outputs yields a stronger lower bound on one side, and an exact inverse and sharp small-gap ceiling on the other.

## 1. Uniform native bound: both endpoints have more than \(2^{18903}\) paths

**Conditional on the existence of a degree-57 Moore graph.** This applies to every complement spine on 3250 vertices. The previous upper bounds \(2^{28860}\) (native) and \(2^{28863}\) (mate) remain. In particular all native \(E_p\pm1\) through \(p=9452\) are excluded, where \(E_p=2^{p-1}(2^p-1)\). No congruence beyond oddness is asserted.

Start with the balanced crossing graph of the displayed spine, with \(k=1625\) original vertices per side. Add \(q=71\) auxiliary vertices on each side; join each to all opposite original vertices, with no auxiliary-to-auxiliary edges. This augmented graph has an \(r=22\) factor, as follows. In the max-flow cut criterion take \(s,t\) original vertices and \(i,j\) auxiliary vertices in the left and right test sets. It is enough that
\[
 e(s,t)+it+js\ge r(s+t+i+j-k-q)
\]
for each such pair of sets. If \(s,t\ge r\), the worst choice is \(i=j=0\), and only \(s+t>k+q\) needs testing. If \(s<r\le t\), the worst choice is \(i=0,j=q\); here \(qs\ge r(s+t-k)\) because \(q\ge r\) and \(t\le k\). The reflected case is identical. If both are below \(r\) the demand is negative, even taking \(i=j=q\).

For \(s,t\ge r\) let \(d=r(s+t-k-q)>0\). The centered spectrum of \(B=2A+I\), where \(B^2=225I+4J\), gives
\[
 2n\, e(s,t)\ \ge\ 115st-15\sqrt{st(n-s)(n-t)},\qquad n=3250.
\]
This certifies \(e(s,t)\ge d\) for 604424 size pairs with \(s\le t\). For each of the 82 remaining pairs use the exact second moment. Suppose \(e(X,Y)\le c=d-1\). Write \(x=|X|,\ y=|Y|\), taking either ordering of the two sets. The nonnegative vector \(g=B\mathbf1_Y\) has total sum \(115y\) and squared norm \(225y+4y^2\). In \(X\) all entries are even and sum to at most \(2c\); outside \(X\) the \(y\) entries on \(Y\) are odd and all others even.

If \((115y-2c)/(n-x)-2c/x>2\), the relaxed minimum of \(\|g\|^2\) under these data takes the maximum \(X\)-sum \(2c\): below that limit transfer two units from a largest entry outside \(X\) to a smallest entry inside, decreasing the square sum and preserving parities and nonnegativity. With grouped sums fixed, start each even entry at 0 and odd entry at 1, then make the cheapest increments of 2 (marginal costs \(4,12,20,\ldots\) and \(8,16,24,\ldots\), respectively). For all 82 size pairs one of the two orderings makes the resulting convex minimum exceed \(225y+4y^2\), by at least 44. This finishes the flow test. This proof keeps one additional true matching edge per matching relative to the simple spectral test using 72 auxiliary vertices.

[Schrijver's lower bound](https://ir.cwi.nl/pub/2154/2154D.pdf) gives at least \((21^{21}/22^{20})^{1696}\) perfect matchings of the regular factor. Upon deleting the auxiliary vertices, each yields a true 1554-edge crossing matching. An original matching can be extended in at most \((22!)^{142/22}\) ways: each completion splits into a perfect matching of each of two \(71\)-by-\(71\) bipartite adjacencies, whose auxiliary rows have degree at most 22. Apply [Brégman's upper bound](https://www.mathnet.ru/php/archive.phtml?jrnid=dan&option_lang=eng&paperid=37755&wshow=paper), using monotonicity of \((a!)^{1/a}\). Divide out these completions, then extend each distinct true matching to paths by [the block-order switching lemma](DEGREE57-FACTORS.md), counting the mate in its converse. The first and last 71 sites are fixed within the unmatched halves for a given matching. The middle ordered pairs read back the actual matching, so distinct matchings and orders do not collide. In the switching lemma there are 1555 objects including the anchor, \(D=59\), and \(E\le52929\). Hence both native counts satisfy
\[
 H\ \ge\ (21^{21}/22^{20})^{1696}\,
 \frac{1554!}{(22!)^{71/11}}
 \left(\frac{1316}{1317}\right)^{52929}>2^{18903}.
\]
The final comparison is an integer comparison after taking eleventh powers, not a floating-point cutoff.

## 2. Exact integer inverse for every one-apex gap

Take a transitive chain, write one bit per site, and let 1 mean the apex beats it. Let \(h(w)\) count its Hamilton paths; with empty word \(h=1\). The [weighted inversion formula](../one_apex_join_representation_mass_pi2_over_12_20261106.md) assigns an inversion of span \(j-i\) weight \(2\) at span 1, and \(2^{j-i-1}\) otherwise, then adds 1. Indeed, if the apex's exit is at \(i<j\), with entry at \(j\), the intermediate chain sites can independently go before or after the apex. The increasing-chain insertions (including apex at an endpoint) contribute 1 plus the number of adjacent inversions. The bit weight budget gives \(1\le h(a)\le a+1\) when \(a\) is the integer encoded by the word without leading zeros: charging each inversion to its left 1 gives at most that bit's place value. Below use \(h(0)=1\).

Fix an edited apex contact. Remove leading zeros **strictly before** it and trailing ones **strictly after** it, since they are common dominating or dominated sites. The remaining word is \(a\,b\,\overline{\operatorname{rev}(c)}\) in binary: \(b\) is the edited bit and \(a,c\ge0\) are written without leading zeros, with 0 represented by an empty flank. Set \(K(u)=\lceil u/2\rceil\). Cross-flank inversions contribute exactly \(2ac\) by the place-value factorization, and the contributions involving \(b\) are \(2K(a)\) for \(b=0\) and \(2K(c)\) for \(b=1\). Reflection preserves each flank's internal count. Thus, up to swapping flanks,
\[
 H_{\rm low}=2ac+h(a)+h(c)-1+2K(a),\qquad
 H_{\rm high}=H_{\rm low}+2(K(c)-K(a)),\quad 0\le a\le c.
\]

For a prescribed pair \((H,H+2d)\), \(d\ge1\), this yields a **complete exact inverse** for this whole one-apex class, not just a modular selector or a search at one fixed word length:

- If \(a=0\), only \(c=2d-1,2d\) can qualify, and the test is \(h(c)=H\).
- Otherwise enumerate four fixed parities \(a=2u-\epsilon,\ c=2(u+d)-\eta\), with \(u\ge1,\ \epsilon,\eta\in\{0,1\}\). For each parity, \(H_{\rm low}\) strictly increases with \(u\): increasing both flanks by two adds \(4(a+c)+8\) to the product term and two to the contact term; the two \(h\)-values together can drop by at most \(a+c\). Binary search each list up to \(u\le H\), test equality on the actual integer, and reflect the surviving collars to recover both orientations. Arbitrary inactive padding adds no more heights.

For example the \(d=1\) family \(a=2u,\ c=a+2\) already yields \(\gg\sqrt X\) distinct pairs. Nothing in the inverse says an arbitrary target has a preimage.

## 3. Sharp gap-window law, beyond \(+2\)

Let \(N(X,D)\) count **distinct unordered height pairs** obtainable by one apex-contact reversal in this class, with both heights at most \(X\) and positive gap at most \(2D\). For \(1\le D\le X\),
\[
 N(X,D)\asymp
 \begin{cases}
 D\sqrt X,&D\le\sqrt X,\\
 X(1+\log(D/\sqrt X)),&D\ge\sqrt X .
 \end{cases}
\]
Here the comparison constants are absolute and \(X\) is sufficiently large. In particular \(D=1\) gives the square-root ceiling even over all chain lengths and padding.

**Upper bound.** Zero flanks produce only \(O(D)\) candidates. Otherwise \(2ac\le X\), \(a\le\sqrt{X/2}\), and \(0\le c-a\le2D+1\). Hence sum \(O(\min(D,X/a))\) over \(1\le a\le\sqrt{X/2}\).

**Lower bound.** Use just positive even flanks \(c=a+2d\), \(a\le\sqrt X/4\), and \(1\le d\le\min(D,\lfloor X/(16a)\rfloor)\). Then \(H_{\rm high}\le2a^2+4ad+3a+4d+1\le X\) for all sufficiently large \(X\). For each fixed \(d\), strict increase in \(a\) rules out collisions; for different \(d\) the actual gaps distinguish the pairs. Sum over even \(a\) in this range. Upper and lower sums have the asserted order: for \(D\le\sqrt X\) a constant fraction of \(\sqrt X\) terms give order \(D\); otherwise split at \(a=X/D\), using a harmonic sum beyond that point (when it lies past the lower range, the constant-order \(X\) term suffices).

The barrier concerns one-apex **arcs**, not all one-apex heights or all single-arc tournament moves. For a density-improving Chan–Pak selector, one has to change the interface, not just add padding or narrow tolerances on one apex. Conversely, a degree-57 matching search can safely use auxiliary sites as a counting relaxation only after bounding how many completions hide the same actual transfer.

Reproduce the 604506 flow-cut size tests, strict factorial bound, and full two-word and inverse checks on every nondecreasing pair of flanks below 256 with `python3 fibre-replay.py`.

# Exact local laws for the two-tail +2 charts

The [two variable-tail families](VARIABLE-TAIL-SELECTOR.md) and [eight-site core](EIGHT-SITE-FIXED-CHANNELS.md) have genuine arc reversals with heights \(H,H+2\). The exponent-31 integer tests show neither of their two absolute target neighbors occurs in these families. This note explains why **no fixed congruence sieve can reproduce that exclusion**.

Write \(G=(H-k)/2\), so
\[
 G=a vu+bv+cu-h(q(v)+q(u)),\qquad
 (h,a,b,c,k)=(29,640,191,158,111)\text{ or }(651,17616,6768,2691,2313).
\]
The transposed chart just exchanges \(v,u\). The [one-variable state mixing theorem](TWO-COORDINATE-RIGIDITY.md) makes each \((v\bmod 2M,q(v)\bmod M)\) uniform in the limit, independent of its fixed high prefix. With independent tail extensions, the two states are independent. Thus the limiting modular law below counts **ordered parameter descriptions**, not distinct heights.

## Prime powers and the biases at 3, 7, 29 and 31

For any prime \(p\nmid h\), conditional on the two flank values the digit sum variables make \(G\) uniform modulo \(p^e\). To translate to physical \(H\) at powers of two keep the fixed odd parity. At odd primes dividing \(h\) (all are simple factors), the \(q\)-term vanishes modulo \(p\). When \(a\) is a unit,
\[
 aG+bc=(av+c)(au+b)\pmod p .
\]
Among the \(p^2\) argument pairs, zero product occurs \(2p-1\) times and each nonzero product \(p-1\) times. Conditional on any fixed pair of flank residues, the state digit sum provides uniform lifts over all higher \(p\)-adic digits because \(v_p(h)=1\). The same lift argument works at 3 for the eight-site chart, where \(a,b,c\) are all zero modulo 3. CRT gives the full independent product law at arbitrary fixed moduli.

| Core | Prime | Physical \(H\) exceptional class | Frequency | Other classes |
|---|---:|---:|---:|---:|
| six-site | 29 | 6 | \(57/841\) | \(28/841\) each |
| eight-site | 3 | 0 | 1 | 0 |
| eight-site | 7 | 1 | \(13/49\) | \(6/49\) each |
| eight-site | 31 | 3 | \(61/961\) | \(30/961\) each |

For example the six-site twin-prime local survival factor at 29 is \(785/841\), since \(H=0,-2\bmod29\) are both ordinary classes. Its ratio to the uniform survival \(27/29\) is \(785/783\). All eight-site lower endpoints are composite, since they are multiples of 3 and at least 2313. These are local factors, not a claim of prime asymptotics.

At every fixed modulus the six-site chart reaches all odd-compatible physical residues. The eight-site chart reaches precisely those also compatible with divisibility by 3. At any odd exponent \(p\), \(E_p=2^{p-1}(2^p-1)\equiv1\bmod3\); among lower queries \(E_p+1,E_p-1,E_p-3\), only \(E_p-1\) passes in the eight-site chart. That is the shared lower query for a potential \(E_p-1\to E_p+1\) pair. In particular **every fixed congruence test admits the six-site exponent-61 queries and the eight-site central query**. There is no claim that any of these integer equalities is attainable.

## Both balanced windows and the full counting range

For binary widths \(L,K\), the leading-order normalized height \(H/2^{L+K}\) is \(2a(v/2^L)(u/2^K)+O(2^{-L}+2^{-K})\) uniformly. Condition on any fixed high prefixes for the two values, and their remaining digits mix independently. Approximate positive-length flank windows by high-prefix intervals; this proves the modular law even with both flanks required in prescribed macroscopic windows. Approximating a macroscopic height interval by rectangles yields the joint limit in such intervals too, since the limiting product of two independent continuous coordinates has no atoms.

The same law holds in the **full height cutoff** \(H\le X\), rather than only when sizes are controlled. For \(u,v\ge1\),
\[
 \#\{(v,u):H(v,u)\le X\}=\frac{X}{2a}\log X+O(X).
\]
Indeed \(2avu\le H\le2avu+O(u+v+1)\), so the difference from the product region is \(O(X)\). For clarity this can be bounded by removing all pairs with either parameter smaller than a fixed constant, costing \(O(X)\). For the rest, at fixed \(v\), comparing
\(\lfloor X/(2av)\rfloor\) with the cutoff intersection loses \(O(1+X/v^2)\): solve the upper linear estimate for \(u\), and summing gives \(O(X)\). The elementary harmonic count of the product region supplies the leading coefficient. The zero axes, when included, contribute just \(O(X)\).

For a fixed modulus and \(\varepsilon>0\), choose an aligned block width \(B=2^r\) so the last \(r\) uniform bits bring the state within \(\varepsilon\) of uniform **from any entering prefix state**. Every full \(B\times B\) box inside the product region then follows the target modular law up to \(O(\varepsilon)\), independently of its high coordinates. A decreasing hyperbola meets \(O_B(X)\) boxes along its boundary, including the axes; each has fixed size \(B^2\). The \(O(X)\) change from the true cutoff also costs lower order. Dividing by \(X\log X\), letting \(X\) increase and then \(\varepsilon\) decrease proves the same local frequencies over the full range.

A growing factor \(X\log X\) in descriptions is not positive density of distinct integers; collisions still have to be controlled. The result explains why fixed modular filters alone cannot decide the exact exponent-61 questions. Run `node tail-local-laws.js` for all exceptional residue frequencies, exact \(p^2\)-lift counts, and the exponent offset check.

# Moving equal-rate cancellation: four consecutive odd counts from one cyclic core

Given \(U=TT_a\to D\to TT_a\), \(a\ge 1\), add a sentinel \(t\) beaten by all of \(U\), and a vertex \(x\) beaten by \(t\) and a chosen mask \(V\subseteq D\), with every other arc from \(x\) into \(U\). This is strong: \(x\) reaches \(U\) through its first guard, and \(U\) reaches \(x\) through \(t\). The single-mask tape formula counts paths by whether \(x\) begins the path, follows \(t\), or follows a mask vertex:
\[
H=h(D)+q(D)4^a+2^a(E(D,V)-L(D,V))-J(D,V).
\]
Here \(q\) counts ordered optional two-path covers, \(E\) constrains the first path's end to lie in \(V\), \(L\) constrains the second path's start to lie in \(V\), and \(J\) imposes both constraints. All four cover counts refer to an ordered partition of \(D\); unconstrained slots may be empty. Thus \(E,L\) may both change drastically in the full graph; it is enough to match their signed difference.

## Physically distinct high reads, four consecutive odds at every height

Use lex-forward code for the arc bits in lex order \(i<j\), with bit 1 meaning \(i\to j\).

| 5-vertex code | mask | h | q | E | L | E−L | J |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 9 | {0,1} | 5 |64|39|15|24|6|
| 9 | {1,2} | 5 |64|41|17|24|8|
| 6 | {0,1,2} | 5 |64|42|18|24|10|
| 6 | {1,2,3} | 5 |64|44|20|24|12|

The two codes each have a single nontransitive SCC of order 4 and Hamilton count 5 (the same isomorphism class), and one transitive singleton. The singleton is a sink for code 9 and a source for code 6. The four masks therefore give four configurations built from *one* fixed cyclic core, not two independent nontransitive factors. In increasing J order the outputs are
\[
64\cdot4^a+24\cdot 2^a-1,\quad
64\cdot4^a+24\cdot 2^a-3,\quad
64\cdot4^a+24\cdot 2^a-5,\quad
64\cdot4^a+24\cdot 2^a-7.
\]
The moving terms \(2^aE\) and \(2^aL\) change at each step by \(2^{a+1},2^a,2^{a+1}\), respectively, canceling exactly. The surviving steps are \(-2,-2,-2\). At \(a=1\) the outputs are 303,301,299,297.

## No nontransitive tape factor is necessary for a three-count window

For \(D=TT_n\) ordered by increasing labels \(0,\ldots,n-1\), direct block assignment gives
\[
E(V)=\sum_{v\in V}2^v,\qquad L(V)=\sum_{v\in V}2^{n-1-v},\qquad
J(V)=\sum_{u<v\in V}\left(2^{v-u-1}+\mathbf 1_{\{v=u+1\}}\right).
\]
Indeed, if the end of the first path is later in the order than the start of the second, every intermediate vertex independently chooses its slot. For the opposite order a cover exists exactly when the named endpoints are adjacent. The two first-order sums agree for reflection-symmetric masks, even though both move when a reflected pair is added or exchanged. If \(n\ge 4\) is even, use the empty mask, the central adjacent pair, or the pair immediately outside it (distance 3). Then \(J=0,2,4\) respectively and the outputs are \(2^{n+2a}+1,2^{n+2a}-1,2^{n+2a}-3\). The outer tournament remains strong for all three choices.

Neither window yet grows in length with the height. These constructions show that high-coordinate equality is not necessary; equal signed projection suffices and supports several exact repairs without a repeated path-count multiplier.

## Both the two-cover count and the guard move: five consecutive odds

With independently chosen nonempty transitive guards of lengths \(a,b\), the formula is
\[
H=q(D)2^{a+b}+E(D,V)2^a-L(D,V)2^b+h(D)-J(D,V).
\]
Write \(a=u+\alpha\), \(b=u+\beta\) for a common growth parameter \(u\). The two high-rate invariants are *only* \(A=q2^{\alpha+\beta}\) and \(B=2^\alpha E-2^\beta L\). The residual is \(h-J\); all orders obey \(H=A4^u+B2^u+(h-J)\). This is also a sharp rate test: if two fixed recipes differ by a bounded amount along unbounded \(u\), their \(A\) and \(B\) must agree, since a nonzero \(4^u\) term dominates any \(2^u\) term and a nonzero \(2^u\) term dominates a fixed residual. In particular the individual \(q,E,L\) need not match.

The following five literal cases test both cancellations at once. Set \(b=u+3\), \(a=u+\alpha\), with \(u\ge0\).

| D order and code | V | h | q | E | L | J | \(\alpha\) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 7, 229304 | {1,3,4,6} | 9 |336|178|124|60|1|
| 6, 4329 | {0,1,2,5} | 9 |168|91|125|62|2|
| 7, 2065790 | {0,2,5} | 9 |336|198|129|64|1|
| 6, 4329 | {0,2,3,5} | 9 |168|95|127|66|2|
| 6, 4329 | {0,3,4,5} | 9 |168|99|129|68|2|

In every row \(q2^{\alpha+3}=5376\), \(E2^\alpha-8L=-636\). The residuals are \(-51,-53,-55,-57,-59\). The five counts are therefore consecutive odds for every \(u\); at \(u=0\) the window in increasing order is \(4681,4683,4685,4687,4689\). The complete tournaments have the same order \(13+2u\): the \(q\)-doubling operation moves a transitive singleton from the outer guard into the masked sector. All have exactly one nontransitive order-five SCC before the wrap; row 1 uses a different SCC isomorphism class, while the remaining rows use the same SCC (the identification from the first five vertices of code 4329 to code 2065790 is \((1,0,5,3,4)\)). Thus this is not a claim of one fixed core for all five rows: four use one fixed core, and five use a two-choice single-SCC bank. There is never multiplication from a second cyclic factor.

## Rigidity even when the guards are allowed to reschedule

Fix any finite bank of profiles \((q,E,L,h,J)\) with \(q>0\). Choose two bank entries and allow all four guard lengths \(a,b,a',b'\) to vary. If their minimum tends to infinity and the two full counts stay a bounded distance apart, then for every sufficiently large member of this sequence:
\[
 q\,2^{a+b}=q'\,2^{a'+b'},\qquad
 E2^a-L2^b=E'2^{a'}-L'2^{b'}.
\]
In particular the surviving signed difference is exactly \((h-J)-(h'-J')\), even though individual high coordinates and their schedules may all move.

For the first equality, the main term \(q2^{a+b}\) dominates both single-guard terms since their ratios to it are bounded by a constant times \(2^{-\min(a,b)}\). Thus the two main terms must have ratio tending to 1. Their ratio belongs to \(\{(q/q')2^d:d\in\mathbb Z,\ q,q'\text{ in the bank}\}\), a discrete set in every compact interval away from zero. Eventually their ratio equals 1. Once the top terms cancel, the second displayed difference is bounded by the initial bound plus the largest residual difference in the bank. It is divisible by \(2^{\min(a,b,a',b')}\) and therefore eventually zero. This reduces the search for a bounded repair within a fixed bank to its two exact high constraints, without imposing equality of \(q,E,L\).

For nonempty masks \(E,L,E',L'>0\), both individual guard-length offsets are bounded as well. The first equality forces \((a'-a)+(b'-b)\) to come from a finite set. If \(a'-a\) becomes very large and positive, then \(E'2^{a'}>E2^a\), and simultaneously \(L2^b>L'2^{b'}\). Their positive signed changes cannot cancel; swapping roles excludes a large negative shift. For empty masks this last conclusion about individual shifts need not apply.

More precisely, write \(c=a'-a\), \(d=b'-b\). On a nonempty-mask subsequence where the two bank entries and \(c,d\) are fixed, the necessary equation is
\[
\big(E-2^cE'\big)2^a=\big(L-2^dL'\big)2^b.
\]
If both parenthesized coefficients vanish, the guards can grow independently (the individually rescaled coordinates already match). Otherwise the coefficients must have the same sign and ratio an integral power of two, and \(a-b\) is uniquely determined. Thus every asymptotically bounded pair from a finite nonempty-mask bank is eventually on one of finitely many exact matched rails: independent growth with both effective readings fixed, or a fixed difference between the guard lengths with the readings genuinely trading off.

## Every raw coordinate can move together

For a small two-count illustration use the following masks with equal guards \(a=b=u+\alpha\), \(u\ge0\).

| D order and code | V | h | q | E | L | J | \(\alpha\) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 4,63 | {0,1} | 1 |16|3|12|2|2|
| 5,4 | {0,4} | 5 |64|19|37|8|1|

Here all of \(h,q,E,L,J\) change. The two exponent-bundled coefficients are nonetheless the same: \(A=256,B=-36\). Thus the full strong counts are \(256\cdot4^u-36\cdot2^u-1\) and \(256\cdot4^u-36\cdot2^u-3\), first \(219\) and \(217\). In normalized height units the moving positive contact readings are \(12\cdot2^u\) and \(38\cdot2^u\) and the moving negative readings are \(48\cdot2^u\) and \(74\cdot2^u\), so each changes by \(26\cdot2^u\) with opposite signs. This is an existence result for coherent growth and does not keep the same cyclic core or order.

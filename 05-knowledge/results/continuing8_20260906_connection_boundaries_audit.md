# Independent referee of Student chart rigidity and modular depth delay

**Full proof and source audit: PASS, no mathematical repair requested.** Both the polynomial scalar-intertwiner theorem and the exact all-prime-power recognition family hold in the stated scopes. The full-source lower bound on the recognition row is proved, not inferred from a particularly successful lift. The saturation distinction in the incoming-work consumer is also sound.

The complete [producer report](continuing8_20260906_connection_boundaries.md) was read at SHA-256 3273facc997e0e977c0fe671960f5c4a97eaffa1eb085419e4b39ac119b772eb. Its frozen source is pinned by the independent referee. The actual source normalization and complete integral diagonal image were checked against the proved incoming planar_jc_long_20260906_depth.md and its routed THM-4308/THM-4369 dependencies. No producer implementation is imported.

## 1. The Student operator and its chart boundary

For characteristic zero and \(m\ne0\), testing the proposed identity on the source functions \(1,x\) first forces \(h=x/\phi\), then
\[
(x^2+6)\phi\phi'=x(\phi^2+6).
\]
If a nonconstant polynomial \(\phi\) has degree \(d\) and leading coefficient \(a\), comparison in degree \(2d+1\) gives \(da^2=a^2\). Thus \(d=1\). Writing \(\phi=ax+b\) gives the three coefficients
\[
-ab,\qquad 6a^2-b^2-6,\qquad 6ab.
\]
Since \(a\ne0\), they force \(b=0,a^2=1\), and hence exactly the two stated charts and scalars. Conversely substitution checks both for every source polynomial. For actual integer indices \(m\ge2\), both test functions already lie in the source degree cap, so this is an obstruction before any additional degree loss from composition.

The polynomial and nonzero-index restrictions are retained. At \(m=0\), the rational scalar \((x^2+6)\phi'/(\phi^2+6)\) works for every nonconstant polynomial chart. The declared \(\phi=x^2\) hostile at nonzero index has residual \(x^3+12x-6/x\). The parameter-covariant scaling identity is exact for nonzero scale \(a\): changing \(x\) to \(ax\) simultaneously changes the quadratic parameter from \(c\) to \(a^2c\). Fixing the parameter at six would lose that coordinate.

No claim is made about arbitrary analytic/algebraic charts or a matrix gauge. Indeed, allowing algebraic charts would change the problem: an equation \(\phi^2+6=A(x^2+6)\) supplies local square-root examples for general constants \(A\). This is consistent with, rather than an exception to, the stated polynomial theorem.

## 2. Full native source and exact first failure

The variables \(x,t\) are independent polynomial indeterminates. Put \(u=x^2t\), \(\pi=t(1+u)\), and \(y=xt\pi\). Direct integer multiplication gives
\[
\pi^3-y^2=t^3(1+u)^2.
\]
For \(P=p^a\), \(L=\lfloor P/2\rfloor\), \(\epsilon=P-2L\), and \(T=3L+\epsilon\), this proves the actual source identity
\[
\pi^\epsilon(\pi^3-y^2)^L=t^T(1+x^2t)^P.
\]
In characteristic \(p\), its only terms are \(t^T\) and \(x^{2P}t^{T+P}\), both with unit coefficients. The exponent \(P\) is a prime power; the coefficient field has characteristic \(p\), not characteristic \(p^a\). No evaluation identity such as \(x^p=x\) is used.

The complete monomial universe is essential to the converse. A monomial \(\pi^c y^e\) has intercept \(2c+3e\). At intercept \(2T\), its nonnegative integer solutions are precisely
\[
c=T-3j,\qquad e=2j,\qquad 0\le j\le L.
\]
There is no extra monomial at this intercept hidden at a higher row. Its unsigned diagonal is \(u^j(1+u)^{T-j}\). Removing \((1+u)^P\) leaves \(u^j(1+u)^{L-j}\); these form an integral unit triangular basis of the full degree-at-most-\(L\) polynomial module. Equivalently, after the producer's signed substitution, the full source is
\[
(1-z)^P k[z]_{\le L}.
\]
This is an equality of full modules over every field, including characteristic two where the displayed signs coincide.

In characteristic \(p\), the inverse of the common factor is
\[
(1-z)^{-P}=(1-z^P)^{-1}=1+z^P+z^{2P}+\cdots .
\]
Since \(L<P\), the unique quotient compatible through order \(L\) is \(q=1\). It agrees through relative order \(P-1\), but cannot have the required degree-\(P\) coefficient. Therefore no source element can match the target through absolute row \(T+P\), whereas the displayed actual lift matches all earlier rows. Projection includes row \(N\); this convention explains the exact first failed row rather than an off-by-one jet length.

Over characteristic zero, the inverse coefficient at relative degree \(L+1\) is the nonzero integer \(\binom{P+L}{L+1}\). Thus the old first failed row is \(T+L+1=\lfloor4T/3\rfloor+1\). Subtracting gives the exact delay \(\lceil P/2\rceil-1\), zero at \(P=2\) and unbounded for every fixed prime as \(a\) grows. The target still never belongs to the full depth-zero algebra; the theorem postpones recognition, not the eventual obstruction.

## 3. Saturation and the one-sided modular consumer

The incoming integral height-two depth image on a diagonal is multiplication by \((1-z)^\rho\) on an integral polynomial space of bounded degree. Every finite row projection has a leading unit triangular block. Before that block is exhausted the projection is the whole ambient integer row space; afterwards it is the graph of an integral map over the first independent coordinates. This proves saturation directly.

Consequently an integral prefix that lies in the rational depth image has an integral lift, so failure of its full projected test modulo a prime certifies rational failure. The argument uses the full integral image, not a generic rank specialization or a nonzero determinant. Modular acceptance is not a rational compatibility theorem, as the family above demonstrates.

The Student target \(2x\) belongs to a different lattice. At index two, the rational source \(-1/2\) maps exactly to \(2x\), while reduction modulo two admits the zero source. No polynomial source over the integers modulo four can map to \(2x\): its \(x\) coefficient is \(12\theta_2-4\theta_0\), always divisible by four. This independently strengthens the displayed degree-one matrix control without changing its stated scope. Thus modulo-four rejection can coexist with rational compatibility when the image is not saturated.

No arithmetic identification with the Mahler reduced denominator or ordinary launch clock is supplied. The report correctly presents that comparison as a restriction on a proposed transfer, not a ring homomorphism or an infinite-orbit consequence.

## 4. Independent verification and freeze

The [referee source](../../04-computation/continuing8_20260906_connection_boundaries_audit.py) uses sparse polynomial coefficient arithmetic, not the producer's symbolic differentiation/expansion engine. It checks surviving charts, five chart obstructions, the \(m=0\) exception, fresh rational scalings, and the unrestricted coefficient obstruction modulo four.

For the depth theorem it enumerates the complete native monomial intercept solutions, expands actual bivariate source products, and computes ranks of the **unfactored** native coefficient matrix augmented by the target. In all eighteen declared cases, membership holds through \(T+P-1\) and fails at \(T+P\). These comprise the producer's fourteen prime-power cases plus \(P=16,32,81,17\). At \(P\le32\), every native source column and the complete integer lift are independently reconstructed; fraction-free determinants give the exact characteristic-zero obstruction \(\binom{P+L}{L+1}\).

There are **5,077 always-active exact gates**. Both normal and optimized outputs were captured as raw subprocess bytes and are identical LF, with no carriage returns.

~~~text
python continuing8_20260906_connection_boundaries_audit.py
python -O continuing8_20260906_connection_boundaries_audit.py
~~~

The full source and proof are accepted as written. The universal claims follow from the complete triangular module and the two-function differential identity, not from the prime or monomial test banks.

Referee source SHA-256: 24c8e8c9fb478af1a12fbf16ee755c06eae79a1caacf45a677b066056956b49a.

Referee output SHA-256: f1479de9a534d51e6ca92241e02833fb0a42a7e1ba1c8a998ae873fb1fdde294.

Pinned producer source SHA-256: bba0c468e6a42142ec623b1a5bd2deb1a3e08d8398629ea13e5ec4f3dd9cae43. Producer output: 69cf4d6b7d8a7686e0544426dfef603bf944070097cc89e750f3bd5caff0320c.

# Exact carry descent and a quantitative bound on multiplicative amplification

**Proved, for the specified physical arithmetic grammars.** This does not close the Chan–Pak conjecture. The inherited least-extra-hole reduction (prime or seven times a prime) and the obstruction to black-box \(H\)-only nonmultiplicative lifts are not reproved here. Instead the descent below decides absolute coverage for a finite common-base affine controller and quantifies exactly what ordinal products can add to a sparse collection of affine and polynomial charts.

## 1. Exact common-base affine controller

Suppose a finite directed graph of port types carries actual tournament constructions. A finite set of starting type/seed pairs is given. Each directed edge acts on a scalar address by
\[
 v\longmapsto B^{k}v+b,\qquad B\ge2,\quad k\ge1,\quad b\in\mathbb Z,
\]
and is physically available for all states specified by its tail type. There are finitely many terminal types, with a final count \(N=Uv+V\), \(U,V\in\mathbb Z\). External guards are allowed only if already encoded in the finite type graph. In particular this statement does not convert a numerical proxy, or a second coordinate that has not been retained, into a physical completion.

**The exact set of nonnegative completed counts is base-\(B\) recognizable.** More concretely, a finite automaton reading least significant digits first accepts exactly the attainable absolute targets, not their residues modulo any restricted collection of moduli.

First fix a terminal read and conjugate every map by \(w=Uv+V\). The edge now has translation \(b'=Ub+(1-B^k)V\), so the desired terminal scalar is the count itself; the new seeds are \(Us+V\). Duplicate this graph for each of the finitely many terminal reads (and use the original permitted terminal states for each copy). A constant read \(U=0\) contributes at most one count if its terminal is reachable, which is a finite graph check.

The automaton starts at a terminal state with carry \(c=0\) and traverses edges **backward**. If the unread numerical prefix is \(X\), the actual scalar at its current type is \(X+c\). On reversing an edge of exponent \(k\) and translation \(b'\), consume \(k\) input digits as \(d=X\bmod B^k\). The predecessor scalar is then
\[
 \lfloor X/B^k\rfloor+c',\qquad
 c'=(d+c-b')/B^k,
\]
and that edge is permitted precisely when the displayed division is integral. A uniform finite bound is
\[
 |c|\le C,\qquad C=\left\lceil 1+\frac{\max_e |b'_e|}{B-1}\right\rceil,
\]
because \(|d|\le B^k-1\) implies \(|c'|\le(C+B^k-1+|b'|)/B^k\le C\). Thus there are finitely many carry/type states. Subdivide a length-\(k\) digit edge into single-digit states; \(k\) is bounded by the fixed graph. At a seed of scalar \(s'\), accept exactly when the remaining digits encode \(s'-c\), with zero high padding allowed. If that integer is negative there is no accepting continuation. A witness preserves all original edge labels; replaying it forward gives an actual construction and equality with the given target. Conversely any physical construction in the grammar reverses through these states and accepts, with sufficient high zero padding. Project back to canonical representations by allowing an arbitrary string of zero digits after the final true input digit: this is just reachability by zero-labelled edges in the finite automaton.

Accordingly each of the following is an exact finite decision for the declared grammar, with no dependence on testing larger and larger sample ranges:
- whether it covers every odd count except any given finite exception list;
- whether it covers all sufficiently large odds, and the complete finite complement when it does;
- the least missed odd count if an odd miss exists;
- the same questions in any fixed arithmetic progression.

For instance, determinize and take the complement within canonical odd input words, then intersect with the fixed exception/progression automata. Emptiness is graph reachability. A productive cycle on an accepting path decides infinitude. Reverse the language into most-significant-first form, determinize again if needed, and use shortest length followed by smallest lexicographic word for the least missed target. These decisions are for positive scalar words themselves; membership in a prime-only lane is not asserted to be a regular language. A polynomial output read can still be coupled to this test one target at a time: list its possible integer preimages by exact polynomial solving, then apply the address automaton to each, rather than approximating a near-square by a root.

## 2. Entropy controls height when the charts grow in the real place

For the same graph, assume every completed word of total exponent length \(L=\sum k\) has positive completed height \(N\gg B^L\), with constant uniform over the finite seeds and terminals. For fixed \(\sigma>0\), form the nonnegative transfer matrix
\[
 M_\sigma(i,j)=\sum_{e:i\to j}B^{-\sigma k_e}.
\]
If \(\rho(M_\sigma)<1\), then the sum of \(N^{-\sigma}\) over all presentations is finite: bound it by the sum of entries of \(\sum_{r\ge0}M_\sigma^r\) times a constant. Therefore the number of distinct output integers at most \(X\) is \(O_\sigma(X^\sigma)\), without any collision assumptions. The abstract growth condition must be checked; an affine graph whose translations permit a bounded recurrent scalar is not covered by this entropy claim. For the exact membership result of §1 no growth condition is needed.

## 3. A multiplicative closure cannot raise a sparse exponent

Let \(\mathcal A\subseteq\{3,4,\ldots\}\) be any actual count bank such that \(A(X)=\#\{a\in\mathcal A:a\le X\}\le C_0X^\delta\) for a fixed \(\delta>0\). Let \(\langle\mathcal A\rangle\) be all finite products of values in the bank. In the tournament application these are actual ordinary Hamilton counts, by ordinal composition; this does **not** preserve an unscaled \(+2\) seam. Then
\[
 \#\{n\in\langle\mathcal A\rangle:n\le X\}
 \ \le\ X^\delta\exp(O_{\mathcal A,\delta}(\sqrt{\log(2X)})).
\]
Indeed the product
\[
 \prod_{a\in\mathcal A}(1-a^{-s})^{-1}
\]
counts all finite multisets of generators with weight equal to product to the power \(-s\). Different multisets giving the same integer can only increase this upper bound. Partial summation gives \(\sum_a a^{-s}=O(1/(s-\delta))\) for \(\delta<s\le\delta+1\), and \(-\log(1-a^{-s})\ll_\delta a^{-s}\). Thus the desired count is at most \(\exp(s\log X+O(1/(s-\delta)))\); optimize just above \(\delta\). If \(A(X)\gg X^\delta\), the multiplicative bank has the **same exact growth exponent** since each generator is itself a permitted product. Adjoining finitely many certified counts leaves this conclusion unchanged. Adjoining finitely many exact polynomial rays \(P(B^j)\) of positive degree and positive leading coefficient also leaves it unchanged, since each contributes \(O(\log X)\) values. An unrestricted positive-variable polynomial ray \(P(j)\) of degree \(d\ge2\) instead replaces \(\delta\) by \(\max(\delta,1/d)\) as an upper bound.

This is a size obstruction, not a residue obstruction: even if an underlying chart has many modular solutions, products of a sublinear-exponent bank still supply only \(o(X)\) ordinary counts. If \(\delta<1\), the prime number theorem also gives a relative density-one set of primes for which **both** ordinary \(p\) and \(7p\) miss the bank. Indeed \(\#\{p\le X:p\ {\rm or}\ 7p\in\langle\mathcal A\rangle\}\) is bounded by the sum of its counts through \(X\) and \(7X\), which is \(X^{\delta+o(1)}\). This concerns the chosen physical constructions and does not rule out a different tournament with either count.

## 4. Apply it to the literal two-digit plus-two family

Use the fixed-gap sieve packet’s raw physical blocks \(AB,BA\), with common scale \(B=490032\) and translations \(360532,362956\). Start at scalar zero, and complete with the fixed inherited suffix and marked seam, giving the two reads \(F(v),F(v)+2\), where
\[
 F(v)=63618252914778017089536\,v+47124685853114258307047.
\]
These are not decoded stand-ins for the Hamilton counts. Both translations are genuine base-\(B\) digits. An exact target test here is especially simple: for each endpoint subtract its intercept, require a nonnegative exact quotient by the common slope, then read that quotient in base \(B\). Every digit must be one of \(360532,362956\), with the zero quotient separately representing the empty word. Thus the accepted integer yields its literal block history.

A depth-\(j\) output has size comparable to \(B^j\), and unique digits give \(2^j\) different lower endpoints at that depth. Therefore the direct bank has \(\Theta(X^\delta)\) distinct counts through \(X\), where \(\delta=\log_B2\). Its entire ordinal-product closure, even with any fixed certified finite bank or finitely many inherited dyadic collision polynomials, has growth exponent exactly \(\delta\). It cannot complete the ordinary prime or sevenfold-prime lane: almost every prime misses both lanes in this construction. This does not conflict with the sieve theorem for this family, which supplies infinitely many rough consecutive odd pairs with unboundedly many address levels, not absolute coverage.

## 5. Lift an exact controller from direct coverage to multiplicative induction

Let \(\mathcal D\) be any finite set of **verified attainable** odd factors greater than one. Let \(\mathcal E=\{7,21\}\), and let \(\mathcal A\) be the recognizable set of direct outputs of physically verified controllers. Form the uncovered obstruction language of positive odd integers outside \(\mathcal E\):
\[
 \mathcal B
 =\{N:\ N\notin\mathcal A,\ 
           \text{every }d\in\mathcal D\text{ either fails to divide }N
           \text{ or has }N/d\in \mathcal E\}.
\]
It is regular in base \(B\): divisibility by each fixed \(d\) is regular, and the quotient exception is simply a check against the finite set \(d\mathcal E\). Thus a test for emptiness or finiteness of \(\mathcal B\) is an actual finite calculation.

- If \(\mathcal B\) is empty, every admissible odd \(N\) is attained. Induct on \(N\): either use a direct output, or take a factor \(d\) and a smaller positive odd quotient outside \(\mathcal E\). Ordinal composition transfers their two certificates to the original \(N\).
- If \(\mathcal B\) is finite, it is enough to construct a witness for every member of this explicit finite list (or add them to the directly accepted outputs). The same induction then applies with no assumption on witness shape for the quotients.
- This can genuinely weaken required direct coverage on multiples of \(\mathcal D\), while it can never rescue a thin chart of exponent below one. To see the latter, restrict to \(N\equiv1\bmod 2\prod_{d\in\mathcal D}d\); no multiplier divides these infinitely many targets, and they have positive linear counting density in that progression. A sublinear-exponent direct bank cannot cover them apart from a finite list.

The criterion is for a certified finite bank; it does not claim that we have found a bank passing it. A proof architecture that simply multiplies finite sparse rays cannot pass. To succeed, an affine controller would need direct capacity on the reduced progression (and on the finitely many quotient exceptions), or some additional physical nonmultiplicative state-preserving operation.

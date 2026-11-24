# From six ordinal generators to a single strong rigid tournament

**PROVED.** Fix \(a,b,c\ge0\), and set \(m=2^a3^b11^c\). Every residue modulo \(m\) admitting an odd integer is realized by a **strong** tournament whose odd-precedence tournament is transitive. All these residues occur at every sufficiently large common order, with exponentially many ordinary isomorphism classes per residue as the order grows (for fixed \(m\)). Unlike ordinal joining alone, the construction can give a prime outside the prime support of the six seeds. It does **not** solve exact integer coverage.

Source: [six ordinal generators](PRECEDENCE-FIBER.md), [pair-precedence identity](ORDER-PARITY.md), and [sharp multiplicative normal form](CRT-SHARP.md). The new map retains one rooted two-piece count at each cap and an ordered two-piece count through the middle. The ordinal selector retains just \(H\).

## One literal reversal, two integer channels

Write \(Q(G)=\sum_{S\subseteq V(G)}H(G[S])H(G[V\setminus S])\), with \(H(\emptyset)=1\). This counts ordered pairs of paths, each allowed to be empty, covering \(G\). For a marked vertex \(x\) let \(U(G,x)\) count such pairs with the second path starting at \(x\). The right-hand analogue uses a path ending at its mark; equivalently, take the converse before reading \(U\).

Make the ordinal join \(A\Rightarrow M\Rightarrow B\), mark \(x\in A,y\in B\), and reverse **only** \(x\to y\). Assume \(M\ne\emptyset\). Old full paths could not use the reversed arc, since they must traverse \(M\). Each new full path uses \(y\to x\). Before the reversal arc it traverses one piece in each block in order; afterwards it traverses the complementary pieces in the same order. The first piece in \(B\) ends at \(y\) and the second piece in \(A\) starts at \(x\). The pieces assemble independently, and every other new boundary uses an unmodified forward arc. This gives the exact integer formula
\[
 H_{\rm new}=H(A)H(M)H(B)+U(A,x)\,Q(M)\,U(B^{\rm op},y).
\]
Under ordinal joining each of the two ordered slots restricts independently to each block, so \(Q\) multiplies as well as \(H\).

If \(|M|\ge3\), the odd-precedence tournament is **unchanged** by the reversal. For each pair \(i,j\), count new directed simple cycles through them. Every such cycle uses \(y\to x\) and a directed \(x\)-to-\(y\) path through the ordered blocks. With its outer pieces fixed, the middle may range over subsets containing the queried vertices of \(M\) (at most two), using the odd number of Hamilton paths on each such subset. The total is even. The only exception is the candidate with empty middle and no outer vertex other than the marks, which is disallowed: it would use the removed arc \(x\to y\). This toggles the cycle count parity for precisely the pair \(x,y\), canceling its arc toggle in the pair-precedence identity. Thus transitive precedence, and hence rigidity, survive. If both caps are strong, every vertex in \(M\) lies on a triangle \(x\to w\to y\to x\), and the whole result is strong.

## Size-tunable caps

Let \(J_{e+3}\) start as an ascending tournament on \(0,\ldots,e+2\), with its only backward arc \(e+2\to0\), and mark \(x=e+1\), where \(e\ge1\). It is strong, has transitive odd precedence, and \(H=1+2^{e+1}\) by the inherited subset formula. Its rooted read is
\[
 U(J_{e+3},e+1)=2+2^{e+1}=2(1+2^e).
\]
If the root slot omits \(0\) and \(e+2\), the complementary slot has \(1+2^e\) paths. If only \(e+2\) joins the root slot, the two slots have one path. If both extremes join it, the path begins \(e+1,e+2,0\), and each of the \(e\) earlier interior vertices chooses freely between the remainder of that slot and the other slot. The remaining case, in which only \(0\) joins it, is impossible.

Choose the left cap with \(e=3^{\max(b-1,0)}\), and the right cap the converse of the cap with \(e=5\cdot11^{\max(c-1,0)}\). By lifting \(2^1+1=3\) and \(2^5+1=33\), their rooted reads are divisible by \(3^b\) and \(11^c\), respectively. Each \(e+1\) is even, so each height is a unit modulo \(3\) and \(11\): modulo \(3\) it is 2, and \(-1\bmod11\) never occurs among even powers of 2. Let \(k\) be the product of the two heights. Use the six-block normal form to select an ordinal middle with \(H(M)\equiv k^{-1}r\bmod m\). Pad it by at least \(\max(a,3)\) singleton vertices. Each singleton multiplies \(Q\) by two and leaves \(H\) alone. The second integer channel now vanishes modulo \(m\), giving \(H_{\rm new}\equiv r\). Padding is adjustable at all higher orders. The extra order above the six-block normal-form bound is at most \(3^{\max(b-1,0)}+5\cdot11^{\max(c-1,0)}+6+\max(a,3)\) (use a projection when necessary).

At fixed modulus take the seed of height 5 to the power \(d=\operatorname{ord}_m(5)\) (use \(d=1\) at modulus 1). In additional chunks of \(4d\) vertices in the middle, independently use either the ordered join of \(d\) copies of this seed or \(4d\) singleton vertices. The \(H\) residue does not change and the second channel still vanishes. The unchanged transitive precedence identifies every vertex by a canonical rank, so different choices of orientation within identical ordered chunks give distinct ordinary classes. At fixed \(m\) this gives at least \(2^{\lfloor(n-N_m)/(4d)\rfloor}\) classes of order \(n\) for each admissible residue, where \(N_m\) is a common base order.

## Integer gain, and its limit

With both caps \(J_4\) (the right one reversed), marked as above, and three singleton vertices in the middle, the literal arc gives \(25+6^2\cdot8=313\), prime and outside the six-block multiplicative bank. This is not an exact inverse for arbitrary integers. Even allowing both cap sizes and singleton padding to vary, for \(H_{\rm new}\le X\) both cap sizes, each of the six seed multiplicities, and the number of padding singletons are each \(O(\log X)\). For padding use \(Q(M)\ge2^{\#\{\text{singletons}\}}\); for other coordinates use the positive product term. Both \(H(M)\) and \(Q(M)\) depend only on multiplicities, not order. Therefore at most \(O((\log X)^9)\) distinct integers up to \(X\) are reached by this entire capped chart. [The square-sum graph has a sharp failure of the same local inference.](SQUARE-WINDOW-CRT.md)

Replay: run strong-and-square-crt-check.py with Python 3.

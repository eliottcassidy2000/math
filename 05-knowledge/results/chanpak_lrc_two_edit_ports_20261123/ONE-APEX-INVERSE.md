# A complete exact inverse for two apex contacts

**PROVED:** deterministic bit-polynomial admission for the entire two-contact transitive-spine one-apex chart; only \(O(\log p)\) Boolean words need be tested for an odd target \(p\ge5\) and upper height \(7p\). **PROVED:** the chart provides at most \(O(\sqrt X\log X)\) starting heights up to \(X\), so it cannot cover the primes. **FINITE-EXACT:** see the replay and examples below. These statements impose no obstruction on two arbitrary arc reversals in a general tournament.

## Setup and all mark patterns

Use the [contact formula and same-side inverse](ONE-APEX-SEVENFOLD.md). Strip pads outside both the lower active span and the marks, and index the lower active word from \(0\) to \(m-1\). Its first contact is O and its last I. Define
\[
 T=2^m,\qquad L=\sum_{\mathrm I_j}2^j,\qquad R=\sum_{\mathrm O_i}2^{m-1-i},\qquad E=m+1.
\]
Thus \(T/2\le L,R<T\). Every word of length \(m\), including an unreduced word, obeys \(|HT-LR|\le ET\). Indeed, \(LR/T\) uses the weight \(2^{j-i-1}\) for all \(\mathrm O_i,\mathrm I_j\). Relative to \(H\), it includes out-of-order pairs of total weight at most \(m/2\), omits the positive adjacency bonuses (at most \(m/2\)), and omits the constant one. As \(2^{m-2}<p<2^m\), the only possible \(m\) are \(\lfloor\log_2 p\rfloor+1,+2\).

The span must grow by one to four, so at least one changed contact is outside the old span. Two exterior marks on the same side have the [unique exact binary inverse](ONE-APEX-SEVENFOLD.md#two-exterior-contacts-on-the-same-side). For the other two patterns we give necessary quadratic bands, then verify actual equality.

**Opposite exterior marks.** Let their extension lengths be \(s,t\ge1\) with \(s+t\le4\). Put \(K_u=2^{u-1}\) and \(\eta_u=0\) for \(u=1\), \(K_u\) otherwise. The exact gain is \(K_sL+K_tR+C\), with \(C=2^{s+t-2}T+\eta_s+\eta_t\): the first two terms and constants are the two exterior gains, and the new exterior-to-exterior pair contributes \(2^{s+t-2}T\). Multiplying the gain equation by \(L\) and using \(LR=pT+e_0\), \(|e_0|\le ET\), gives the band
\[
 |K_sL^2+(C-6p)L+K_t pT|\le K_t ET .
\]

**One exterior and one interior mark.** Converse puts the exterior mark on the left, at \(-s\), \(s\ge1\). If the interior mark does not erase the last I then \(s\le4\). Otherwise the right edge retreats by \(k\ge1\). The orientation after the interior flip cannot be transitive: from a transitive tournament one flip yields only 1 or \(1+2^j\), never a multiple of seven. Hence it has an active span, with the old first O unchanged. The growth bound now gives \(1\le s-k\le4\). From the intermediate orientation, both flips to the endpoints are extensions, by \(k\) on the right and \(s\) on the left. If \(k\ge3\) both gains are divisible by four, contradicting their difference \(6p\). Thus \(k\le2\) and \(s\le6\).

For each interior index \(j\), set \(\delta=+1\) for O to I and \(-1\) for I to O. Then \(L'=L+l,\ R'=R+r\) with \(l=\delta 2^j,\ r=-\delta 2^{m-1-j}\). Put \(K=2^{s-1}\) and \(A=KT+r\ge T/2\). The gain after the interior flip from adding the left exterior O is \(K L'+\eta\), where \(\eta=K\) for \(s\ge2\), and for \(s=1\) it is one if index zero was switched to I and zero otherwise. This last correction counts the adjacency bonus. The product approximation gives
\[
 (L+l)(R+A)=(7p-\eta)T+e_1,\qquad |e_1|\le ET.
\]
Eliminate \(R\) against \(LR=pT+e_0\), multiplying by \(L\), to obtain
\[
 |F(L)|\le 2ET^2,\qquad F(L)=AL^2+(Al+(\eta-6p)T)L+lpT .
\]
The error is \(Le_1-(L+l)e_0\), and both binary readings are less than \(T\).

## Decode with integer bounds

For \(a>0\), \(|ax^2+bx+c|\le e\) is equivalent to requiring \((2ax+b)^2\) between \(b^2-4ac-4ae\) and \(b^2-4ac+4ae\). Integer square roots and signed floor and ceiling division give at most two integer intervals. Their total width is \(O(\sqrt{e/a})=O(\sqrt{mT})\) in both bands.

The mixed case has an additional exact filter. At a true candidate
\[
 F'(L)=AL-lR+e_1-e_0.
\]
If \(l<0\), \(AL-lR\ge T^2/4\). If \(l>0,j\ge1\), the old contact is O, so \(j\le m-2\) and \(l/T\le1/4\), while \(A/T\ge3/4\); hence \(AL-lR\ge T^2/8\). For \(l>0,j=0\), \(l=1\) and \(A\ge T/2\), giving the same \(T^2/8\) bound when \(m\ge3\). For \(m\ge9\) we have \(T\ge32E\), so \(F'(L)\ge T^2/16\) at any actual solution. Intersect the band with this derivative half-line. On it \(F\) is increasing at that rate, so the retained interval has width at most \(64E\). Small \(m\) use the original band.

Split \(L\) into buckets by its top \(\lfloor m/2\rfloor\) bits. They determine the last \(\lfloor m/2\rfloor\) contacts (I bits). Within one bucket the product bound
\[
 (p-E)T/L\le R\le(p+E)T/L
\]
has union of width \(O(\sqrt T)\), because \(L\ge T/2,\ p<T\), and the bucket has width \(O(\sqrt T)\). Thus only \(O(1)\) buckets of the top \(\lceil m/2\rceil\) bits of \(R\) can occur. Each such bucket determines the first \(\lceil m/2\rceil\) contacts (O bits in reverse binary order). The two buckets recover the entire word; check the marked directions, actual heights, and then include converses. There are \(O(1)\) opposite-side patterns with \(O(\sqrt m)\) left buckets each, and \(O(m)\) mixed patterns with \(O(1)\) left buckets each by the derivative test. Hence only \(O(m)\) literal words are tested overall. The arithmetic operates on \(O(m)\)-bit integers.

## Sparsity without fixing a target

Write \(a=L/T,\ b=R/T\). For opposite-side marks, any valid word lies within \(O(m/T)\) in value of a zero of
\[
 6ab-K_sa-K_tb-2^{s+t-2}.
\]
At least one of \(s,t\) is at most two, so the derivative in the corresponding variable is at least one on \([1/2,1]^2\). For mixed marks, put \(\alpha=l/T,\ \beta=r/T\). The analogous residual is
\[
 7ab-(a+\alpha)(b+\beta+K).
\]
Its derivative in \(b\) is \(6a-\alpha\ge5/2\). All partial derivatives are uniformly bounded in the allowed patterns since \(s\le6\). Fix the word half that pins the other variable to a bucket of width \(O(T^{-1/2})\). Any two solutions within that bucket have their remaining variable within \(O(T^{-1/2})\) of one another by the derivative bound and the \(O(m/T)\) residuals. This allows only \(O(1)\) buckets of the opposing half. There are \(O(2^{m/2})\) words per pattern and \(O(m)\) patterns, plus the smaller same-side contribution. Sum over \(m\le\log_2 X+2\) to obtain \(O(\sqrt X\log X)\) distinct possible starting heights. This is \(o(X/\log X)\), below even the Chebyshev lower bound for the number of primes.

## Literal checks

In [one-apex-inverse.js](one-apex-inverse.js), run the exact BigInt inverse against an independent enumeration of active words and reduced pad placements. It agrees on all 2,046 odd targets 5 through 4,095, with 40 labelled marked packets (6 same-side, 2 opposite-side and 32 mixed, counting both converses). The stated smallest hosts can be checked from this scan: for any smaller host with at most ten contacts, all heights are \(<2^{10}\), so none are omitted by the height cutoff.

- **Opposite-side minimum:** on 11 vertices, contact word \(\mathtt{0110100011}\) (index zero first, 1 means O) and marks 0,9 give \(97\), middle reads \(213,307\), and \(679=7\cdot97\).
- **Mixed minimum:** on 12 vertices, word \(\mathtt{00010100110}\) and marks 0,5 give \(103\), middle reads \(723,85\), and \(721=7\cdot103\). Both middle reads are outside \([103,721]\). The complete inverse gives just this square and its converse at \(p=103\) after removing pads, so **there is no monotone route in this chart at any order** for that target. The same all-order nonmonotonicity test passes for \(p=227,1399\). This does not constrain general two-edit witnesses at those primes.

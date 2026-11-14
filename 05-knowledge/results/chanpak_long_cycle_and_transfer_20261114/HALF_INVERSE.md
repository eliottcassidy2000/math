# Split-axis inverse: recovering a literal source or certifying the one-apex obstruction

The width-renewal theorem proves that one apex with a transitive spine cannot cover Chan–Pak. This note gives a direct target inverse for the whole obstructed family, rather than merely a proof that some target fails.

## Exact boundary product

Cut a binary contact word into nonempty consecutive halves $A|B$. Put $x=[A]_2$, $y=[\operatorname{reverse}(\operatorname{complement}(B))]_2$, and $p(a)=a\bmod 2$. For normalized words the halves each start with 1 in their own reading directions, so their lengths are the binary lengths of $x,y$. Let $h$ be the one-bit count recurrence from APEX_RENEWAL.md. Then
$$
 H(A B)=x y +h(x)+h(y)-1+p(x)p(y).
$$
Indeed the cross-inversions contribute $xy$ by splitting $2^{j-i-1}$ across the seam. The sole extra adjacency correction is $p(x)p(y)$. Internal contributions are additive after removing their duplicate constants; converse gives $h(B)=h(y)$. This formula can be used in *either* direction: every pair of positive integers yields a literal normalized contact-word matrix.

For a prescribed target $N$ this is a product close to $N$ with controlled but genuinely signed corrections to naive product descent. In fact $1\le h(a)\le a+1$ by the recurrence; in the even-input step $a=2p$, if $p$ is odd its last 1 can be stripped and $h(p)\le(p+1)/2\le p$. Hence $xy < H\le (x+1)(y+1)+1$. There is no license to round the product: the corrections retain the two entire half words and their boundary bits.

## Exact monotonicity at a balanced cut

Let $x$ have at least as many bits as $y$. On all $y$ of a fixed bit length $\ell$, the formula is nondecreasing. At an even $y=2q$ it changes by precisely $x-q-p(q)+p(x)$ on incrementing $y$; this is nonnegative since $q\le 2^{\ell-1}-1$ and $x\ge 2^{\ell-1}$. For odd $y$ and $y+1$ remaining at the same bit length, write $y=P0\,1^r$ in bits, so $y+1=P1\,0^r$, with $r\ge 1$. Then $h(y)=h(P)+P+p(P)$, since the trailing ones add nothing. On the other side the very first appended zero incurs $2P+2\ge P+p(P)+2$, the others adding nonnegative amounts. Thus $h(y+1)-h(y)\ge 2$, more than the maximum loss of one in the parity term. This handles the actual carry cases; no smoothness assumption is needed.

## Target-driven implementation

If $N>1$ is at all representable, its normalized width obeys $2\le L\le \lfloor\log_2(N-1)\rfloor+2$, by the end-to-end inversion. For each such $L$, take the left cut of length $\lceil L/2\rceil$ and the right cut of length $\lfloor L/2\rfloor$. Try each $x$ in the corresponding dyadic bin, binary search for the first $y$ whose exact formula is at least $N$, and accept only equality. The sum of the numbers of $x$'s tested over the widths is $O(\sqrt N)$. Each binary search uses $O(\log N)$ exact evaluations and $h$ costs $O(\log N)$ bit-level recurrence steps (these steps are on $O(\log N)$-bit integers). The overall bound is $O(\sqrt N\log^2 N)$ arithmetic operations. Absence in this algorithm certifies absence from every unrestricted one-apex/transitive-spine order, not just all short words.

A monotone plateau is real: $(x,y)=(4,6),(4,7)$ both give $33$, with words `100100` and `100000`. The product may change by four while the physical Hamiltonian count does not move. Even so, the entire width/no-width obstruction is exact. Any unbounded proof by movable sources must leave this family or add a complementary realization mechanism.

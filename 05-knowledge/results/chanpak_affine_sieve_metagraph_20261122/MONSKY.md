# A three-state Monsky character of arbitrary one-apex contact words

**PROVED**, for all one-apex words. The arithmetic interpretation is **CITED**; no elliptic-curve rank, Tate–Shafarevich finiteness, general BSD, or Chan–Pak completeness is asserted.

## Statement

Use a transitive spine $1\to\cdots\to m$ and one apex $a$, with contact bit $b_i=1$ iff $i\to a$. Let $A$ have off-diagonal entries $A_{ij}=1_{i\to j}$ over $\mathbb F_2$, with diagonal chosen to give row sum zero. Put $\nu=\dim_{\mathbb F_2}\ker(A^2+A+I)$ and $D=\sum_{i=1}^m(-1)^{m-i}b_i$. Then
$
\nu=2\quad\Longleftrightarrow\quad
D\equiv \begin{cases}
1\pmod 3 & m\text{ even},\\
2\pmod 3 & m\text{ odd};
\end{cases}
\qquad \nu=0\text{ otherwise}.
$
For $m\ge 1$, the number of words with $\nu=2$ is $(2^m-1)/3$ for even $m$, and $(2^m-2)/3$ for odd $m$. Among the $2^{m-2}$ strong words for $m\ge2$, the respective counts are $(2^{m-2}+2)/3$ and $(2^{m-2}-2)/3$.

## Proof: the hidden three-cycle

Work over $\mathbb F_4$, with $\lambda^2+\lambda+1=0$. The spine block of $A-\lambda I$ is triangular with nonzero diagonal, so an eigenvector is determined by its apex coordinate. For a nonzero eigenvector normalize this coordinate to $1$. With $d_i=m-i+b_i\pmod2$ and $t_i=\sum_{j\ge i}x_j$, its spine equations are exactly
$
 t_{m+1}=0,\qquad
 t_i+b_i=\left(1+\frac{1}{\lambda+d_i}\right)(t_{i+1}+b_i).
$
Every column of $A$ sums to $m\bmod 2$, since $A+A^\top=J+I$ and the row sums vanish. When all spine equations hold, the single apex residual is therefore $(\lambda+m)(1+t_1)$; it vanishes iff $t_1=1$. Frobenius conjugates the two eigenspaces for the distinct roots of $X^2+X+1$. Thus $\nu$ is twice the dimension of the $\lambda$ eigenspace and can only be $0$ or $2$.

Scan from right to left. On an even-distance position, bit $0$ acts as $F_0(t)=\lambda t$ and bit $1$ as $F_1(t)=\lambda^2t+\lambda$. On an odd-distance position the maps are $G_0(t)=\lambda^2t$ and $G_1(t)=\lambda t+\lambda^2$. Pair each even-distance position with the odd-distance position immediately left of it. Then
$
 G_0F_0=G_1F_1=\mathrm{id},\qquad
 G_0F_1=R,\qquad G_1F_0=R^{-1},\qquad R(t)=\lambda t+1.
$
The map $R$ fixes $\lambda$ and cycles $0,1,\lambda^2$ in that order. Starting at $0$, a paired word gives $R^D(0)$ for even $m$, hence acceptance exactly at $D\equiv1$. For odd $m$, its unpaired leftmost bit $b_1$ contributes $+b_1$ to $D$ and then applies $F_{b_1}$. Since $F_0^{-1}(1)=\lambda^2$ and $F_1^{-1}(1)=1$, acceptance requires paired exponent $2-b_1$, exactly $D\equiv 2$.

Each pair contributes $0,+1,-1$ with multiplicities $2,1,1$. The number of $k$-pair words with sum zero modulo 3 is $(4^k+2)/3$, and for either nonzero sum it is $(4^k-1)/3$: the three Fourier eigenvalues are $4,1,1$. With odd $m$ each leftmost bit asks for one of the two nonzero sums. For strongness, $b_1=0$ and $b_m=1$ are necessary (otherwise the first spine vertex is a source or the last a sink) and sufficient (the apex closes the whole spine into a directed cycle). Their contribution to $D$ is $1$. The remaining alternating sum has $m-2$ bits, with rightmost sign negative. In the even case it must be zero, giving the zero-residue formula. In the odd case leave one negative-sign bit unpaired: its two possibilities ask from the pairs for residues $1$ and $2$. This gives the stated strong counts.

## Two literal consequences

* Reversing the extreme arc of a transitive tournament of order $q+2$ is the word $0^q1$, $q\ge1$. By the [Grinberg–Stanley packing formula](https://arxiv.org/abs/2307.05569), $H=1+2^q$: all directed odd cycles use the apex and the last spine vertex, and the odd subsets of the $q$ intervening vertices give $2^{q-1}$ cycles. Here $D=1$, and hence $\nu=2$ iff $q$ is odd. On this rail alone, $\nu=0$ is exactly $3\nmid H$.
* At the same order $5$, the strong words $0001$ and $0101$ both give $H=9$, but $\nu=2$ and $0$, respectively. In either word there are exactly four directed odd cycles through the apex: for $0001$ it is three choices of two vertices ending at vertex $4$ plus the entire four-vertex spine; for $0101$ the three pairs are $(1,2),(1,4),(3,4)$, plus the whole spine. No disjoint pair can coexist. Thus equality of absolute count does not fix even this arithmetic coordinate within one-apex strong sources.

## What transfers, and what does not

The [finite Selmer receipts](../bsd_selmer_growth_20261013a.md) give the arithmetic setting: by CRT and Dirichlet one can realize every tournament as the quadratic-reciprocity tournament of distinct primes $p_i\equiv3\pmod8$. Monsky's matrix then identifies $\nu$ with the pure 2-Selmer dimension of $y^2=x^3-N^2x$, $N=\prod_i p_i$. The [selective-port Schur complement](../session321_selective_returns_and_fano_wiring_20261113/BSD_PARITY.md) already proves parity dependence of a **single constant-contact block**. [THM-5156](../../../01-canon/theorems/THM-5156-the-monsky-nullity-zero-sector-attains-every-odd-count-below-125-except-3-7-13-19-21.md) already proves additivity of $\nu$ under ordinal sums and gives a finite Selmer-refined spectrum. The new step is an exact three-state character and strong-label census for **all** interleavings of the two one-apex contacts, plus the scoped dyadic rail result. The retained arithmetic coordinate is genuinely independent of absolute height. It is not a method for constructing an arbitrary height from its residue.

In the RH and Hodge analogies the transport is limited and explicit: finite character diagonalization makes the residue-count cancellation exact, and the literal tournament words supply the physical witnesses. The character is not a zeta zero, and simultaneous attainment of two scalar invariants is not an algebraic cycle-lifting theorem.

Finite check: run `python3 check_monsky.py`.

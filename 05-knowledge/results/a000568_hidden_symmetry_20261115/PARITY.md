# A000568 splits into exact dyadic parity recurrences

## The all-order parity cutoff

Fix $s\ge2$ and an odd residue $r\bmod 2^s$. For every $n\ge 2^s$, **the number of ordinary tournament isomorphism classes with $H\equiv r\bmod 2^s$ is even**. By converse pairing, the same holds for the number of self-converse classes in that stratum. This is a statement about multiplicities, not an obstruction to realization.

Fix the directed path $1\to2\to\cdots\to n$. The remaining $B=\binom{n-1}{2}$ orientations are independent Boolean coordinates. Each class supplies $H/|\operatorname{Aut}(T)|$ elements of this cube: marking a Hamiltonian path kills the automorphism group. Both numerator and denominator are odd; therefore the cube and class counts have the same parity, including on each residue stratum.

In the full orientation cube the path polynomial is multilinear of degree at most $n-1$. Using centered signs, reversing every arc preserves $H$ and negates every sign, so the true degree is at most the largest even integer $d\le n-1$. Fixing the path arcs does not increase it. The function $Q=(H-1)/2$ is integral on the Boolean cube by Rédei. Its multilinear interpolation has integer coefficients by Boolean Möbius inversion and degree at most $d$.

Put $q=2^{s-1}$. Lucas's theorem gives bit $i$ of $Q$ modulo two as $\binom{Q}{2^i}$. The indicator of a prescribed residue of $Q\bmod q$ is a product of these bits or their complements, so it is represented modulo two by an integer-valued function of degree at most $(q-1)d$ after multilinearizing. Every cube monomial of degree smaller than $B$ sums to an even integer. The strict inequality $(q-1)d<B$ follows at $n\ge2q$, both for even $n$, where $d=n-2$, and for odd $n$, where $d=n-1$. This proves the cutoff. Chiral classes outside the self-converse ones occur in pairs with the same $H$, proving the self-converse refinement.

## Strong tournament residues follow an eventual cycle

Work in $\mathbf F_2[G]$, $G=(\mathbf Z/2^s)^\times$, with basis elements $[r]$. Let $A_s(x)$ be the parity series for ordinary classes, each weighted by $[H\bmod2^s]$, and $S_s(x)$ the corresponding nonempty strong series. The cutoff says $\deg A_s\le2^s-1$. The inherited ordered SCC law gives $A_s=(1-S_s)^{-1}$. Its multiplication is multiplication of residues. The sum of the residue coefficients is A000568 parity, the series $1+x+x^2$: the path cube has even cardinality at every $n\ge3$.

Let $K$ be the exponent of $G$ ($K=2$ for $s=2,3$, and $K=2^{s-2}$ for $s\ge3$). Frobenius gives $A_s(x)^K=(1+x+x^2)^K$ since each group element raised to $K$ is the identity. Thus
$$
 S_s(x)=1+\frac{(1+x^K)A_s(x)^{K-1}}{1+x^{3K}} .
$$
The strong residue parities have eventual period dividing $3K$. If the numerator has degree $D$, then $S_{n+3K}=S_n$ at positive $n>D-3K$. No large-order census appears in this recurrence. The same series gives the self-converse strong parities, by pairing.

Small fixed-path cube counts (see the replay) give $A_2=1+x+x^2+(1+[3])x^3$ and
$$
 A_3=1+x+x^2+(1+[3])x^3+(1+[5])(x^4+x^5).
$$

| $n\bmod6$ | Residues mod 4 with odd strong multiplicity | Residues mod 8 with odd strong multiplicity |
|---:|:---|:---|
| 0 | 1 | 5 |
| 1 | 1 | 5, except order 1 has 1 |
| 2 | none | none |
| 3 | 3 | 3 |
| 4 | 1 | 5 |
| 5 | 1, 3 | 3, 5 |

Every entry of odd multiplicity guarantees a strong self-converse realization in that residue, not just an ordinary one. An even multiplicity makes no claim about support. These are infinite existence statements in fixed residue classes, not absolute-count coverage.

## Higher valuation for the path-marked cube

For $a\ge1$, let $f(z)$ be the indicator of the required residue of $z\bmod q$. In the Newton expansion $f(Q)=\sum_j\Delta^jf(0)\binom Qj$ for $Q\ge0$, coefficients with $j\ge(a+1)q/2$ vanish modulo $2^a$. Indeed the difference operator $D$ on $q$-periodic sequences satisfies
$$
D^q=-\sum_{j=1}^{q-1}\binom qj D^j.
$$
In reducing $2^bD^m$ at $m\ge q$, the weight $m+bq/2$ never decreases: the binomial coefficients are even at $j\ge q/2$ and divisible by four below $q/2$. At $m<q$, every term with $b<a$ has weight at most $(a+1)q/2-1$. This proves the truncation. Each surviving $\binom Qj$ has integer Boolean interpolation of degree at most $jd$; summing monomials now gives divisibility by $2^{B-jd}$ whenever $B\ge jd$. Thus the **cube count**, i.e. class count weighted by $H/|\operatorname{Aut}|$, is divisible by $2^a$ whenever
$$
 B-\big((a+1)q/2-1\big)d\ge a.
$$
Only the first power of two was transferred to ordinary class multiplicities.

## More than one Hamiltonian mark

For a tournament $T$, the isomorphism-orbit generating polynomial for unordered subsets of Hamiltonian paths is
$$
 \frac1{|\operatorname{Aut}(T)|}\sum_{g\in\operatorname{Aut}(T)}
 (1+z^{d(g)})^{H(T)/d(g)}.
$$
An element of order $d(g)$ acts in free cycles of length $d(g)$ on the paths. Summing over classes interpolates A000568 at $z=0$ and the fixed-path cube at $[z]$. Coefficients at positive powers of two are identity-only, hence are exact labelled binomial moments divided by $n!$.

For converse merging, use Burnside over $S_n$ and global converse. In the converse sector a relabeling $g$ fixing the labelled object as an anti-automorphism has order $2d$ with odd $d$, since $g^2$ belongs to the odd ordinary automorphism group. The action $P\mapsto g(P^{\mathrm{rev}})$ has cycles of length $d$ or $2d$: its square has cycles of length $d$. If $R$ counts paths reversed by $g^d$, the fixed-subset polynomial is $(1+z^d)^{R/d}(1+z^{2d})^{(H-R)/(2d)}$. This counts higher marked subsets; it does not reconstruct unmarked coefficients from $H/|\operatorname{Aut}|$ alone.

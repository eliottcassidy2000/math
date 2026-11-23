# All odd-prime continuation states from three guarded physical apex tiles

For a tournament $T$, let $Q_s(T)$ count its ordered $s$-slot directed path covers, empty slots allowed, and $Q_1=H$. Fix an odd prime $p$, put $d=(p-1)/2$, and prescribe either parity for the order of $T$.

**Theorem.** Every vector in $\mathbb F_p^d$ is $(Q_1(T),\ldots,Q_d(T))\bmod p$ for an actual strong tournament of the prescribed parity. One can take $T$ to be a single apex over a transitive ordered chain whose apex contacts are a concatenation of the three fixed words $OI,OOI,OII$. Moreover, for every finite set of odd primes, the continuation vectors can be prescribed simultaneously and independently, with common chosen parity, by one actual such tournament. This is a new construction, not a claim of transitivity for the two older reentry moves $P,W$.

## Literal path-cover read as a matrix word

Write $O$ at a chain vertex beaten by the apex $z$, and $I$ at one beating $z$. Every nonempty concatenation of the three tiles starts with $O$ and ends with $I$; its tournament is strong, since $z$ reaches the whole chain through its first vertex and the chain reaches $z$ through its last.

For a contact word $w$ of expanded length $n$, a disjoint odd-cycle packing either contains no cycle or contains just one cycle through $z$, since the deletion of $z$ is transitive. The latter is specified by $i<j$ with $w_i=O,w_j=I$, and an even-cardinality subset of the $j-i-1$ intervening chain vertices. A cycle contributes weight two and one common path-slot color; the other vertices contribute one color each. The [Grinberg–Stanley odd-cycle-packing formula](https://arxiv.org/abs/2307.05569) therefore gives, for any nonzero $s\bmod p$,
$$
\frac{Q_s(T)}{s^{n+1}}
=1+s^{-2}\sum_{i<j:\,w_i=O,w_j=I}
\left[(1+1/s)^{j-i-1}+(1-1/s)^{j-i-1}\right].
$$
Indeed choosing an intervening vertex for the cycle removes one independent slot color, and averaging the two binomial expansions implements the even-subset condition with the weight two. Adopt $0^0=1$.

For each nonzero $\lambda\ne1$ in $\mathbb F_p$, use
$$
O_\lambda=\begin{pmatrix}1&1&0\\0&\lambda&0\\0&0&1\end{pmatrix},\qquad
I_\lambda=\begin{pmatrix}1&0&0\\0&\lambda&1\\0&0&1\end{pmatrix}.
$$
The top-right coordinate $g_\lambda$ of the product for $w$ is precisely $\sum_{i<j:O_i,I_j}\lambda^{j-i-1}$. Each $\lambda=1\pm1/s$, for $1\le s\le d$, occurs once; as $s$ runs over the representatives, the two signs of $1/s$ run over all nonzero residues. At $\lambda=0$, the sum is the number $g_0$ of adjacent $OI$'s. It adds across guarded words, since junctions have adjacency $IO$.

Thus the three macro words map multiplicatively into the product of the invertible matrix channels, the additive coordinate $g_0$, and expanded length modulo $p-1$. For several primes take the product of all channels and track a single common length modulo the least common multiple of the $p-1$'s. The image of this finite monoid is a group: every inverse is a positive power. In fact an exponent for this ambient group is $E=\operatorname{lcm}_p p(p-1)$: raising a matrix channel to $p-1$ makes it unipotent of degree three, hence of exponent dividing $p$ because $p$ is odd. The scalar channels and length also divide $E$.

## Generate a basis of literal central read changes

Let $C,P,S$ denote the group elements of the physical tiles $OI,OOI,OII$. Set
$$
\mathcal O=PC^{-1},\qquad\mathcal I=C^{-1}S,\qquad u=\mathcal O\mathcal I^{-1}.
$$
At each nonzero channel these are the usual single-letter matrices and their displayed product. At the zero channel, all three physical tiles have $g_0=1$, so $\mathcal O,\mathcal I$ have zero and $C(\mathcal O\mathcal I)^{-1}$ is a pure unit $g_0$ translation. Its length is zero. It affects only read coordinate $s=1$, whose two channels are $2,0$.

Write a unipotent channel matrix as $(a,b,g)$, with $a,b,g$ the upper $12,23,13$ entries and diagonal one. Its additive quotient retains just $(a,b)$; the commutator of two such matrices is central with top-right coordinate $a b'-a'b$. The element $u$ is unipotent in all nonzero channels, with quotient pair
$$
(a_\lambda,b_\lambda)=(\lambda^{-1},-1).
$$
Conjugating by $\mathcal O$ multiplies these two quotient coordinates by $\lambda^{-1},\lambda$, respectively. For each $\tau\ne0$, the polynomial $f_\tau(X)=1-(X-\tau)^{p-1}$ is the spectral indicator of eigenvalue $\tau$ over the chosen prime. An actual word lifting the projection of $u$ is obtained by writing $f_\tau(X)=\sum_i c_i X^i$ over $\mathbb F_p$, choosing $c_i\in\{0,\ldots,p-1\}$, and multiplying the words $(\mathcal O^i u\mathcal O^{-i})^{c_i}$; the order is irrelevant in the additive quotient. Each such lifted word is unipotent in every prime channel and has zero expanded length and zero $g_0$ coordinates.

If $\tau\ne\pm1$, the commutator of the lifts for $\tau,\tau^{-1}$ is therefore central. At the selected prime it is zero except in channels $\lambda=\tau,\tau^{-1}$, and both these entries are nonzero: the projection for $\tau$ takes $b_\tau=-1$ and $a_{\tau^{-1}}=\tau$, and the other takes $a_\tau=\tau^{-1}$ and $b_{\tau^{-1}}=-1$. With the commutator convention $XYX^{-1}Y^{-1}$, the two coordinates are $+\tau^{-1}$ at $\tau$ and $-\tau$ at $\tau^{-1}$.

For each $j=1,\ldots,d-1$, use $\tau=1+1/j$. It is neither $0$, $1$, nor $-1$, and $\tau^{-1}=1-1/(j+1)$. Thus its two nonzero channel shifts affect exactly read coordinates $j,j+1$, both with nonzero coefficients after multiplication by $s^{-2}$ in the normalized formula. The pure zero-channel shift and these $d-1$ adjacent-pair shifts span the full read space: eliminate from the last coordinate downward and finish with the $s=1$ shift. At $p=3$, there are no edges, and the pure shift already spans the single coordinate. Each central correction has length zero.

For simultaneous primes, each constructed shift is central in **all** channels at all primes, even though its values elsewhere were not projected. Raise it to a nonnegative exponent congruent to $1$ at the intended prime and $0$ at all other primes, by CRT. Central additives have those prime exponents. This isolates the shift without changing its intended effect or length, proving independent simultaneous control.

## Constructive output, chosen parity and a size bound

For odd tournament order start at the identity of channel groups (expanded length zero). For even order start with $\mathcal O$ (expanded length one, all top-right and zero-channel coordinates zero). These are instructions in the finite group, not claims that their signed-word spellings are literal tournaments. In either case the base vector of normalized reads is all ones. Given desired $Q_s$ residues, divide them by $s$ in the odd case or $s^2$ in the even case, subtract ones, and solve the edge-basis corrections described above, with CRT isolation when relevant. Append the corresponding powers of the central words. The final length classes are unchanged, so the normalization is exact.

Finally expand into signed occurrences of the three macro letters, replace each inverse macro letter by its $E-1$ positive repetitions, and if the result is empty append $C^E$. The resulting nonempty positive word is a literal strong tournament. Since every $p-1$ is even, its parity of expanded length is exactly the one specified. This is a deterministic constructor using only finite-field polynomial arithmetic, a back substitution on a path, and optionally elementary CRT.

For one prime, combine consecutive conjugations before expanding: the signed product is
$$u^{c_0}\mathcal O u^{c_1}\mathcal O\cdots\mathcal O u^{c_{p-1}}\mathcal O^{-(p-1)}.$$
Here $\mathcal O$ has signed macro length 2 and $u$ length 4, so the shared-conjugation expression is at most $4(p^2-1)$ signed macro letters. A commutator is at most $16(p^2-1)$; its outer scalar power is at most $p-1$. The zero-channel generator is at most five signed macro letters. With $E=p(p-1)$, the order is at most
$$1+2E+3E\bigl(2+5(p-1)+16(d-1)(p-1)(p^2-1)\bigr)=O(p^6).$$
This generous bound avoids any finite orbit search. It gives no bound on, or exact specification of, the *integer Hamilton height* beyond what the actual constructed tournament provides.

## Prime-power lift of the first $(p-1)/2$ readings

The same construction prescribes arbitrary residues for the *first* $d=(p-1)/2$ readings modulo $p^e$ for every $e\ge1$, with both parities and arbitrary finite simultaneous odd prime-power moduli. It does not purport to prescribe all of the larger continuation vector that may be needed modulo $p^e$.

Replace each field by $\mathbb Z/p^e\mathbb Z$, still with $s=1,\ldots,d$; all these $s$ are units. The $\lambda=1\pm1/s$ channels that were nonzero modulo $p$ remain units; the only missing one is precisely $\lambda=0$ for $s=1$ and the minus sign, so its adjacent count remains exactly additive on guarded words. Track length modulo $\varphi(p^e)$ (jointly, the lcm of such lengths). The signed words constructed above for the projected commutators are still exactly unipotent before commutation and exactly central afterward, since having each channel diagonal one follows from the word definition, not from the field projection. Their central read-change matrix modulo $p^e$ reduces to the spanning edge basis modulo $p$. Hence its determinant is a unit, and Gaussian elimination over $\mathbb Z/p^e\mathbb Z$ (always pivot at a unit) solves for independent shifts. The pure zero-channel generator remains exact. The same CRT isolation works by taking exponents one modulo the selected prime power and zero modulo the others.

The inverse replacement uses $E=\operatorname{lcm}_{p^e}p^e\varphi(p^e)$: the $\varphi(p^e)$ power makes every channel unitriangular, and $(1+N)^{p^e}=1$ since $N^3=0$ and $p^e\mid\binom{p^e}{2}$ for odd $p$. Thus signed words again convert to literal nonempty physical words of chosen parity. In particular every pair $(H,Q_2)$ of residues can be simultaneously prescribed modulo any odd $m$ coprime to $3$, even when either residue is zero, unlike a reagent construction restricted to unit pairs.

For a single prime power $q=p^e$, the same signed-size argument uses coefficients at most $q-1$ and gives the explicit order bound
$$
1+2E+3E\bigl(2+5(q-1)+16(d-1)(q-1)(p^2-1)\bigr)=O(q^3p^3)\subseteq O(q^6),
\quad E=q\varphi(q).
$$
Use $\texttt{chanpak-all-prime-apex-powers.py --prime P --power E [--even] TARGET...}$ for an explicit positive constructor. The script uses a straight-line representation: each printed block is an ordered sequence of positive tile runs with an outer repetition count. It does not expand the resulting huge literal tournament. Its modular checker certifies examples at squares and cubes of $3,5,7,11,17,23$, and both parities of a joint target modulo $27,25,343,121$.

## Complete two-reading image including powers of three

For fixed order parity, the only constraint on $(H,Q_2)$ modulo $3^e$ is
$$
Q_2\equiv(-1)^{|T|}H\pmod 3.
$$
Necessity follows from the parity polynomial $Q(-t)=(-1)^{|T|}Q(t)$. To prove sufficiency, use the longer fixed tiles $C=O^eI^e,\ P=OC,\ S=CI$. For $s=1,2$, the two nonunit channels are $\lambda=0$ and $\lambda=3/2$, where $\lambda^e=0$. Across guarded macro junctions any $O,I$ pair has at least $2e$ intervening vertices, so $g_\lambda$ is additive at these nilpotent channels. In both, $g_\lambda(C)=(1+\lambda+\cdots+\lambda^{e-1})^2$; prepending $O$ or appending $I$ contributes zero because the new geometric sum starts at $\lambda^e$. Thus the augmented $\mathcal O,\mathcal I$ have zero in both additive channels and retain their usual matrices in unit channels.

The centered word $C(\mathcal O^e\mathcal I^e)^{-1}$ is pure in the additive channels. Their values are $1$ and $4=(1-3/2)^{-2}$, so with the $s^{-2}$ coefficient it gives a shift $(1,1)$ in the pair of normalized readings. The unit channels for $s=1,2$ are $\lambda=2,1/2$, respectively. The commutator $[u,\mathcal O u\mathcal O^{-1}]$ has top-right coordinate $\lambda^{-2}-1$ there, from the determinant of the two quotient pairs. Its normalized shift is $(-3/4,3/4)$. Its length and additive channel increments are zero. Over $\mathbb Z/3^e$, these two vectors span exactly $\{(x,y):x=y\bmod3\}$, since the difference of the second vector's coordinates is $3/2$. Start at the identity for odd order or at $\mathcal O$ for even order, and apply the two central shifts as before. The normalization factors at the two coordinates are respectively $(1,2)$ and $(1,4)$; the required relation is exactly the condition for the normalized target difference to be in this span. This includes $e=1$.

The same argument permits joint control with the other odd prime-power components: use these same longer guards on all components. At other primes the unit channels retain the earlier $\mathcal O,\mathcal I$ matrices, and the pure center word becomes $C(\mathcal O^e\mathcal I^e)^{-1}$, giving exactly the old pure unit shift in their sole zero channel. All correction words are central in every component, so isolate them by prime-power CRT exponents. Therefore for any odd modulus and either fixed order parity, the full pair image is restricted only by the displayed relation when the modulus is a multiple of three.

Run `python3 chanpak-all-prime-apex-three-check.py` for the long-guard generator replay at powers of three.

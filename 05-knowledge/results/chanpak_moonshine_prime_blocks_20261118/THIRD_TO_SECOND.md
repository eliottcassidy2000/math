# A quantitative prime-bank obstruction using just the third-to-second read

This strengthens the [fixed-core prime-support theorem](PRIME_SUPPORT.md) along a different axis. The core may vary, and no restriction is placed on the prime support of the ordinary count or two-slot count separately. Only the ratio $z=f(3)/f(2)$ is restricted. The main result requires equality of marked **two-slot** end and start totals, not equality at one slot. An unequal-response refinement quantifies what changing the core must buy. All results concern this neutral-switch architecture, not general Chan–Pak coverage.

## 1. Height-neutral pairs must have unusually close ratios

Use the [two-pad transfer](PADDING_OBSTRUCTION.md) on a nonempty marked core $U,I$ of order $m$, with marked end/start totals $\mathcal E_j,\mathcal S_j$. Fix nonempty pads, and for a pad write $x=f(1)/f(2)$ and $z=f(3)/f(2)$. The neutral equation is
$$
 \mathcal E_2 z_L-\mathcal S_2 z_R
       =\mathcal S_1 x_L-\mathcal E_1 x_R,
$$
and the shift is nonzero exactly when the right-hand side is nonzero. Each $x$ is at most $1/2$ by using either single nonempty slot. The [marked slot split](MIRROR_PADDING.md) gives $\mathcal E_1/\mathcal E_2,\mathcal S_1/\mathcal S_2\le 2/(m+1)$. Thus for $\rho=\mathcal E_2/\mathcal S_2$,
$$
0<|\rho z_L-z_R|\le\frac{\max(\rho,1)}{m+1}.
$$
The factor 1, rather than 2, in the numerator follows because the absolute difference of two nonnegative numbers is at most their maximum. In particular if $\rho=1$, then $0<|z_L-z_R|\le1/(m+1)$ regardless of whether $\mathcal E_1=\mathcal S_1$.

## 2. Effective scarcity for balanced two-slot cores

Let $S$ be any fixed finite prime set. Suppose $z_L,z_R$ are positive $S$-units, with no support restriction on $x_L,x_R$. Let $N$ be the larger order of the pads and $M$ the larger number of their strong components (including singleton components). For nonzero neutral shifts at $\rho=1$ there is an effectively computable constant $C_S$ such that
$$
M\log(3/2)+\log(m+1)\le C_S\log(N+2).
$$

Indeed, for a pad of order $n$, $f(q)\le n!\binom{n+q-1}{q-1}$: a concatenation order and weak composition determine an ordered slot cover, though some fail the adjacency tests. Consequently every exponent of a prime in $z=f(3)/f(2)$ is $O(N\log(N+2))$ in absolute value. Distinct positive $S$-units of this exponent height satisfy
$$
 |z_L-z_R|/\max(z_L,z_R)\ge (N+2)^{-C_S}
$$
by the real rational specialization of [Matveev's lower bound for logarithmic forms](https://www.mathnet.ru/eng/im314). Specifically apply it to the nonzero form $\sum_{\ell\in S}b_\ell\log\ell=\log(z_L/z_R)$; for $|\log(z_L/z_R)|\le1$ the relative difference is at least $1-e^{-1}$ times its magnitude, and otherwise is bounded below by $1-e^{-1}$. The effective constant absorbs the polynomial bound on exponent height.

The [ordinary slot split](PRIME_SUPPORT.md) gives $f(3)/f(2)\ge(n+2)/2\ge3/2$ for every nonempty component, and ordinal joins multiply ratios. Hence $\max(z_L,z_R)\ge(3/2)^M$. Combining this with §1 gives the bound.

**Consequences.**
- With fixed support, the larger pad must have a strong component on $\Omega_S(N/\log(N+2))$ vertices. A bounded-size alphabet of strong components cannot escape by changing ordinal words.
- Arbitrary ordinal words over a fixed finite palette of pad shapes give only finitely many nonzero neutral switches, **even when the marked core varies** within the class $\mathcal E_2=\mathcal S_2$. The palette bounds each component order, so $N\le kM$ for a fixed $k$, giving an effective upper bound on $N$, then on $m$. The conclusion is finite up to isomorphism, not a useful numerical census.
- For a palette whose $z$ values are powers of a single real rational $g\ge3/2$, there are no such switches at all. Nonempty pads have positive powers; distinct powers differ by at least $g(g-1)\ge3/4$, contradicting §1. This includes arbitrary powers of one fixed tournament, and in particular transitive pads of any nonempty lengths. This special case needs no logarithmic-form theorem.

## 3. Unequal marked response, controlled arithmetically

If $\rho$ is also a positive $S$-unit, put $B_*=2+\max_{\ell\in S}|\nu_\ell(\rho z_L/z_R)|$. Applying Matveev to the adjusted ratio and using $\max(\rho z_L,z_R)\ge\min(\rho,1)\max(z_L,z_R)$ gives
$$
 M\log(3/2)+\log(m+1)
 \le C_S'\log B_*+|\log\rho|.
$$
Nonzero shift is crucial here: it excludes the case $\rho z_L=z_R$ in which the logarithmic form vanishes.

A finite set of allowed values of $\rho$ therefore gives the same finite-palette obstruction as equality, after adjoining their prime support. More generally, if $\rho$ stays between two fixed positive constants but remains an $S$-unit, the exponent size $B_*$ is $O((N+m)\log(N+m+2))$. To bound its core part, the marked totals at 2 are at most $f_U(3)$, because they count three-slot covers subject to a restriction on one marked slot. Thus with a fixed pad palette, unbounded pad order would force core order at least exponential in $N$: $N=O_S(\log(N+m+2))$. In a fixed reader $p=\alpha f_T(2)\pm H(T)$ with $\alpha\ge1$ on that full source, $p\ge(\alpha-\tfrac12)2^{|T|}$, so those pads would have at most $O(\log\log(p+2))$ vertices. These conclusions do not cover unbounded asymmetry with unrestricted prime support.

## 4. Literal control with balanced two-slot but unequal one-slot totals

The two kinds of balance really differ. On $U=(24,109,121,80,2,25,49)$, mark mask 3. The exact totals are $(\mathcal E_1,\mathcal S_1)=(259,211)$ and $\mathcal E_2=\mathcal S_2=1619$. The two matched pads have profiles
$$
 (H,B,f(3))_L=(1659,58560,777555),\qquad
 (H,B,f(3))_R=(207,7920,105165).
$$
The one cut at $v$ gives a literal same-lattice move
$$
 2788528491\ \longrightarrow\ 3155717691,\qquad B=20639790800640
$$
at order 27, with exact shift $367189200$. The [replay](marked_balance_replay.py) uses subset DP on the pads and counts outer templates independently, rather than deducing the whole counts from the transfer formula. Both ordinary counts and the shared height are divisible by 3; this fixture is a **structural control**, not a prime pair.

There is an exact filter before any prime search. Each pad has at most three runs over two slots: at most one per path before accounting for the single visit to the switching vertex. If a pad has $3\mid H$, then $3$ also divides $f(2)$ by slot parity, and divides $f(3)$ because each nonzero monomial has positive degree. All its ordered nonempty cover readings up through three are consequently divisible by 3, by inclusion–exclusion. Factoring each template through its pad readings shows that both host reads are then divisible by 3, not just the signed offset. Thus prime reads of this form above 3 must have $3\nmid L_1R_1$. Reducing the neutrality law modulo 3 then gives the cheap necessary condition
$$
\mathcal S_1(-1)^{|R|}\equiv \mathcal E_1(-1)^{|L|}\pmod3.
$$
This explains why the displayed fixture, with same-residue one-slot marked totals and opposite pad parities, cannot be repaired by feeding it to an integer reader.

The new mathematical leverage is the separation of **marked response** from **prime support of a single slot ratio**. Growing the core does not rescue a finite pad alphabet when its two-slot response is balanced or ranges over a fixed finite list of biases. Any successful route in this architecture has to change one of these hypotheses and still satisfy the exact absolute target.

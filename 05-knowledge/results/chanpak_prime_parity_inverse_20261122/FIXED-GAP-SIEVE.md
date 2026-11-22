# Sifting both endpoints of a literal +2 tournament switch

**PROVED**, using elementary Fourier analysis and the classical fundamental lemma of the combinatorial sieve. No RH, prime-gap conjecture, or unproved distribution estimate is used.

For infinitely many **actual** pairs $N,N+2$ of ordinary Hamiltonian path counts of strong tournaments, every prime divisor of $N(N+2)$ exceeds $\exp(c\sqrt{\log N})$ for some $c>0$. In particular $\Omega(N(N+2))=O(\sqrt{\log N})$, with multiplicity. The proof gives two-sided order of magnitude $\asymp 2^n/n$ for the number of sifted pairs at each sufficiently large depth $n$. It also gives joint estimates for cutoff Liouville signs. Neither estimate proves primality or bounds $\Omega$ independently of $N$.

## 1. Two positive moves become a two-digit alphabet

Use the inherited [(2,1) strong +2 switch](../chanpak_affine_sieve_metagraph_20261122/UNIT-SWITCHES.md). Its lower count for left and right increments $u,v$ is
$
20311uv+14572u+8654v+6207,
$
and the second literal center contributes exactly 2 more. Fix the left history $P$, giving $u=120$, so that the lower output for a final right increment $v$ is $2445974v+1754847$.

The inherited right maps are $A(v)=82v+60$, $B(v)=5976v+4396$, $C(v)=96264v+71856$. Starting at $v=0$, choose at every step either $AB$ or $BA$, read chronologically. Their scalar maps are $av+b$ and $av+b-2424$, with $a=490032$, $b=362956$. After $n$ such blocks append the fixed chronological suffix $AACAAB$ and apply the two legal switches. Writing $v_n$ for the increment before this suffix, the outputs are $F(v_n)$ and $F(v_n)+2$, where
$
 F(v)=63618252914778017089536v+47124685853114258307047 .
$
These are two actual independently strong hosts, by inheritance. The length-$n$ choices are injective: differences are $2424$ times differences between binary-digit expansions in base $a$ (the largest differing power exceeds the sum of all smaller ones). Also $v_n\asymp a^n$ for $n\ge1$ by the geometric sum formula with positive digits, so $N\asymp a^n$, and the explicit host orders are $O(n)$.

## 2. No exceptional prime divides either output

The exceptional set $E$ of prime divisors of $2a\cdot2424\cdot F'(v)$ is
$
 \{2,3,7,41,47,83,101,191,26021\}.
$
For the eight primes other than 101 the slope $F'$ vanishes modulo the prime; the residue pairs $(F(0),F(0)+2)$, in increasing prime order, are $(1,1),(2,1),(3,5),(33,35),(8,10),(3,5),(2,4),(11440,11442)$. Neither coordinate vanishes. At 101 both blocks agree with $v\mapsto81v+63$, whose orbit from zero returns after 25 steps and avoids both forbidden residues throughout. These finite data are checkable by the companion replay. Thus every depth is allowed. For any $p\notin E$, $F(v)$ and $F(v)+2$ vanish on two distinct residues.

## 3. All-modulus mixing, then the two-dimensional sieve

Take independent fair block choices and let $\mu_{n,q}$ be the law of $v_n\bmod q$, where $q\ge2$ is coprime to the primes in $E$. For a nonzero character $t\bmod q$ its Fourier coefficient has absolute value
$
 \prod_{j=0}^{n-1} |\cos(\pi\cdot2424\,t a^j/q)|.
$
Set $L=1+\lceil\log_a q\rceil$. There is an index in each consecutive interval of $L$ factors whose phase's distance to the nearest integer is at least $1/(2a)$: if all distances were smaller, signed centered representatives would multiply exactly by $a$ at each step, yet the first distance is at least $1/q$ since $2424a$ is invertible modulo $q$. This would give a distance at least 1 after $L-1$ steps, a contradiction. Each complete interval contributes at most $\cos(\pi/(2a))<1$. Thus for suitable absolute effective constants,
$
 \max_{t\ne0}|\widehat\mu_{n,q}(t)|
 \le C\exp(-c_0 n/\log(2q)).
$
This also bounds the error per residue after Fourier inversion. Total variation from uniform is at most $Cq\exp(-c_0 n/\log(2q))$.

For squarefree $d$ coprime to $E$, CRT therefore gives
$
 \Pr[d\mid F(v_n)(F(v_n)+2)]
 =\beta(d)+O\!\left(2^{\omega(d)}e^{-c_0n/\log(2d)}\right),
 \qquad \beta(d)=2^{\omega(d)}/d .
$
Set $D=\exp(\delta\sqrt n)$ with fixed $\delta>0$ sufficiently small. As $2^{\omega(d)}\le d$, the sum of normalized errors up to $D$ is $O(e^{-c_1\sqrt n})$, choosing $\delta$ so that the factor $D^2$ is absorbed. For example for large $n$, use $\log(2D)\le2\delta\sqrt n$ and any $\delta^2<c_0/4$.

Apply the classical [fundamental lemma, in the events formulation of Tao, Corollary 19](https://terrytao.wordpress.com/2015/01/21/254a-notes-4-some-sieve-theory/), with sieve dimension 2 and $z=D^{1/s}$ for a sufficiently large **fixed** $s$. Its main factor is $\prod_{p\le z,\ p\notin E}(1-2/p)\asymp(\log z)^{-2}$ by Mertens and $\sum p^{-2}<\infty$; the dimension condition follows from $\sum_{w\le p\le z}(\log p)/p=\log (z/w)+O(1)$. Its normalized remainder is bounded by the above sum of errors, which is $o((\log z)^{-2})$. This use allows $D$ to be subpower in the sample size $2^n$; a theorem requiring power-saving errors in that size is not being invoked. Thus a proportion $\asymp(\log z)^{-2}\asymp n^{-1}$ of words survive, and the injectivity gives the count estimate. The prime cutoff translates to $\exp(c\sqrt{\log N})$; comparing the log of the product with the log of this cutoff bounds $\Omega$.

In addition, for any fixed $m$ coprime to $E$ and any $r\bmod m$ satisfying $\gcd(r(r+2),m)=1$, one may impose $N\equiv r\bmod m$. It selects a single residue for $v_n\bmod m$, since $F'$ is a unit, and the modulus for CRT in the argument becomes $md$ with $\gcd(m,d)=1$. The initial mass is asymptotic to $1/m$, the new sieve product omits primes dividing $m$, and the aggregate error remains exponentially small in $\sqrt n$. Hence the same two-sided order $\asymp_m 2^n/n$ holds in every such class.

## 4. Joint cutoff Liouville signs on the physical pairs

Write $\lambda_z(k)=(-1)^{\sum_{p\le z}v_p(k)}$. With $z=\eta\sqrt n$ for a sufficiently small fixed $\eta>0$, for every nonempty $J\subseteq\{0,2\}$ we have
$
 {\mathbb E}\prod_{j\in J}\lambda_z(F(v_n)+j)
 \ \sim\ \frac{C_{|J|}}{(\log z)^{2|J|}},\qquad C_{|J|}>0 .
$
Expectation is over all physical length-$n$ words, not over the survivors of the roughness sieve.

To see this, let $e_p=\lceil10\log n/\log p\rceil$ and
$
 q=\prod_{p\le z,\ p\notin E}p^{e_p}.
$
For large $n$, $\log q\le10\pi(z)\log n+\sum_{p\le z}\log p=O(z)$, with the implied constant independent of fixed $\eta$ once $n$ is sufficiently large depending on $\eta$. Taking $\eta$ small enough in the total-variation bound gives $O(e^{-c_2\sqrt n})$ distance to the uniform residue law modulo $q$.

Capping each valuation at $e_p$ gives a statistic of the residue class modulo $q$. Under uniform residues, the difference from the uncapped $p$-adic local average has error at most $O(|J|\sum_{p\le z}p^{-e_p})=O(n^{-9})$. The same bound plus the total-variation error compares the capped and true signs on the literal words, since $p^{e_p}\mid F(v_n)+j$ is a single residue and its union over $p,j$ has that small probability. For a good $p$, $F'$ is a unit, so a single linear form has odd valuation with $p$-adic probability $\sum_{k\ge0}(p^{-(2k+1)}-p^{-(2k+2)})=1/(p+1)$. The odd-valuation events for the two endpoints are disjoint at every odd $p$, because the forms differ by 2. The resulting local sign expectation for a set of $s=|J|$ endpoints is $1-2s/(p+1)$. CRT gives independent local factors over $p\le z$ outside $E$; the excluded primes contribute signs $+1$. All local factors are positive, since $p\notin E$ implies $p\ge5$, and the product is asymptotic to $C_s(\log z)^{-2s}$ by Mertens (the ratio to $\prod_{p\le z}(1-1/p)^{2s}$ converges to a finite positive constant). This dominates both errors above and proves the estimate.

In particular the two cutoff signs each approach an unbiased sign and their joint law approaches independent fair signs. This is genuinely parity-sensitive but only at the stated cutoff. It does not describe the signs of the large remaining factors or circumvent the prime-sieve parity barrier.

## Scope of the connection

The [Allahyar–Yianni preprint, §2](https://github.com/sahibzada-allahyar/math/blob/main/preprint/prime-gaps-180/prime_gaps_at_most_180.pdf) uses a marked endpoint pair to distinguish pair detection from a marginal prime count. The marked-pair distinction is relevant here, but no analytic estimate from that preprint is used: its bound of 180 itself depends on the unestablished Type III hypothesis A4. Our +2 is exact because two literal tournament switches already give it; the sieve does not create the gap and certifies only roughness and truncated parity, not primality.

## 5. A power-scale obstruction for this particular Fourier method

In fact no estimate claiming uniform convergence to the uniform distribution for **all** good moduli up to a fixed exponential in $n$ can hold for this alphabet. Let $k>\max E$ be prime and put $q=(a^k-1)/(a-1)$. This modulus is good: if $p\in E$ divides $a$ then $q=1\bmod p$; if $a=1\bmod p$ then $q=k\bmod p$; otherwise $p\mid q$ would force the order of $a\bmod p$ to equal $k>p-1$.

Take frequency $t=2424^{-1}\bmod q$. Since $a^k=1\bmod q$, the magnitudes in its Fourier product repeat in blocks of $k$. For $0\le j\le k-2$, the centered phase distance is at most $a^{j-k+1}$, by $q\ge a^{k-1}$; for $j=k-1$, the phase $a^{k-1}/q$ is within $1/a$ of 1. The sum of squared distances over a period is therefore less than $3/a^2$. Using $\cos x\ge 1-x^2/2$ and $\prod_j(1-u_j)\ge1-\sum_j u_j$ for $0\le u_j\le1$, the absolute Fourier product per period is at least $1-3\pi^2/(2a^2)>0$. A partial period cannot lower the product more than a whole period. At $n=\lceil Bk\rceil$ for any fixed $B>0$ the corresponding coefficient stays bounded below by a positive constant depending only on $a,B$. The total-variation distance is at least half this coefficient. Taking $B$ sufficiently large puts these witnesses below $\exp(\varepsilon n)$ for any prescribed fixed $\varepsilon>0$. This obstructs uniform power-scale **full-residue mixing**, not average or specially aligned bounds for the two forbidden residues used by the sieve.

Exact finite replay: `python3 check_fixed_gap_sieve.py`.

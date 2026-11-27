# Additive closure inside an Euler product: the radical defect of \(a+b=c\)

This revisits the [old abc Cayley experiment](https://github.com/eliottcassidy2000/math/blob/7a36391d197a164914f464a1814f6c8b00853695/05-knowledge/results/famous_problems_s116.out) and its warning that prime-address coordinates do not linearize addition. The precise abc inequality has a constant \(K_\epsilon\): \(c\le K_\epsilon\operatorname{rad}(abc)^{1+\epsilon}\) for positive primitive triples, or equivalently finiteness of exceptions to \(c\le\operatorname{rad}(abc)^{1+\epsilon}\) for each \(\epsilon>0\). See the formulation in [Congruence ABC implies ABC](https://arxiv.org/abs/math/9909098). Here we compute a limiting law and a logarithmic mean at the additive slice, not its uniform extreme bound. Robert, Stewart and Tenenbaum's [refinement](https://tenenb.perso.math.cnrs.fr/PPP/abc.pdf) studies related kernel statistics and a heuristic for extremes; the explicit local limit below does not invoke their independence heuristic.

## 1. Three changes of coordinates, and two different meanings of phase

For \(\Re s>1\), the Euler formula and its logarithm give
\[
\sum_{n\ge1} n^{-s}
 =\prod_p(1-p^{-s})^{-1},\qquad
\log\zeta(s)=\sum_p\sum_{j\ge1}\frac{p^{-js}}{j}.
\]
Here the branch of the logarithm comes from the convergent power series. The first equality turns **naturals** into **prime exponent vectors**. The logarithm turns products into sums, and the identity \(n^{-s}=n^{-\sigma}e^{-it\log n}\) for \(s=\sigma+it\) adds trigonometric phase in the *logarithm* of \(n\). The harmonic/convergent boundary is visible from the same formula: as real \(s\downarrow1\), the \(j\ge2\) terms remain bounded, while the integral comparison for the harmonic sum gives \(\zeta(s)=(s-1)^{-1}+O(1)\). Thus \(\sum_p p^{-s}=\log(1/(s-1))+O(1)\). This is an absolute-convergence statement about the product, not a claim about its boundary on the line \(\Re s=1\).

The phase that detects \(a+b=c\) is of a different kind. In a finite box,
\[
{\bf1}_{a+b=c}=\int_0^1 e^{2\pi i\theta(a+b-c)}\,d\theta .
\]
This uses ordinary \(n\), not \(\log n\). For example the **unfiltered** sum over all pairwise coprime \(a,b,c\), weighted by \((abc)^{-s}\operatorname{rad}(abc)^{-t}\), has Euler factor \(1+3p^{-t}/(p^s-1)\) for real \(s>0,\ s+t>1\). Each prime can be absent or assigned positively to exactly one leg. The factor \(e^{2\pi i\theta a}\) for an ordinary additive filter is not multiplicative in \(a\), so the same product cannot be kept while inserting that filter. In real logarithmic coordinates the obstruction is explicit:
\[
\log c=\frac{\log a+\log b}{2}+
 \log\!\left(2\cosh\frac{\log a-\log b}{2}\right).
\]
Equivalently \((\sqrt{a/c},\sqrt{b/c})=(\cos\alpha,\sin\alpha)\). The archimedean balance angle, the additive Fourier phase, and the multiplicative \(\log n\) phase are three different observables. In particular, none is a prime-independent replacement for the exact equation.

The earlier [Farey prefix work](../../../07-reflections/the-cut-prefix-is-a-continued-fraction-times-an-euler-product.md) used \(1/\zeta(2)\) as a coprime density next to an additive constraint. Here we can keep \(a+b=c\) exact throughout and still recover an Euler product for one statistic. The [tournament zeta](../../../01-canon/theorems/THM-1926-tournament-zeta-euler-product-over-strong-core.md) is the corresponding caution in another setting: it factors over primitive closed walks and strong components, but Hamiltonian paths need the extra constraint that a walk does not reuse vertices.

## 2. An exact limiting Euler product on the additive slice

Sample \((a,b)\) uniformly from \(1\le a,b\le X\) with \(\gcd(a,b)=1\), and set \(c=a+b\). Put
\[
R=\operatorname{rad}(abc),\qquad D=abc/R .
\]
A prime can divide at most one of the three legs. For any prime \(p\) and \(k\ge1\), its limiting local probability is
\[
\Pr(v_p(abc)\ge k)=\frac{3p^{-k}(1-p^{-1})}{1-p^{-2}}
                  =\frac{3}{p^{k-1}(p+1)}.
\]
In the primitive local domain \(p\nmid\gcd(a,b)\), the events \(p^k\mid a,b,a+b\) (each read separately) are disjoint, and each has unconditional measure \(p^{-k}(1-p^{-1})\) in the full residue square. Distinct primes are independent in this fixed-modulus limit by CRT. Primitivity at other primes multiplies by \(\prod(1-p^{-2})\); this can be justified by first imposing finitely many primes and bounding the fraction with a common prime \(p>Y\) by \(\sum_{p>Y}p^{-2}\) before dividing by the limiting primitive density \(1/\zeta(2)\). Independence is **between primes**, not between the three legs at a fixed prime.

The limit law of \(D\) therefore has independent local exponents. To justify passing from finite primes to the entire integer \(D\), use the uniform bound \(O(1/p^j)\) on the primitive probability that \(p^j\) divides any leg for \(p^j\le 2X\). For the sum leg, after fixing \(a\) there are at most \(X/p^j+1\le3X/p^j\) choices for \(b\); for the other legs it is simpler. Hence the probability that **any** prime \(p>Y\) divides \(D\) is \(O(\sum_{p>Y}p^{-2})\), uniformly. Also each fixed exponent has the displayed geometric tail, so no mass escapes. Writing \(D_\infty\) for the limit gives the explicit law
\[
\Pr(D_\infty=d)=\frac{C}{d}\prod_{p\mid d}\frac{3(p-1)}{p^2+p-3},
\qquad C=\prod_p\left(1-\frac{3}{p(p+1)}\right).
\]
The atom \(D=1\) is the squarefree-density case (compare [Poonen's squarefree-value framework](https://math.mit.edu/~poonen/papers/sqfree.pdf)); the local calculation here tracks *all* exponents beyond the first. The product for \(C\) is positive by the summability of its local losses. The replay estimate is \(C\approx0.29262\).

### Two ways to read the boundary

- **Logarithmic mean.** The local mean of \((v_p(abc)-1)_+\) is the sum of its tail \(3/(p^2-1)\). The uniform bound above lets us sum \(\log p\) times the tails over all \(p,j\ge2\), since \(\sum_p \log p\sum_{j\ge2}p^{-j}<\infty\). Thus the finite-box mean itself converges:
  \[
  \lim_{X\to\infty}\mathbb E\log D
        =3\sum_p\frac{\log p}{p^2-1}=-3\frac{\zeta'(2)}{\zeta(2)}
        \ \approx1.7099 .
  \]
  The primitive pair, scaled by \(X\), tends to uniform area on \([0,1]^2\) by the same CRT tail bound; the endpoint logarithms are uniformly integrable by bounding the fraction in \(a/X\le\delta\) or \(b/X\le\delta\) by \(O(\delta)\). Integrating \(\log u+\log v+\log(u+v)\) gives \(2\log2-7/2\). Hence
  \[
  \mathbb E\log R=3\log X+2\log2-\tfrac72+3\frac{\zeta'(2)}{\zeta(2)}+o(1).
  \]
  This is a literal decomposition of the mean into an additive area contribution and a prime contribution from a logarithmic derivative of an Euler product.
- **Harmonic tail.** In the limit \(\mathbb E D_\infty^u<\infty\) for \(0\le u<1\), and diverges for \(u\ge1\) already at the prime 2. For \(u<1\), each local moment differs from 1 by \(O_u(p^{u-2})\), so this is a convergent Euler product. The limiting tail is \(t^{-1+o(1)}\), but it is not bounded by \(t^{-1}\) times a fixed power of \(\log t\). Indeed fix 2 and any \(r\) other primes, vary each of the \(r\) positive odd-prime exponents in an interval of size proportional to \(\log t\) small enough that their product is at most \(\sqrt t\), and choose the power of 2 putting \(d\) in \([t,2t)\). Each of the \(\gg_r(\log t)^r\) distinct integers has probability \(\gg_r1/t\) in the explicit limit law. The upper exponent follows by taking moments \(u\) arbitrarily close to 1.

## 3. What this does and does not say about abc

In the uniform box the quality \(q=\log c/\log R\) tends in probability to \(1/3\): scaled \(a,b\) avoid tiny corners with arbitrarily high probability, and \(\log D=o(\log X)\) in probability by the bounded mean. This says nothing pointwise about exceptionally large \(q\), or about one small summand next to a very large one.

There is a modest **uniform** estimate that does not exchange the cutoff and tail. For \(0<u<1/3\), let \(d(n)=n/\operatorname{rad}(n)\). The one-variable density bound \(\sum_{n\le Y}d(n)^v\ll_v Y\) for \(0<v<1\) follows by expanding the nonnegative divisor sum \(d(n)^v=\sum_{e\mid n}g_v(e)\). Here \(g_v(p)=0\), \(g_v(p^k)=p^{v(k-1)}-p^{v(k-2)}\) for \(k\ge2\), and \(\sum_e g_v(e)/e\) converges. Hölder with exponent 3 applied to \(d(a)^u,d(b)^u,d(a+b)^u\) now gives \(\mathbb E_X D^u\ll_u1\); the multiplicity of a sum \(a+b\) is at most \(X\), and the primitive denominator is \(\asymp X^2\). For fixed \(0<\eta<1\), among \(\eta X\le a,b\le X\) the number with \(q\ge q_0>1/3\) is therefore
\[
 O_{\eta,q_0,\epsilon}\!\left(X^{1+1/(3q_0)+\epsilon}\right).
\]
Indeed that quality forces \(D\gg_{\eta,q_0}X^{3-1/q_0}\), and we take \(u\) just below \(1/3\). The sharper \(u<1\) moment of the limiting distribution cannot simply be plugged into this growing-threshold problem. This is the same distinction as the repo's [conditional prime-gluing work](../conditional_determinant_vitali_20261116/CONDITIONAL.md): a multiplicative baseline or a density does not provide the missing uniform control on a particular additive fiber.

`node replay.js` checks the local counts for five primes and three powers, and compares nine atomic masses and the two logarithmic averages in finite coprime boxes.

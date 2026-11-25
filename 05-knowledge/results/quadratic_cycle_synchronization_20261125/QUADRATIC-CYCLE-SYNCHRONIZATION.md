# Quadratic cycles with one parameter: the conductor 7 norm and the period-six gate

25 November 2026. A transfer from the attached [matrix multiplication preprint](https://mathcommunity.github.io/oai-math-revised/papers/The-matrix-multiplication-exponent-is-nine-fourths-November-24-2026/paper.pdf) into the repo's [Farey/Mandelbrot work](../farey_mandelbrot_20261019_collatz_angles_and_cross_conjecture_dictionary.md), [parabolic clocks](../parabolic_scaling_session_20261019.md), and [local-versus-exact tournament checks](../moore_fixed_order_spectra_20261124/RADIUS-TWO-PRIME-BALANCE.md). The two results are a sharp normalized gap at the period-three real root, and a finite-prime reduction of Poonen's question when the parameter denominator has at most five prime factors. No uniform polynomial point bound is claimed.

## 1. Read the attachment as a lesson about synchronization, not an exponent to import

Sections 1–2 and Figure 2 of the preprint combine two distinct estimates: under a hypothetical exponent strictly below \(9/4\), tensor powers compress, yet independently tagged distinctness instances demand more volume. The critical tripod in Proposition 8.1 compares three edge labels with one common rank label; in the equal-arm normalization the comparison is \(6/(8/3)=9/4\). Theorem 12.1 is about a different but related synchronization: under that hypothesis, the generic assisted price for polynomial tensor families with independent parameter vectors agrees with the price on a common parameter. Parts II–IV then have to preserve the identities of the records through stable sorting, lookup and concentration of the network wiring. None of those theorems computes rational periodic points or supplies an arithmetic point bound.

What transfers is a question: **when is it safe to forget which parameter a local witness belongs to?** For this family it is generally not safe. Write \(f_z(x)=x^2-z\), \(z=A/d^2\) reduced. All rational finite-orbit points have reduced denominator \(d\): at a bad prime an orbit valuation must match half the parameter valuation, or squaring drives it off to infinity; at a good prime a negative starting valuation does the same. Thus on numerator sites \(X,Y\), the edge records the exact common tag \(A=X^2-dY\). In the real Mandelbrot slice \(z\in[-1/4,2]\), every finite orbit lies in \([-2,2]\), so a given parameter has a finite labelled verification graph. The transfer matrix \(M_A(X,Y)=[X^2-dY=A]\) is the matrix of a partial function in this graph. It can be pruned by deleting any vertex with no surviving successor, or powered by following pointers without generic matrix multiplication.

For a finite parameter bank, \(\operatorname{tr}((\bigoplus_A M_A)^n)=\sum_A\operatorname{tr}(M_A^n)\) counts fixed points of the \(n\)-th iterate. Möbius inversion gives exact periods. Replacing the direct sum by a superposition before taking powers instead allows the parameter to change at each visit. On denominator \(4\), for example, the numerator walk
\[
-5\to-1\to-3\to1\to-5
\]
uses four **different** tags \((29,13,5,21)\). All four parameters lie on the real Mandelbrot slice; no one of their maps has period four. The error resembles keeping only local shell or phase readings in a tournament, then treating those readings as an exact compatible construction. An exact CRT check for an edge may reconstruct its tag (the absolute residual is at most \(25d^2/4\) in the box), but every edge must still carry the *same* reconstructed tag.

## 2. Complete integer layer: square sparsity, not absence at the next three periods

For distinct periodic points \(x,y\), with a common return time \(L\),
\[
1=\prod_{i=0}^{L-1}\big(f_z^i(x)+f_z^i(y)\big).
\]
Indeed, divide \(f_z(u)-f_z(v)=(u-v)(u+v)\) along the two orbits and telescope. This is the pairwise form of the [Eliahou–Fares Vandermonde argument](https://arxiv.org/abs/1903.08865). For integers, every factor must be \(+1\) or \(-1\). Every pair of distinct periodic points therefore has opposite parity; at most two can exist. Equivalently, on the periodic set \(f_z\) is an isometry of the integer line. For integral \(z\), the denominator rule shows the integer and rational affine finite-orbit sets coincide.

- **Fixed points:** precisely \(z=m(m+1)\), \(m\ge0\), with fixed sites \(-m,m+1\).
- **Exact 2-cycles:** precisely \(z=m(m+1)+1\), \(m\ge0\), with sites \(m,-m-1\). Use \(f_z^2(x)-x=(x^2-x-z)(x^2+x+1-z)\); the factors cannot share an integral root.
- Their sign mates are the only other finite sites except \(0\) at \(z=2\). For the fixed case, a further preimage would require \(m(m+2)=(m+1)^2-1\) or \(m^2-1\) to be a square; only the named degeneracies yield one. For the 2-cycle case the further squares would be \(m^2+1\) or \((m+1)^2+1\); the sole zero-\(m\) exception just repeats the existing cycle. All other integer parameters have no finite rational orbit.

| \(z\) | Finite cycles | Tails |
|---:|---|---|
| 0 | \(0,1\) fixed | \(-1\to1\) |
| 1 | \(0\leftrightarrow-1\) | \(1\to0\) |
| 2 | \(-1,2\) fixed | \(1\to-1,\ 0\to-2\to2\) |
| 3 | \(1\leftrightarrow-2\) | \(-1\to-2,\ 2\to1\) |
| 4,5 | none | none |
| 6 | \(-2,3\) fixed | \(2\to-2,\ -3\to3\) |
| 7 | \(2\leftrightarrow-3\) | \(-2\to-3,\ 3\to2\) |

This explains the stated parameter gaps at \(3,4,5\) for fixed points and \(4,5,6\) for the 2-cycle. The Mandelbrot condition is a different question: whether the *critical orbit at zero* is bounded for complex \(c=-z\). The real slice is \(z\in[-1/4,2]\); \(z=3\) is outside even though its rational 2-cycle exists.

## 3. Three rational gaps, one shared parameter, and the field of conductor 7

Let \(\alpha<\beta<\gamma\) be a real 3-cycle. Its cyclic order is \(\alpha\to\gamma\to\beta\to\alpha\): in the reverse order the interior point \(\beta\) would need the largest absolute value, since \(f_z(\beta)\) would be largest. Put \(t=(\beta-\alpha)/(\gamma-\beta)>0\). Taking differences of the three equations with **the same** \(z\) gives
\[
 \alpha+\beta=-1-1/t,\qquad
 \beta+\gamma=t,\qquad
 \alpha+\gamma=-1/(t+1).
\]
The factors multiply to 1, as required. The three sum values pass through the Möbius rotation \(t\mapsto-1/(1+t)\mapsto-(1+t)/t\). Conversely solving these sums gives a rational 3-cycle for each rational \(t>0\): the two adjacent differences are positive and proportional to \(t:1\), and all three reconstructed edge parameters agree. This is an ordered-gap presentation of the known [Morton–Raianu parametrization](https://arxiv.org/abs/2105.07435); see also the [consecutive-sum and S-unit formulation](https://borowiecki.dev/pdf/2610.39820) of Käenmäki–Kosunen.

Write \(t=a/b\) with coprime positive integers. Calculating \(z=\alpha^2-\gamma\) in this chart yields
\[
d=2ab(a+b),\qquad F=a^3+a^2b-2ab^2-b^3,\qquad
\boxed{z=\frac74+\frac{F^2}{d^2}} .
\]
The denominator is exactly \(d^2\): \(4\mid d\), \(F\) is odd, and reducing \(F\) modulo a prime dividing \(a,b,a+b\) gives a unit. This is a norm identity for the *detuning* from \(7/4\), not the distinct norms for the parameter numerator tabulated by Morton–Raianu.

Now \(F=N_{\mathbb Q(2\cos(2\pi/7))/\mathbb Q}(a-b\,2\cos(2\pi/7))\). At the parabolic root \(z=7/4\), the gap coordinate \(t_0=2\cos(2\pi/7)\) obeys
\[
t_0^2-2=-1/(1+t_0).
\]
So \(1\mapsto2\mapsto4\bmod 7\) on cyclotomic angles is exactly the order-three rotation of gap coordinates in the parabolic fiber, via \(\theta\mapsto 2\cos(2\pi\theta)\). This refines the prior Farey/Collatz connection without equating external ray angles with rational orbit points. In particular the real parabolic period-three root has no rational 3-cycle: the cubic for \(t_0\) has no rational root. For \(z=29/16\) the rational cycle is **repelling**, with multiplier \(35/8\); it is not the critical orbit.

For primitive \((a,b)\), a prime \(p\ne7\) dividing \(F\) must be \(\pm1\bmod7\): the roots modulo \(p\) are \(\xi^j+\xi^{-j}\) for a primitive seventh root \(\xi\), and Frobenius fixes one exactly when \(pj=\pm j\bmod7\). At 7 the cubic is \((t-2)^3\bmod7\); for \(a=2b+7k\) the norm is \(7b^3\bmod49\), so \(v_7(F)=1\) whenever \(7\mid F\).

## 4. Two uniqueness statements on the real slice

A real 3-cycle with \(z\le2\) must satisfy \(1/2<t<2\). If \(r=1/t\ge2\), then \(-2\alpha=1+r+1/r+r/(1+r)>4\), outside the escape bound. For \(t\ge2\), the numerator of \(z-2\), over denominator \(4t^2(t+1)^2\), is
\[
(t^2-4)t^3(t+2)+t^2+4t+1>0.
\]
Consequently if \(z\) has at most two denominator primes, the pairwise coprime numbers \(a,b,a+b\) cannot all exceed 1. If one of \(a,b\) is 1, the interval \(1/2<t<2\) forces both to be 1. **Thus \(29/16\) is the only real Mandelbrot parameter with at most two denominator primes and an exact rational 3-cycle.** By the Eliahou–Fares theorem it is also the only such parameter supporting any rational period greater than two. The previously known power-of-two uniqueness is a special case. Every other real parameter with a rational 3-cycle needs at least three distinct denominator primes.

There is a different uniqueness relative to the parabolic root. With \(d^2=\operatorname{den}(z)\), the smallest possible normalized detuning \(d^2(z-7/4)\) is **1**, attained at exactly two parameters. The positive primitive solutions of \(F=\pm1\) from [Thomas's simplest-cubic classification, or Hoshi's Table 1](https://arxiv.org/pdf/0810.3374), are \((1,1),(5,4)\). The primitive solutions to \(F=\pm7\) in the same table give just \((2,1)\) positive, which fails the real-slice bound \(t<2\). Since the next possible prime divisor is 13, the next possible normalized detuning is at least \(13^2\), attained at \((3,2)\).

| \(t\) | \(z\) | cycle numerators in directed order | \(d\) | \(d^2(z-7/4)\) | affine rational finite-orbit points |
|---|---|---|---:|---:|---:|
| \(1\) | \(29/16\) | \(-7,5,-1\) | 4 | 1 | 8 |
| \(5/4\) | \(226801/129600\) | \(-629,469,-19\) | 360 | 1 | 6 |
| \(3/2\) | \(6469/3600\) | \(-107,83,7\) | 60 | 169 | 6 |

The middle parameter is in the Morton–Raianu table; the observation here is its exact parabolic gap and uniqueness alongside \(29/16\) at normalized gap one. These are not the two closest rational 3-cycle parameters by absolute distance. Rational approximations \(a/b\) to \(t_0\) yield parameters arbitrarily close from above; for convergents, the \(O(b^{-2})\) approximation gives \(z-7/4=O(d^{-4/3})\), because \(d\) grows as \(b^3\) and the detuning is quadratic in the ratio error.

## 5. A period-six gate and what a 12-point bound would actually imply

The absence of rational cycles of periods [4 (Morton)](https://matwbn.icm.edu.pl/ksiazki/aa/aa62/aa6243.pdf) and [5 (Flynn–Poonen–Schaefer)](https://arxiv.org/abs/math/9508211) is uniform in the parameter; it does **not** follow from the integer gaps or one example at \(29/16\). The 3-cycle assertion alone also does not prove the classification conjecture. [Poonen's classification](https://arxiv.org/abs/math/9512217), under the absence of all higher periods, gives **12 possible graphs** and **at most nine points on \(\mathbb P^1(\mathbb Q)\)**, counting infinity. Twelve points for arbitrary quadratic *rational functions* is a separate issue: [Blanc–Canci–Elkies](https://arxiv.org/abs/1305.1054) exhibit examples with at least fourteen.

Here \(\omega(d)\) counts distinct prime divisors of \(d\). The following is a narrower unconditional reduction beyond the \(\omega(d)\le2\) theorem. **If \(\omega(d)\le5\), the only possible rational period above three is 6, and period 6 requires \(5\mid d\).** If \(\omega(d)\le6\), the same filter leaves only 6 and 8 above three; period 8 forces exactly the prime support \(\{2,3,7,11,19,23\}\) in this six-prime class.

1. For \(n\ge3\), Eliahou–Fares give \(4\mid d\) and \(n\le2^{\omega(d)}+2\). If \(3\nmid d\), their Lemma 3.13 forces the whole cycle to one nonzero residue mod 3. Zieve's lifting theorem (their Theorem 2.11) gives \(n\in\{1,2,3,6\}\) in that case, because the reduced period is one, the multiplier order divides two and the extra 3-power has exponent at most one.
2. If \(p\ge5\) is a good prime, the same lifting theorem allows only \(m\) or \(mr\), where \(m\) is the reduced period and \(r\) is the order of the reduced multiplier, with no \(p\)-power enlargement. If the multiplier is zero, only \(m\) is allowed. Exhaust the \(p\) maps \(X^2-Z\) and their cycles over \(\mathbb F_p\), keeping each multiplier with its *own* cycle and parameter. For \(7\le n\le34\), the replay finds at least four blocking primes in \(\{5,7,11,13,17,19,23,29,31,37,41,43,47\}\) for each \(n\); period 8 has exactly four, \(7,11,19,23\). Together with 2 and 3, this requires at least six primes. With the upper endpoint extended to \(66\), at least five blocking primes occur for every \(7\le n\le66\) except 8, so only 8 survives at support six. For 6 the table shows \(p=5\) blocks all parameters of good reduction.
3. A rational 6-cycle already gives **13 projective preperiodic points**: its six negatives are distinct nonperiodic mates (a permutation cannot send two different periodic points to the same image), and infinity is fixed. Thus a 12-point bound including infinity would exclude the remaining branch for the five-prime class. Conversely, ruling out 6 there gives Poonen's nine-point conclusion in that class. This leaves the same period-six problem studied conditionally by [Stoll](https://arxiv.org/abs/0803.2836), now restricted to a denominator budget if the target is that class.

Replay: **python period3_sync_check.py**. The good-prime table certifies only incompatibility with good reduction, not lifts from a parameter that passes the screen.

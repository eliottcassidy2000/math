# The common-center capacity obstruction for AP-paired square-sum constructions

**Status: PROVED in the stated research scope; not a canon promotion.** The maximal-order sharpness uses the **CITED** prime number theorem in the fixed progression `1 mod 4`. All finite controls are **FINITE-EXACT**. No contact-chain upper bound or novelty claim is asserted.

## 1. Statement: polynomial numerical labels are impossible

Let `k` distinct positive AP square-sum triangles share one center `M`, and let `N` be their largest vertex label. This includes every common-center AP-paired chain, fan, and contact tree from the inherited construction, whether or not additional inter-petal edges exist.

**Theorem 1.** Necessarily `M=2H²` for an integer `H>=1`. Put

\[
 D(H)=\prod_{p\equiv1\pmod4}(2v_p(H)+1).
\]

The number of possible primitive unordered Pythagorean petals before the positivity condition is imposed is exactly

\[
 U(H)=\frac{D(H)-1}{2}.                                      \tag{1}
\]

Consequently

\[
 k\le U(H),\qquad 2H^2<N<4H^2,                              \tag{2}
\]

and, uniformly as `k` tends to infinity,

\[
 \boxed{\quad
 \log N\ \ge\ \left(\frac2{\log3}-o(1)\right)
                    \log k\,\log\log k.
 \quad}                                                     \tag{3}
\]

In particular, no fixed constants `A,B` can bound these labels by `N<=B k^A` for arbitrarily large `k`. This is an obstruction to **numerical** height; it does not prohibit polynomial output-bit complexity.

The coefficient `2/log 3` is sharp for the capacity problem alone: there are common centers supporting `k_j -> infinity` positive petals for which equality in (3) holds asymptotically. One may also require any fixed finite collection of admissible primitive petals to occur. This assertion makes no claim that the available petals form a connected contact graph, much less a compatible chain.

## 2. Inheritance and exact denominator accounting

The closest proved mechanism is the complete AP/PPT classification in [ap_square_pythagorean_20261019.md](ap_square_pythagorean_20261019.md), Theorem 1. Every positive AP square-sum triangle is uniquely

\[
 2g^2(C^2-4ab,\ C^2,\ C^2+4ab),\qquad C^2>4ab,
                                                               \tag{4}
\]

for a primitive PPT with unordered positive legs and integer `g>=1`. Hence a shared center forces `H=gC`, so `C|H`; conversely every admissible primitive PPT with `C|H` gives one such triangle. There are no additional scale choices after `H` and the primitive petal are fixed.

The canonical hostile example is that a finite atlas bounds each successive denominator in terms of its predecessor, but does not control their common multiple. The corrected near miss is to seek a polynomial numerical bound in the growing number of petals; (3) excludes it even before contact compatibility is tested. The least-used relevant sidecar is the prime-support/angle separation in [endpoint_pythagorean_fibers_20261019.md](endpoint_pythagorean_fibers_20261019.md), Theorems 1–2.

The live board is: common-center scaling; Gaussian prime allocation; angular positivity; signed-contact paths; rational-map height. Its new distinction is between **capacity** (which petals can fit one denominator) and **connectivity** (which of those petals can be joined through the required ports).

The cited internal Gaussian-fiber theorem gives exactly `2^(omega(C)-1)` primitive unordered positive PPTs for a hypotenuse `C>1` supported only on primes `1 mod 4`, and none otherwise. Summing over divisors of `H` gives

\[
\begin{aligned}
 2U(H)+1
 &=\sum_{\substack{C|H\\p|C\Rightarrow p\equiv1(4)}}2^{\omega(C)}\\
 &=\prod_{p\equiv1(4)}(1+2v_p(H))=D(H).
\end{aligned}
\]

This proves (1). Positive petals form a subset, which proves the first part of (2). Formula (4), with `0<4ab<C²`, proves the two strict label inequalities.

## 3. An elementary maximal-order upper bound

No external divisor-function estimate is needed. Write `L=log H`. For any real cutoff `y>=2`, split the product defining `D(H)` at `y`. Every exponent is at most `L/log 5`, and there are at most `floor(y)` small primes. For every integer `e>=1`,

\[
 2e+1\le3^e.
\]

Therefore the following explicit bound holds:

\[
 \log D(H)\le
 \lfloor y\rfloor\log\left(1+\frac{2L}{\log5}\right)
       +\frac{\log3}{\log y}L.                            \tag{5}
\]

Indeed the sum of the exponents on primes greater than `y` is at most `L/log y`. For sufficiently large `L`, choose

\[
 y=\frac{L}{(\log L)^3}.
\]

The first term in (5) is `O(L/(log L)²)`, and the second is

\[
 (\log3+o(1))\frac{L}{\log L}.
\]

Thus

\[
 \log(2k+1)\le\log D(H)
       \le(\log3+o(1))\frac{\log H}{\log\log H}.          \tag{6}
\]

Also `D(H)<=H`, since `2e+1<=p^e` for each participating prime `p>=5`. Thus `log H>=log(2k+1)`, allowing `log log H>=log log(2k+1)` in the inversion of (6). It follows that

\[
 \log H\ge(1/\log3-o(1))\log k\log\log k.
\]

Since `N>2H²`, this proves (3). The `o(1)` is uniform over every admissible family; the argument used only `H` and `k`.

## 4. Sharpness for capacity, with positivity retained

The angular part of this proof is an extension of the inherited single-prime equidistribution argument, with uniformity in the other-prime phase supplied explicitly here. The only analytic input beyond that argument is the prime number theorem for `p=1 mod 4`: see J.-H. Evertse, [Chapter 6, Theorem 6.2.1](https://pub.math.leidenuniv.nl/~evertsejh/ant20-6.pdf), with fixed modulus `4` and residue `1`. It implies that the `j`th such prime is asymptotic to `2j log j`, and hence that the logarithm of the product of the first `s` of them is asymptotic to `s log s`.

Fix a split prime, say `p=5`, and a Gaussian prime `pi` above it. Let `p_1,p_2,...` be the other primes `1 mod 4`, in increasing order. Set

\[
 m_s=\lfloor\sqrt s\rfloor,\qquad
 H_s=5^{m_s}\prod_{j=1}^s p_j.
                                                               \tag{7}
\]

Then

\[
 D(H_s)=(2m_s+1)3^s,\qquad \log H_s\sim s\log s.        \tag{8}
\]

### Gaussian exponent vectors

For a general split-prime-supported `H=prod p^e`, choose a Gaussian prime `pi_p` above each `p`. Exponent vectors `j_p in {-e,...,e}` give rational unit-circle points

\[
 z_{\boldsymbol j}
       =\prod_{p|H}\left(\frac{\pi_p}{\overline{\pi_p}}\right)^{j_p}.
                                                               \tag{9}
\]

Their primitive hypotenuses divide `H`. Unique factorization makes the vectors distinct; nonzero vectors have neither coordinate zero. After taking absolute values of the odd and even legs, precisely the pair `j,-j` gives one unordered primitive petal. There are no other identifications: an equality up to conjugation and a Gaussian unit forces the same prime valuations, hence the same vector up to sign. The zero vector is the excluded degenerate petal.

For `z=exp(i theta)`, the positivity condition is

\[
 |\sin(2\theta)|<1/2.                                  \tag{10}
\]

Its accepted set on the circle is a finite union of intervals of normalized measure `1/3`. A nonzero exponent vector cannot lie at the equality boundary, since it would give `C²=4ab` with `C` odd.

### Uniformity in an arbitrary phase

Let `alpha=arg(pi/bar(pi))`. As proved in the inherited note, `alpha/(2pi)` is irrational: otherwise `pi/bar(pi)` would be a root of unity in `Q(i)`, inconsistent with its nonzero Gaussian-prime valuations. For each nonzero integer frequency `l`,

\[
 \left|\frac1{2m+1}\sum_{j=-m}^m
       e^{il(\phi+j\alpha)}\right|
 \le\frac{2}{(2m+1)|1-e^{il\alpha}|}\longrightarrow0,    \tag{11}
\]

uniformly in the phase `phi`. Approximate the indicator of (10) above and below by fixed trigonometric polynomials with arbitrarily small integral gap. Equation (11) proves

\[
 \sup_\phi\left|
 \frac{\#\{-m\le j\le m:\phi+j\alpha\text{ satisfies (10)}\}}
      {2m+1}-\frac13\right|\longrightarrow0.            \tag{12}
\]

This is uniform in the number, size, and phase distribution of all the other primes. Apply (12) to each of the `3^s` choices of their exponent coordinates in (7). After deleting the zero vector and quotienting by conjugation, the number `K(H_s)` of positive petals satisfies

\[
 K(H_s)=(1/6+o(1))(2m_s+1)3^s.                         \tag{13}
\]

In particular,

\[
 \log K(H_s)\sim s\log3,\qquad
 \frac{\log K(H_s)\log\log H_s}{\log H_s}\longrightarrow\log3.
                                                               \tag{14}
\]

Take all these positive petals at common center `2H_s²`; their maximum label `N_s` obeys `2H_s²<N_s<4H_s²`. Equations (8), (13), and (14) show

\[
 \frac{\log N_s}{\log K(H_s)\log\log K(H_s)}
       \longrightarrow\frac2{\log3}.                  \tag{15}
\]

This proves the claimed capacity sharpness. To include any fixed finite collection of admissible primitive petals, multiply `H_s` by their fixed common hypotenuse multiple. The prime support and exponent changes outside `5` affect (8) by only a fixed multiplicative factor and an additive constant in `log H_s`; the `5` exponent still tends to infinity. The same uniform phase argument and asymptotic coefficient apply.

**Boundary:** the construction only supplies a large stock of petals with shared denominator. It supplies no inter-petal square contacts. In particular, (15) is not an upper bound for a chain or Hamiltonian cycle.

## 5. Exact finite controls and a reusable optimization

The companion script uses two independent coordinate routes: direct square-root enumeration of AP labels and Euclidean PPT enumeration. It also checks the Gaussian product capacity against the complete Euclidean divisor sieve for every `H<=100000`. No existing contact census is rerun.

Examples of `(H, unfiltered capacity, positive capacity)` are:

- `(41,1,1)`;
- `(145,4,2)`;
- `(1105,13,5)`;
- `(32045,40,14)`.

The largest positive capacity through `H=100000` is `15`, first attained at `H=69745`. This is a finite fact, not an asymptotic extrapolation and not a contact-path count.

The script also exactly minimizes `H` under the necessary constraint `D(H)>=2k+1` for selected `k` through `1024`. The optimizer uses increasing split primes with nonincreasing positive exponents. This search is complete because removing other prime factors and rearranging exponents in that order cannot increase `H` and does not reduce `D(H)`. A finite product of the first split primes gives an initial upper bound; recursive branches with product already at that bound are discarded. Small outputs are independently checked by exhaustive search.

For example `k=1024` forces `H>=2576450045` and therefore `N>13276189668761004050`, even when positivity and every contact requirement are relaxed. This is a necessary lower bound, not a construction at that scale.

Reproduce from the repository root after integration:

```bash
python3 04-computation/experiments/square_sum_common_center_capacity_20261019.py
python3 -O 04-computation/experiments/square_sum_common_center_capacity_20261019.py
```

The universal proof is independent of this finite computation.

## 6. Transfer, distinctions, and the revised target

Source: integer square-sum triangles with a shared AP center. Target: a divisor box of Gaussian prime-allocation vectors. Map: classify each triangle, retain its primitive hypotenuse divisor, and encode that divisor's Gaussian sign choices. Preserved: distinct petals and their exact common-center scale. Lost: angular positivity and contact edges. Sidecars: the angle inequality (10) and the signed contact predicate `1+epsilon f_T+eta f_U in Q²`. The transfer is the lower bound (3); the hostile test is that equal capacities can have different positive counts.

For arbitrary tree or path square-sum labeling without the AP/common-center constraint, this obstruction does not apply. It also does not address graceful differences, consecutive interval coverage, Collatz trajectories, or a claimed equivalence with any of them.

The correct long-chain target is now a **bit-length** theorem, or a numerical bound of at least superpolynomial size in `k`. At the capacity scale the optimal lower order is `log k log log k`; actual contact compatibility may force substantially larger labels. The remaining task is to construct many compatible contacts while controlling the prime support and exponent growth of their common denominator. Repeated one-step height bounds do not provide that control.

# Prime powers, biased races, and the limits of finite evidence

Research synthesis, 14 October 2026. Existing repository theorem files were not edited; no claim of novelty is made for the classical identities. The identities below are derived directly; historical and analytic results are separately attributed. Repository comparisons use the available local snapshot `7435450cc4d104002588215451e5cdd21544b3d3`, not an asserted latest remote state.

## 1. The strongest defensible synthesis

The suggested mechanism is real: removing repeated prime contributions from the natural Euler-product counting function creates a deterministic square-scale correction. In the mod-4 prime race this correction favors primes congruent to 3 over primes congruent to 1 by an amount asymptotic to `sqrt(x)/log(x)`.

Three qualifications are essential:

1. This describes a term in an exact decomposition, not an asymptotic formula for the entire prime-count difference.
2. The residual oscillation comes from zeros of the relevant zeta or L-function and can defeat that drift at arbitrarily large scales. There is no final transition to a permanently reversed regime.
3. The counting weights matter: `log(p)` at `p^k` for Chebyshev's `psi`, or `1/k` at `p^k` for Riemann's `J`/`Pi`. Simply counting every prime power once does not remove the drift.

A sound replacement for “a finite microcosmic region followed by an infinite macrocosmic region” is:

> A long finite prefix can display a highly persistent bias. The same exact formula continues to govern arbitrarily large inputs: a structural drift and an oscillatory component coexist, and rare phase alignments overturn the empirical pattern. A first failure is not a change in the governing law.

This is a promising research principle for the corpus: preserve multiplicity, labels, source constraints, and the actual observable before making an asymptotic inference. It does not by itself prove a correspondence between number-theoretic zeros and the corpus's geometric or dynamical operators.

There are two different *scales*, rather than two successive regions: the prime number theorem in progressions gives `pi_1(x)/pi_3(x) -> 1` at the main scale, while the smaller signed difference can remain statistically biased and change sign infinitely often. Relative equality and a biased fluctuation law are compatible.

## 2. Exactly what 26,861 and 1,473 mean

Write `pi_a(x) = #{p <= x : p prime, p = a mod 4}`. John Leech's original 1957 paper reports the first strict reversal:

`pi_1(26861) = 1473`, while `pi_3(26861) = 1472`.

It also discusses the connection to zero-phase alignment near later reversals. [Leech, *Note on the Distribution of Prime Numbers*, pp. 56–58, especially p. 57](https://personal.math.ubc.ca/~gerg/teaching/592-Fall2018/papers/1957.Leech.pdf), [publisher record](https://doi.org/10.1112/jlms/s1-32.1.56).

An independent Eratosthenes sieve run for this note confirms:

| x | pi_1(x) | pi_3(x) |
|---:|---:|---:|
| 26,833 | 1,471 | 1,471 |
| 26,849 | 1,472 | 1,472 |
| 26,861 | 1,473 | 1,472 |
| 26,863 | 1,473 | 1,473 |
| 26,879 | 1,473 | 1,474 |

Thus the first strict lead lasts for real `x` in `[26861,26863)`, or the two integer cutoffs 26,861 and 26,862. The 1,473rd prime of the first class precedes the 1,473rd of the second class, 26,863. It is not a new permanent ordering.

## 3. Why prime powers enter naturally

For `Re(s)>1`, Euler's product gives the absolutely convergent identities

`zeta(s) = product_p (1-p^(-s))^(-1)`,

`log zeta(s) = sum_p sum_(k>=1) p^(-ks)/k`,

`-zeta'(s)/zeta(s) = sum_p sum_(k>=1) log(p) p^(-ks)`.

Define

`J(x) = sum_(p^k<=x) 1/k`,

`theta(x) = sum_(p<=x) log p`,

`psi(x) = sum_(p^k<=x) log p = sum_(n<=x) Lambda(n)`.

The exact finite identities are

`J(x) = sum_(k>=1) pi(x^(1/k))/k`,

`psi(x) = sum_(k>=1) theta(x^(1/k))`.

All these sums stop once `2^k>x`. Möbius inversion gives, also exactly,

`pi(x) = sum_(k>=1) mu(k) J(x^(1/k))/k`,

`theta(x) = sum_(k>=1) mu(k) psi(x^(1/k))`.

For example, substituting the first formula for `J` into the inversion gives coefficient `(1/m) sum_(k|m) mu(k)`, which vanishes unless `m=1`. This is a proof, not an asymptotic analogy.

The Stieltjes integral identities

`J(x) = integral_(2-)^x (1/log t) d psi(t)`,

`pi(x) = integral_(2-)^x (1/log t) d theta(t)`

explain the weight `1/k`: the jump at `t=p^k` becomes `log p / log(p^k) = 1/k`.

Prime powers carry no independent information beyond the primes that generate them. Their advantage is that these weights linearize the Euler product.

## 4. The exact mod-4 correction

Set

`Delta_pi(x) = pi_3(x)-pi_1(x)`,

`J_a(x) = sum_(p^k<=x, p^k=a mod4) 1/k`,

`Delta_J(x) = J_3(x)-J_1(x)`,

`pi_odd(y) = #{odd primes <= y}`.

For odd `k`, an odd prime and its `k`th power have the same class modulo 4. For even `k`, every odd prime power is 1 modulo 4. Therefore

`Delta_J(x)`

`= sum_(odd k>=1) Delta_pi(x^(1/k))/k`

`  - sum_(even k>=2) pi_odd(x^(1/k))/k`.

Isolating `k=1,2` yields the exact identity

`Delta_pi(x) = Delta_J(x) + (1/2) pi_odd(sqrt x) + R(x)`,

where

`R(x) = -sum_(odd k>=3) Delta_pi(x^(1/k))/k`

`       +sum_(even k>=4) pi_odd(x^(1/k))/k`.

A crude bound `|R(x)| = O(x^(1/3) log x)` suffices here. Indeed, every count is at most `x^(1/3)`, there are `O(log x)` nonzero exponents, and all weights are at most one. Consequently `R(x)=o(sqrt(x)/log x)`.

By the prime number theorem,

`(1/2) pi_odd(sqrt x) ~ sqrt(x)/log x`.

Hence the unconditional decomposition is

`Delta_pi(x) = Delta_J(x) + (1+o(1)) sqrt(x)/log x`.        (A)

Equation (A) does **not** say `Delta_pi(x) ~ sqrt(x)/log x`: that would require `Delta_J(x)=o(sqrt(x)/log x)`, which is not supplied and is contrary to the oscillatory behavior under discussion.

The shortage of the 1-class relative to the 3-class has square drift `sqrt(x)/log x`. Its shortage relative to half the total odd-prime count has half that drift. These are different reference baselines.

### An exact numerical decomposition at the first reversal

At `x=26861`, there are 37 odd primes at most `sqrt(x)`. Direct rational counting gives

`Delta_J(26861) = -50831/2520`,

`(1/2) pi_odd(sqrt(26861)) = 37/2`,

`R(26861) = 1691/2520`.

Thus

`-50831/2520 + 37/2 + 1691/2520 = -1`.

The residual prime-power imbalance has overcome the square correction at this point. The leading approximation `sqrt(x)/log x` is about 16.0704, while the exact square correction is 18.5; an asymptotic scale need not be numerically exact this early.

### A common but important weighting error

Let `Q_a(x)` count each prime power in class `a` with weight one. Then

`Q_3(x)-Q_1(x) = Delta_pi(x)-pi_odd(sqrt x)+O(x^(1/3) log x)`

`= Delta_J(x)-(1/2)pi_odd(sqrt x)+o(sqrt(x)/log x)`.

Thus the unit-counted prime-power race has the opposite square drift. The phrase “prime powers are balanced” must name the weighting and the sense of balance.

## 5. The same statement with theta and psi

Define `theta_a`, `psi_a` by restricting their defining jumps to the class `a mod4`, and let `Delta_theta=theta_3-theta_1`, `Delta_psi=psi_3-psi_1`. Exactly the same parity decomposition gives

`Delta_theta(x) = Delta_psi(x) + theta_odd(sqrt x) + R_theta(x)`,

`R_theta(x) = -sum_(odd k>=3) Delta_theta(x^(1/k))`

`             +sum_(even k>=4) theta_odd(x^(1/k))`.

The Chebyshev bound `theta(y)=O(y)` gives `R_theta=O(x^(1/3)log x)`, and the prime number theorem gives `theta_odd(sqrt x)~sqrt x`. This is the logarithmically weighted square drift.

The passage back to prime counts is not division by `log x` alone; exact partial summation is

`Delta_pi(x) = Delta_theta(x)/log x`

`              + integral_2^x Delta_theta(t)/(t log^2 t) dt`.

This matters whenever lower-order terms or a rigorous error estimate are required.

## 6. Where the oscillation enters

For `x>1`, use the half-jump convention `psi_0` at prime powers. The von Mangoldt explicit formula is

`psi_0(x) = x - sum_rho x^rho/rho - log(2 pi)`

`           - (1/2) log(1-x^(-2))`,

where nontrivial zeros of `zeta` are summed symmetrically in height, not absolutely. A precise exposition including endpoint and convergence conventions is [Strömbergsson's notes, Theorem 13.1](https://www2.math.uu.se/~ast10761/analtalt08/www_notes.pdf).

For the mod-4 race the relevant character is `chi_4(1)=1`, `chi_4(3)=-1`, and `chi_4(even)=0`. Its Euler product gives

`-L'(s,chi_4)/L(s,chi_4) = sum_n Lambda(n)chi_4(n)n^(-s)`.

The twisted prime-power sum is `psi_1-psi_3=-Delta_psi`. Its explicit formula has zero terms `-sum_(rho of L) x^rho/rho`, plus constant/trivial-zero terms; there is no principal-character pole producing an `x` term. Reversing the sign produces the plus sign for the zero terms in `Delta_psi`.

Under GRH, zeros have `rho=1/2+i gamma`, so after normalization the oscillatory factors are `exp(i gamma log x)`. In the standard truncated/limiting-distribution interpretation, the normalized prime race is represented by

`(log x)/sqrt x * Delta_pi(x)`

`= 1 + 2 Re sum_(gamma>0) exp(i gamma log x)/(1/2+i gamma)`.

The equality here denotes the explicit-formula limiting model, not an absolutely convergent pointwise series with its remainder silently discarded. For a rigorous finite evaluation, truncate and bound the remainder; for distributional statements, use the theorem's prescribed limiting procedure. The constant 1 is the square correction. The frequencies and amplitudes come from the zeros of `L(s,chi_4)`, not from zeta's zeros.

Under GRH plus the relevant linear-independence hypothesis on positive zero ordinates, Rubinstein–Sarnak prove a limiting logarithmic distribution. For the mod-4 race the logarithmic density of `Delta_pi(x)>0` is approximately `0.9959`; the opposite lead has positive density approximately `0.0041`. The density uses

`lim_(X->infinity) [1/log X] integral_2^X 1_(Delta_pi(t)>0) dt/t`,

not uniform sampling of integers up to `X`. The centered prime-power model is symmetric under these assumptions; “centered” does not mean equal at every cutoff. [Rubinstein–Sarnak, *Chebyshev's Bias* (1994)](https://doi.org/10.1080/10586458.1994.10504289), [author-hosted paper](https://www.math.uwaterloo.ca/~mrubinst/publications/Chebyshev.pdf). No unconditional 99.59% claim is intended.

## 7. The general root-multiplicity law

For a fixed modulus `q` and reduced residue `a`, define

`r_q(a) = #{b mod q : b^2=a mod q}`.

Exact exponent bookkeeping gives

`J_a(x;q) = pi_a(x;q)`

`           +(1/2) sum_(b^2=a modq) pi_b(sqrt x;q)`

`           +O_q(x^(1/3)log x)`.

For two reduced residues `a,b`, the prime number theorem in progressions consequently gives

`pi_b(x;q)-pi_a(x;q)`

`= J_b(x;q)-J_a(x;q)`

`  +[(r_q(a)-r_q(b))/phi(q)] sqrt(x)/log x`

`  +o_q(sqrt(x)/log x)`.

For `q=4`, `r_4(1)=2`, `r_4(3)=0`, and `phi(4)=2`, recovering coefficient 1. This exposes the actual multiplicity: the fiber size of the squaring map on the unit group. An unsupported numerical constant has become a structural root count.

Higher powers have their own maps `b -> b^k` and root multiplicities. Squares dominate this correction scale because the next exponent contributes only on the `x^(1/3)` scale before logarithms.

## 8. Skewes: the same square offset, a different zero spectrum

From the total-count identity in Section 3,

`pi(x) = J(x)-(1/2)pi(sqrt x)+O(x^(1/3)log x)`.

Thus

`li(x)-pi(x) = [li(x)-J(x)]`

`              +(1/2)pi(sqrt x)+O(x^(1/3)log x)`.

Again the deterministic term is asymptotic to `sqrt(x)/log x`. The term `li-J` is governed by the nontrivial zeros of `zeta`, rather than `L(s,chi_4)`. Under RH the analogous normalized limiting expression has the same constant 1 but a different collection of frequencies and amplitudes.

Skewes's classical work supplied enormous upper bounds ensuring that a reversal exists below those bounds. These are not numbers after which the reverse inequality must hold forever, and the historical bounds are not the first crossing itself. Littlewood's sign-change theorem gives arbitrarily large values of both signs. [Skewes, *On the Difference pi(x)−li(x), I* (1933)](https://doi.org/10.1112/jlms/s1-8.4.277), [Ford–Hudson's research paper on sign changes](https://www.ford126.web.illinois.edu/wwwpapers/lehman.pdf).

Changing `li` to an offset version differing by a fixed constant affects numerical crossing conventions but not this square-scale analysis or infinite oscillation. Any quoted exact crossover must specify the convention.

## 9. Mertens: related cancellation, not the same square-drift proof

The Mertens function is `M(x)=sum_(n<=x)mu(n)`. Here `mu(n)=0` whenever a prime square divides `n`, and on squarefree numbers its sign is parity of the number of prime factors. Its generating series is

`sum_n mu(n)n^(-s) = 1/zeta(s)` for `Re(s)>1`.

The conjecture `|M(x)|<sqrt(x)` for `x>1` was disproved by Odlyzko and te Riele. Their proof establishes `limsup M(x)/sqrt(x)>1.06` and `liminf M(x)/sqrt(x)<-1.009`, without exhibiting an individual counterexample. This already rules out an eventual regime of only positive or only negative scaled excess. The disproof does not disprove RH. [Original 1985 paper, especially pp. 138–144](https://ir.cwi.nl/pub/1823/1823D.pdf).

There is an elementary stronger consequence relevant to the proposed two-region picture. Since `M(n+1)-M(n)` is one of `-1,0,1`, a passage from positive to negative values must hit zero. Arbitrarily large values of both signs therefore force arbitrarily large zeros as well. The conjectured inequality is violated infinitely often and satisfied infinitely often; even its absolute-value form does not become permanently false past one boundary.

For simple zeta zeros, the zero-residue terms in a properly truncated/smoothed explicit formula are

`x^rho/[rho zeta'(rho)]`.

Thus the oscillation has a different coefficient from the prime-count formula's `x^rho/rho`. The proof exploits simultaneous phase alignment and rigorous smoothing bounds. Removing prime squares from a balanced race is not a proof of the Mertens counterexample.

There is nevertheless an exact algebraic connection: logarithmic differentiation and Dirichlet inversion give `Lambda = mu * log`, namely

`Lambda(n) = sum_(d|n) mu(d) log(n/d)`.

This connects the objects without conflating their summatory inequalities. Higher powers are built into the Euler product and its inverse, but the cancellation statistic has changed.

Later work substantially strengthened the bounds and finite computation; for example Hurst's paper reports complete evaluation up to `10^16` and stronger liminf/limsup bounds. These remain evidence and theorems of different types. [Hurst, *Computations of the Mertens Function and Improved Bounds on the Mertens Conjecture*](https://arxiv.org/abs/1610.08551). No claim that this is the latest record is needed here.

## 10. A concrete, falsifiable bridge to the repository

The relevant common structure is **primitive data plus repeated/multiple contributions, observed through a specified functional**. The following local source comparisons were read directly:

- **THM-1740**: its correction replaces a guessed universal factor 2 with radial multiplicity. Its stated scope distinguishes tested single-straddle cases from the conjectural general bound. The prime-race analogue is `r_q(a)`, the actual fiber multiplicity of the squaring map. This is a useful methodological analogy; these two multiplicities are not identified by a proved map.
- **THM-2839**: prime-power unit mass gives full convolution spectrum. Its own statement distinguishes signed algebraic reconstruction from a positive physical action.
- **THM-2870**: convolution by a mask and physical multiplication by the same mask need not be intertwinable; the kernel/rank obstruction is explicit. This directly cautions against substituting the prime-power observable for the prime observable without the correction identity and its weights.
- **THM-3163**: every finite terminal law has a state-dependent sequential realization, while value/response/current-compatible transport remains additional structure. Matching a long finite distribution is consequently insufficient evidence of the required physical mechanism.

An actionable investigation should retain rows `(primitive object, exponent, source label, destination label, weight)` and test:

1. An exact finite pushforward identity before discarding exponents or labels.
2. The first nontrivial repeated-object correction and its fiber multiplicity.
3. Whether that correction is small at the main scale but leading at the fluctuation scale.
4. The remaining signed/oscillatory term; never infer its smallness from a large finite prefix.
5. A source- and observable-preserving intertwiner if a result is to be transported between problems.

For signed accelerated Collatz maps, valuation parity can be indispensable to output residues. That makes a valuation-labelled audit sensible. It does not import the 99.59% logarithmic prime-race density, independent phase models, or zeta zeros into Collatz.

There is a precise small common algebraic component. The mod-4 prime calculation uses `chi_4(p^k)=chi_4(p)^k` in the two-element group `{+1,-1}`. For `U_epsilon(n)=(3n+epsilon)/2^r`, where `epsilon` is `+1` or `-1` and `r=v_2(3n+epsilon)`, reduction modulo 3 gives `U_epsilon(n)=epsilon*(-1)^r mod3`. Both identities retain exponent parity through a two-element character. The objects, admissible-source fibers, and sampling measures differ; identifying this character does not identify the distributions or global dynamics.

The companion [signed Collatz audit](signed_collatz_pairing_and_exponent_bias_20261014.md) develops this source-labelled valuation calculation separately; its density statements use their explicitly defined arithmetic source measures.

## 11. A stronger formal bridge: primitive periodic orbits

There is an exact analogue of the Euler-product mechanism for a finite directed graph with adjacency matrix `A` (nonnegative integral entries, with parallel edges distinguished). Let `P_d` count cyclic-equivalence classes of primitive closed **edge-walks** of length `d`: a primitive edge-word is not a repetition of a shorter closed edge-word. These walks may revisit vertices and edges; they are not necessarily simple cycles.

Each primitive orbit of length `d|n` contributes exactly `d` based length-`n` closed walks. Therefore

`tr(A^n) = sum_(d|n) d P_d`,

`P_n = (1/n) sum_(d|n) mu(d) tr(A^(n/d))`.

As formal power series,

`-log det(I-zA) = sum_(n>=1) tr(A^n) z^n/n`

`                 = sum_(primitive gamma) sum_(k>=1) z^(k length(gamma))/k`.

Equivalently,

`det(I-zA)^(-1) = product_(primitive gamma) (1-z^(length(gamma)))^(-1)`.

This is a genuine identity, with the same repeat weight `1/k`, and can be proved just by the orbit decomposition above and the formal trace-log identity. It provides a usable bridge if a corpus observable actually counts closed walks or primitive orbits. It does not immediately apply to Hamiltonian paths, simple cycles, first-hit/stopped trajectories, or an arbitrary response quotient. Those require separate source-preserving constructions. A finite graph's eigenvalues also need not reproduce the zeta/L zero spectrum or its limiting distribution.

This is standard dynamical-zeta theory, not a new theorem of this note. The finite-type shift determinant formula appears in [Bowen–Lanford, *Zeta Functions of Restrictions of the Shift Transformation* (1970), pp. 43–49](https://people.math.harvard.edu/~knill/history/lanford/papers/BowenLanford.pdf). The present role is to identify a correctly typed candidate bridge to the corpus.

There is also a weight-preserving formal version. Assign each directed edge a commuting weight `w_e` in a commutative algebra over the rationals, set `W_ij=sum_(e:i->j)w_e`, and put `w_gamma=product_(e in gamma)w_e`, retaining repetitions. Then, as formal power series in `z`,

`-log det(I-zW) = sum_(n>=1) tr(W^n) z^n/n`

`               = sum_(primitive gamma) sum_(k>=1) w_gamma^k z^(k length(gamma))/k`,

`det(I-zW)^(-1) = product_(primitive gamma) (1-w_gamma z^(length(gamma)))^(-1)`.

No analytic convergence or positivity is needed for the formal identity. The `k`th repeat carries `w_gamma^k`, not `w_gamma`. The unweighted scalar Möbius formula for `P_n` must therefore not be naively applied to weighted traces while discarding this exponent information. Signed weights additionally prevent an automatic probabilistic interpretation.

## 12. Verification and status boundaries

The companion `04-computation/experiments/prime_power_bias_exact_check_20261014.py` uses only the Python standard library. It checks the exact mod-4 decomposition at every integer from 1 through 100,000 using integer-scaled rational weights, and independently locates the first strict reversal. Checks use explicit exceptions and remain active under Python optimization.

- **Proved directly here:** finite prime-power decompositions, Möbius inversion identities, root-multiplicity correction (using the stated PNT input), the unit-weight warning, and formal primitive-orbit identities.
- **Established literature:** Littlewood oscillation, classical Skewes bounds, Mertens disproof.
- **Conditional literature:** GRH/LI logarithmic race densities and the stated centered random-phase interpretation.
- **Finite checked:** first reversal and 100,000 exact decomposition instances.
- **Research analogy only:** transfer to radial multiplicity, sidecars, Collatz, or geometric current problems absent a specified compatible map.

The note establishes a precise common analytical lesson and a concrete route to test it. It does not establish a universal finite/infinite phase boundary or solve any repository frontier problem.


Companion exact output: `05-knowledge/results/prime_power_bias_exact_check_20261014.json`. Reproduce from the repository root with `python3 -B 04-computation/experiments/prime_power_bias_exact_check_20261014.py`; the same command with `-O` produces identical output.

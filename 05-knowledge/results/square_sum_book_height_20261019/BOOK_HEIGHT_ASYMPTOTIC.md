# Sharp height of induced books in the square-sum graph

**Status:** proved, with an explicit construction and a separate exact checker. All asymptotic estimates below are established elementarily. Classical sources are recorded for attribution and context, rather than used as unproved dependencies.

## 1. Definitions and results

The square-sum graph on the positive integers has an edge between **distinct** labels `a,b` exactly when `a+b` is an integer square. Its restriction to `[1,N]` is denoted `G_N`.

For an integer `m >= 1`, the book `B_m = K_2 ∨ I_m` consists of an adjacent pair of anchors and `m` independent vertices, each adjacent to both anchors. It has `m+2` vertices and is 2-degenerate: remove the degree-two pages first. Let

- `H(m)` be the least `N` such that `G_N` contains an **induced** copy of `B_m`;
- `b(N)` be the maximum number of pages in an induced book in `G_N` (zero if there is none).

All logarithms in asymptotic statements are natural.

**Theorem.**

\[
 \log H(m)=\left(\frac1{\log2}+o(1)\right)\log m\,\log\log m
 \qquad(m\to\infty).
\]

Equivalently,

\[
 \log b(N)=(\log2+o(1))\frac{\log N}{\log\log N}
 \qquad(N\to\infty).
\]

In particular `H(m)` grows faster than every fixed power of `m`. Consequently there is **no uniform polynomial height bound for all finite 2-degenerate graphs**, even if inducedness is dropped from the proposed polynomial-height assertion. This does not contradict qualitative induced realizability of every finite 2-degenerate graph, the near-quadratic induced-forest result, or the quartic simultaneous one-subdivision result. It does not identify the optimal worst-case height for the whole 2-degenerate class.

The proof connects common-neighbor capacity with the divisor function; an explicit residue class makes a constant fraction of the relevant divisors into mutually nonadjacent pages.

## 2. Exact common-neighbor parametrization and parity

Let `A<B` be positive integer labels, put `d=B-A`, and consider a common neighbor `X` satisfying `1 <= X <= N` and `X` different from both anchors. Write

\[
 X+A=r^2,\qquad X+B=s^2,
 \qquad 0<r<s.
\]

Then

\[
 u=s-r,\quad v=s+r,\quad uv=d,\quad 0<u<v,\quad u\equiv v\pmod2,
\]

and conversely every such factor pair gives

\[
 X=\left(\frac{v-u}{2}\right)^2-A.
\]

Keep precisely the pairs for which `1 <= X <= N` and `X` is not an anchor. This is a bijection: the positive roots recover the factor pair uniquely. Thus it is an exact formula, including positivity, the label ceiling, and the no-loop convention.

Writing `tau(t)` for the number of positive divisors of `t`, the number of common neighbors is at most

\[
 C(d)=
 \begin{cases}
 \lfloor\tau(d)/2\rfloor,&d\text{ odd},\\
 0,&d\equiv2\pmod4,\\
 \lfloor\tau(d/4)/2\rfloor,&4\mid d.
 \end{cases}
\]

For odd `d`, all factor pairs have matching odd parity. For `4|d`, matching parity forces `u=2u'`, `v=2v'`, so `u'v'=d/4`. The strict inequality discards the diagonal factorization when the integer being factored is a square. If `d=2 mod 4`, matching parity is impossible. Positivity and all other restrictions can only reduce these upper bounds.

Let `D(x)=max_{1<=t<=floor(x)} tau(t)` for `x>=1`. Since `d<=N-1`, every book with labels at most `N`, induced or otherwise, satisfies

\[
 m\le C(d)\le\tfrac12D(N-1).                 \tag{2.1}
\]

The extra requirement that the anchors themselves are adjacent does not weaken this necessary bound.

## 3. An elementary uniform divisor upper bound

We include the direction of Wigert's maximal-order estimate needed here. Set `L=log N` and, for sufficiently large `N`,

\[
 y=\frac{L}{(\log L)^3}>1.
\]

If `t<=N` and `t=prod_p p^{a_p}`, split `log tau(t)=sum_p log(a_p+1)` at `p=y`.

For small primes, there are at most `y` possible primes and each exponent is at most `L/log2`. Their total contribution is at most

\[
 y\log(1+L/\log2)
 =O\left(\frac{L}{(\log L)^2}\right)
 =o\left(\frac{L}{\log L}\right).
\]

For large primes, `a+1<=2^a` for every nonnegative integer `a`, and

\[
 \sum_{p>y}a_p\le\frac{\log t}{\log y}\le\frac L{\log y}.
\]

Their contribution is therefore at most

\[
 (\log2)\frac L{\log y}
 =(\log2+o(1))\frac L{\log L}.
\]

Uniformly for `t<=N`,

\[
 \log D(N)\le(\log2+o(1))\frac{\log N}{\log\log N}.       \tag{3.1}
\]

No prime-distribution theorem enters this upper bound.

To invert (2.1), put `X=log H(m)` and `Y=log m`. Distinctness gives `H(m)>=m+2`, hence `X>=Y`. Equations (2.1) and (3.1) give

\[
 Y\le(\log2+o(1))\frac X{\log X},
\]

so

\[
 X\ge\left(\frac1{\log2}-o(1)\right)Y\log Y.             \tag{3.2}
\]

This already rules out every fixed polynomial height bound. In particular, for every fixed `C>0`, `H(m)/m^C` tends to infinity.

## 4. An explicit induced construction for every number of pages

Let `k>=4`, and let `p_1=5,p_2=7,p_3=11,...,p_k=P` be the first `k` primes at least 5. Put

\[
 Q=\prod_{i=1}^k p_i,\qquad
 h=\min\{j\in\mathbb Z_{>0}:j^4\ge k^3\}
   =\lceil k^{3/4}\rceil.
\]

Thus `Q` is squarefree, is divisible by 5, and is coprime to 6; also `Q>=5005`.

Let `S` be the least odd positive integer not divisible by 3 for which `S>3 sqrt(Q)`. Define the two anchors

\[
 U=\frac{S^2-9Q}{2},\qquad V=U+9Q.                       \tag{4.1}
\]

They are positive distinct integers and `U+V=S^2`. The admissible integers for `S` have residues 1 or 5 modulo 6, so a deliberately loose bound is

\[
 3\sqrt Q<S\le3\sqrt Q+6.
\]

Consequently

\[
 0<U\le18\sqrt Q+18,\qquad
 V\le9Q+18\sqrt Q+18\le10Q.                              \tag{4.2}
\]

The last inequality holds for every `Q>=5005`.

Retain each divisor `a|Q` satisfying the two **exact integer inequalities**

\[
 4a^2<Q,\qquad a^2P^{2h}\ge Q.                           \tag{4.3}
\]

For such a divisor, let `c=Q/a` and define

\[
 r_a=\frac{3(c-a)}2,\qquad s_a=\frac{3(c+a)}2,\qquad
 X_a=r_a^2-U.                                            \tag{4.4}
\]

These are integers because `a,c` are odd, and

\[
 X_a+U=r_a^2,\qquad X_a+V=s_a^2.                         \tag{4.5}
\]

### Positivity and injectivity

The first inequality in (4.3) gives `a<sqrt(Q)/2` and `c>2sqrt(Q)`, hence

\[
 r_a>\frac94\sqrt Q,\qquad r_a^2>\frac{81}{16}Q.
\]

For `Q>=5005`, the latter exceeds `18sqrt(Q)+18`, so `X_a>0` by (4.2).

For positive `a<sqrt(Q)`, the function `Q/a-a` is strictly decreasing. Thus distinct retained divisors give distinct positive `r_a`, and hence distinct page labels.

### Inducedness and anchor-page separation

Because `3` divides `9Q`, whereas `S` is not divisible by 3,

\[
 U\equiv V\equiv2\pmod3.
\]

Both roots in (4.4) are divisible by 3, so every page satisfies

\[
 X_a\equiv1\pmod3.
\]

Therefore pages cannot coincide with anchors, and the sum of two pages is 2 modulo 3, which is not a square. Together with (4.5) and the anchor edge, these facts prove that **all retained pages simultaneously induce a book**.

### Height

The second inequality in (4.3) gives `c<=sqrt(Q)P^h`. Hence

\[
 X_a<r_a^2<\frac94c^2\le\frac94QP^{2h}.
\]

By (4.2), all anchors and pages lie below the explicit integer ceiling

\[
 N_k=10QP^{2h}.                                          \tag{4.6}
\]

## 5. Counting retained divisors

There are `2^k` divisors of `Q`.

First pair each divisor `d` not divisible by 5 with the divisor `5d`. The central closed interval

\[
 [\sqrt Q/2,\ 2\sqrt Q]
\]

has endpoint ratio 4, while the ratio in each pair is 5. At most one member of each pair lies in the interval. Thus at least half of all divisors lie outside it. Reciprocal pairing `d -> Q/d` is a bijection between the lower and upper outside portions. At least `2^(k-2)` divisors therefore satisfy `4d^2<Q`.

Now choose a divisor `A` uniformly at random. Equivalently,

\[
 \log A=\sum_{i=1}^k B_i\log p_i,
\]

where the `B_i` are independent Bernoulli variables of mean `1/2`. Thus

\[
 \mathbb E\log A=\tfrac12\log Q,\qquad
 \operatorname{Var}(\log A)=\tfrac14\sum_i(\log p_i)^2
 \le\tfrac k4(\log P)^2.
\]

Divisors failing the second inequality of (4.3) satisfy

\[
 \log A<\tfrac12\log Q-h\log P.
\]

By Chebyshev's inequality their fraction is at most `k/(4h^2)`. Since `h^2>=k^(3/2)>=2k` for `k>=4`, at most `2^k/8` divisors are discarded by that test. Therefore the number satisfying both tests is at least

\[
 \left(\frac14-\frac{k}{4h^2}\right)2^k
 \ge2^{k-3}.                                             \tag{5.1}
\]

This finite bound is valid for every `k>=4`; no unspecified asymptotic threshold is needed.

Given any `m>=1`, choose

\[
 k=\max\{4,\ \lceil\log_2m\rceil+3\}.                    \tag{5.2}
\]

At least `m` retained divisors exist. Selecting any `m` of their page labels, together with `U,V`, yields an induced `B_m` below `N_k`. Formula (5.2) can be evaluated without floating point as `k=max(4,(m-1).bit_length()+3)`.

## 6. Elementary prime growth and the matching asymptotic

Only the elementary bound `p_j=O(j log j)` is needed. Here is a short proof sufficient for this purpose.

For every positive integer `t`,

\[
 \binom{2t}{t}\ge\frac{4^t}{2t+1}.
\]

For each prime `p<=2t`, the usual factorial valuation formula gives

\[
 v_p\binom{2t}{t}
 =\sum_{j\ge1}\left(\left\lfloor\frac{2t}{p^j}\right\rfloor
                  -2\left\lfloor\frac{t}{p^j}\right\rfloor\right)
 \le\frac{\log(2t)}{\log p},
\]

since each summand is 0 or 1 and vanishes for `p^j>2t`. Thus

\[
 2t\log2-\log(2t+1)
 \le\log\binom{2t}{t}
 \le\pi(2t)\log(2t).
\]

Taking the largest even integer at most a real `x` yields `pi(x)>=c x/log x` for some absolute `c>0` and all sufficiently large `x`. Substituting `x=Cj log j` with a fixed sufficiently large `C` gives at least `j` primes below `x`, proving `p_j=O(j log j)`.

Our `P` is the `(k+2)`nd prime overall. Therefore `P=O(k log k)`. Also `p_i>=i`, so `Q>=k!`. The elementary integral estimate `log(k!)=k log k+O(k)` now gives

\[
 k\log k+O(k)\le\log Q\le k\log P
             \le k\log k+O(k\log\log k),
\]

and hence `log Q=(1+o(1))k log k`. Since `h=O(k^(3/4))`,

\[
 \log N_k=\log10+\log Q+2h\log P
          =(1+o(1))k\log k.                              \tag{6.1}
\]

The choice (5.2) has `k=log m/log2+O(1)`. Consequently

\[
 \log H(m)\le\log N_k
   =\left(\frac1{\log2}+o(1)\right)\log m\,\log\log m.
\]

Together with (3.2), this proves the asserted sharp asymptotic.

## 7. The equivalent local graph-capacity statement

Equation (2.1) and the elementary divisor bound give

\[
 \log b(N)\le(\log2+o(1))\frac{\log N}{\log\log N}.
\]

For any fixed `epsilon` with `0<epsilon<log2`, set

\[
 m_N=\left\lfloor\exp\left((\log2-\epsilon)
                           \frac{\log N}{\log\log N}\right)\right\rfloor.
\]

Then `log log m_N ~ log log N`, and the proved height asymptotic implies

\[
 \log H(m_N)=(1-\epsilon/\log2+o(1))\log N<\log N
\]

for all sufficiently large `N`. Thus `b(N)>=m_N`. Letting `epsilon` decrease to zero proves the matching lower bound, and hence

\[
 \boxed{\log b(N)\sim(\log2)\frac{\log N}{\log\log N}}.
\]

This is a statement about actual induced books, not merely an upper bound on arbitrary codegrees.

## 8. Scope, provenance, sources, and verification

- The construction combines divisor concentration with a modulus-3 independence guard. The finite cutoff makes the construction valid for every number of pages. The arithmetic estimates are proved here, making the sharp asymptotic self-contained.
- The existing pair-extension packet proves qualitative induced realizability for every finite 2-degenerate graph. This addendum answers its quantitative obstruction separately. Previously frozen files are unchanged.
- The closest primary background source is A. Dujella and C. Elsholtz, *Sumsets being squares*, Acta Mathematica Hungarica 141 (2013), 353–357, introduction on pp. 353–354: [DOI](https://doi.org/10.1007/s10474-013-0334-8), [author-hosted full paper](https://www.math.tugraz.at/~elsholtz/WWW/papers/papers44.pdf). It explicitly gives the two-anchor difference-of-squares factorization and the same maximal-divisor upper bound for the number of common neighbors. Those mechanisms are classical background, not claimed as new here. The present conclusion additionally imposes the anchor edge and simultaneous independence of all pages while attaining the matching leading asymptotic.
- The maximal-order statement is classical, due to Wigert. A modern research source explicitly states the full equality `max_{n<=x} log tau(n)=(log2+o(1))log x/log log x`: Christian Elsholtz, Marc Technau, and Niclas Technau, *The maximal order of iterated multiplicative functions*, abstract and §1.1, [arXiv:1709.04799](https://arxiv.org/abs/1709.04799), [full PDF](https://arxiv.org/pdf/1709.04799). We inspected the full paper, not only search snippets. Our proof needs only the elementary upper direction, reproduced in §3.
- A sharper classical prime estimate is in J. Barkley Rosser and Lowell Schoenfeld, *Approximate formulas for some functions of prime numbers*, Illinois Journal of Mathematics 6 (1962), 64–94, Theorem 3 and its corollary, printed p. 69: [DOI](https://doi.org/10.1215/ijm/1255631807), [scan of the original paper](https://denisevellachemla.eu/Rosser-Schoenfeld-1962.pdf). The original paper was inspected. The weaker estimate needed here is proved in §6.
- These citations establish classical arithmetic background, not a literature-priority claim about the graph-theoretic conclusion. No exhaustive novelty claim is made.
- `verify_book_height.py` independently enumerates common neighbors and factor pairs on small domains, checks parity capacities, validates the exact construction and divisor counts, and exhaustively verifies induced adjacency for a finite family of constructed books. Its checks are explicit exceptions, so they remain active under Python `-O`. Exact computational tests support implementation correctness; the arguments above establish the unbounded results.
